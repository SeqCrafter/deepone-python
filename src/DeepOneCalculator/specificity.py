"""Genome-wide CRISPRspec using RIsearch candidates and our energy model."""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field
from importlib import import_module
from os import PathLike, replace
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from types import MappingProxyType, ModuleType
from typing import TYPE_CHECKING, Literal

import numpy as np

from .features import calculate_crisproff

if TYPE_CHECKING:
    import polars as pl


@dataclass(frozen=True, slots=True)
class CRISPRspecConfig:
    """Candidate and partition-function settings for :func:`calculate_crisprspec`.

    Only ungapped 20 nt spacers with at most ``max_mismatches`` substitutions
    are included. ``pam_suffixes`` selects a subset of GG, AG and GA, ignoring
    the first PAM base. ``temperature_kelvin`` sets the Boltzmann factor;
    thermodynamic tables and RNA folding remain evaluated at 37 degrees C.

    ``pam_model='crisproff'`` uses the original NGG=1, NAG=.9, NGA=.8
    corrections. ``'deepone'`` instead uses the training feature's PAM table.
    Missing on-targets raise by default; ``require_on_target=False`` adds one
    theoretical on-target, as the original pipeline's ``--evaluate_all`` does.
    """

    max_mismatches: int = 6
    pam_suffixes: tuple[str, ...] = ("GG", "AG", "GA")
    ignored_chromosomes: frozenset[str] = frozenset()
    temperature_kelvin: float = 310.15
    threads: int | None = None
    require_on_target: bool = True
    pam_model: Literal["deepone", "crisproff"] = "crisproff"

    def __post_init__(self) -> None:
        if (isinstance(self.max_mismatches, bool)
                or not isinstance(self.max_mismatches, int)
                or not 0 <= self.max_mismatches <= 6):
            raise ValueError("max_mismatches must be an integer between 0 and 6")
        if (not self.pam_suffixes or isinstance(self.pam_suffixes, str)
                or any(pam not in {"GG", "AG", "GA"} for pam in self.pam_suffixes)):
            raise ValueError("pam_suffixes must be a nonempty subset of ('GG', 'AG', 'GA')")
        object.__setattr__(self, "pam_suffixes", tuple(dict.fromkeys(self.pam_suffixes)))
        if (isinstance(self.ignored_chromosomes, str)
                or any(not isinstance(name, str) for name in self.ignored_chromosomes)):
            raise ValueError("ignored_chromosomes must contain chromosome names")
        object.__setattr__(self, "ignored_chromosomes", frozenset(self.ignored_chromosomes))
        if (isinstance(self.temperature_kelvin, bool)
                or not math.isfinite(self.temperature_kelvin) or self.temperature_kelvin <= 0):
            raise ValueError("temperature_kelvin must be finite and positive")
        thermal_energy = 0.001987 * self.temperature_kelvin
        if thermal_energy == 0 or not math.isfinite(1.0 / thermal_energy):
            raise ValueError("temperature_kelvin is too small for a finite Boltzmann factor")
        if self.threads is not None and (
            isinstance(self.threads, bool) or not isinstance(self.threads, int) or self.threads < 1
        ):
            raise ValueError("threads must be a positive integer or None")
        if not isinstance(self.require_on_target, bool):
            raise ValueError("require_on_target must be a bool")
        if self.pam_model not in ("deepone", "crisproff"):
            raise ValueError("pam_model must be 'deepone' or 'crisproff'")


def _risearch() -> ModuleType:
    if sys.version_info < (3, 11):
        raise ImportError("CRISPRspec requires Python >=3.11 and DeepOneCalculator[specificity]")
    try:
        module = import_module("risearch")
    except ImportError as exc:
        raise ImportError(
            "Install DeepOneCalculator[specificity] to calculate CRISPRspec"
        ) from exc
    if not all(hasattr(module, name) for name in ("search", "index", "TargetRegistry")):
        raise ImportError("Install DeepOneCalculator[specificity]; RIsearch bindings are missing")
    return module


def _read_fasta(path: Path) -> Mapping[str, str]:
    sequences: dict[str, str] = {}
    alphabet = frozenset("ACGTRYSWKMBDHVN")
    parts: list[str] = []
    name: str | None = None
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if not parts:
                        raise ValueError(f"Empty reference sequence: {name}")
                    sequences[name] = "".join(parts)
                identifiers = line[1:].split()
                if not identifiers or identifiers[0] in sequences:
                    raise ValueError("Reference FASTA must have nonempty, unique sequence IDs")
                name, parts = identifiers[0], []
            else:
                sequence = line.upper()
                if name is None or set(sequence) - alphabet:
                    raise ValueError("Reference must be DNA FASTA (ACGT and IUPAC ambiguity bases)")
                parts.append(sequence)
    if name is None or not parts:
        raise ValueError("Reference FASTA must contain nonempty sequences")
    sequences[name] = "".join(parts)
    return MappingProxyType(sequences)


def _open_persistent_index(
    native: ModuleType, fasta: Path, destination: Path,
    sequences: Mapping[str, str], threads: int | None,
) -> object:
    import fcntl

    metadata_path = destination.with_name(destination.name + ".deepone.json")
    lock_path = destination.with_name(destination.name + ".lock")
    if fasta in {destination, metadata_path.resolve(), lock_path.resolve()}:
        raise ValueError("index_path and its sidecar files must not overwrite the reference FASTA")
    digest = hashlib.sha256()
    for name, sequence in sequences.items():
        digest.update(name.encode("utf-8") + b"\0")
        for offset in range(0, len(sequence), 1024 * 1024):
            digest.update(sequence[offset:offset + 1024 * 1024].encode("ascii"))
        digest.update(b"\0")
    fingerprint = digest.hexdigest()

    # Keep the native <index>.tmp staging file away from user files, and
    # serialize publication of the index and its source fingerprint.
    with lock_path.open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if destination.exists():
            status = destination.stat()
            if status.st_mtime_ns < fasta.stat().st_mtime_ns:
                raise ValueError("RIsearch index is older than the reference FASTA; rebuild it")
            try:
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                raise ValueError(
                    "Cannot verify RIsearch index provenance; use a new index_path to rebuild it"
                ) from exc
            expected = {
                "reference_sha256": fingerprint,
                "index_size": status.st_size,
                "index_mtime_ns": status.st_mtime_ns,
            }
            if metadata != expected:
                raise ValueError("RIsearch index fingerprint does not match; rebuild at a new path")
        else:
            with TemporaryDirectory(prefix=".deepone-index-", dir=destination.parent) as directory:
                staged = Path(directory) / "reference.idx"
                native.index(fasta, staged, threads=threads)
                status = staged.stat()
                metadata = {
                    "reference_sha256": fingerprint,
                    "index_size": status.st_size,
                    "index_mtime_ns": status.st_mtime_ns,
                }
                staged_metadata = Path(directory) / "reference.json"
                staged_metadata.write_text(json.dumps(metadata), encoding="utf-8")
                replace(staged, destination)
                replace(staged_metadata, metadata_path)
        targets: object = native.TargetRegistry.open(destination)
    return targets


@dataclass(frozen=True, slots=True)
class ReferenceGenome:
    """Reusable reference sequences and RIsearch's mmap-backed target index.

    Construct with :meth:`from_fasta`. DNA FASTA can be plain or gzip-compressed.
    Sequences are loaded once to extract PAMs from DataFrame coordinates; the
    current RIsearch Python bindings do not expose target bases. Hold this
    object across guide calls to reuse both sequences and the opened index.
    """

    fasta_path: Path
    index_path: Path
    _sequences: Mapping[str, str] = field(repr=False, compare=False)
    _targets: object = field(repr=False, compare=False)
    _temporary_index: TemporaryDirectory[str] | None = field(
        default=None, repr=False, compare=False,
    )

    @classmethod
    def from_fasta(
        cls,
        fasta: str | PathLike[str],
        *,
        index_path: str | PathLike[str] | None = None,
        threads: int | None = None,
    ) -> ReferenceGenome:
        """Open a reference, building its index if needed.

        With ``index_path=None``, build a temporary index owned by this object.
        Supply an index path to persist it and reuse it in future processes.
        Persisted indexes have a source fingerprint sidecar and a build lock.
        Existing indexes must have been built through this method from the
        same reference; unverified, changed or stale indexes are rejected.
        """
        native = _risearch()
        if threads is not None and (
            isinstance(threads, bool) or not isinstance(threads, int) or threads < 1
        ):
            raise ValueError("threads must be a positive integer or None")
        fasta_path = Path(fasta).resolve(strict=True)
        sequences = _read_fasta(fasta_path)
        temporary = TemporaryDirectory(prefix="deepone-reference-") if index_path is None else None
        if temporary is not None:
            destination = Path(temporary.name) / "reference.idx"
        else:
            assert index_path is not None
            destination = Path(index_path).resolve()
        try:
            if destination == fasta_path:
                raise ValueError("index_path must not overwrite the reference FASTA")
            if temporary is None:
                targets = _open_persistent_index(
                    native, fasta_path, destination, sequences, threads,
                )
            else:
                native.index(fasta_path, destination, threads=threads)
                targets = native.TargetRegistry.open(destination)
        except Exception:
            if temporary is not None:
                temporary.cleanup()
            raise
        return cls(fasta_path, destination, sequences, targets, temporary)


def _candidate_counts(
    hits: pl.DataFrame,
    reference: ReferenceGenome,
    guide: str,
    config: CRISPRspecConfig,
) -> Counter[str]:
    import polars as pl

    columns = {"target_name", "q_start", "q_end", "t_start", "t_end", "strand"}
    if not columns <= set(hits.columns):
        raise ValueError("RIsearch result is missing required coordinate columns")
    candidates = hits.with_columns(
        pl.col("target_name").str.strip_chars().str.extract(r"^(\S+)", 1).alias("target_name"),
    ).filter(
        (pl.col("q_start") == 0)
        & (pl.col("q_end") == 19)
        & (pl.col("t_end") == pl.col("t_start") + 19)
        & ~pl.col("target_name").is_in(list(config.ignored_chromosomes))
    ).select("target_name", "t_start", "t_end", "strand").unique()

    counts: Counter[str] = Counter()
    complement = str.maketrans("ACGT", "TGCA")
    spacer = guide[:20]
    alphabet = frozenset("ACGT")
    for chromosome, start, end, strand in candidates.iter_rows():
        if chromosome not in reference._sequences:
            raise ValueError(f"RIsearch index has an unknown reference sequence: {chromosome}")
        sequence = reference._sequences[chromosome]
        # RIsearch reports the physical interacting strand, opposite the
        # protospacer strand. Its coordinates are zero-based and inclusive.
        if strand == "-":
            if start < 0 or end + 4 > len(sequence):
                continue
            target = sequence[start:end + 4]
        elif strand == "+":
            if start < 3 or end >= len(sequence):
                continue
            target = sequence[start - 3:end + 1].translate(complement)[::-1]
        else:
            raise ValueError(f"Unexpected RIsearch strand: {strand}")
        if (len(target) != 23 or target[-2:] not in config.pam_suffixes
                or set(target) - alphabet):
            continue
        if sum(a != b for a, b in zip(spacer, target[:20], strict=True)) <= (
            config.max_mismatches
        ):
            counts[target] += 1
    return counts


def calculate_crisprspec(
    grna_seq: str,
    reference: ReferenceGenome | str | PathLike[str],
    *,
    config: CRISPRspecConfig | None = None,
    index_path: str | PathLike[str] | None = None,
) -> float:
    """Calculate CRISPRspec for a 23 nt guide (20 nt spacer + PAM).

    ``reference`` is a FASTA path or a reusable :class:`ReferenceGenome`.
    ``index_path`` optionally persists/reuses an index when passing a FASTA
    path; do not supply it with an already opened reference. Installing the
    ``specificity`` extra and Python >=3.11 are required for RIsearch.

    Search uses a full-length, non-wobble seed with no extension and a
    permissive energy cutoff. Candidates remain in a Polars DataFrame;
    RIsearch energies are ignored and guide-target pairs are rescored with
    :func:`calculate_crisproff`. No interaction TSV/gzip files are written
    or parsed. The current RIsearch API still needs one tiny query FASTA.

    Returns ``-log10(Z_off / (Z_off + Z_on))`` as a float, evaluated in the
    log domain. Higher values mean greater predicted specificity; no retained
    off-targets gives ``math.inf``. Exactly one matching genomic on-target
    is excluded; other copies still compete. Ignored chromosomes are excluded
    before selecting the on-target. Missing on-targets raise ``ValueError``
    unless disabled in the config. Distinct genomic positions retain their
    multiplicity, while duplicate search rows for one position are collapsed.

    Examples
    --------
    >>> reference = ReferenceGenome.from_fasta("genome.fa", index_path="genome.idx")
    >>> score = calculate_crisprspec("GGTGGACAAGCGGCAGATAGCGG", reference)
    """
    if config is None:
        config = CRISPRspecConfig()
    elif not isinstance(config, CRISPRspecConfig):
        raise TypeError("config must be a CRISPRspecConfig")
    on_energy = calculate_crisproff(grna_seq, pam_model=config.pam_model)
    guide = grna_seq.upper()
    if guide[-2:] not in config.pam_suffixes:
        raise ValueError("The guide PAM must be included in pam_suffixes")
    if isinstance(reference, ReferenceGenome):
        if index_path is not None:
            raise ValueError("index_path cannot be supplied with a ReferenceGenome")
    else:
        reference = ReferenceGenome.from_fasta(
            reference, index_path=index_path, threads=config.threads,
        )

    native = _risearch()
    with NamedTemporaryFile(mode="w", suffix=".fa", encoding="ascii") as query:
        query.write(f">guide\n{guide[:20]}\n")
        query.flush()
        hits = native.search(
            query.name, reference._targets,
            seed_start=1, seed_end=20, seed_length=20,
            mismatches=config.max_mismatches, mismatch_prefix=0, mismatch_suffix=0,
            seed_wobble=False, max_extension=0, energy_threshold=10000.0,
            no_max_prune=True, alignment=False, threads=config.threads,
        )
    counts = _candidate_counts(hits, reference, guide, config)
    if counts[guide]:
        counts[guide] -= 1
        if not counts[guide]:
            del counts[guide]
    elif config.require_on_target:
        raise ValueError("No exact on-target found on the included reference chromosomes")
    if not counts:
        return math.inf

    beta = 1.0 / (0.001987 * config.temperature_kelvin)
    log_weights = np.fromiter(
        (math.log(count) + beta * (
            calculate_crisproff(guide, target, pam_model=config.pam_model) - on_energy
        ) for target, count in counts.items()), dtype=np.float64,
    )
    log_ratio = np.logaddexp.reduce(log_weights, keepdims=True)[0]
    return float(np.logaddexp(0.0, -log_ratio) / math.log(10.0))
