"""Check the native candidate search against an independent genome scan."""

import gzip
import math
import os
import subprocess
import sys
import threading
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import FrozenInstanceError
from pathlib import Path
from types import ModuleType

import pytest

from DeepOneCalculator import CRISPRspecConfig, ReferenceGenome, calculate_crisprspec
from DeepOneCalculator.features import calculate_free_energy
from DeepOneCalculator.specificity import _candidate_counts, _read_fasta

pl = pytest.importorskip("polars")
risearch = pytest.importorskip("risearch")
if not hasattr(risearch, "TargetRegistry"):
    pytest.skip("Install the specificity extra for native RIsearch tests", allow_module_level=True)

pytestmark = pytest.mark.specificity

GUIDE = "GGTGGACAAGCGGCAGATAGCGG"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(sequence: str, positions: tuple[int, ...]) -> str:
    values = list(sequence)
    for position in positions:
        values[position] = "ACGT"[("ACGT".index(values[position]) + 1) % 4]
    return "".join(values)


def write_reference(path: Path, sequences: dict[str, str]) -> Path:
    text = "".join(f">{name} description\n{sequence}\n" for name, sequence in sequences.items())
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


def brute_force_counts(
    sequences: dict[str, str], guide: str, config: CRISPRspecConfig,
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for name, sequence in sequences.items():
        if name in config.ignored_chromosomes:
            continue
        sequence = sequence.upper()
        for oriented in (sequence, reverse_complement(sequence)):
            for start in range(len(oriented) - 22):
                target = oriented[start:start + 23]
                if set(target) - set("ACGT") or target[-2:] not in config.pam_suffixes:
                    continue
                if sum(a != b for a, b in zip(guide[:20], target[:20], strict=True)) <= (
                    config.max_mismatches
                ):
                    counts[target] += 1
    return counts


def reference_score(guide: str, counts: Counter[str], config: CRISPRspecConfig) -> float:
    counts = counts.copy()
    if counts[guide]:
        counts[guide] -= 1
    beta = 1.0 / (0.001987 * config.temperature_kelvin)
    off = sum(
        count * math.exp(beta * calculate_free_energy(guide, target, pam_model=config.pam_model))
        for target, count in counts.items() if count
    )
    on = math.exp(beta * calculate_free_energy(guide, pam_model=config.pam_model))
    return -math.log10(off / (off + on)) if off else math.inf


@pytest.fixture
def genome_sequences() -> dict[str, str]:
    return {
        "chr1": "N" * 8 + ("N" * 7).join([
            GUIDE,
            mutate(GUIDE, (0,)),
            mutate(GUIDE, (19,)),
            mutate(GUIDE[:20], (0, 3, 7, 11, 16, 19)) + "AAG",
            mutate(GUIDE[:20], (3, 4, 5, 6, 7, 8)) + "TGA",
            mutate(GUIDE[:20], (0, 1, 2, 3, 4, 5, 6)) + "TGG",
            GUIDE[:20] + "AAA",
            GUIDE,
        ]) + "N" * 8,
        "chr2": reverse_complement(
            "N" * 5 + GUIDE + "N" * 9 + mutate(GUIDE, (8, 9)) + "N" * 5,
        ),
        "chr_lowercase": mutate(GUIDE, (5,)).lower(),
        "chr_no_pam": GUIDE[:20],
    }


@pytest.mark.parametrize("max_mismatches", [0, 1, 3, 6])
def test_native_candidates_match_brute_force(
    tmp_path: Path, genome_sequences: dict[str, str], max_mismatches: int,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", genome_sequences)
    reference = ReferenceGenome.from_fasta(fasta)
    config = CRISPRspecConfig(max_mismatches=max_mismatches, threads=1)
    expected = brute_force_counts(genome_sequences, GUIDE, config)
    query = tmp_path / "query.fa"
    query.write_text(f">guide\n{GUIDE[:20]}\n")
    hits = risearch.search(
        query, reference._targets, seed_start=1, seed_end=20, seed_length=20,
        mismatches=max_mismatches, mismatch_prefix=0, mismatch_suffix=0,
        max_extension=0, energy_threshold=10000, seed_wobble=False,
        no_max_prune=True, alignment=False, threads=1,
    )
    assert _candidate_counts(hits, reference, GUIDE, config) == expected


@pytest.mark.parametrize("max_mismatches", [0, 1, 3, 6])
def test_partition_score_matches_brute_force(
    tmp_path: Path, genome_sequences: dict[str, str], max_mismatches: int,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", genome_sequences)
    config = CRISPRspecConfig(max_mismatches=max_mismatches, threads=1)
    expected = brute_force_counts(genome_sequences, GUIDE, config)
    assert calculate_crisprspec(GUIDE.lower(), fasta, config=config) == pytest.approx(
        reference_score(GUIDE, expected, config), abs=1e-12,
    )


def test_duplicate_on_targets_retain_locus_multiplicity(tmp_path: Path) -> None:
    sequences = {"a": GUIDE + "N" * 10 + GUIDE, "b": reverse_complement(GUIDE)}
    fasta = write_reference(tmp_path / "repeated.fa", sequences)
    score = calculate_crisprspec(GUIDE, fasta, config=CRISPRspecConfig(max_mismatches=0))
    assert score == pytest.approx(-math.log10(2 / 3))


@pytest.mark.parametrize("compressed", [False, True])
def test_single_on_target_has_infinite_specificity(tmp_path: Path, compressed: bool) -> None:
    fasta = write_reference(tmp_path / ("single.fa.gz" if compressed else "single.fa"),
                            {"a": GUIDE})
    assert calculate_crisprspec(GUIDE, fasta) == math.inf


@pytest.mark.parametrize("pam_model", ["crisproff", "deepone"])
def test_pam_filters_ignored_chromosomes_and_energy_modes(
    tmp_path: Path, genome_sequences: dict[str, str], pam_model: str,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", genome_sequences)
    config = CRISPRspecConfig(
        pam_suffixes=("GG", "GA"), ignored_chromosomes=frozenset({"chr2"}),
        pam_model=pam_model,  # type: ignore[arg-type]
    )
    expected = brute_force_counts(genome_sequences, GUIDE, config)
    assert calculate_crisprspec(GUIDE, fasta, config=config) == pytest.approx(
        reference_score(GUIDE, expected, config), abs=1e-12,
    )


def test_missing_on_target_is_explicit(tmp_path: Path) -> None:
    sequences = {"off": mutate(GUIDE, (0, 19))}
    fasta = write_reference(tmp_path / "absent.fa", sequences)
    with pytest.raises(ValueError, match="No exact on-target"):
        calculate_crisprspec(GUIDE, fasta)
    config = CRISPRspecConfig(require_on_target=False)
    expected = brute_force_counts(sequences, GUIDE, config)
    assert calculate_crisprspec(GUIDE, fasta, config=config) == pytest.approx(
        reference_score(GUIDE, expected, config), abs=1e-12,
    )


def test_ignored_on_target_does_not_satisfy_requirement(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"on": GUIDE, "off": mutate(GUIDE, (4,))})
    with pytest.raises(ValueError, match="No exact on-target"):
        calculate_crisprspec(
            GUIDE, fasta, config=CRISPRspecConfig(ignored_chromosomes=frozenset({"on"})),
        )


def test_dataframe_energy_is_ignored_and_query_file_is_cleaned(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"a": GUIDE + "N" * 10 + mutate(GUIDE, (5,))})
    reference = ReferenceGenome.from_fasta(fasta)
    baseline = calculate_crisprspec(GUIDE, reference)
    native_search = risearch.search
    query_paths: list[Path] = []

    def search(query: str, targets: object, **kwargs: object) -> object:
        query_paths.append(Path(query))
        assert kwargs["mismatch_prefix"] == kwargs["mismatch_suffix"] == 0
        assert kwargs["max_extension"] == 0
        assert kwargs["alignment"] is False
        hits = native_search(query, targets, **kwargs)
        assert isinstance(hits, pl.DataFrame)
        return pl.concat([hits, hits]).with_columns(pl.lit(-1e12).alias("energy"))

    monkeypatch.setattr(risearch, "search", search)
    assert calculate_crisprspec(GUIDE, reference) == pytest.approx(baseline)
    assert query_paths and all(not path.exists() for path in query_paths)
    assert list(tmp_path.iterdir()) == [fasta]


def test_guide_folding_is_shared_across_off_targets(
    tmp_path: Path, genome_sequences: dict[str, str], monkeypatch: pytest.MonkeyPatch,
) -> None:
    from DeepOneCalculator import features

    cached = [value for value in vars(features).values() if hasattr(value, "cache_clear")]
    for function in cached:
        function.cache_clear()
    original_fold = features.fold_energy
    calls: list[str] = []

    def fold(sequence: str) -> float:
        calls.append(sequence)
        return original_fold(sequence)

    monkeypatch.setattr(features, "fold_energy", fold)
    fasta = write_reference(tmp_path / "genome.fa", genome_sequences)
    try:
        assert math.isfinite(calculate_crisprspec(GUIDE, fasta))
        assert calls == [GUIDE[:20]]
    finally:
        for function in cached:
            function.cache_clear()


def test_reuses_open_registry_and_persistent_index(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"a": GUIDE})
    index = tmp_path / "genome.idx"
    native_index = risearch.index
    calls = 0

    def build(*args: object, **kwargs: object) -> None:
        nonlocal calls
        calls += 1
        native_index(*args, **kwargs)

    monkeypatch.setattr(risearch, "index", build)
    reference = ReferenceGenome.from_fasta(fasta, index_path=index)
    assert calculate_crisprspec(GUIDE, reference) == math.inf
    assert calculate_crisprspec(GUIDE, reference) == math.inf
    assert calculate_crisprspec(GUIDE, fasta, index_path=index) == math.inf
    assert calls == 1
    with pytest.raises(ValueError, match="cannot be supplied"):
        calculate_crisprspec(GUIDE, reference, index_path=index)


def test_stale_index_and_overwrite_are_rejected(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"a": GUIDE})
    index = tmp_path / "genome.idx"
    ReferenceGenome.from_fasta(fasta, index_path=index)
    os.utime(index, ns=(1, 1))
    with pytest.raises(ValueError, match="older than"):
        ReferenceGenome.from_fasta(fasta, index_path=index)
    with pytest.raises(ValueError, match="overwrite"):
        ReferenceGenome.from_fasta(fasta, index_path=fasta)


def test_native_index_staging_does_not_overwrite_reference(tmp_path: Path) -> None:
    index = tmp_path / "genome.idx"
    fasta = write_reference(tmp_path / "genome.idx.tmp", {"a": GUIDE})
    original = fasta.read_bytes()
    reference = ReferenceGenome.from_fasta(fasta, index_path=index)
    assert fasta.read_bytes() == original
    assert calculate_crisprspec(GUIDE, reference) == math.inf


@pytest.mark.parametrize("suffix", [".deepone.json", ".lock"])
def test_index_sidecars_cannot_overwrite_reference(tmp_path: Path, suffix: str) -> None:
    index = tmp_path / "genome.idx"
    fasta = write_reference(tmp_path / (index.name + suffix), {"a": GUIDE})
    original = fasta.read_bytes()
    with pytest.raises(ValueError, match="overwrite"):
        ReferenceGenome.from_fasta(fasta, index_path=index)
    assert fasta.read_bytes() == original


def test_reference_fingerprint_detects_same_mtime_edit(tmp_path: Path) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"a": GUIDE})
    index = tmp_path / "genome.idx"
    ReferenceGenome.from_fasta(fasta, index_path=index)
    original = fasta.stat()
    write_reference(fasta, {"a": mutate(GUIDE, (0,))})
    os.utime(fasta, ns=(original.st_atime_ns, original.st_mtime_ns))
    with pytest.raises(ValueError, match="fingerprint"):
        ReferenceGenome.from_fasta(fasta, index_path=index)


def test_concurrent_index_builds_are_serialized(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"a": GUIDE})
    index = tmp_path / "genome.idx"
    native_index = risearch.index
    calls: list[Path] = []
    ready = threading.Barrier(2)

    def build(source: Path, destination: Path, **kwargs: object) -> None:
        calls.append(destination)
        native_index(source, destination, **kwargs)

    def open_reference() -> ReferenceGenome:
        ready.wait(timeout=10)
        return ReferenceGenome.from_fasta(fasta, index_path=index)

    monkeypatch.setattr(risearch, "index", build)
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(open_reference) for _ in range(2)]
        references = [future.result(timeout=30) for future in futures]
    assert len(calls) == 1
    assert all(reference.index_path == index for reference in references)


@pytest.mark.parametrize("contents", [
    "", "ACGT\n", ">a\n", ">a\nACGT\n>a\nTTTT\n", ">\nACGT\n", ">a\nAC-T\n",
])
def test_invalid_reference_fasta(tmp_path: Path, contents: str) -> None:
    fasta = tmp_path / "invalid.fa"
    fasta.write_text(contents)
    with pytest.raises(ValueError):
        _read_fasta(fasta)


@pytest.mark.parametrize("kwargs", [
    {"max_mismatches": -1}, {"max_mismatches": 7}, {"max_mismatches": True},
    {"max_mismatches": 1.5}, {"threads": 0}, {"threads": True}, {"threads": -1},
    {"temperature_kelvin": 0}, {"temperature_kelvin": math.inf},
    {"temperature_kelvin": math.nan}, {"temperature_kelvin": 1e-309},
    {"temperature_kelvin": 5e-324}, {"pam_suffixes": ()}, {"pam_suffixes": ("GC",)},
    {"pam_suffixes": "GG"}, {"ignored_chromosomes": "chr1"},
    {"require_on_target": 1}, {"pam_model": "unknown"},
])
def test_invalid_config(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        CRISPRspecConfig(**kwargs)  # type: ignore[arg-type]


def test_configuration_is_immutable() -> None:
    config = CRISPRspecConfig()
    with pytest.raises(FrozenInstanceError):
        config.max_mismatches = 1  # type: ignore[misc]


@pytest.mark.parametrize("guide", ["", "A" * 20, "N" * 20 + "AGG", GUIDE[:20] + "AAA"])
def test_invalid_guide_fails_before_loading_reference(guide: str) -> None:
    with pytest.raises(ValueError):
        calculate_crisprspec(guide, "does-not-exist.fa")


def test_guide_pam_must_be_included() -> None:
    with pytest.raises(ValueError, match="guide PAM"):
        calculate_crisprspec(
            GUIDE, "does-not-exist.fa", config=CRISPRspecConfig(pam_suffixes=("AG",)),
        )


def test_rejects_incorrect_config_type() -> None:
    with pytest.raises(TypeError, match="CRISPRspecConfig"):
        calculate_crisprspec(GUIDE, "does-not-exist.fa", config={})  # type: ignore[arg-type]


def test_same_span_opposite_strands_and_first_pam_base(tmp_path: Path) -> None:
    guide = "AT" * 10 + "CGG"
    other_pam = guide[:20] + "GGG"
    fasta = write_reference(tmp_path / "palindrome.fa", {
        "a": "CCG" + guide + "N" * 10 + other_pam,
    })
    config = CRISPRspecConfig(max_mismatches=0, temperature_kelvin=298.15, pam_model="deepone")
    beta = 1.0 / (0.001987 * config.temperature_kelvin)
    on = math.exp(beta * calculate_free_energy(guide))
    off = on + math.exp(beta * calculate_free_energy(guide, other_pam))
    assert calculate_crisprspec(guide, fasta, config=config) == pytest.approx(
        -math.log10(off / (on + off)), abs=1e-12,
    )


@pytest.mark.parametrize("off_energy", [-1000.0, 1000.0])
def test_log_domain_avoids_exponent_overflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, off_energy: float,
) -> None:
    fasta = write_reference(tmp_path / "genome.fa", {"a": GUIDE, "b": mutate(GUIDE, (4,))})

    def energy(guide: str, target: str | None = None, **kwargs: object) -> float:
        return off_energy if target is not None else 0.0

    monkeypatch.setattr("DeepOneCalculator.specificity.calculate_free_energy", energy)
    score = calculate_crisprspec(GUIDE, fasta)
    assert math.isfinite(score)
    if off_energy < 0:
        assert score == pytest.approx(1000 / (0.001987 * 310.15 * math.log(10)))
    else:
        assert score == 0.0


def test_energy_import_does_not_load_optional_dependencies() -> None:
    result = subprocess.run([
        sys.executable, "-c",
        "import sys; from DeepOneCalculator import calculate_free_energy; "
        "calculate_free_energy('AT' * 10 + 'AGG'); "
        "assert not {'risearch', 'polars', 'tensorflow'} & sys.modules.keys()",
    ], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_missing_native_bindings_has_actionable_error(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        "DeepOneCalculator.specificity.import_module", lambda name: ModuleType(name),
    )
    with pytest.raises(ImportError, match="specificity"):
        ReferenceGenome.from_fasta("unused.fa")
