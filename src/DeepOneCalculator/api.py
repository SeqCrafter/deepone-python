"""Public scoring interface."""

import numpy as np
import pandas as pd

from .features import _encode, _energy, _find_guides
from .inference import Device, _predict, _resolve_device
from .resources import _model_specs, _normalize_cas9_types, _normalize_cell_lines

_COLUMNS = ["sequence_id", "sequence", "gRNA_id", "gRNA_seq", "cell_line", "cas9_type",
            "deepone_score"]


def calculate_deepone(
    sequences: list[str] | str,
    cell_lines: list[str] | str = "all",
    batch_size: int = 256,
    device: Device = "auto",
    *,
    cas9_types: list[str] | str = "SpCas9",
) -> pd.DataFrame:
    """Score SpCas9, SpCas9-NG or SpG guides on both strands of DNA inputs.

    Parameters
    ----------
    sequences : list of str or str
        DNA sequences of at least 30 bp (4 bp prefix, 20 nt spacer, 3 bp
        PAM, 3 bp suffix). Only A, C, G and T are accepted, case insensitive.
        Longer inputs are scanned for every complete 30 bp context: NGG
        for SpCas9, NGN for SpCas9-NG and SpG, matching the upstream CLI.
    cell_lines : list of str or str, default "all"
        One or more supported human cell lines; names are case insensitive.
        ``"all"`` selects HEK293, K562, HAP1, iPSC and RPE-1. This selection
        applies only to SpCas9. Variant models have no cell-specific weights
        and each is evaluated once regardless of the selected cell lines.
        Cell line names are validated even for variant-only requests.
    batch_size : int, default 256
        Maximum number of guides per TensorFlow inference batch.
    device : {"auto", "cpu", "cuda"}, default "auto"
        Auto selects CUDA when available, otherwise CPU. Explicit CUDA
        raises RuntimeError when CUDA is unavailable.
    cas9_types : list of str or str, default "SpCas9"
        Keyword-only Cas9 selection: SpCas9, SpCas9-NG, SpG or ``"all"``.
        Names are case insensitive; SpCas9_NG is accepted as an alias.
        Repeated selections are deduplicated, preserving their first order.

    Returns
    -------
    pandas.DataFrame
        Columns: sequence_id, sequence, gRNA_id, gRNA_seq, cell_line,
        cas9_type, deepone_score. Variant rows have a missing cell_line,
        indicating unspecified cell context, not validation in all cells.
        IDs are stable and one-based; gRNA_id contains p/m
        for strand and the oriented 30 bp window start. gRNA_seq is the
        20 nt spacer plus 3 bp PAM. Scores are the mean of five pretrained
        models, predicting templated +1 bp insertion frequency on the
        original training scale (not general editing efficiency).
        Inputs with no valid guide produce no rows. Empty input returns
        an empty frame with the same schema.

    Raises
    ------
    ValueError
        For invalid sequences, cell lines, batch size or device.
    RuntimeError
        When CUDA is explicitly requested and unavailable.

    Notes
    -----
    Output is grouped by requested Cas9 type and, for SpCas9, requested cell
    line. Each group preserves input order, followed by forward/reverse
    strand and ascending oriented window position. Input order is never
    shuffled; groups can repeat it. Duplicate sequences retain distinct IDs.
    """
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("batch_size must be a positive integer")
    if device not in ("auto", "cpu", "cuda"):
        raise ValueError("device must be 'auto', 'cpu', or 'cuda'")
    selected = _normalize_cell_lines(cell_lines)
    selected_cas9 = _normalize_cas9_types(cas9_types)
    allow_ngn = any(cas9 != "SpCas9" for cas9 in selected_cas9)
    values = [sequences] if isinstance(sequences, str) else sequences
    if not isinstance(values, list):
        raise ValueError("sequences must be a string or list of strings")
    guides = []
    for index, sequence in enumerate(values, 1):
        if not isinstance(sequence, str) or len(sequence) < 30:
            raise ValueError(f"sequence {index} must be a DNA string of at least 30 bp")
        normalized = sequence.upper()
        if set(normalized) - set("ATCG"):
            raise ValueError(f"sequence {index} contains non-ACGT characters")
        guides.extend(_find_guides(normalized, f"sequence_{index}", allow_ngn=allow_ngn))
    if not guides:
        empty: pd.DataFrame = pd.DataFrame({
            name: pd.Series(dtype="float64" if name == "deepone_score" else "object")
            for name in _COLUMNS
        })
        empty["cell_line"] = empty["cell_line"].astype("string")
        return empty
    resolved = _resolve_device(device)
    encoded = _encode([guide.context for guide in guides])
    energies = np.array([[_energy(g.grna_seq)] for g in guides], dtype=np.float32)
    rows = []
    for spec in _model_specs(selected_cas9, selected):
        indices = np.array([
            index for index, guide in enumerate(guides)
            if spec.cas9_type != "SpCas9" or guide.context[25:27] == "GG"
        ], dtype=np.intp)
        if not len(indices):
            continue
        predictions = _predict(encoded[indices], energies[indices], spec.weight_key,
                               resolved, batch_size)
        for index, score in zip(indices, predictions, strict=True):
            guide = guides[index]
            rows.append((guide.sequence_id, guide.sequence, guide.grna_id, guide.grna_seq,
                         spec.cell_line, spec.cas9_type, float(score)))
    scores: pd.DataFrame = pd.DataFrame(rows, columns=_COLUMNS)
    # Use pandas' nullable string dtype consistently in mixed and variant-only calls.
    scores["cell_line"] = scores["cell_line"].astype("string")
    return scores
