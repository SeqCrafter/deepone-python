"""Public scoring interface."""

from typing import cast

import numpy as np
import pandas as pd

from .features import _encode, _energy, _find_guides
from .inference import Device, _predict, _resolve_device
from .resources import _normalize_cell_lines

_COLUMNS = ["sequence_id", "sequence", "gRNA_id", "gRNA_seq", "cell_line", "deepone_score"]


def score_grnas(
    sequences: list[str] | str,
    cell_lines: list[str] | str = "all",
    batch_size: int = 256,
    device: Device = "auto",
) -> pd.DataFrame:
    """Score all SpCas9 NGG guides found on both strands of DNA inputs.

    Parameters
    ----------
    sequences : list of str or str
        DNA sequences of at least 30 bp (4 bp prefix, 20 nt spacer, 3 bp
        PAM, 3 bp suffix). Only A, C, G and T are accepted, case insensitive.
        Longer inputs are scanned for every complete 30 bp NGG context.
    cell_lines : list of str or str, default "all"
        One or more supported human cell lines; names are case insensitive.
        ``"all"`` selects HEK293, K562, HAP1, iPSC and RPE-1.
    batch_size : int, default 256
        Maximum number of guides per TensorFlow inference batch.
    device : {"auto", "cpu", "cuda"}, default "auto"
        Auto selects CUDA when available, otherwise CPU. Explicit CUDA
        raises RuntimeError when CUDA is unavailable.

    Returns
    -------
    pandas.DataFrame
        Columns: sequence_id, sequence, gRNA_id, gRNA_seq, cell_line,
        deepone_score. IDs are stable and one-based; gRNA_id contains p/m
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
    """
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("batch_size must be a positive integer")
    if device not in ("auto", "cpu", "cuda"):
        raise ValueError("device must be 'auto', 'cpu', or 'cuda'")
    selected = _normalize_cell_lines(cell_lines)
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
        guides.extend(_find_guides(normalized, f"sequence_{index}"))
    if not guides:
        return cast(pd.DataFrame, pd.DataFrame({
            name: pd.Series(dtype="float64" if name == "deepone_score" else "object")
            for name in _COLUMNS
        }))
    resolved = _resolve_device(device)
    encoded = _encode([guide.context for guide in guides])
    energies = np.array([[_energy(g.grna_seq)] for g in guides], dtype=np.float32)
    rows = []
    for cell_line in selected:
        predictions = _predict(encoded, energies, cell_line, resolved, batch_size)
        for guide, score in zip(guides, predictions, strict=True):
            rows.append((guide.sequence_id, guide.sequence, guide.grna_id, guide.grna_seq,
                         cell_line, float(score)))
    return cast(pd.DataFrame, pd.DataFrame(rows, columns=_COLUMNS))
