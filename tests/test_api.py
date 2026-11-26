from unittest.mock import patch

import numpy as np
import pytest

from DeepOneCalculator import SUPPORTED_CELL_LINES, score_grnas
from DeepOneCalculator.features import _encode, _find_guides
from DeepOneCalculator.inference import _resolve_device

TARGET = "AAAA" + "ATATATATATATATATATAT" + "AGG" + "AAA"


@pytest.mark.parametrize("sequences", [TARGET, [TARGET, TARGET.lower()]])
@pytest.mark.parametrize("cells", ["K562", ["HEK293", "K562"], "all"])
def test_shapes(sequences: str | list[str], cells: str | list[str]) -> None:
    count = 1 if isinstance(sequences, str) else len(sequences)
    ncell = 5 if cells == "all" else 1 if isinstance(cells, str) else len(cells)
    with patch("DeepOneCalculator.api._predict", side_effect=lambda x, *args: np.arange(len(x))), \
         patch("DeepOneCalculator.api._resolve_device", return_value="/CPU:0"):
        frame = score_grnas(sequences, cells, batch_size=1)
    assert len(frame) == count * ncell
    assert list(frame.columns) == ["sequence_id", "sequence", "gRNA_id", "gRNA_seq",
                                   "cell_line", "deepone_score"]
    assert frame.gRNA_seq.eq(TARGET[4:27]).all()
    assert frame.sequence.eq(TARGET).all()
    assert frame.gRNA_id.nunique() == count


def test_reverse_complement_and_window_edges() -> None:
    reverse = TARGET.translate(str.maketrans("ATCG", "TAGC"))[::-1]
    guides = _find_guides(reverse, "sequence_1")
    assert len(guides) == 1
    assert guides[0].grna_id == "sequence_1_m1"
    assert guides[0].context == TARGET
    assert _encode([TARGET]).shape == (1, 30, 4)
    assert np.array_equal(_encode(["ATCG" * 7 + "AT"])[0, :4], np.eye(4))
    assert len(_find_guides(TARGET + "A", "s")) == 1


@pytest.mark.parametrize("sequences", ["A" * 29, "N" * 30, " " * 30, [TARGET, 2], None])
def test_bad_sequences(sequences: object) -> None:
    with pytest.raises(ValueError):
        score_grnas(sequences)  # type: ignore[arg-type]


@pytest.mark.parametrize("kwargs", [
    {"batch_size": 0}, {"batch_size": True}, {"batch_size": 1.5},
    {"device": "gpu"}, {"cell_lines": "HeLa"}, {"cell_lines": "HEK293T"},
    {"cell_lines": []}, {"cell_lines": ["all", "K562"]},
])
def test_bad_options(kwargs: dict) -> None:
    with pytest.raises(ValueError):
        score_grnas(TARGET, **kwargs)


def test_empty_and_no_pam() -> None:
    for sequences in ([], "A" * 30):
        result = score_grnas(sequences)
        assert result.empty
        assert len(result.columns) == 6
        assert result.deepone_score.dtype == np.float64


def test_cell_line_normalization() -> None:
    from DeepOneCalculator.resources import _normalize_cell_lines

    assert _normalize_cell_lines("ALL") == list(SUPPORTED_CELL_LINES)
    assert _normalize_cell_lines(["k562", "K562"]) == ["K562"]


def test_device_fallback() -> None:
    with patch("tensorflow.test.is_built_with_cuda", return_value=True), \
         patch("tensorflow.config.list_physical_devices", return_value=[]):
        assert _resolve_device("auto") == "/CPU:0"
        assert _resolve_device("cpu") == "/CPU:0"
        with pytest.raises(RuntimeError, match="CUDA"):
            _resolve_device("cuda")


def test_cuda_selection_and_non_cuda_plugin() -> None:
    with patch("tensorflow.config.list_physical_devices", return_value=[object()]):
        with patch("tensorflow.test.is_built_with_cuda", return_value=True):
            assert _resolve_device("auto") == "/GPU:0"
            assert _resolve_device("cuda") == "/GPU:0"
        with patch("tensorflow.test.is_built_with_cuda", return_value=False):
            assert _resolve_device("auto") == "/CPU:0"
