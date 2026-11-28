import hashlib
import json
from importlib.resources import files
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

from DeepOneCalculator import calculate_crispron
from DeepOneCalculator.crispron import _build_model, _encode, _manifest, _predict
from DeepOneCalculator.features import _find_guides, calculate_crisproff

TARGET = "AAAA" + "ATATATATATATATATATAT" + "AGG" + "AAA"
COLUMNS = ["sequence_id", "sequence", "gRNA_id", "gRNA_seq", "crispron_score"]


@pytest.mark.parametrize("sequences", [TARGET, [TARGET, TARGET.lower()]])
def test_interface(sequences: str | list[str]) -> None:
    count = 1 if isinstance(sequences, str) else len(sequences)
    with patch("DeepOneCalculator.crispron._predict", return_value=np.arange(count)), \
         patch("DeepOneCalculator.crispron._resolve_device", return_value="/CPU:0"):
        result = calculate_crispron(sequences, batch_size=1)
    assert list(result.columns) == COLUMNS
    assert len(result) == count
    assert result.sequence.eq(TARGET).all()
    assert result.gRNA_seq.eq(TARGET[4:27]).all()
    assert result.sequence_id.tolist() == [f"sequence_{i}" for i in range(1, count + 1)]
    assert result.crispron_score.dtype == np.float64


def test_encoding_and_features() -> None:
    np.testing.assert_array_equal(_encode(["ACGT" * 7 + "AC"])[0, :4], np.eye(4))
    assert _encode([TARGET]).dtype == np.float32
    contexts = [TARGET[:24] + pam + TARGET[27:] for pam in ("AGG", "CGG", "GGG", "TGG")]
    with patch("DeepOneCalculator.crispron._predict", return_value=np.zeros(4)) as predict, \
         patch("DeepOneCalculator.crispron._resolve_device", return_value="/CPU:0"):
        calculate_crispron(contexts, device="cpu", batch_size=2)
    encoded, energies, device, batch_size = predict.call_args.args
    assert encoded.shape == (4, 30, 4)
    assert energies.shape == (4, 1)
    assert energies.dtype == np.float32
    assert device == "/CPU:0" and batch_size == 2
    expected = calculate_crisproff(TARGET[4:27], pam_model="crisproff")
    np.testing.assert_allclose(energies, expected)
    assert not np.isclose(calculate_crisproff(contexts[2][4:27]), expected)


def test_reverse_and_multiple_windows() -> None:
    reverse = TARGET.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    sequence = TARGET + reverse
    guides = _find_guides(sequence, "sequence_1")
    with patch("DeepOneCalculator.crispron._predict", return_value=np.arange(len(guides))), \
         patch("DeepOneCalculator.crispron._resolve_device", return_value="/CPU:0"):
        result = calculate_crispron(sequence)
    assert result.gRNA_id.tolist() == [g.grna_id for g in guides]
    assert result.gRNA_id.str.contains("_p").any()
    assert result.gRNA_id.str.contains("_m").any()
    assert "sequence_1_p1" in set(result.gRNA_id)
    assert "sequence_1_m1" in set(result.gRNA_id)


@pytest.mark.parametrize("sequences", [[], "A" * 30, "AAAA" + "AT" * 10 + "AGA" + "AAA"])
def test_empty_without_inference(sequences: str | list[str]) -> None:
    with patch("DeepOneCalculator.crispron._predict") as predict, \
         patch("DeepOneCalculator.crispron._resolve_device") as resolve:
        result = calculate_crispron(sequences)
    assert result.empty
    assert list(result.columns) == COLUMNS
    assert result.crispron_score.dtype == np.float64
    predict.assert_not_called()
    resolve.assert_not_called()


@pytest.mark.parametrize("sequences", [None, (), "A" * 29, "N" * 30, "U" * 30, [TARGET, 2]])
def test_invalid_sequences(sequences: object) -> None:
    with pytest.raises(ValueError):
        calculate_crispron(sequences)  # type: ignore[arg-type]


@pytest.mark.parametrize("kwargs", [
    {"batch_size": 0}, {"batch_size": -1}, {"batch_size": True}, {"batch_size": 1.5},
    {"device": "gpu"}, {"device": None},
])
def test_invalid_options(kwargs: dict) -> None:
    with pytest.raises(ValueError):
        calculate_crispron(TARGET, **kwargs)


def test_explicit_cuda_error() -> None:
    with patch("DeepOneCalculator.crispron._resolve_device", side_effect=RuntimeError("CUDA")):
        with pytest.raises(RuntimeError, match="CUDA"):
            calculate_crispron(TARGET, device="cuda")


def test_packaged_weights() -> None:
    manifest = _manifest()
    assert len(manifest) == 6
    assert len({entry["file"] for entry in manifest}) == 6
    resource = files("DeepOneCalculator").joinpath("models/crispron")
    for entry in manifest:
        data = resource.joinpath(entry["file"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == entry["sha256"]
    for name in ("LICENSE", "COPYING", "AUTHORS", "NOTICE", "provenance.json"):
        assert resource.joinpath(name).read_text()


@pytest.mark.integration
def test_original_models_and_ensemble() -> None:
    reference = json.loads(Path(__file__).with_name("crispron_reference.json").read_text())
    encoded = _encode(reference["contexts"])
    energies = np.array(reference["energies"], dtype=np.float32).reshape(-1, 1)
    model = _build_model()
    for entry, expected in zip(_manifest(), reference["scores"], strict=True):
        resource = files("DeepOneCalculator").joinpath("models/crispron", entry["file"])
        model.load_weights(str(resource))
        actual = np.asarray(model([encoded, energies], training=False)).ravel()
        np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-5)
    expected_mean = np.mean(reference["scores"], axis=0)
    for batch_size in (1, 256):
        np.testing.assert_allclose(_predict(encoded, energies, "/CPU:0", batch_size),
                                   expected_mean, rtol=1e-5, atol=1e-5)
    computed = [calculate_crisproff(c[4:27], pam_model="crisproff")
                for c in reference["contexts"]]
    np.testing.assert_allclose(computed, reference["energies"], rtol=1e-6, atol=1e-5)


@pytest.mark.integration
def test_end_to_end_batch_and_strands() -> None:
    reverse = TARGET.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    result = calculate_crispron([TARGET, reverse, TARGET.lower()], device="cpu", batch_size=1)
    assert result.gRNA_id.tolist() == ["sequence_1_p1", "sequence_2_m1", "sequence_3_p1"]
    assert np.isfinite(result.crispron_score).all()
    np.testing.assert_allclose(result.crispron_score, result.crispron_score.iloc[0])
