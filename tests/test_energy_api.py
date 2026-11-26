import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from DeepOneCalculator import calculate_free_energy, fold_energy


def test_public_energy_against_upstream_reference() -> None:
    reference = json.loads(Path(__file__).with_name("variant_reference.json").read_text())
    actual = [calculate_free_energy(context[4:27]) for context in reference["contexts"]]
    np.testing.assert_allclose(actual, reference["energies"], rtol=1e-7, atol=1e-7)
    guide = reference["contexts"][0][4:27]
    assert calculate_free_energy(guide.lower()) == calculate_free_energy(guide)
    assert isinstance(actual[0], float)


def test_fold_energy_remains_separate() -> None:
    spacer = "AT" * 10
    assert fold_energy(spacer) == pytest.approx(-2.5)
    expected = 11.078522641620026 - 2.5
    assert calculate_free_energy(spacer + "AGG") == pytest.approx(expected)
    assert calculate_free_energy(spacer + "AAG") == pytest.approx(expected * 0.9)


@pytest.mark.parametrize("sequence", ["", "A" * 20, "A" * 30, "N" * 20 + "AGG",
                                       "U" * 20 + "AGG", "AT" * 10 + "AAA"])
def test_invalid_energy_input(sequence: str) -> None:
    with pytest.raises(ValueError):
        calculate_free_energy(sequence)


@pytest.mark.parametrize("value", [None, ["AT" * 10 + "AGG"], 23])
def test_invalid_energy_type(value: object) -> None:
    with pytest.raises(TypeError):
        calculate_free_energy(value)  # type: ignore[arg-type]


def test_energy_does_not_load_tensorflow() -> None:
    result = subprocess.run(
        [sys.executable, "-c", "import sys; from DeepOneCalculator import "
         "calculate_free_energy; calculate_free_energy('AT' * 10 + 'AGG'); "
         "assert 'tensorflow' not in sys.modules"], capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
