import json
import os
import random
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from DeepOneCalculator._viennarna import fold_energy
from DeepOneCalculator.features import _energy

SEQUENCES = ["A" * 20, "GGGAAACCC", "GCGCUUCGCG", "GCGCGCGCGCGCGCGCGCGC",
             "ATGCTAGCTAGCTAGCATGC", "TTATCTTCGCTATCACCTCC"]


def test_recorded_cli_reference() -> None:
    reference = json.loads(Path(__file__).with_name("folding_reference.json").read_text())
    for case in reference["cases"]:
        assert fold_energy(case["sequence"]) == pytest.approx(case["energy"], abs=0.011)


@pytest.mark.parametrize("sequence", SEQUENCES)
def test_reference_python(sequence: str) -> None:
    RNA = pytest.importorskip("RNA")
    assert fold_energy(sequence) == pytest.approx(RNA.fold(sequence)[1], abs=1e-5)


@pytest.mark.cli
@pytest.mark.parametrize("sequence", SEQUENCES)
def test_reference_cli(sequence: str) -> None:
    cli = shutil.which("RNAfold")
    if cli is None:
        if os.environ.get("DEEPONE_REQUIRE_RNAFOLD") == "1":
            pytest.fail("RNAfold reference CLI is required in this test environment")
        pytest.skip("reference RNAfold is not installed (not needed at runtime)")
    result = subprocess.run([cli, "--noPS"], input=sequence + "\n", text=True,
                            capture_output=True, check=True)
    matches = re.findall(r"\(\s*(-?\d+\.\d+)\s*\)", result.stdout)
    assert len(matches) == 1, result.stdout
    assert fold_energy(sequence) == pytest.approx(float(matches[0]), abs=0.011)


@pytest.mark.cli
def test_reference_cli_varied_sequences() -> None:
    cli = shutil.which("RNAfold")
    if cli is None:
        if os.environ.get("DEEPONE_REQUIRE_RNAFOLD") == "1":
            pytest.fail("RNAfold reference CLI is required in this test environment")
        pytest.skip("reference RNAfold is not installed")
    rng = random.Random(2345)
    sequences = ["".join(rng.choices("ACGU", k=length))
                 for length in (1, 4, 10, 20, 30, 60, 100, 200) for _ in range(8)]
    result = subprocess.run([cli, "--noPS"], input="\n".join(sequences) + "\n",
                            text=True, capture_output=True, check=True)
    energies = re.findall(r"\(\s*(-?\d+\.\d+)\s*\)", result.stdout)
    assert len(energies) == len(sequences)
    for sequence, energy in zip(sequences, energies, strict=True):
        assert fold_energy(sequence) == pytest.approx(float(energy), abs=0.011), sequence


def test_input_and_thread_safety() -> None:
    assert fold_energy("gggaaaccc") == fold_energy("GGGAAACCC")
    assert fold_energy("GCGCTTCGCG") == fold_energy("GCGCUUCGCG")
    for invalid in ("", "ACGN", "A C", "A\x00C"):
        with pytest.raises(ValueError):
            fold_energy(invalid)
    expected = [fold_energy(sequence) for sequence in SEQUENCES * 4]
    with ThreadPoolExecutor(max_workers=4) as executor:
        assert list(executor.map(fold_energy, SEQUENCES * 4)) == expected


@pytest.mark.parametrize(("spacer", "reference"), [
    ("ATATATATATATATATATAT", 11.078522641620026),
    ("GCGCGCGCGCGCGCGCGCGC", 34.74726040124747),
    ("CTTCGCTATCACCTCCGCCG", 23.688161588368516),
])
def test_original_nonfolding_energy(spacer: str, reference: float) -> None:
    # Golden values from utils/Energy_cal.py with grna_folding=False.
    expected = reference + fold_energy(spacer)
    assert _energy(spacer + "AGG") == pytest.approx(expected, abs=1e-10)
    assert _energy(spacer + "GGG") == pytest.approx(expected * 1.1, abs=1e-10)
