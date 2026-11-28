"""Pairwise guide-target energies checked against the original CRISPRoff code.

The REFERENCE values below were produced by running
crisproff-1.1.2/CRISPRspec_CRISPRoff_pipeline.py (calcRNADNAenergy/get_eng
with pos_weight, pam_corr, grna_folding and dna_opening enabled) against
energy_dics.pkl, with DeepOneCalculator.fold_energy substituted for the
RNAfold subprocess. The runtime has no dependency on that checkout.

Each entry is (name, guide, target, core, crisproff, deepone): ``core`` is the
PAM-unaffected value (pam_corr=False), ``crisproff`` the original pam_corr=True
value (None where the original maps the PAM to its "OTHERS" 0.0 bucket, which
this API rejects with ValueError instead), and ``deepone`` core times the
DeepOne PAM ratio.
"""
import json
from pathlib import Path

import pytest

from DeepOneCalculator import calculate_crisproff, fold_energy
from DeepOneCalculator.features import _crisproff_parameters, _dna_opening, _pair_hybrid

GUIDE = "ATGCTAGCTAGCTAGCATGCAGG"
MM9 = "ATGCTAGCTCGCTAGCATGCAGG"
ISOLATED = "CCCGCCCGCACGCCCGCCCGAGG"
NOMATCH = "CCCGCCCGCCCGCCCGCCCGAGG"

REFERENCE = [
    ("self", GUIDE, GUIDE, 18.662363089630308, 18.662363089630308, 18.662363089630308),
    ("mm_p9", GUIDE, MM9, 11.06749521344138, 11.06749521344138, 11.06749521344138),
    ("mm_p0", GUIDE, "CTGCTAGCTAGCTAGCATGCAGG", 16.602760273557244, 16.602760273557244,
     16.602760273557244),
    ("mm_p19", GUIDE, "ATGCTAGCTAGCTAGCATGGAGG", 12.989752766448309, 12.989752766448309,
     12.989752766448309),
    ("mm_p1", GUIDE, "ACGCTAGCTAGCTAGCATGCAGG", 16.05895035872012, 16.05895035872012,
     16.05895035872012),
    ("mm_p18", GUIDE, "ATGCTAGCTAGCTAGCATCCAGG", 10.70823679671156, 10.70823679671156,
     10.70823679671156),
    ("bulge2_p910", GUIDE, "ATGCTAGCTCCCTAGCATGCAGG", 4.911635185697278, 4.911635185697278,
     4.911635185697278),
    ("bulge2_p01", GUIDE, "CAGCTAGCTAGCTAGCATGCAGG", 13.961260250157167, 13.961260250157167,
     13.961260250157167),
    ("loop3", GUIDE, "ATGCTAGCCCCCTAGCATGCAGG", 3.9337531716814738, 3.9337531716814738,
     3.9337531716814738),
    ("loop5", GUIDE, "ATGCTAGGCCCGTAGCATGCAGG", -9.230473124087126, -9.230473124087126,
     -9.230473124087126),
    ("loop11", GUIDE, "ATGCCCCGCCCGCCCCATGCAGG", -23.99647065029669, -23.99647065029669,
     -23.99647065029669),
    ("term3", GUIDE, "ATGCTAGCTAGCTAGCACCGAGG", 6.024225116620478, 6.024225116620478,
     6.024225116620478),
    ("lead3", GUIDE, "CCCCTAGCTAGCTAGCATGCAGG", 8.268516131325498, 8.268516131325498,
     8.268516131325498),
    ("mix", GUIDE, "ATGCTCCCCCGCTAGCATGCAGG", -1.207170207404758, -1.207170207404758,
     -1.207170207404758),
    ("isolated", GUIDE, ISOLATED, -45.220000190734865, -45.220000190734865,
     -45.220000190734865),
    ("nomatch", GUIDE, NOMATCH, -46.010000190734864, -46.010000190734864,
     -46.010000190734864),
    ("pamAAG", GUIDE, "ATGCTAGCTCGCTAGCATGCAAG", 11.06749521344138, 9.960745692097241,
     9.960745692097241),
    ("pamAGA", GUIDE, "ATGCTAGCTCGCTAGCATGCAGA", 11.06749521344138, 8.853996170753105,
     11.06749521344138),
    ("role", MM9, GUIDE, 11.943630588206144, 11.943630588206144, 11.943630588206144),
    ("self_g2", "GCGAAAGCCGTTAAGATCGTAAG", "GCGAAAGCCGTTAAGATCGTAAG", 18.60274517897955,
     16.742470661081594, 16.742470661081594),
    ("self_g3", "CCTGATGCTACAGCGGCTTCAGA", "CCTGATGCTACAGCGGCTTCAGA", 24.244923292579784,
     19.39593863406383, 24.244923292579784),
    ("self_g4", "ACCGAGTTCACTAAAAGTATGGA", "ACCGAGTTCACTAAAAGTATGGA", 18.945874696415498,
     15.1566997571324, 20.840462166057048),
]
CRISPROFF_SELF = {"self": 18.662363089630308, "self_g2": 16.742470661081594,
                  "self_g3": 19.39593863406383, "self_g4": 15.1566997571324}
DEEPONE_SELF = {"self": 18.662363089630304, "self_g2": 16.74247066108159,
                "self_g3": 24.244923292579788, "self_g4": 20.84046216605705}


@pytest.mark.parametrize(("name", "guide", "target", "core", "crisproff", "deepone"), REFERENCE)
def test_pairs_match_original_algorithm(name, guide, target, core, crisproff, deepone) -> None:
    assert calculate_crisproff(guide, target, pam_model="crisproff") == pytest.approx(
        crisproff, rel=1e-12, abs=1e-12)
    assert calculate_crisproff(guide, target) == pytest.approx(deepone, rel=1e-12, abs=1e-12)
    assert calculate_crisproff(guide, target, pam_model="deepone") == pytest.approx(
        deepone, rel=1e-12, abs=1e-12)
    assert isinstance(calculate_crisproff(guide, target), float)


@pytest.mark.parametrize(("pam", "deepone_value"), [
    ("AGC", 8.853996170753105),
    ("GGT", 13.280994256129654),
])
def test_crisproff_mode_rejects_other_pams(pam: str, deepone_value: float) -> None:
    target = "ATGCTAGCTCGCTAGCATGC" + pam
    # The original get_eng multiplies these PAMs by its "OTHERS" 0.0 ratio;
    # this API rejects them instead.
    with pytest.raises(ValueError, match="CRISPRoff"):
        calculate_crisproff(GUIDE, target, pam_model="crisproff")
    with pytest.raises(ValueError, match="CRISPRoff"):
        calculate_crisproff(target, target, pam_model="crisproff")
    assert calculate_crisproff(GUIDE, target) == pytest.approx(deepone_value, rel=1e-12,
                                                                 abs=1e-12)


@pytest.mark.parametrize(("name", "guide", "target", "core", "crisproff", "deepone"),
                         [entry for entry in REFERENCE if entry[0] in CRISPROFF_SELF])
def test_self_score_semantics(name, guide, target, core, crisproff, deepone) -> None:
    # Legacy one-argument behavior, exactly.
    assert calculate_crisproff(guide) == calculate_crisproff(guide, None)
    assert calculate_crisproff(guide, None) == calculate_crisproff(
        guide, None, pam_model="deepone")
    assert calculate_crisproff(guide) == DEEPONE_SELF[name]
    # 'crisproff' self: the original on-target score; None means target == guide.
    assert calculate_crisproff(guide, None, pam_model="crisproff") == pytest.approx(
        crisproff, rel=1e-12, abs=1e-12)
    assert calculate_crisproff(guide, guide, pam_model="crisproff") == pytest.approx(
        crisproff, rel=1e-12, abs=1e-12)
    # An explicit self target in 'deepone' mode uses the mismatch-aware hybrid
    # (per-base weights), not the historical dinucleotide-weighted feature.
    assert calculate_crisproff(guide, guide) == pytest.approx(deepone, rel=1e-12, abs=1e-12)


def test_mismatch_positions_clusters_and_terminal_handling() -> None:
    def value(target: str) -> float:
        return calculate_crisproff(GUIDE, target, pam_model="crisproff")

    # Terminal mismatches (outside the first-to-last matched base) carry less
    # penalty than internal ones; clusters deepen monotonically.
    assert value(GUIDE) > value("CTGCTAGCTAGCTAGCATGCAGG") > value(MM9)
    assert value(MM9) > value("ATGCTAGCTCCCTAGCATGCAGG")
    assert value("ATGCTAGCCCCCTAGCATGCAGG") > value("ATGCTAGGCCCGTAGCATGCAGG")
    assert value("ATGCTAGGCCCGTAGCATGCAGG") > value("ATGCCCCGCCCGCCCCATGCAGG")
    # A single isolated match (or none) leaves no hybrid energy at all: the
    # score is exactly folding plus target opening for both cases.
    for target in (ISOLATED, NOMATCH):
        assert _pair_hybrid(GUIDE[:20], target[:20]) == (0.0,) * 20
        assert calculate_crisproff(GUIDE, target, pam_model="crisproff") == (
            fold_energy(GUIDE[:20]) + _dna_opening(target[:20]))


def test_folding_uses_guide_and_opening_uses_target() -> None:
    target = "ATGCTAGCTCCCTAGCATGCAGG"
    core = 4.911635185697278  # reference value; PAM AGG has ratio 1.0
    assert fold_energy(GUIDE[:20]) != fold_energy(target[:20])
    value = calculate_crisproff(GUIDE, target, pam_model="crisproff")
    assert value == pytest.approx(core, rel=1e-12, abs=1e-12)
    # The original folds the guide spacer (never the off-target), and opens
    # the target DNA (never the guide).
    wrong_fold = core - fold_energy(GUIDE[:20]) + fold_energy(target[:20])
    wrong_opening = core - _dna_opening(target[:20]) + _dna_opening(GUIDE[:20])
    assert value != pytest.approx(wrong_fold, rel=1e-9)
    assert value != pytest.approx(wrong_opening, rel=1e-9)


def test_guide_and_target_roles_are_directed() -> None:
    assert calculate_crisproff(GUIDE, MM9, pam_model="crisproff") == pytest.approx(
        11.06749521344138, rel=1e-12, abs=1e-12)
    assert calculate_crisproff(MM9, GUIDE, pam_model="crisproff") == pytest.approx(
        11.943630588206144, rel=1e-12, abs=1e-12)


def test_lowercase_inputs_are_case_insensitive() -> None:
    for pam_model in ("deepone", "crisproff"):
        assert calculate_crisproff(
            GUIDE.lower(), MM9.lower(), pam_model=pam_model
        ) == calculate_crisproff(GUIDE, MM9, pam_model=pam_model)
        assert calculate_crisproff(GUIDE.lower()) == calculate_crisproff(GUIDE)


def test_internal_loop_beyond_parameter_range() -> None:
    # The original pipeline aborts with KeyError(12) on a 12 nt internal loop;
    # this API reports it as a ValueError.
    with pytest.raises(ValueError, match="internal loop of 12 nt"):
        calculate_crisproff(GUIDE, "ATGCCCCGCCCGCCCGATGCAGG", pam_model="crisproff")
    with pytest.raises(ValueError, match="internal loop of 12 nt"):
        calculate_crisproff(GUIDE, "ATGCCCCGCCCGCCCGATGCAGG")


@pytest.mark.parametrize("target", ["", "A" * 20, "A" * 22, "A" * 24, "N" * 20 + "AGG",
                                    "U" * 20 + "AGG", "ATGCTAGCTAGCTAGCATGCAAA"])
def test_invalid_off_target(target: str) -> None:
    with pytest.raises(ValueError):
        calculate_crisproff(GUIDE, target)
    with pytest.raises(ValueError):
        calculate_crisproff(GUIDE, target, pam_model="crisproff")


@pytest.mark.parametrize("target", [23, 1.0, ["A" * 23]])
def test_invalid_off_target_type(target: object) -> None:
    with pytest.raises(TypeError):
        calculate_crisproff(GUIDE, target)  # type: ignore[arg-type]


@pytest.mark.parametrize("pam_model", [None, 5, "off"])
def test_invalid_pam_model(pam_model: object) -> None:
    if pam_model is None or isinstance(pam_model, int):
        with pytest.raises(TypeError):
            calculate_crisproff(GUIDE, MM9, pam_model=pam_model)  # type: ignore[arg-type]
    else:
        with pytest.raises(ValueError, match="pam_model"):
            calculate_crisproff(GUIDE, MM9, pam_model=pam_model)  # type: ignore[arg-type]


def test_pam_validation_uses_target_and_mode() -> None:
    # Target PAM outside the selected table.
    with pytest.raises(ValueError, match="DeepOne"):
        calculate_crisproff(GUIDE, "ATGCTAGCTAGCTAGCATGCAAA")
    with pytest.raises(ValueError, match="CRISPRoff"):
        calculate_crisproff(GUIDE, "ATGCTAGCTAGCTAGCATGCAAA", pam_model="crisproff")
    # Self PAM (None target) is validated against the selected table as well.
    ngc_guide = "ATGCTAGCTAGCTAGCATGCGGC"
    assert calculate_crisproff(ngc_guide)  # NGC is scored by DeepOne
    with pytest.raises(ValueError, match="CRISPRoff"):
        calculate_crisproff(ngc_guide, None, pam_model="crisproff")


def test_legacy_one_argument_calls_unchanged() -> None:
    reference = json.loads(Path(__file__).with_name("variant_reference.json").read_text())
    for context, energy in zip(reference["contexts"][:8], reference["energies"][:8], strict=True):
        guide = context[4:27]
        assert calculate_crisproff(guide) == pytest.approx(energy, rel=1e-7, abs=1e-7)
        assert calculate_crisproff(guide) == calculate_crisproff(guide, None)
        assert calculate_crisproff(guide, None) == calculate_crisproff(
            guide, None, pam_model="deepone")


def test_bundled_crisproff_parameters() -> None:
    internal_loop, dimer, bulge1, bulge2 = _crisproff_parameters()
    assert set(internal_loop) == {str(size) for size in range(3, 12)}
    assert len(dimer) == 16 and len(bulge1) == 64 and len(bulge2) == 256
    assert all(len(sub) == 3 for sub in bulge1.values())
    assert all(len(sub) == 9 for sub in bulge2.values())
    complement = str.maketrans("ACGT", "TGCA")
    for guide_dimer in dimer:
        # Fully matched stacks are favorable; G/T double wobbles are not.
        assert dimer[guide_dimer][guide_dimer.translate(complement)] < 0.0
        if set(guide_dimer) <= {"G", "T"}:
            wobble = guide_dimer.translate(str.maketrans("GT", "TG"))
            assert wobble in dimer[guide_dimer]
