"""Guide discovery and the training model's thermodynamic feature."""

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from ._viennarna import fold_energy

_COMPLEMENT = str.maketrans("ATCG", "TAGC")
_POSITION_WEIGHTS = (
    1.80067099242007, 1.95666668400006, 1.90472004401173, 2.13047270152512,
    1.37853848098249, 1.46460783730408, 1.0, 1.387220146823, 1.51401000729362,
    1.98058344620751, 1.87939168587699, 1.7222593588838, 2.02228445489326,
    1.92692086621503, 2.08041972716723, 1.94496755678903, 2.14539112893591,
    2.04277109036766, 2.24911493451185, 2.25,
)
_DNA_NN = dict(zip(
    ("AA", "TT", "AT", "TA", "CA", "TG", "GT", "AC", "CT", "AG", "GA", "TC",
     "CG", "GC", "GG", "CC"),
    (-1., -1., -.88, -.58, -1.45, -1.45, -1.44, -1.44, -1.28, -1.28, -1.30,
     -1.30, -2.17, -2.24, -1.84, -1.84), strict=True,
))
# Effective values from upstream Energy_cal.py after resolving duplicate keys.
_PAM_RATIOS = {
    "AGA": 1.0, "AGC": 0.8, "AGG": 1.0, "AGT": 1.2,
    "CGA": 1.0, "CGC": 0.9, "CGG": 1.0, "CGT": 1.1,
    "GGA": 1.1, "GGC": 0.9, "GGG": 1.1, "GGT": 1.2,
    "TGA": 0.9, "TGC": 0.8, "TGG": 1.0, "TGT": 1.1,
    "AAG": 0.9, "CAG": 0.9, "GAG": 0.9, "TAG": 0.9,
}
# Original CRISPRoff 1.1.2 pam_ratios from the pam_corr branch of get_eng.
_CRISPROFF_PAM_RATIOS = {
    "AGG": 1.0, "CGG": 1.0, "GGG": 1.0, "TGG": 1.0,
    "AAG": 0.9, "CAG": 0.9, "GAG": 0.9, "TAG": 0.9,
    "AGA": 0.8, "CGA": 0.8, "GGA": 0.8, "TGA": 0.8,
}
# Watson-Crick pairs only, as in CRISPRoff's RI_MATCH_noGU (no GU pairs).
_WATSON_CRICK = frozenset({"AT", "TA", "CG", "GC"})
_PAM_MODELS = ("deepone", "crisproff")
_StrDict = dict[str, float]
_StrTable = dict[str, dict[str, float]]


@dataclass(frozen=True, slots=True)
class _Guide:
    sequence_id: str
    sequence: str
    grna_id: str
    context: str

    @property
    def grna_seq(self) -> str:
        return self.context[4:27]


def _find_guides(sequence: str, sequence_id: str, *, allow_ngn: bool = False) -> list[_Guide]:
    guides = []
    for strand, oriented in (("p", sequence), ("m", sequence.translate(_COMPLEMENT)[::-1])):
        for start in range(len(oriented) - 29):
            context = oriented[start:start + 30]
            if context[25] == "G" and (allow_ngn or context[26] == "G"):
                guides.append(_Guide(sequence_id, sequence, f"{sequence_id}_{strand}{start + 1}",
                                     context))
    return guides


@lru_cache(maxsize=1)
def _rna_parameters() -> dict[str, float]:
    parameters: dict[str, float] = json.loads(
        files("DeepOneCalculator").joinpath("energy_parameters.json").read_text()
    )
    return parameters


@lru_cache(maxsize=8192)
def _energy(guide: str) -> float:
    spacer = guide[:20]
    parameters = _rna_parameters()
    hybrid = 0.0
    opening = 0.0
    for i in range(19):
        pair = spacer[i:i + 2]
        contribution = parameters[pair]
        if (i == 0 and spacer[i] in "AT") or (i == 18 and spacer[i + 1] in "AT"):
            contribution += 0.25
        hybrid += contribution * _POSITION_WEIGHTS[i]
        opening += _DNA_NN[pair]
    pam_ratio = _PAM_RATIOS.get(guide[-3:], 0.0)
    return (-hybrid + fold_energy(spacer) + opening) * pam_ratio


@lru_cache(maxsize=1)
def _crisproff_parameters() -> tuple[_StrDict, _StrTable, _StrTable, _StrTable]:
    """Load the bundled CRISPRoff 1.1.2 RNA-DNA duplex parameters.

    The JSON is a safe conversion of the local crisproff-1.1.2 checkout:
    ``dimer``, ``bulge1`` and ``bulge2`` are the RNA_DNA tuple from
    energy_dics.pkl (loaded with a restricted no-globals unpickler), and
    ``internal_loop`` is RNA_DNA_internal_loop from
    CRISPRspec_CRISPRoff_pipeline.py.  Runtime code never touches the pickle
    or the checkout.
    """
    raw = json.loads(
        files("DeepOneCalculator").joinpath("crisproff_parameters.json").read_text()
    )
    return raw["internal_loop"], raw["dimer"], raw["bulge1"], raw["bulge2"]


def _is_match(guide_base: str, target_base: str) -> bool:
    return guide_base + target_base in _WATSON_CRICK


@lru_cache(maxsize=8192)
def _pair_hybrid(spacer: str, target_spacer: str) -> tuple[float, ...]:
    """Per-base RNA-DNA hybrid energies of the aligned 20 nt guide-target pair.

    Reimplements CRISPRoff 1.1.2 calcRNADNAenergy: ``spacer`` is the guide
    (RNA strand) and ``target_spacer`` the target DNA strand.  Stacked
    matches, 1-2 nt bulges and 3-11 nt internal loops take their parameters
    from the bundled tables; an AU pair at the first or last matched
    position of the hybrid receives the 0.25 kcal/mol terminal penalty.
    Positions outside the first-to-last matched base (terminal mismatches)
    contribute no energy, matching the original walk.
    """
    internal_loop, dimer, bulge1, bulge2 = _crisproff_parameters()
    tables = (dimer, bulge1, bulge2)
    seq = target_spacer.translate(_COMPLEMENT)
    energy = [0.0] * 20
    spos = -1
    epos = -1
    for i in range(20):
        if _is_match(spacer[i], seq[i]):
            if spos == -1:
                spos = i
            epos = i
    i = spos
    while i < epos:
        j = i + 1
        while not _is_match(seq[j], spacer[j]):
            j += 1
            if j > epos:
                break
        if j > epos:
            break
        loop_size = j - i - 1
        if loop_size < 3:
            contribution = float(tables[loop_size][spacer[i:j + 1]][seq[i:j + 1]])
            if loop_size == 0 and (
                (i == spos and (spacer[i] == "T" or seq[i] == "T"))
                or (j == epos and (spacer[j] == "T" or seq[j] == "T"))
            ):
                contribution += 0.25
        else:
            if loop_size > 11:
                raise ValueError(
                    f"internal loop of {loop_size} nt exceeds the "
                    "11 nt CRISPRoff parameter range"
                )
            contribution = (
                float(internal_loop[str(loop_size)])
                + float(dimer[spacer[i:i + 2]][seq[i:i + 2]])
                + float(dimer[spacer[j - 1:j + 1]][seq[j - 1:j + 1]])
            )
        share = contribution / (loop_size + 1)
        for k in range(loop_size + 1):
            energy[i + k] += share
        i = j
    return tuple(energy)


@lru_cache(maxsize=8192)
def _spacer_fold(spacer: str) -> float:
    return fold_energy(spacer)


def _dna_opening(spacer: str) -> float:
    """Target DNA-DNA nearest-neighbor opening energy (first base contributes 0)."""
    return sum([0.0] + [_DNA_NN[spacer[i - 1:i + 1]] for i in range(1, 20)])


def _pam_ratio(pam_model: str, pam: str) -> float:
    if pam_model == "deepone":
        ratio = _PAM_RATIOS.get(pam)
        if ratio is None:
            raise ValueError("PAM has no DeepOne correction factor; expected NGN or NAG")
        return ratio
    ratio = _CRISPROFF_PAM_RATIOS.get(pam)
    if ratio is None:
        raise ValueError("PAM has no CRISPRoff correction factor; expected NGG, NAG or NGA")
    return ratio


def _pair_energy(guide: str, target: str, ratio: float) -> float:
    """CRISPRoff 1.1.2 get_eng with pos_weight, pam_corr, grna_folding and
    dna_opening enabled (dna_pos_wgh disabled), for normalized 23 nt inputs.

    Folding always uses the guide spacer and opening always uses the target
    DNA; the PAM ratio is applied to their combined sum, last.
    """
    hybrid = _pair_hybrid(guide[:20], target[:20])
    scores = [weight * value for weight, value in zip(_POSITION_WEIGHTS, hybrid, strict=True)]
    off = sum(scores) + 0.0
    off = (-1.0) * off
    off += _spacer_fold(guide[:20])
    off += _dna_opening(target[:20])
    return off * ratio


def calculate_crisproff(
    grna_seq: str,
    off_target_seq: str | None = None,
    *,
    pam_model: Literal["deepone", "crisproff"] = "deepone",
) -> float:
    """Calculate the CRISPRoff weighted energy score of a guide against a target site.

    The default PAM correction remains ``'deepone'``; select
    ``pam_model='crisproff'`` for the original CRISPRoff PAM corrections.

    Parameters
    ----------
    grna_seq : str
        Exactly 23 DNA bases: the oriented 20 nt spacer followed by its
        3 bp PAM, without flanking context. Lowercase DNA is accepted.
    off_target_seq : str or None, optional
        23 bp target DNA (20 nt spacer + 3 bp PAM) to score against the
        guide; lowercase accepted. ``None`` scores the guide against
        itself.
    pam_model : {'deepone', 'crisproff'}, optional, keyword-only
        PAM correction table applied to the target PAM. ``'deepone'``
        (default) uses the DeepOne NGN/NAG table; ``'crisproff'`` uses the
        original CRISPRoff table (NGG = 1.0, NAG = 0.9, NGA = 0.8) and
        rejects every other PAM.

    Returns
    -------
    float
        With ``off_target_seq=None`` and ``pam_model='deepone'`` this is the
        historical DeepOne         self-binding feature
        ``(-weighted_hybridization + spacer_fold_energy + dna_opening)
        * pam_ratio`` computed on the guide itself, unchanged from the
        single-argument function.  Otherwise it is the original CRISPRoff
        pair score
        ``(-sum_p POS_WGH[p] * hybrid_p + fold_energy(guide[:20])
        + dna_opening(target)) * pam_ratio(target[-3:])``, where position
        ``p`` counts 5' to 3' and ``hybrid_p`` are the per-base aligned 20
        nt RNA-DNA hybrid energies (matched stacks, exact 1-2 nt bulges,
        3-11 nt internal loops, terminal AU penalty); bases outside the
        first-to-last matched position carry no energy.  The PAM
        correction always uses
        the target PAM and folding always uses the guide spacer, never the
        off-target.  Components use kcal/mol parameters, but the returned
        value is the model's weighted composite feature, not an unweighted
        physical binding free energy or the DeepOne score.

    Raises
    ------
    TypeError
        If grna_seq or off_target_seq is not a string (or None), or if
        pam_model is not a string.
    ValueError
        If a sequence is not exactly 23 A/C/G/T bases, if pam_model is not
        'deepone' or 'crisproff', if the target PAM has no correction
        factor in the selected table, or if the mismatch pattern spans an
        internal loop of more than 11 nt (beyond the CRISPRoff parameter
        range; the original pipeline aborts on such inputs).

    Notes
    -----
    This function does not discover guides or select Cas9 models.  With an
    explicit off-target, both modes share the mismatch-aware CRISPRoff
    hybrid (the DeepOne table has no mismatch parameters); only the PAM
    correction table differs.  PAM correction factors for NAG exist in the
    upstream feature table even though calculate_deepone scans only NGG (SpCas9)
    or NGN (the two variants).  Use fold_energy for RNA secondary-structure
    minimum free energy alone.  TensorFlow, H5 weights and a GPU are not
    required for this calculation.

    Examples
    --------
    >>> isinstance(calculate_crisproff("AT" * 10 + "AGG"), float)
    True
    """
    if not isinstance(grna_seq, str):
        raise TypeError("grna_seq must be a DNA string")
    if off_target_seq is not None and not isinstance(off_target_seq, str):
        raise TypeError("off_target_seq must be a DNA string or None")
    if not isinstance(pam_model, str):
        raise TypeError("pam_model must be 'deepone' or 'crisproff'")
    if pam_model not in _PAM_MODELS:
        raise ValueError("pam_model must be 'deepone' or 'crisproff'")
    normalized = grna_seq.upper()
    if len(normalized) != 23 or set(normalized) - set("ACGT"):
        raise ValueError("grna_seq must contain exactly 23 A/C/G/T bases (20 nt spacer + PAM)")
    if off_target_seq is None and pam_model == "deepone":
        if normalized[-3:] not in _PAM_RATIOS:
            raise ValueError("PAM has no DeepOne correction factor; expected NGN or NAG")
        return _energy(normalized)
    target = normalized if off_target_seq is None else off_target_seq.upper()
    if len(target) != 23 or set(target) - set("ACGT"):
        raise ValueError(
            "off_target_seq must contain exactly 23 A/C/G/T bases (20 nt spacer + PAM)"
        )
    return _pair_energy(normalized, target, _pam_ratio(pam_model, target[-3:]))


def _encode(contexts: list[str]) -> NDArray[np.float32]:
    mapping = {base: i for i, base in enumerate("ATCG")}
    indices = np.array([[mapping[b] for b in s] for s in contexts], dtype=np.intp)
    return np.eye(4, dtype=np.float32)[indices]
