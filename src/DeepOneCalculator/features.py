"""Guide discovery and the training model's thermodynamic feature."""

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

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


def calculate_free_energy(grna_seq: str) -> float:
    """Calculate the complete thermodynamic feature used by DeepOne.

    Parameters
    ----------
    grna_seq : str
        Exactly 23 DNA bases: the oriented 20 nt spacer followed by its
        3 bp PAM, without flanking context. Lowercase DNA is accepted.
        The PAM must have an upstream correction factor (NGN or NAG).

    Returns
    -------
    float
        ``(-weighted_hybridization + spacer_fold_energy + dna_opening)
        * pam_ratio``. Components use kcal/mol energy parameters, but the
        returned value is the model's weighted composite feature, not an
        unweighted physical binding free energy or the DeepOne score.
        It is the same scalar that score_grnas supplies to the network.

    Raises
    ------
    TypeError
        If grna_seq is not a string.
    ValueError
        If the input is not exactly 23 A/C/G/T bases, or its PAM has no
        DeepOne correction factor.

    Notes
    -----
    This function does not discover guides or select Cas9 models. PAM
    correction factors for NAG exist in the upstream feature table even
    though score_grnas scans only NGG (SpCas9) or NGN (the two variants).
    Use fold_energy for RNA secondary-structure minimum free energy alone.
    TensorFlow, H5 weights and a GPU are not required for this calculation.

    Examples
    --------
    >>> isinstance(calculate_free_energy("AT" * 10 + "AGG"), float)
    True
    """
    if not isinstance(grna_seq, str):
        raise TypeError("grna_seq must be a DNA string")
    normalized = grna_seq.upper()
    if len(normalized) != 23 or set(normalized) - set("ACGT"):
        raise ValueError("grna_seq must contain exactly 23 A/C/G/T bases (20 nt spacer + PAM)")
    if normalized[-3:] not in _PAM_RATIOS:
        raise ValueError("PAM has no DeepOne correction factor; expected NGN or NAG")
    return _energy(normalized)


def _encode(contexts: list[str]) -> NDArray[np.float32]:
    mapping = {base: i for i, base in enumerate("ATCG")}
    indices = np.array([[mapping[b] for b in s] for s in contexts], dtype=np.intp)
    return np.eye(4, dtype=np.float32)[indices]
