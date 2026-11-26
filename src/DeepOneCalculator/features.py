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


@dataclass(frozen=True, slots=True)
class _Guide:
    sequence_id: str
    sequence: str
    grna_id: str
    context: str

    @property
    def grna_seq(self) -> str:
        return self.context[4:27]


def _find_guides(sequence: str, sequence_id: str) -> list[_Guide]:
    guides = []
    for strand, oriented in (("p", sequence), ("m", sequence.translate(_COMPLEMENT)[::-1])):
        for start in range(len(oriented) - 29):
            context = oriented[start:start + 30]
            if context[25:27] == "GG":
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
    pam_ratio = 1.1 if guide[-3:] == "GGG" else 1.0
    return (-hybrid + fold_energy(spacer) + opening) * pam_ratio


def _encode(contexts: list[str]) -> NDArray[np.float32]:
    mapping = {base: i for i, base in enumerate("ATCG")}
    indices = np.array([[mapping[b] for b in s] for s in contexts], dtype=np.intp)
    return np.eye(4, dtype=np.float32)[indices]
