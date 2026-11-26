"""Metadata and selection for the cell line and Cas9 variant ensembles."""

import json
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files

SUPPORTED_CELL_LINES: tuple[str, ...] = ("HEK293", "K562", "HAP1", "iPSC", "RPE-1")
SUPPORTED_CAS9_TYPES: tuple[str, ...] = ("SpCas9", "SpCas9-NG", "SpG")


@dataclass(frozen=True, slots=True)
class _ModelSpec:
    cas9_type: str
    cell_line: str | None
    weight_key: str


@lru_cache(maxsize=1)
def _manifest() -> dict[str, list[dict[str, str]]]:
    data: dict[str, list[dict[str, str]]] = json.loads(
        files("DeepOneCalculator").joinpath("models/manifest.json").read_text()
    )
    return data


def _normalize_cell_lines(cell_lines: list[str] | str) -> list[str]:
    if isinstance(cell_lines, str) and cell_lines.lower() == "all":
        return list(SUPPORTED_CELL_LINES)
    values = [cell_lines] if isinstance(cell_lines, str) else cell_lines
    if not isinstance(values, list) or not values:
        raise ValueError("cell_lines must be a nonempty string or list of strings")
    names = {name.lower(): name for name in SUPPORTED_CELL_LINES}
    result = []
    for value in values:
        if not isinstance(value, str) or value.lower() not in names:
            raise ValueError(f"Unsupported cell line {value!r}; choose {SUPPORTED_CELL_LINES}")
        canonical = names[value.lower()]
        if canonical not in result:
            result.append(canonical)
    return result


def _normalize_cas9_types(cas9_types: list[str] | str) -> list[str]:
    if isinstance(cas9_types, str) and cas9_types.lower() == "all":
        return list(SUPPORTED_CAS9_TYPES)
    values = [cas9_types] if isinstance(cas9_types, str) else cas9_types
    if not isinstance(values, list) or not values:
        raise ValueError("cas9_types must be a nonempty string or list of strings")
    names = {name.replace("-", "").lower(): name for name in SUPPORTED_CAS9_TYPES}
    result = []
    for value in values:
        key = value.replace("_", "").replace("-", "").lower() if isinstance(value, str) else ""
        if key not in names:
            raise ValueError(f"Unsupported Cas9 type {value!r}; choose {SUPPORTED_CAS9_TYPES}")
        canonical = names[key]
        if canonical not in result:
            result.append(canonical)
    return result


def _model_specs(cas9_types: list[str], cell_lines: list[str]) -> list[_ModelSpec]:
    return [
        spec
        for cas9_type in cas9_types
        for spec in (
            [_ModelSpec(cas9_type, cell_line, cell_line) for cell_line in cell_lines]
            if cas9_type == "SpCas9" else [_ModelSpec(cas9_type, None, cas9_type)]
        )
    ]
