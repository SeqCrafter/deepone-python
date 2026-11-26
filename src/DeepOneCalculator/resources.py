"""Immutable metadata for the pretrained human cell line ensembles."""

import json
from functools import lru_cache
from importlib.resources import files

SUPPORTED_CELL_LINES: tuple[str, ...] = ("HEK293", "K562", "HAP1", "iPSC", "RPE-1")


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
