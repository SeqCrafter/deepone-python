"""CRISPRon indel efficiency predictions with the six pretrained models."""

import json
from contextlib import ExitStack
from functools import lru_cache
from importlib.resources import as_file, files
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .features import _find_guides, calculate_crisproff
from .inference import Device, _resolve_device

_COLUMNS = ["sequence_id", "sequence", "gRNA_id", "gRNA_seq", "crispron_score"]


@lru_cache(maxsize=1)
def _manifest() -> list[dict[str, str]]:
    entries: list[dict[str, str]] = json.loads(
        files("DeepOneCalculator").joinpath("models/crispron/manifest.json").read_text()
    )
    return entries


def _encode(contexts: list[str]) -> NDArray[np.float32]:
    indices = np.array([["ACGT".index(base) for base in context] for context in contexts])
    encoded: NDArray[np.float32] = np.eye(4, dtype=np.float32)[indices]
    return encoded


def _build_model() -> Any:
    from tensorflow import keras

    layers = keras.layers
    sequence = layers.Input((30, 4), name="input_onehot")
    energy = layers.Input((1,), name="input_dGB")
    branches = []
    for width, filters in ((3, 100), (5, 70), (7, 40)):
        x = layers.Conv1D(filters, width, activation="relu", name=f"conv_{width}")(sequence)
        x = layers.Dropout(0.3, name=f"drop_{width}")(x)
        x = layers.AveragePooling1D(2, padding="same", name=f"pool_{width}")(x)
        branches.append(layers.Flatten(name=f"flatten_{width}")(x))
    x = layers.Concatenate()(branches)
    x = layers.Dense(80, activation="relu", name="dense_0")(x)
    x = layers.Dropout(0.3, name="drop_d0")(x)
    x = layers.Concatenate()([x, energy])
    for index, units in ((1, 80), (2, 60)):
        x = layers.Dense(units, activation="relu", name=f"dense_{index}")(x)
        x = layers.Dropout(0.3, name=f"drop_d{index}")(x)
    output = layers.Dense(1, name="output")(x)
    return keras.Model([sequence, energy], output)


def _predict(
    encoded: NDArray[np.float32], energies: NDArray[np.float32],
    device: str, batch_size: int,
) -> NDArray[np.float64]:
    import tensorflow as tf

    weights = _manifest()
    predictions = np.zeros(len(encoded), dtype=np.float64)
    with tf.device(device), ExitStack() as stack:
        model = _build_model()
        for entry in weights:
            resource = (
                files("DeepOneCalculator").joinpath("models/crispron").joinpath(entry["file"])
            )
            path = stack.enter_context(as_file(resource))
            model.load_weights(str(path))
            for start in range(0, len(encoded), batch_size):
                end = start + batch_size
                values = model([encoded[start:end], energies[start:end]], training=False)
                predictions[start:end] += np.asarray(values).reshape(-1)
    return predictions / len(weights)


def calculate_crispron(
    sequences: list[str] | str,
    batch_size: int = 256,
    device: Device = "auto",
) -> pd.DataFrame:
    """Predict indel efficiency for both-strand NGG sites in ACGT DNA of at least 30 bp."""
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size <= 0:
        raise ValueError("batch_size must be a positive integer")
    if device not in ("auto", "cpu", "cuda"):
        raise ValueError("device must be 'auto', 'cpu', or 'cuda'")
    values = [sequences] if isinstance(sequences, str) else sequences
    if not isinstance(values, list):
        raise ValueError("sequences must be a string or list of strings")
    guides = []
    for index, sequence in enumerate(values, 1):
        if not isinstance(sequence, str) or len(sequence) < 30:
            raise ValueError(f"sequence {index} must be a DNA string of at least 30 bp")
        normalized = sequence.upper()
        if set(normalized) - set("ACGT"):
            raise ValueError(f"sequence {index} contains non-ACGT characters")
        guides.extend(_find_guides(normalized, f"sequence_{index}"))
    if not guides:
        empty: pd.DataFrame = pd.DataFrame({
            name: pd.Series(dtype="float64" if name == "crispron_score" else "object")
            for name in _COLUMNS
        })
        return empty
    resolved = _resolve_device(device)
    encoded = _encode([guide.context for guide in guides])
    energies = np.array([
        [calculate_crisproff(guide.grna_seq, pam_model="crisproff")] for guide in guides
    ], dtype=np.float32)
    predictions = _predict(encoded, energies, resolved, batch_size)
    scores: pd.DataFrame = pd.DataFrame([
        (guide.sequence_id, guide.sequence, guide.grna_id, guide.grna_seq, float(score))
        for guide, score in zip(guides, predictions, strict=True)
    ], columns=_COLUMNS)
    return scores
