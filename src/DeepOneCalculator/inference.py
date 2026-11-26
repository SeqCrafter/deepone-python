"""Lazy TensorFlow loading and bounded ensemble inference."""

from contextlib import ExitStack
from importlib.resources import as_file, files
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from .resources import _manifest

Device = Literal["auto", "cpu", "cuda"]


def _resolve_device(device: Device) -> str:
    if device not in ("auto", "cpu", "cuda"):
        raise ValueError("device must be 'auto', 'cpu', or 'cuda'")
    if device == "cpu":
        return "/CPU:0"
    import tensorflow as tf

    # TensorFlow's GPU enumeration can include non-CUDA plugins (e.g. Metal).
    cuda = tf.test.is_built_with_cuda() and bool(tf.config.list_physical_devices("GPU"))
    if device == "cuda" and not cuda:
        raise RuntimeError("CUDA was requested but no CUDA-capable TensorFlow GPU is available")
    return "/GPU:0" if cuda else "/CPU:0"


def _build_model() -> Any:
    from tensorflow import keras

    layers = keras.layers
    sequence = layers.Input((30, 4), name="input_onehot")
    energy = layers.Input((1,), name="input_energy")
    branches = []
    for width, filters in ((3, 100), (5, 70), (7, 40)):
        x = layers.Conv1D(filters, width, activation="relu", name=f"conv_{width}")(sequence)
        x = layers.Dropout(0.3, name=f"drop_{width}")(x)
        x = layers.AveragePooling1D(2, name=f"pool_{width}")(x)
        branches.append(layers.Flatten(name=f"flatten_{width}")(x))
    # Cropping equals the original Lambda(x[:,18:28,:]), without unsafe deserialization.
    focus = layers.Cropping1D((18, 2), name="focus_seq")(sequence)
    focus = layers.Conv1D(10, 2, activation="relu", name="focus_conv")(focus)
    focus = layers.Dropout(0.3, name="focus_drop")(focus)
    focus = layers.AveragePooling1D(2, name="focus_pool")(focus)
    branches.append(layers.Flatten(name="focus_flatten")(focus))
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
    cell_line: str, device: str, batch_size: int,
) -> NDArray[np.float64]:
    import tensorflow as tf

    weights = _manifest()[cell_line]
    predictions = np.zeros(len(encoded), dtype=np.float64)
    with tf.device(device), ExitStack() as stack:
        model = _build_model()
        for entry in weights:
            resource = files("DeepOneCalculator").joinpath("models").joinpath(entry["file"])
            path = stack.enter_context(as_file(resource))
            model.load_weights(str(path))
            for start in range(0, len(encoded), batch_size):
                end = start + batch_size
                values = model([encoded[start:end], energies[start:end]], training=False)
                predictions[start:end] += np.asarray(values).reshape(-1)
    return predictions / len(weights)
