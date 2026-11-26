# DeepOneCalculator

`DeepOneCalculator` packages pretrained DeepOne ensembles for predicting SpCas9-induced,
templated **+1 bp insertion frequencies**. This is the original model's outcome
score, not a general gRNA editing-efficiency probability. Scores retain the training
scale and are neither clipped nor converted to percentages.

## Install

```bash
pip install DeepOneCalculator
```

This command applies once a release has been published to your package index.
For this standalone source checkout (Python 3.10 or newer):

```bash
python -m pip install .
# Build a distributable wheel, including the H5 models:
pip wheel . -w dist/
```

Source builds need a C compiler and Python development headers. Build isolation
installs Cython, NumPy and setuptools automatically. Binary wheels require no C
compiler. Neither installation nor inference requires a system ViennaRNA library,
the ViennaRNA Python package, or an `RNAfold` executable. The private Cython
extension statically includes a portable subset of ViennaRNA 2.7.2 RNAlib.

## Quick start

```python
from DeepOneCalculator import SUPPORTED_CELL_LINES, score_grnas

scores = score_grnas(
    sequences=["TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT"],
    cell_lines=["HEK293", "K562"],
    batch_size=256,
    device="auto",
)
print(scores)
# A single sequence and all human cell lines:
all_scores = score_grnas("AAAA" + "ATATATATATATATATATAT" + "AGG" + "AAA")
```

Inputs contain at least **30 bp = 4 bp prefix + 20 nt spacer + 3 bp PAM + 3 bp
suffix**. The 23 bp gRNA field includes the PAM. Longer sequences are scanned for
every complete NGG context on both strands. Lowercase DNA is normalized to uppercase;
ambiguous bases, whitespace and sequences shorter than 30 bp raise `ValueError`.
An empty list or a sequence without complete NGG sites produces an empty DataFrame.

The output has one row per discovered guide and requested cell line:

| Column | Meaning |
| --- | --- |
| `sequence_id` | One-based input identifier, e.g. `sequence_1` |
| `sequence` | Complete normalized input DNA |
| `gRNA_id` | Input ID, p/m strand, one-based oriented context start; e.g. `sequence_1_p1` |
| `gRNA_seq` | Oriented 20 nt spacer followed by 3 bp PAM |
| `cell_line` | Canonical model identifier |
| `deepone_score` | Unrounded arithmetic mean of the five pretrained models |

Rows follow the requested cell line order, then input order, then forward and
reverse context order. Repeated cell-line names are deduplicated; repeated input
sequences keep distinct IDs. Names are case insensitive.

## Supported human cell lines

| Identifier | Models |
| --- | --- |
| `HEK293` | 5 pretrained H5 weights |
| `K562` | 5 pretrained H5 weights |
| `HAP1` | 5 pretrained H5 weights |
| `iPSC` | 5 pretrained H5 weights |
| `RPE-1` | 5 pretrained H5 weights |

`cell_lines="all"` selects these five. HEK293T and HeLa weights are absent from this
repository and are not substituted with HEK293. Nonhuman CHO/mESCs and Cas9-variant models are not part of this human SpCas9 API.
Each packaged H5 file is approximately 1 MB. `models/manifest.json` records their
original filenames and SHA-256 checksums.

## Devices and folding

`device="auto"` selects TensorFlow's first CUDA GPU if TensorFlow was built with
CUDA and detects a GPU; otherwise it uses CPU. `device="cpu"` always uses CPU.
`device="cuda"` fails explicitly if CUDA is unavailable. On Linux a compatible
CUDA TensorFlow environment can be installed with `pip install 'tensorflow[and-cuda]'`.
Metal plugins are not treated as CUDA devices.

Pip does not detect GPU hardware or choose CUDA dependencies automatically. The
base install requests `tensorflow`, without its optional CUDA runtime libraries.
On a supported Linux x86_64 or Windows WSL2 Linux environment with an NVIDIA GPU:

```bash
pip install 'DeepOneCalculator[cuda]'
# For this local source checkout:
pip install '.[cuda]'
# Verify that TensorFlow detects a GPU:
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

The `cuda` extra adds `tensorflow[and-cuda]`, which installs the matching CUDA/cuDNN
user-space dependencies. It does not install the NVIDIA driver. If compatible
libraries and a driver already exist, the base TensorFlow install may already
use the GPU; a separate `tensorflow-gpu` package is unnecessary. macOS has no CUDA
support, and current TensorFlow CUDA support on Windows requires WSL2. The package's
C extension currently targets POSIX Linux/macOS. See the
[official TensorFlow installation guide](https://www.tensorflow.org/install/pip).


```python
from DeepOneCalculator._viennarna import fold_energy

energy = fold_energy("GGGAAACCC")  # kcal/mol, default 37 °C Turner 2004 parameters
```

The inference feature combines the original RNA/DNA nearest-neighbor energies,
position weighting, DNA opening and PAM correction with the actual 20 nt spacer
fold energy. The legacy script accidentally sent individual bases to RNAfold;
this corrected term can change scores relative to that script. ViennaRNA version
changes can also affect predictions. See [IMPLEMENTATION_NOTES.md](IMPLEMENTATION_NOTES.md).

## Development and validation

```bash
pip install '.[dev]' ViennaRNA==2.7.2
ruff check .
mypy
pytest
python -m build
```

Tests cover real H5 inference across all five ensembles, batching, input validation,
both-strand discovery, CUDA selection/fallback, packaged checksums, folding and
concurrent folding calls. The optional ViennaRNA Python package is an independent
folding reference. CLI comparisons use a separately installed **RNAfold 2.7.2**;
without it those tests skip. Set `DEEPONE_REQUIRE_RNAFOLD=1` to make absence a failure.
CI installs the reference CLI, requires these comparisons, runs Ruff/mypy, and
provides a cibuildwheel workflow for manylinux x86_64 CPython 3.10–3.12 wheels.
The reference tool is used exclusively in tests.

## License and citation

The repository is GPL-3.0-only; bundled ViennaRNA retains its upstream license in
`vendor/viennarna/COPYING`. See the vendor provenance file for the exact source
archive and modifications.

Original model: [HaoDK12/DeepOne_model](https://github.com/HaoDK12/DeepOne_model).
Citation: Hao Yuan, Xiaoguang Pan et al., *Template-free single nucleotide insertion
with deep-learning CRISPR gRNAs for highly efficient gene disruption and systemic
correction of DMD mutations*, 2026 (manuscript under revision).

## Standalone repository

This directory is an independent Git repository and can be cloned, built and tested
without its parent DeepOne repository. All 25 H5 models and the required C sources
are included here; no imports, symlinks or build paths refer to the parent checkout.
The parent repository excludes this directory using `/deepone-python/`.

The vendored implementation supports only the energy of a single linear RNA using
fixed default RNAfold settings. It excludes structure traceback, parameter-file
loading, alternative parameter sets, alignment interfaces, plotting, partition
function entry points, and other standalone ViennaRNA tools. Shared internal
allocation, energy evaluation and folding helpers remain as required dependencies.
The extension exports only its Python initialization symbol; RNAlib internals are
hidden. No ViennaRNA command-line programs are built or bundled.
