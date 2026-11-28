# DeepOneCalculator

[![Test](https://github.com/SeqCrafter/deepone-python/actions/workflows/ci.yml/badge.svg)](https://github.com/SeqCrafter/deepone-python/actions/workflows/ci.yml)
[![Manylinux wheels](https://github.com/SeqCrafter/deepone-python/actions/workflows/wheels.yml/badge.svg)](https://github.com/SeqCrafter/deepone-python/actions/workflows/wheels.yml)
[![Publish to PyPI](https://github.com/SeqCrafter/deepone-python/actions/workflows/publish.yml/badge.svg)](https://github.com/SeqCrafter/deepone-python/actions/workflows/publish.yml)
![PyPI Version](https://img.shields.io/pypi/v/DeepOneCalculator)

## Installation

Python 3.10 or later is required. Installations with RIsearch require Python 3.11
or later.

### Standard

```bash
python -m pip install DeepOneCalculator
```

### With CUDA

For NVIDIA GPU inference on a supported Linux system with a compatible NVIDIA
driver:

```bash
python -m pip install "DeepOneCalculator[cuda]"
```

### With RIsearch

Install the `specificity` extra to enable genome-wide CRISPRspec scoring:

```bash
python -m pip install "DeepOneCalculator[specificity]"
```

### With CUDA and RIsearch

```bash
python -m pip install "DeepOneCalculator[cuda,specificity]"
```

## API Overview

```python
from DeepOneCalculator import (
    calculate_deepone,
    calculate_crispron,
    calculate_crisprspec,
    calculate_crisproff,
    fold_energy,
    CRISPRspecConfig,
    ReferenceGenome,
)
```

| API | Purpose | Returns |
| --- | --- | --- |
| `calculate_deepone` | Predict templated +1 bp insertion frequency | DataFrame with `deepone_score` |
| `calculate_crispron` | Predict indel editing efficiency | DataFrame with `crispron_score` |
| `calculate_crisprspec` | Calculate genome-wide specificity using RIsearch | `float` |
| `calculate_crisproff` | Calculate a guide-target weighted energy score | `float` |
| `fold_energy` | Calculate RNA secondary-structure minimum free energy | `float` |
| `CRISPRspecConfig` | Configure CRISPRspec search and scoring | Configuration object |
| `ReferenceGenome.from_fasta` | Load a reference genome and build or reuse its search index | `ReferenceGenome` |

## Usage

### DeepOne

```python
from DeepOneCalculator import calculate_deepone

sequences = ["TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT"]
scores = calculate_deepone(
    sequences,
    cell_lines=["HEK293", "K562"],
    cas9_types="SpCas9",
    batch_size=256,
    device="auto",
)
print(scores)

# Include SpCas9, SpCas9-NG and SpG.
variant_scores = calculate_deepone(sequences, cas9_types="all")
```

- **Input:** a DNA string or list of strings, each containing at least 30 A/C/G/T
  bases, case insensitive. Both strands are scanned for complete contexts:
  4 bp prefix + 20 nt spacer + 3 bp PAM + 3 bp suffix.
- **Cas9 types:** `SpCas9` (NGG PAM), `SpCas9-NG` and `SpG` (NGN PAM).
- **Cell lines:** `HEK293`, `K562`, `HAP1`, `iPSC`, `RPE-1`, or `"all"` (default).
  Cell-line selection applies only to SpCas9.
- **Device:** `"auto"` (default), `"cpu"`, or `"cuda"`. Auto uses CUDA when
  available; explicitly requesting unavailable CUDA raises an error.
- **Output columns:** `sequence_id`, `sequence`, `gRNA_id`, `gRNA_seq`,
  `cell_line`, `cas9_type`, `deepone_score`. Variant rows have no cell-line value.

### CRISPRon

```python
from DeepOneCalculator import calculate_crispron

sequence = "ACTGAACTTGAAAAGCAAAAAGAAACTGGCCATACTTTCGAAGAAATGCTACTGACTG"
scores = calculate_crispron(sequence, batch_size=256, device="auto")
print(scores)
```

Accepts a DNA string or list of strings with the same 30 bp minimum context as
DeepOne, scanning NGG sites on both strands. `batch_size` defaults to `256`;
`device` accepts `"auto"`, `"cpu"`, or `"cuda"`.

Returns a DataFrame with columns `sequence_id`, `sequence`, `gRNA_id`,
`gRNA_seq`, and `crispron_score`. The score predicts indel editing efficiency,
not DeepOne's +1 bp insertion frequency. Both APIs return an empty DataFrame
with the same columns when no eligible guides are found.

### CRISPRspec

Requires the `specificity` extra. The guide must contain exactly 23 DNA bases:
20 nt spacer + 3 bp PAM. The reference can be a plain or gzip-compressed DNA
FASTA file.

```python
from DeepOneCalculator import (
    CRISPRspecConfig,
    ReferenceGenome,
    calculate_crisprspec,
)

guide = "GGTGGACAAGCGGCAGATAGCGG"
score = calculate_crisprspec(guide, "genome.fa", index_path="genome.idx")

# Reuse the reference and index when scoring multiple guides.
reference = ReferenceGenome.from_fasta("genome.fa", index_path="genome.idx")
config = CRISPRspecConfig(
    max_mismatches=6,
    pam_suffixes=("GG", "AG", "GA"),
    ignored_chromosomes=frozenset(),
    threads=4,
)
score = calculate_crisprspec(guide, reference, config=config)
```

Higher scores indicate greater predicted specificity; no retained off-targets
returns `math.inf`. Defaults search both strands with up to six spacer
substitutions and NGG/NAG/NGA PAMs. An exact on-target must be present unless
`require_on_target=False` is set in `CRISPRspecConfig`.

Omit `index_path` to use a temporary index. Persistent indexes must match the
reference; keep their `.deepone.json` and `.lock` sidecar files alongside them.
CRISPRspec defaults to `pam_model="crisproff"`.

### CRISPRoff

```python
from DeepOneCalculator import calculate_crisproff

guide = "GGTGGACAAGCGGCAGATAGCGG"
off_target = "AGTGGACAAGCGGCAGATAGCGG"

# Original CRISPRoff PAM corrections.
on_target_score = calculate_crisproff(guide, pam_model="crisproff")
off_target_score = calculate_crisproff(guide, off_target, pam_model="crisproff")

# Default: DeepOne PAM corrections.
deepone_energy = calculate_crisproff(guide)
```

Both sequences must contain exactly 23 A/C/G/T bases (20 nt spacer + 3 bp PAM),
case insensitive. Omitting `off_target_seq` scores the guide against itself.
The default is `pam_model="deepone"`; explicitly select `"crisproff"` for the
original CRISPRoff PAM table. Returns a weighted energy score, not a DeepOne
prediction or an unweighted physical binding free energy.

### RNA Folding Energy

```python
from DeepOneCalculator import fold_energy

energy = fold_energy("ATCTTCGCTATCACCTCCGCC")
```

Returns RNA secondary-structure minimum free energy in kcal/mol. No external
RNAfold executable is required.

## References

### DeepOne

Hao Yuan, Xiaoguang Pan, Yonglun Luo et al.
*Predicting CRISPR-SpCas9-induced single nucleotide insertion by deep learning.*
**Nature Communications**, 2026 (manuscript under revision).

[Original model](https://github.com/HaoDK12/DeepOne_model)

### CRISPRon

Xiang, X., Corsi, G. I., Anthon, C. et al.
*Enhancing CRISPR-Cas9 gRNA efficiency prediction by data integration and deep learning.*
**Nature Communications** 12, 3238 (2021).

[Paper](https://doi.org/10.1038/s41467-021-23576-0) ·
[CRISPRon](https://rth.dk/resources/crispr/crispron/)
