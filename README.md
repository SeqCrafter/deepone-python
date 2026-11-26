# DeepOneCalculator
[![Test](https://github.com/SeqCrafter/deepone-python/actions/workflows/ci.yml/badge.svg)](https://github.com/SeqCrafter/deepone-python/actions/workflows/ci.yml)
[![Manylinux wheels](https://github.com/SeqCrafter/deepone-python/actions/workflows/wheels.yml/badge.svg)](https://github.com/SeqCrafter/deepone-python/actions/workflows/wheels.yml)
[![Publish to PyPI](https://github.com/SeqCrafter/deepone-python/actions/workflows/publish.yml/badge.svg)](https://github.com/SeqCrafter/deepone-python/actions/workflows/publish.yml)
![PyPI Version](https://img.shields.io/pypi/v/DeepOneCalculator)

## Installation

Requires Python 3.10 or later.

```bash
python -m pip install DeepOneCalculator

# Optional CUDA dependencies
python -m pip install "DeepOneCalculator[cuda]"
```

To install from a local checkout, run `python -m pip install .`.

## Usage

```python
from DeepOneCalculator import score_grnas

sequences = ["TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT"]

# SpCas9 predictions for selected cell lines
scores = score_grnas(
    sequences,
    cell_lines=["HEK293", "K562"],
    device="auto",
)
print(scores)

# Include SpCas9-NG and SpG
variant_scores = score_grnas(sequences, cas9_types="all")
```

you can also calculate the total energy which also called *CRISPRoff score* by the following codes:

```python
from DeepOneCalculator import calculate_free_energy, fold_energy

guide = "ATCTTCGCTATCACCTCCGCCGG"  # 20 nt spacer + 3 bp PAM
energy = calculate_free_energy(guide)
```
Pass a DNA string or a list of strings containing only A/C/G/T (case insensitive).
Each input must contain at least 30 bp: 4 bp prefix + 20 nt spacer + 3 bp PAM + 3 bp
suffix. Both strands are scanned for NGG (SpCas9) or NGN (SpCas9-NG and SpG) sites.

Supported SpCas9 cell lines are `HEK293`, `K562`, `HAP1`, `iPSC`, and `RPE-1`;
`cell_lines="all"` selects all five by default. `device="auto"` uses CUDA when
available, otherwise CPU.

## Output

Returns a `pandas.DataFrame` with one row per eligible guide and model group.

| Column | Meaning |
| --- | --- |
| `sequence_id` | Input identifier, e.g. `sequence_1` |
| `sequence` | Uppercase input DNA sequence |
| `gRNA_id` | Input ID, strand (`p`/`m`), and one-based oriented context start |
| `gRNA_seq` | 20 nt spacer followed by the 3 bp PAM |
| `cell_line` | SpCas9 cell line; missing for SpCas9-NG and SpG |
| `cas9_type` | `SpCas9`, `SpCas9-NG` or `SpG` |
| `deepone_score` | Mean prediction from five pretrained models |

Scores predict templated **+1 bp insertion frequency** on the original training
scale, without conversion to percentages. SpCas9-NG and SpG use models without
cell-specific weights; `cell_lines` applies only to SpCas9. Inputs with no matching
guides produce no rows.

## Authors and Citation

If you use this package, please cite:

Hao Yuan, *Xiaoguang Pan*,..,Yonglun Luo et al. *Predicting CRISPR-SpCas9-induced single nucleotide insertion by deep learning*. **Nature Communications**, 2026 (manuscript under revision).

Original model: [DeepOne_model](https://github.com/HaoDK12/DeepOne_model).
