# DeepOneCalculator

`DeepOneCalculator` packages pretrained DeepOne ensembles for predicting SpCas9-, SpCas9-NG- and SpG-induced,
templated **+1 bp insertion frequencies**. This is the original model's outcome
score, not a general gRNA editing-efficiency probability. Scores retain the training
scale and are neither clipped nor converted to percentages.

## Install

For no GPU device :
```bash
python -m pip install .
```

For device with GPU:

```bash
pip install '.[cuda]'
```

## Quick start

```python
from DeepOneCalculator import SUPPORTED_CAS9_TYPES, SUPPORTED_CELL_LINES, score_grnas

scores = score_grnas(
    sequences=["TTATCTTCGCTATCACCTCCGCCGGGGTCACCCATTAT"],
    cell_lines=["HEK293", "K562"],
    batch_size=256,
    device="auto",
)
print(scores)
# A single sequence and all human cell lines:
all_scores = score_grnas("AAAA" + "ATATATATATATATATATAT" + "AGG" + "AAA")
# All three Cas9 types (variants use one shared ensemble each):
variant_scores = score_grnas(
    "AAAA" + "ATATATATATATATATATAT" + "AGA" + "AAA",
    cas9_types=["SpCas9", "SpCas9-NG", "SpG"],
    device="auto",
)
```

Inputs contain at least **30 bp = 4 bp prefix + 20 nt spacer + 3 bp PAM + 3 bp
suffix**. The 23 bp gRNA field includes the PAM. Longer sequences are scanned for
every complete NGG context for SpCas9, or NGN context for SpCas9-NG/SpG, on both
strands. These PAM rules match the original DeepOne CLI. Lowercase DNA is normalized to uppercase;
ambiguous bases, whitespace and sequences shorter than 30 bp raise `ValueError`.
An empty list or a sequence without a complete PAM site for the selected Cas9 types produces an
empty DataFrame.

The output has one row per eligible guide and model ensemble. It includes the
new `cas9_type` column even when the default SpCas9 selection is used:

| Column | Meaning |
| --- | --- |
| `sequence_id` | One-based input identifier, e.g. `sequence_1` |
| `sequence` | Complete normalized input DNA |
| `gRNA_id` | Input ID, p/m strand, one-based oriented context start; e.g. `sequence_1_p1` |
| `gRNA_seq` | Oriented 20 nt spacer followed by 3 bp PAM |
| `cell_line` | Canonical SpCas9 cell line; missing for SpCas9-NG and SpG |
| `cas9_type` | `SpCas9`, `SpCas9-NG` or `SpG` |
| `deepone_score` | Unrounded arithmetic mean of the five pretrained models |

Rows are grouped by the requested Cas9 order and, for SpCas9, the requested cell
line order. **Within each group, input sequences retain their original order**,
then forward/reverse strand and ascending oriented window position. No score-based
sorting or shuffling occurs. With multiple models the input order repeats for each
model group, so the entire frame is not globally grouped by input sequence.
Inputs without eligible sites have no rows; repeated sequences retain distinct
one-based `sequence_id` values. IDs beyond 9 are kept in numeric input order.
Repeated cell-line and Cas9 selections are deduplicated, preserving first order.

`cas9_types` is keyword-only and defaults to `"SpCas9"`. Use `"all"` for all three;
case-insensitive names and the alias `"SpCas9_NG"` are accepted.

## Supported human cell lines

| Identifier | Models |
| --- | --- |
| `HEK293` | 5 pretrained H5 weights |
| `K562` | 5 pretrained H5 weights |
| `HAP1` | 5 pretrained H5 weights |
| `iPSC` | 5 pretrained H5 weights |
| `RPE-1` | 5 pretrained H5 weights |

`cell_lines="all"` selects these five. HEK293T and HeLa weights are absent from this
repository and are not substituted with HEK293. Nonhuman CHO/mESCs are not packaged.

## Cas9 variant models and cell context

| Cas9 type | PAM | Ensemble | Cell context |
| --- | --- | --- | --- |
| `SpCas9` | NGG | 5 H5 files per selected cell line | Five named human cell lines above |
| `SpCas9-NG` | NGN | 5 H5 files from upstream `variants/ng` | No cell-specific model selection |
| `SpG` | NGN | 5 H5 files from upstream `variants/spg` | No cell-specific model selection |

The upstream variant directories contain **one ensemble per variant**, without
separate cell-type subdirectories. `cell_lines` controls only SpCas9 prediction;
each requested variant is computed once, regardless of that list. Variant rows
have a missing `cell_line` (nullable string dtype; `pandas.isna` is true). This means cell context is
unspecified by the model metadata; it does not claim validation in all cell types.
Cell line names are still validated for variant-only requests.

For example, `cas9_types="all", cell_lines=["HEK293", "K562"]` selects four model
groups: SpCas9/HEK293, SpCas9/K562, SpCas9-NG and SpG. An NGG guide receives scores
from all four; an NGN guide without NGG receives only the two variant scores.

Each packaged H5 file is approximately 1 MB. `models/manifest.json` records their
original filenames and SHA-256 checksums. There are 35 packaged H5 files in total.
`models/provenance.json` pins the variant weights to the upstream Git commit and
records source paths, Git blob hashes and SHA-256 checksums.

## License and citation

The repository is GPL-3.0-only; bundled ViennaRNA retains its upstream license in
`vendor/viennarna/COPYING`. See the vendor provenance file for the exact source
archive and modifications.

Original model: [SeqCrafter/DeepOne_model](https://github.com/HaoDK12/DeepOne_model).
Citation: Hao Yuan, Xiaoguang Pan et al., *Template-free single nucleotide insertion
with deep-learning CRISPR gRNAs for highly efficient gene disruption and systemic
correction of DMD mutations*, 2026 (manuscript under revision).


