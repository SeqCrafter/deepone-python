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

## CRISPRon 编辑效率评分

```python
from DeepOneCalculator import calculate_crispron

sequence = "ACTGAACTTGAAAAGCAAAAAGAAACTGGCCATACTTTCGAAGAAATGCTACTGACTG"
scores = calculate_crispron(sequence, device="cpu")
batch_scores = calculate_crispron([sequence, sequence.lower()], batch_size=256)
print(scores)
```

`calculate_crispron(sequences, batch_size=256, device="auto")` 接受单条 DNA
或字符串列表；仅允许 A/C/G/T（不区分大小写），每条至少 30 bp。
它扫描两条链上具有完整上下文的 NGG 位点：4 bp 前缀 + 20 nt spacer +
3 bp PAM + 3 bp 后缀，不支持仅有 20 nt 或 23 nt 的输入。
`device` 可选 `auto`、`cpu` 或 `cuda`，显式请求不可用的 CUDA 会报错。
六套预训练权重随包提供，运行时无需参考仓库、额外下载或 RNAfold 命令行。

返回 `pandas.DataFrame`，列为 `sequence_id, sequence, gRNA_id, gRNA_seq,
crispron_score`。顺序为输入顺序、正链/反链、定向窗口起点升序，重复输入保留
独立 ID；ID 沿用本包约定，例如 `sequence_1_p1` 表示正链第一个 30 bp
上下文，不等同于原 CRISPRon CLI 的坐标命名。`gRNA_seq` 为 20 nt spacer +
3 bp PAM；空列表或没有完整 NGG 位点时返回相同列的空表。

`crispron_score` 是六套模型输出的等权平均，预测原始训练标度上的
**indel 编辑效率**，不是 DeepOne 的 +1 插入频率，也不是 CRISPRspec 特异性。
不再除以 100、不裁剪、不四舍五入；无细胞系选择参数。输入采用 CRISPRon 的
A/C/G/T 编码及原始 CRISPRoff 能量（NGG 的 PAM 校正均为 1，包括 GGG）。
包内 RNAfold 版本与上游测试的 ViennaRNA 2.2.5 不同时，折叠能和最终分数可能
略有不同；转换后的网络已使用上游固定能量逐模型核对。

CRISPRon 组件及权重采用 **AGPL-3.0-or-later**，其他组件保留原有许可证。
组合分发和网络服务须遵循适用的 AGPL 条款，包括对应源代码提供义务；不能将
整个组合包视为仅有原 GPL 声明。上游许可、作者、修改声明及权重来源记录
随包保存在 `DeepOneCalculator/models/crispron/`。

使用时请引用：Xiang, X., Corsi, G.I., Anthon, C. et al.
*Enhancing CRISPR-Cas9 gRNA efficiency prediction by data integration and deep learning.*
Nature Communications 12, 3238 (2021).
[论文](https://doi.org/10.1038/s41467-021-23576-0) ·
[上游 CRISPRon](https://rth.dk/resources/crispr/crispron/)。

## CRISPRspec Specificity

Genome-wide specificity requires Python 3.11 or later and the optional RIsearch
Python bindings:

```bash
python -m pip install "DeepOneCalculator[specificity]"
```

```python
from DeepOneCalculator import (
    CRISPRspecConfig,
    ReferenceGenome,
    calculate_crisprspec,
)

guide = "GGTGGACAAGCGGCAGATAGCGG"  # 20 nt spacer + 3 bp PAM

# A single call accepts a plain or gzip-compressed DNA FASTA.
score = calculate_crisprspec(guide, "genome.fa", index_path="genome.idx")

# For many guides, reuse the loaded sequences and mmap-backed search index.
reference = ReferenceGenome.from_fasta("genome.fa", index_path="genome.idx")
config = CRISPRspecConfig(
    max_mismatches=6,
    pam_suffixes=("GG", "AG", "GA"),
    ignored_chromosomes=frozenset(),
    threads=4,
)
score = calculate_crisprspec(guide, reference, config=config)
```

The returned `float` is `-log10(Z_off / (Z_off + Z_on))`, where each site
contributes `exp(E / (0.001987 * 310.15))`. Higher scores indicate greater
predicted specificity. Defaults include both strands, NGG/NAG/NGA PAMs and
up to six spacer substitutions, without bulges. An exact on-target must be
present on the included chromosomes; set `require_on_target=False` to add a
theoretical one when absent. Exactly one on-target copy is removed from the
off-target sum; other genomic copies retain their multiplicity. No retained
off-targets gives `math.inf`. Ambiguous reference windows are excluded.

Candidates pass directly from RIsearch's Polars DataFrame to the package's
guide-target energy calculation. RIsearch energies are not used. No result
TSV/gzip files or subprocesses are involved; the current binding still needs
one small temporary query FASTA per call. Reference sequences are held in
memory for PAM extraction, so a reusable `ReferenceGenome` is recommended.
Without `index_path`, the reference owns a temporary index; providing a path
persists it. An existing index must come from the same FASTA and must not be
older than it. Persistent indexes created by `ReferenceGenome.from_fasta` have
a `.deepone.json` source fingerprint and a `.lock` file; keep these beside the
index. Unverified external indexes and changed references are rejected; use a
new index path to rebuild. Concurrent initial builds share the same lock.

The original CRISPRoff PAM corrections (NGG=1, NAG=0.9, NGA=0.8) differ from
DeepOne's training feature table. `CRISPRspecConfig` defaults to
`pam_model="crisproff"`; select `pam_model="deepone"` explicitly to use the
training corrections instead. Existing `calculate_free_energy(guide)` calls
retain their DeepOne behavior. Pairwise energies can also be calculated:

```python
from DeepOneCalculator import calculate_free_energy

off_target = "AGTGGACAAGCGGCAGATAGCGG"
energy = calculate_free_energy(guide, off_target, pam_model="crisproff")
```

## Guide Inputs

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
