---
title: Quick start
layout: default
nav_order: 3
description: From an .h5ad file to an interactive cell-cell communication result in three commands, plus a comparison against a healthy tissue reference panel.
---

# Quick start

You need an `.h5ad` file with a cell-type annotation in `adata.obs` and gene symbols as `var_names`. Counts, normalised or log-normalised matrices are all fine. If you want something to try, the [CellphoneDB tutorial dataset](https://github.com/ventolab/CellphoneDB/raw/v5.0.1/notebooks/data_tutorial.zip) (placenta/decidua, 3,312 cells, column `cell_labels`) and scanpy's `pbmc3k` are what the test suite uses.

## 1. Look at the file

```bash
fastccc inspect data.h5ad
```

```text
data.h5ad: 3,312 cells x 30,800 genes
adata.X looks like: log1p (confidence high: per-cell totals are constant after expm1 (normalize_total + log1p))
raw counts available: True (from raw)
candidate cell-type columns (score, unique values):
  cell_labels                      6.3    42  name mentions cell type; 33% of values look like cell types
suggested celltype_key: cell_labels (confidence high)
standard workflow input: ok
reference workflow input: ok
```

`inspect` tells you what the matrix is, whether raw counts are available, which column holds the cell types, and whether anything needs attention. Nothing is run.

## 2. Cell–cell communication of the dataset

```bash
fastccc run data.h5ad -c cell_labels -o results/placenta
```

Raw counts are normalised automatically (`normalize_total(1e4)` + `log1p`); log-normalised input is used as is. A few minutes later:

```text
RunResult(name='placenta', task_id='8f3a1c', n_celltypes=42, n_significant=8768, result_dir='results/placenta')
explorer: results/placenta/explorer.html   (open in a browser)
```

`--db liana_consensus` switches the ligand–receptor database; 22 ship with the package, including a mouse one ([Databases]({{ site.baseurl }}/guides/databases.html)).

Open `explorer.html`. Start with the **Highlights** button: it keeps the interactions that are specific (significant in few cell-type pairs), well supported by expression and robust across the scoring modules, with thresholds chosen from this dataset's distributions. Click any dot or row to see its evidence.

## 3. Compare against healthy tissue

```bash
fastccc map data.h5ad -c cell_labels -r blood -o map.json       # how do my labels match the panel's?
fastccc ref data.h5ad -c cell_labels -r blood -m map.json -o results/vs_blood
```

Built-in panels use Cell Ontology names, so `fastccc map` proposes a mapping (review it; the ambiguous ones are listed) and `ref` uses it. For a quick look, `--auto-map` applies the suggestion directly. The result labels every interaction:

```text
ReferenceResult(query='placenta', reference='blood', Up=57, Down=610, BothSig=494, BothNS=232)
cell types: 8/8 query cell types covered by the reference (100.0% of cells); mapping: user
explorer: results/vs_blood/explorer.html
```

- **Up**: significant in your data, not in the reference.
- **Down**: significant in the reference, not in your data.
- **Both Sig / Both NS**: the same in both.

The explorer's *Comparability* panel compares the depth of your data with the panel's, and the *Evidence* panel shows, for any interaction, the reference-derived null distributions and where your data falls.

A control group inside the same dataset works the same way, without any mapping:

```bash
fastccc ref data.h5ad -c cell_labels -r condition:control:disease -o results/disease_vs_control
```

## 4. Conditions or samples

```bash
fastccc run data.h5ad -c cell_labels --split-by condition -o results/by_condition
fastccc run data.h5ad -c cell_labels --split-by sample -o results/by_sample
```

Every group is analysed separately (plus the full dataset). With two groups the explorer adds a *difference between conditions* heatmap with gained and lost pairs; with three or more it adds *replication*: how many samples each interaction is significant in, and a samples × interactions matrix.

## 5. The same from Python

```python
import fastccc

info = fastccc.inspect("data.h5ad")
res  = fastccc.run("data.h5ad", "cell_labels", out_dir="results/placenta")
res.summary(); res.significant.head(); res.explorer_path

sug = fastccc.suggest_celltype_mapping("data.h5ad", "cell_labels", reference="blood")
print(sug.format())
q = fastccc.compare_to_reference("data.h5ad", "cell_labels", reference="blood",
                                 celltype_mapping=sug.mapping, out_dir="results/vs_blood")
q.trend_counts; q.interaction_summary("Down"); q.celltype_coverage

later = fastccc.load_result("results/placenta")    # any result folder, any time
```

## Next

- [Inputs]({{ site.baseurl }}/guides/inputs.html): matrix types, cell-type columns, overriding the detection.
- [Reference comparison]({{ site.baseurl }}/guides/reference.html): panels, mapping review, reading Down calls.
- [Result explorer]({{ site.baseurl }}/guides/explorer.html): every panel explained.
- [Command line]({{ site.baseurl }}/reference/cli.html) and [Python API]({{ site.baseurl }}/reference/api.html).
