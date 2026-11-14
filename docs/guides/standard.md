---
title: One dataset
layout: default
parent: Guides
nav_order: 2
description: The standard FastCCC analysis of one dataset - what it computes, which files it writes, and how to pick the interactions that matter.
---

# One dataset

```bash
fastccc run data.h5ad -c cell_type -o results/sample --name "Sample 1"
```

```python
res = fastccc.run("data.h5ad", "cell_type", out_dir="results/sample", name="Sample 1")
```

## What is computed

For every ordered pair of cell types (sender, receiver; autocrine pairs included) and every ligand–receptor pair of the database whose ligand is expressed in at least 10% of the sender cells and whose receptor is expressed in at least 10% of the receiver cells (`--min-percentile`), FastCCC computes a communication score and its p-value under random cell-type labels. Sixteen scoring modules are run (expression summarised by mean, median, upper quartile or 90th percentile; complexes aggregated by minimum or average; the two sides combined arithmetically or geometrically) and their p-values are combined with the Cauchy combination test. Interactions with a combined p below 0.05 are reported as significant. The [methods page]({{ site.baseurl }}/methods/how-it-works.html) explains how the null is obtained without permutations.

## Outputs

In `results/sample/`:

| File | Content |
|---|---|
| `explorer.html` | the interactive result page |
| `fastccc_run.json` | the run manifest: inputs, detected matrix type, parameters, timing, output paths |
| `<task>_significant_results.tsv` | one row per significant sender, receiver and ligand–receptor pair with its Cauchy p-value |
| `<task>_Cauchy_pvals.tsv` | the full p-value matrix (rows `sender|receiver`, columns interaction ids) |
| `<task>_average_interactions_strength.tsv` | the communication-score matrix averaged over the modules |
| `<task>_<timestamp>_<module>_{pvals,interactions_strength,percents_analysis}.tsv` | per-module intermediates (16 × 3 files; large) |

`<task>` is a six-character id, also stored in the manifest. The per-module tables are what makes a run several gigabytes for large datasets; keep them if you want the explorer's module-based *signal type* and bootstrap features, delete them otherwise.

## Reading the result

A dataset usually has thousands of significant interactions. The explorer gives you four handles to pick the ones that matter, all also exported in the table:

- **Specificity** (`Pairs`): in how many cell-type pairs the ligand–receptor pair is significant. Ubiquitous pairs (APP–CD74, HLA-*) score in dozens of pairs; the interesting ones in a few.
- **Robustness** (`Modules`): how many of the 16 modules call the interaction on their own. The pattern also classifies the **signal type**: *broad* when the median-based modules agree (a typical cell expresses the ligand), *sub-population* when only the upper-quantile modules do (a fraction of the cells carries the signal), plus *complex-limited* (one subunit of a complex limits it) and *one-sided* (ligand and receptor levels are unbalanced).
- **Support**: the fraction of expressing cells on both sides.
- **Effect size**: for the mean-based module, how many null standard deviations the score sits above its expectation. This is the ranking that p-values at their floor cannot give.

The **Highlights** preset combines these with thresholds taken from the dataset's own distributions and relaxes them until at least 20 interactions remain. Add `--bootstrap 10` to measure, for every interaction, how often it stays significant when 80% of the cells of every type are resampled (`bootstrap_support.tsv`, shown in the *Stability* panel and used by Highlights).

## Python result object

```python
res.summary()          # dict: n_celltypes, n_significant, top senders / receivers / LR pairs
res.significant        # DataFrame of the significant table
res.pvals              # Cauchy p-value matrix
res.strength           # average communication-score matrix
res.celltypes          # list
res.explorer_path      # Path to explorer.html
res.explore()          # rebuild the explorer (options: cells=False, max_cells=2000, genesets=False, n_patterns=4)
res.report()           # the older static HTML report
fastccc.load_result("results/sample")   # the same object from a folder
```

## Options

| Option | Default | Purpose |
|---|---|---|
| `--min-percentile` | 0.1 | expression filter (fraction of expressing cells on each side) |
| `--input-type` | detected | state what `adata.X` holds |
| `--bootstrap N` | 0 | bootstrap support over N resamples |
| `--no-cells` | | smaller explorer without per-cell violins and embeddings |
| `--legacy-report` | | also render the older static HTML report |
| `--db` | `CPDBv5.0.0` | one of the 22 packaged databases or an imported one ([Databases]({{ site.baseurl }}/guides/databases.html)) |

Python: `run(..., single_unit_summary_list=..., complex_aggregation_list=..., LR_combination_list=..., use_DEG=..., select_list=..., filter_=...)` exposes the module selection of the original function.
