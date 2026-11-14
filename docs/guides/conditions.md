---
title: Conditions and samples
layout: default
parent: Guides
nav_order: 3
description: Run FastCCC per condition or per sample with --split-by, compare conditions, measure replication across samples and bootstrap stability.
---

# Conditions and samples

```bash
fastccc run cohort.h5ad -c cell_type --split-by condition -o results/by_condition
fastccc run cohort.h5ad -c cell_type --split-by sample    -o results/by_sample
fastccc run cohort.h5ad -c cell_type --split-by condition -g treated -g control --no-all -o results/pair
```

Every value of the column is analysed as its own dataset (`results/by_condition/<group>/`), plus the full dataset in `all/` unless `--no-all`. Groups with fewer than 20 cells are skipped with a warning and listed under `skipped_groups` in the manifest (`--min-cells` changes the threshold); a handful of cells cannot support per-cell-type summaries and would only add noise to the replication fractions. The matrix is normalised once on the full dataset, then subset. One `explorer.html` covers all groups.

## Two conditions

The explorer gets a *Difference between conditions* heatmap: for every sender–receiver pair, the number of significant interactions in A minus the number in B. Clicking a cell lists the ligand–receptor pairs present only in A, only in B, and shared. This is a descriptive comparison of two independent runs; for a statistical Up/Down call, make the control condition the reference instead:

```bash
fastccc ref cohort.h5ad -c cell_type -r condition:control:treated -o results/treated_vs_control
```

See [Reference comparison]({{ site.baseurl }}/guides/reference.html).

## Three or more groups: replication

With samples or donors as groups, the explorer adds, for every interaction of the full dataset, the fraction of samples in which it is significant on its own (`Replication` column and filter), and a *Samples × interactions* panel: the most recurrent interactions as columns, samples as rows ordered by the first principal component of their communication profile, plus a PCA of the samples. Interactions that replicate across most samples are the ones to build on; those found only in the pooled data may be driven by one sample.

<p align="center">
  <img src="{{ site.baseurl }}/images/explorer/samples.png" width="860" alt="Samples x interactions matrix and sample PCA">
</p>

Samples need enough cells of each type to produce significant interactions on their own; with fewer than a hundred cells per sample the per-sample runs are sparse and replication fractions are low by construction.

## Bootstrap stability

```bash
fastccc run data.h5ad -c cell_type --bootstrap 10 -o results/sample
```

Re-runs the analysis on 10 subsamples that keep 80% of the cells of every cell type (`--bootstrap-fraction`), and writes `bootstrap_support.tsv`: for every significant interaction of the main run, the fraction of subsamples in which it stayed significant. Run time is about N + 1 times a single run; the per-module tables of the resamples are deleted, only their significant sets are kept. The explorer shows a *Stability* histogram, a support filter, and the Highlights preset requires support of at least 80% when bootstrap data are present.

## Python

```python
multi = fastccc.run("cohort.h5ad", "cell_type", split_by="condition", out_dir="results/by_condition")
multi.groups            # ['control', 'treated']
multi["treated"].summary()
multi.all.summary()
multi.explorer_path
res = fastccc.run("data.h5ad", "cell_type", out_dir="results/sample", bootstrap=10)
res.manifest["bootstrap"]        # n, fraction, n_stable_80, ...
```
