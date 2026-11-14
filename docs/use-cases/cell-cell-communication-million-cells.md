---
title: Cell–cell communication for a million cells
layout: default
parent: Use cases
nav_order: 2
description: Cell-cell communication analysis for single-cell datasets of a million cells or more, and a faster alternative when CellPhoneDB's permutation test takes too long - FastCCC computes the same score with an analytic p-value.
---

# Cell–cell communication for a million cells, or when CellPhoneDB is too slow

**Short answer.** FastCCC computes the null distribution of the communication score analytically instead of permuting cell-type labels, so there is no permutation loop to wait for and no number of permutations to choose. One of its scoring modules is the score of CellPhoneDB's statistical method (mean expression, minimum over complex subunits, arithmetic mean of ligand and receptor), and the CellPhoneDB v5.0.0 database ships inside the package.

```bash
pip install fastccc

# the CellPhoneDB score only, with an analytic p-value
fastccc run atlas.h5ad -c cell_type -o results/atlas \
    --summary Mean --aggregation Minimum --combination Arithmetic

# the default: sixteen scoring modules combined with the Cauchy test
fastccc run atlas.h5ad -c cell_type -o results/atlas16
```

## Measured in the paper

All numbers below are from the benchmark of [Hou, Ma & Zhou, *Nat Commun* 16, 11428 (2025)](https://doi.org/10.1038/s41467-025-66272-z), Figs. 2 and 3.

| Question | Result |
|---|---|
| Does the analytic p-value agree with CellPhoneDB's permutation p-value? | Pearson correlation 0.996 to 0.999 across ten datasets. With 10,000 permutations in CellPhoneDB, precision above 99% and intersection over union above 97% for the significant calls. |
| How much faster is the same score? | At least 40 times faster than CellPhoneDB above 100,000 cells; nearly 500 times faster on a two-million-cell dataset, which finished in under 20 minutes. |
| Memory | On the eleven smaller datasets, 69% of CellPhoneDB's on average. About 60% of FastCCC's memory and 25% of its time go into reading the data rather than computing. |
| Datasets above one million cells (AIDA, HLCA) | In the comparison with CellPhoneDB v4 and v5, CellChat, SingleCellSignalR, ICELLNET, NicheNet and CellCall, FastCCC with its default sixteen modules was the only method that completed them. |

The agreement improves as CellPhoneDB uses more permutations, as expected when a permutation estimate converges towards the analytic value.

## What to expect on your machine

A few thousand cells take under a minute, a hundred thousand a few minutes and a million under an hour on a workstation with the default sixteen modules. Memory follows the size of the expression matrix.

## Why sixteen modules instead of one

A single score answers one question. The sixteen modules summarise expression by the mean, the median, the upper quartile or the 90th percentile, aggregate complex subunits by minimum or average, and combine ligand and receptor arithmetically or geometrically; the Cauchy combination gives one p-value. Their pattern is informative: when the median-based modules are significant a typical cell expresses the genes, when only the upper-quantile modules are, a sub-population does. See [How FastCCC works]({{ site.baseurl }}/methods/how-it-works.html).

See [FastCCC compared with other tools]({{ site.baseurl }}/use-cases/fastccc-vs-cellphonedb-cellchat-liana.html) for the side-by-side table.
