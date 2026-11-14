---
title: Toolkits
layout: default
nav_order: 3
parent: Legacy API
description: The building blocks behind FastCCC that are usable on their own - ligand-receptor database tools, gene identifier handling, the analytic distribution class, input inspection, order statistics, the DEG filter and the result explorer.
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# FastCCC Toolkits

FastCCC is built from a few components that are useful on their own: tools for ligand–receptor databases, the analytic distribution arithmetic that replaces permutations, input inspection, and the result explorer. This page shows how to call each of them directly. Everything here is part of the installed package; nothing needs the repository checkout.

## Build your own LRIs database
A ligand–receptor database is a folder of five CSV tables in CellPhoneDB's layout (`interaction_table.csv`, `gene_table.csv`, `protein_table.csv`, `complex_table.csv`, `complex_composition_table.csv`) plus a `database.json` with species and provenance. Twenty-two such databases ship with the package (`fastccc db list`), and any table of ligand–receptor pairs can be turned into one:

```bash
fastccc db import pairs.csv --name mydb --species human --classification-col pathway \
    --source "my curation" --citation "Doe et al. 2026"
fastccc db validate ~/.fastccc/databases/mydb
fastccc run data.h5ad -c cell_type --db mydb
```

`pairs.csv` needs a `ligand` and a `receptor` column (other names via `--ligand-col` / `--receptor-col`); complex subunits are joined with `_` or `,` (`ITGA1_ITGB1`). An optional pathway column becomes the explorer's classification, `--id-col` keeps your own interaction ids, `--source-col` an evidence column. The same from Python:

```python
import pandas as pd
from fastccc import dbtools

pairs = pd.DataFrame({'ligand': ['TGFB1', 'IL6', 'CXCL12'],
                      'receptor': ['TGFBR1_TGFBR2', 'IL6R_IL6ST', 'CXCR4'],
                      'classification': ['TGFb', 'IL6', 'CXCL']})
path = dbtools.build_database(pairs, 'mydb', out_dir='databases', species='human')
dbtools.validate_database(path)          # {'ok': True, 'errors': [], 'warnings': [...], 'summary': {...}}
dbtools.import_database(path)            # copies it into ~/.fastccc/databases so that --db mydb finds it
```

A folder already in the five-table layout (for example one exported from CellPhoneDB's own database generator) is imported the same way. Details, the list of packaged databases and their citations: [Databases guide]({{ site.baseurl }}/guides/databases.html), [databases reference]({{ site.baseurl }}/reference/databases.html).

## Ensembl and HGNC symbols converter
The human databases carry both HGNC symbols and Ensembl gene ids (`gene_table.csv`, columns `hgnc_symbol` and `ensembl`), and the package ships the mapping they were built with, `fastccc/data/db/ensembl_hgnc_Homo_sapiens.GRCh38.112.csv` (42,670 Ensembl ids, GRCh38 release 112).

Data whose `var_names` are Ensembl ids does not need to be renamed: tell FastCCC what the identifiers are and it matches them against the Ensembl column of the database.

```bash
fastccc run data.h5ad -c cell_type --gene-id ensembl
fastccc db list --data data.h5ad        # detects Ensembl ids and reports the coverage of every database
```

```python
res = fastccc.run("data.h5ad", "cell_type", convert_type="ensembl")
```

To translate identifiers yourself, in either direction:

```python
from fastccc.dbtools import load_ensembl_map

sym2ens = load_ensembl_map()                                   # {'TGFB1': 'ENSG00000105329', ...}
ens2sym = {v: k for k, v in sym2ens.items()}
adata.var['symbol'] = [ens2sym.get(e.split('.')[0]) for e in adata.var_names]   # version suffixes stripped
adata = adata[:, adata.var['symbol'].notna()].copy()
adata.var_names = adata.var['symbol'].values
adata.var_names_make_unique()
```

The mapping is one symbol per Ensembl id; genes absent from it (about 40% of a whole-transcriptome matrix, mostly non-coding) are dropped by the snippet above and simply not scored when `--gene-id ensembl` is used. Mouse databases carry MGI symbols only.

## How to use FastCCC's probability tools for calculating distributions of random variables after algebraic operations
`fastccc.distrib.Distribution` adds, scales, multiplies and convolves distributions (sums and means of $$n$$ i.i.d. draws, complex minimum and average, arithmetic and geometric ligand–receptor combination) and `get_pvalue_from_pmf` reads exact tail probabilities. The operators, with examples, are documented in the [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html).

## How to use FastCCC's tools to quantify and normalize scRNA-seq data
`fastccc inspect` reports what an `.h5ad` holds; the same functions are available in Python and are what `fastccc.run` applies before scoring:

```python
from fastccc.datacheck import detect_matrix_type, find_raw_counts, prepare_for_standard, prepare_for_reference, suggest_celltype_key

info = detect_matrix_type(adata.X)          # kind: raw | normalized | log1p | scaled, with confidence and evidence
raw = find_raw_counts(adata)                # (AnnData, MatrixInfo) from .raw or a counts layer, or None
prepared, info, applied = prepare_for_standard(adata)    # normalize_total(1e4) + log1p for raw counts; log1p data used as is
prepared, info, applied = prepare_for_reference(adata)   # rank-based workflows: counts, normalised and log1p give the same result
key, confidence, candidates = suggest_celltype_key(adata.obs)
```

```bash
fastccc inspect data.h5ad                   # the same, as a report; --json for machine-readable output
```

`--input-type raw|normalized|log1p` overrides the detection when you know what the matrix is. See [Inputs]({{ site.baseurl }}/guides/inputs.html).

## How to use FastCCC's tools to compute the distribution of the $$k^{th}$$ order statistic for an arbitrary discrete distribution of $$n$$ random variables
The textbook order-statistic formula is for continuous variables and fails on discrete data, where ties among the $$n$$ draws are the rule. FastCCC computes the exact distribution through the cdf (for the values 0, 1, 2 with probabilities 0.5, 0.3, 0.2 the 5th of 10 draws has pmf 0.6230, 0.3706, 0.0064, which the continuous formula cannot produce). `fastccc.distrib.order_statistic_pmf(values, probs, n, k)` gives it for any support and any $$n$$; the derivation, the worked example and the routine used inside FastCCC are in the [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html).

## How to use FastCCC to get DEGs
FastCCC does not produce a differential-expression table. What `use_DEG=True` does is tighten the test: an interaction that is significant by the Cauchy combination is kept only if, in addition, the ligand in the sender and the receptor in the receiver are each individually significant against their own null distributions (the cell-type summary of the gene compared with the distribution of that summary under random labels), both at $$p < 0.05$$. It is a "both sides must stand out" filter rather than a gene-level DEG analysis.

```python
res = fastccc.run("data.h5ad", "cell_type", use_DEG=True)
```

For per-gene differential expression use scanpy's `rank_genes_groups`; its results can be used to restrict FastCCC's test to a chosen set of interactions through `select_list` / `filter_`.

## Plot tools
Plots are no longer made with R scripts. Every run writes `explorer.html`, an interactive page with all figures (heatmaps, chord diagram, dot plots, flows, communication patterns, pathway families, evidence panels, violins and embeddings); each chart exports SVG, and the figure builder composes panels into one lettered figure. See [Result explorer]({{ site.baseurl }}/guides/explorer.html). The older static report is still available with `fastccc run --legacy-report` or `fastccc report <result_dir>`.
