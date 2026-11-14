---
title: Inputs
layout: default
parent: Guides
nav_order: 1
description: What FastCCC needs from an .h5ad file, how it detects raw counts versus log-normalised data and the cell-type column, and how to override the detection.
---

# Inputs

## The file

One `.h5ad` (AnnData) file, or an in-memory `AnnData` in Python, with:

- an expression matrix in `adata.X` (counts, normalised or log-normalised; see below),
- a column of `adata.obs` with the cell-type labels, or a two-column TSV (cell id, cell type) passed with `--celltype-file`,
- gene symbols as `adata.var_names` (HGNC symbols for the human databases, MGI symbols for the mouse one); for Ensembl gene ids pass `--gene-id ensembl` (human databases).

Whole-transcriptome matrices (all genes, not only highly variable genes) give better results, especially for reference panels.

## `fastccc inspect`

```bash
fastccc inspect data.h5ad            # human readable
fastccc inspect data.h5ad --json     # for scripts and agents
```

It reports:

| Field | Meaning |
|---|---|
| `matrix.kind` | `raw_counts`, `normalized`, `log1p`, `scaled` or `empty`, with `confidence` and the evidence (integer fraction, per-cell total variability before and after `expm1`) |
| `raw_counts_available` | whether raw counts are in `X`, in `layers['counts']` (also `raw_counts`, `raw`, `X_counts`, `umi`) or in `adata.raw` |
| `celltype_candidates` | categorical columns scored by name (cell_type, annotation, lineage, ...), by their values (T cell, macrophage, fibroblast, ...) and by cardinality; metadata-like names (sample, donor, batch, condition, sex, ...) are penalised |
| `suggested_celltype_key` | the top candidate, with `high` or `low` confidence |
| `recommended_workflow_input` | whether the standard and the reference workflows can use the matrix as is |

`run` and `ref` always require `-c <column>` (or `--celltype-file`): the suggestion is shown, never applied silently.

## Matrix types and what each workflow does

| `adata.X` | Standard analysis (`run`) | Reference workflows (`ref`, `build-ref`) |
|---|---|---|
| raw counts, integer or fractional (SoupX, alevin, RSEM) | `normalize_total(1e4)` + `log1p` applied to a copy | used as is |
| normalised, not logged | `log1p` applied | used as is |
| log-normalised (`log1p`) | used as is | used as is |
| scaled / centred (negative values) | error: pass the unscaled matrix | error |

The reference workflows rank genes within each cell, so per-cell normalisation and `log1p` do not change their results; raw counts are preferred only because they usually carry the full gene set. When `adata.X` is log-normalised and raw counts exist in a layer or in `.raw`, the reference workflows take the raw counts.

If the detection is reported with low confidence and you know what the matrix is, tell `run` explicitly:

```bash
fastccc run data.h5ad -c cell_type --input-type raw        # raw | normalized | log1p
```

Everything that was detected and applied is recorded in `fastccc_run.json` under `input`.

## Cell numbers

For cell types with few cells, the explorer shows cell numbers everywhere (dot sizes, tooltips) and `--bootstrap` reports how stable each interaction is across cell resamples. Tiny clusters can also be merged into their parent type before running.

## Other databases

FastCCC ships 22 ligand–receptor databases: CellPhoneDB v4/v5, CellChatDB, NicheNet, the LIANA+ consensus and the 15 resources it is built from, and a mouse consensus. `fastccc db list` prints them with species, size and whether they carry a pathway classification; add `--data your.h5ad` to see what fraction of each database's genes your data contains (a low value means a species or gene-identifier mismatch). Choose one with `--db <name>` (`database=` in Python). The full table with citations is in the [databases reference](../reference/databases.html).

Your own pair table or a directory in the five-table layout is installed with `fastccc db import` (species, pathway column and provenance are recorded) and then used by name. A reference panel and the query it is compared with must use the same database.
