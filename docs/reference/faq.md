---
title: FAQ
layout: default
parent: Reference
nav_order: 6
description: Frequently asked questions about FastCCC - how it compares with CellPhoneDB, CellChat and LIANA, how to compare cell-cell communication with healthy tissue, spatial transcriptomics with Squidpy, mouse data, inputs, mapping, run time and memory.
---

# FAQ

**How is FastCCC different from CellPhoneDB, CellChat or LIANA?**
It answers the same question (which ligand–receptor pairs connect which cell types) with the same kind of ligand–receptor database (CellPhoneDB v5 by default; CellChatDB, NicheNet and the LIANA+ resources also ship with the package), but computes the null distribution analytically instead of by permutation, so the p-values are exact and a million cells run in minutes; it scores every interaction with sixteen modules instead of one; and it can compare a dataset against built-in healthy tissue reference panels or a control group, labelling each interaction Up or Down. The paper's benchmark (ten datasets, false-positive rate and F1 against known interactions) is in Fig. 3. A side-by-side comparison is in [FastCCC compared with other tools]({{ site.baseurl }}/use-cases/fastccc-vs-cellphonedb-cellchat-liana.html).

**How do I compare cell–cell communication in my disease sample with healthy tissue?**
`fastccc ref data.h5ad -c cell_type -r liver --auto-map` (replace `liver` by the tissue; `fastccc panels` lists the 19). Review the cell-type mapping with `fastccc map` first for a publication. The result labels every interaction Up or Down versus the healthy panel and the explorer shows the evidence. See [Reference comparison]({{ site.baseurl }}/guides/reference.html).

**Can I use my own control group instead of a built-in panel?**
Yes: `fastccc ref data.h5ad -c cell_type -r condition:control:disease`. The control cells become the panel; no mapping is needed because the labels are the same.

**Do I need raw counts?**
No. The standard analysis normalises raw counts itself and uses log-normalised data as is; the reference workflows rank genes within cells and give identical results for counts, normalised or log-normalised input. Scaled (centred) matrices are rejected; pass the unscaled one.

**Which cell-type column should I use?**
`fastccc inspect` scores the candidates and suggests one; `run` and `ref` still need `-c` explicitly. Use the finest annotation you trust that still has a reasonable number of cells per type.

**My labels do not match the panel's names.**
`fastccc map data.h5ad -c cell_type -r liver -o map.json` proposes a mapping (CD4 T cells → the panel's CD4-positive subtypes, FCGR3A+ monocytes → non-classical monocytes, ...), lists what is ambiguous and what it could not match. Edit the JSON and pass it with `-m`.

**Which of the thousands of significant interactions should I look at?**
Use the explorer's Highlights preset or its sliders: specific (significant in few cell-type pairs), supported (enough expressing cells), robust (most of the 16 modules agree; high bootstrap support with `--bootstrap`), with a large effect size. Specificity alone removes the ubiquitous pairs (APP–CD74, HLA–*) that dominate raw lists.

**How long does it take and how much memory?**
A few thousand cells take under a minute, a hundred thousand a few minutes, a million under an hour on a workstation; memory follows the expression matrix. A standard run keeps its 16 per-module tables for the explorer's module-based features; they can be deleted afterwards to save disk space.

**Can FastCCC run offline?**
Yes, once the panels you need are downloaded (`fastccc panels -d liver`) and, optionally, the enrichment libraries are cached (built the first time with network). The explorer page itself never needs the internet.

**Does it work with mouse data?**
Yes: keep the MGI gene symbols and run with `--db liana_mouseconsensus` (the LIANA+ mouse consensus, 3,989 pairs), or import another mouse database with `fastccc db import ... --species mouse`. The built-in healthy tissue panels are human, so a mouse reference comparison uses a control group of the same dataset (`-r condition:control:treated`) or a panel built from a mouse atlas with `fastccc build-ref --db liana_mouseconsensus`. The explorer skips the human-only gene-set libraries. See [Databases]({{ site.baseurl }}/guides/databases.html).

**Can I use another ligand–receptor database, or my own?**
Twenty-two ship with the package: `fastccc db list --data data.h5ad` shows them with species, size and the fraction of their genes present in your data, and `--db <name>` selects one. A CSV of pairs (complex subunits joined by `_`, optional pathway column) or a directory in CellPhoneDB's five-table layout is installed with `fastccc db import` and then used by name. A reference panel and the query compared with it must use the same database.

**Spatial transcriptomics?**
Run FastCCC on the expression and the cell-type labels, then keep the interactions between cell types that are neighbours in the tissue, with Squidpy's neighbourhood enrichment test:

```python
import scanpy as sc, squidpy as sq, fastccc

adata = sc.read_h5ad("spatial.h5ad")       # expression, obs["cell_type"], coordinates in obsm["spatial"]
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

res = fastccc.run(adata, "cell_type", out_dir="results/spatial")
sig = res.significant                       # which ligand-receptor pairs connect which cell types

sq.gr.spatial_neighbors(adata)              # which cell types are neighbours in the tissue
sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
z = adata.uns["cell_type_nhood_enrichment"]["zscore"]
types = list(adata.obs["cell_type"].cat.categories)
neighbours = {(a, b) for i, a in enumerate(types) for j, b in enumerate(types) if z[i, j] > 2}

near = sig[[pair in neighbours for pair in zip(sig["sender_celltype"], sig["receiver_celltype"])]]
```

`z > 2` keeps the cell-type pairs whose neighbourhood enrichment is about two standard deviations above the permutation expectation; move the cut-off to be stricter or more inclusive. `near` is the table of significant interactions between neighbouring cell types (checked with Squidpy 1.8).

**How do I cite FastCCC?**
Hou, S., Ma, W. & Zhou, X. *Nat Commun* 16, 11428 (2025), [doi:10.1038/s41467-025-66272-z](https://doi.org/10.1038/s41467-025-66272-z). The repository has a `CITATION.cff`.
