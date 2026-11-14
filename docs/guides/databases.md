---
title: Databases
layout: default
parent: Guides
nav_order: 5
description: "Choosing among the 22 packaged ligand-receptor databases, checking which one fits your data, analysing mouse data, and importing your own pairs."
---

# Ligand–receptor databases

Every FastCCC analysis scores the ligand–receptor pairs of one database. Twenty-two ship inside the package, so nothing has to be downloaded or converted. CellPhoneDB v5 is the default: curated complexes, signalling direction and a pathway classification that the explorer uses. The others are there for comparison with published work, for broader coverage, or for mouse data. The full table with sizes, citations and licenses is in the [databases reference]({{ site.baseurl }}/reference/databases.html).

## Which one

| You want | Use |
|---|---|
| the curated human default with pathway families | `CPDBv5.0.0` |
| results comparable to a CellChat analysis | `CellChatv1.6.1` (CellChatDB v1) or `liana_cellchatdb` |
| the broadest human coverage | `liana_consensus` (pairs supported by several resources) or `NicheNetv2.1.5` |
| mouse data | `liana_mouseconsensus` |
| the resource a paper used | the matching `liana_*` database (CellTalkDB, connectomeDB2020, ICELLNET, ...) |

See them all, with the fraction of each database's genes present in your file:

```bash
fastccc db list --data data.h5ad
```

```text
  name                     species  interactions  genes complexes  pathways origin   in data
  CPDBv5.0.0               human           2,911  1,354       362  yes      packaged    97%
  liana_consensus          human           4,620  1,876       182  no       packaged    96%
  liana_mouseconsensus     mouse           3,989  1,615       164  no       packaged     1%
  ...
```

A low *in data* value means the species or the gene identifiers do not match (here a human dataset against the mouse database). Every workflow checks this before scoring and stops below 5%.

`fastccc db show liana_consensus` prints the provenance, citation and license of one database.

## Running with another database

```bash
fastccc run data.h5ad -c cell_type --db liana_consensus -o results/consensus
fastccc run data.h5ad -c cell_type --db liana_consensus --split-by condition -o results/by_condition
```

```python
res = fastccc.run("data.h5ad", "cell_type", database="liana_consensus", out_dir="results/consensus")
```

Reference comparisons use the database of the panel. The built-in healthy tissue panels were built with `CPDBv5.0.0`, so they are compared with that database; a panel you build yourself takes any database, and the query must then use the same one (FastCCC refuses a mismatch):

```bash
fastccc build-ref atlas.h5ad -c cell_type -n healthy_pancreas --db liana_consensus -o panels
fastccc ref query.h5ad -c cell_type -r panels/healthy_pancreas --db liana_consensus -o results/vs_pancreas
```

Choose one database per project and state it in the methods (the explorer's methods paragraph names it): a larger database tests more pairs and reports more significant interactions.

## Mouse data

Mouse datasets use MGI gene symbols (`Il2`, `Notch1`) and the mouse database:

```bash
fastccc run mouse.h5ad -c cell_type --db liana_mouseconsensus -o results/mouse
fastccc ref mouse.h5ad -c cell_type -r condition:control:treated --db liana_mouseconsensus -o results/treated_vs_control
```

For mouse data, the reference is a control group of the same dataset (`-r condition:control:treated`) or a mouse panel built with `fastccc build-ref --db liana_mouseconsensus`. The explorer works the same; its gene-set groupings (Hallmark, KEGG, GO) apply to the human databases.

## What the database changes in the explorer

- **Pathway families** groups interactions by the database's own classification when it has one (CellPhoneDB v5), otherwise by Hallmark / KEGG / GO gene sets (human databases).
- **Key findings** and the **methods paragraph** name the database, its species, size and citation; the **About** tab shows its full provenance.
- Interaction ids follow the database (`CPI-...` for CellPhoneDB, `LIGAND__RECEPTOR` for the converted ones).

## Your own database

Any table with one ligand–receptor pair per row can be imported. Complex subunits are joined with `_` (or `,`); an optional pathway column becomes the explorer's classification:

```text
ligand,receptor,pathway
TGFB1,TGFBR1_TGFBR2,TGFb
IL6,IL6R_IL6ST,IL6
CXCL12,CXCR4,CXCL
```

```bash
fastccc db import pairs.csv --name mydb --species human --classification-col pathway \
    --source "my curation 2026" --citation "Doe et al. 2026"
fastccc run data.h5ad -c cell_type --db mydb -o results/mydb
```

The database is written to `~/.fastccc/databases/mydb` (`FASTCCC_HOME` moves it) in the five-table layout FastCCC reads (`interaction_table.csv`, `gene_table.csv`, `protein_table.csv`, `complex_table.csv`, `complex_composition_table.csv`) plus a `database.json` with species and provenance. Column names other than `ligand` / `receptor` are given with `--ligand-col` / `--receptor-col`; `--id-col` keeps your own interaction ids and `--source-col` an evidence column. A directory already in that layout, for example one generated by CellPhoneDB's database tools, is imported the same way, and `fastccc db validate <dir>` checks it first.

In Python:

```python
from fastccc import dbtools
dbtools.build_database(pairs_df, "mydb", out_dir, species="mouse")        # columns ligand, receptor, [classification]
dbtools.import_database("pairs.csv", name="mydb", species="mouse",
                        pairs_kwargs={"classification_col": "pathway"})
dbtools.validate_database(path)                                           # {"ok", "errors", "warnings", "summary"}
fastccc.list_databases()                                                  # packaged and imported names
```

## Citing

Cite the database together with FastCCC: `fastccc db show <name>` prints the reference, and the [databases reference]({{ site.baseurl }}/reference/databases.html) lists them all. The `liana_*` databases are the resources distributed with LIANA+ (Dimitrov et al. 2022) and curated through OmniPath (Türei et al. 2021); their original publications are listed there too.
