---
title: Databases
layout: default
parent: Reference
nav_order: 5
description: The 22 ligand-receptor databases shipped with FastCCC (species, size, pathway annotation, source) and how to add your own.
---

# Ligand–receptor databases

FastCCC ships 22 ligand–receptor databases (21 human, 1 other species). The default is `CPDBv5.0.0`; choose another with `--db <name>` (`database=` in Python). Every workflow (standard, split, reference) accepts every database; a reference panel and the query compared with it must use the same one.

`fastccc db list` prints this table (add `--data your.h5ad` to see the fraction of each database's genes present in your data: a low value means a species or gene-identifier mismatch), `fastccc db show <name>` prints provenance and citation.

| Name | Species | Interactions | Genes | Complexes | Pathway classification | Source |
|---|---|---:|---:|---:|---|---|
| `CPDBv4.1.0` | human | 2,923 | 1,354 | 363 | no | [CellPhoneDB v4.1.0 (cellphonedb-data)](https://github.com/ventolab/cellphonedb-data) |
| `CPDBv5.0.0` | human | 2,911 | 1,354 | 362 | yes | [CellPhoneDB v5.0.0 (cellphonedb-data)](https://github.com/ventolab/cellphonedb-data) |
| `CellChatv1.6.1` | human | 1,103 | 342 | 65 | no | [CellChatDB (human) from CellChat v1.6.1](https://github.com/sqjin/CellChat) |
| `NicheNetv1.1.1` | human | 12,019 | 1,430 | 0 | no | [NicheNet ligand–receptor network, nichenetr v1.1.1](https://github.com/saeyslab/nichenetr) |
| `NicheNetv2.1.5` | human | 4,986 | 2,034 | 0 | no | [NicheNet ligand–receptor network, nichenetr v2.1.5](https://github.com/saeyslab/nichenetr) |
| `liana_baccin2019` | human | 1,648 | 1,037 | 23 | no | [Baccin et al. 2019 ligand–receptor pairs (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_cellcall` | human | 1,123 | 457 | 24 | no | [CellCall (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_cellchatdb` | human | 1,912 | 944 | 151 | no | [CellChatDB (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_cellinker` | human | 3,700 | 1,799 | 135 | no | [Cellinker (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_cellphonedb` | human | 1,223 | 924 | 109 | no | [CellPhoneDB (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_celltalkdb` | human | 3,392 | 1,592 | 1 | no | [CellTalkDB (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_connectomedb2020` | human | 2,266 | 1,412 | 1 | no | [connectomeDB2020 of NATMI (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_consensus` | human | 4,620 | 1,876 | 182 | no | [LIANA+ consensus of its human resources](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_embrace` | human | 1,650 | 935 | 0 | no | [EMBRACE (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_guide2pharma` | human | 665 | 536 | 0 | no | [IUPHAR/BPS Guide to PHARMACOLOGY (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_hpmr` | human | 595 | 661 | 0 | no | [Human Plasma Membrane Receptome (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_icellnet` | human | 738 | 567 | 117 | no | [ICELLNET (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_italk` | human | 2,566 | 1,407 | 1 | no | [iTALK (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_kirouac2010` | human | 152 | 222 | 0 | no | [Kirouac et al. 2010 ligand–receptor pairs (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_lrdb` | human | 3,228 | 1,525 | 1 | no | [LRdb of SingleCellSignalR (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_mouseconsensus` | mouse | 3,989 | 1,615 | 164 | no | [LIANA+ consensus resource translated to mouse orthologs](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |
| `liana_ramilowski2015` | human | 1,889 | 1,227 | 1 | no | [Ramilowski et al. 2015 ligand–receptor pairs (via OmniPath)](https://github.com/scverse/liana/blob/main/src/liana/resource/omni_resource.csv) |

Only CellPhoneDB v5 carries a pathway classification of its interactions; the explorer's *Pathway families* panel falls back to the Hallmark/KEGG/GO gene sets for the others (human only). `liana_*` databases are the resources distributed with [LIANA+](https://github.com/scverse/liana) (curated through OmniPath), converted one to one; their complexes are the `A_B` subunit lists of that table. Please cite the original resource (printed by `fastccc db show`) together with FastCCC.

## Citations

| Name | Citation | License |
|---|---|---|
| `CPDBv4.1.0` | Garcia-Alonso L., Lorenzi V. et al. Single-cell roadmap of human gonadal development. Nature 607, 540–547 (2022); Efremova M. et al. Nat Protoc 15, 1484–1506 (2020). | see source |
| `CPDBv5.0.0` | Troulé K. et al. CellPhoneDB v5: inferring cell–cell communication from single-cell multiomics data. Nat Protoc 20, 3412–3440 (2025); Efremova M. et al. CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes. Nat Protoc 15, 1484–1506 (2020). | see source |
| `CellChatv1.6.1` | Jin S. et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun 12, 1088 (2021). | GPL-3.0 |
| `NicheNetv1.1.1` | Browaeys R., Saelens W., Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods 17, 159–162 (2020). | GPL-3.0 |
| `NicheNetv2.1.5` | Browaeys R., Saelens W., Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods 17, 159–162 (2020); Sang-aram C., Browaeys R., Seurinck R., Saeys Y. Unraveling cell–cell communication with NicheNet by inferring active ligands from transcriptomics data. Nat Protoc 20, 1439–1467 (2025). | GPL-3.0 |
| `liana_baccin2019` | Baccin C. et al. Combined single-cell and spatial transcriptomics reveals the molecular, cellular and spatial bone marrow niche organization. Nat Cell Biol 22, 38–48 (2020). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_cellcall` | Zhang Y. et al. Nucleic Acids Res 49, 8520–8534 (2021). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_cellchatdb` | Jin S. et al. Nat Commun 12, 1088 (2021). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_cellinker` | Zhang Y. et al. Bioinformatics 37, 2025–2032 (2021). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_cellphonedb` | Efremova M. et al. Nat Protoc 15, 1484–1506 (2020). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_celltalkdb` | Shao X. et al. Brief Bioinform 22, bbaa269 (2021). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_connectomedb2020` | Hou R. et al. Nat Commun 11, 5011 (2020). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_consensus` | LIANA+ consensus resource (Dimitrov et al. 2022). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_embrace` | Sheikh B.N. et al. Systematic identification of cell-cell communication networks in the developing brain. iScience 21, 273–287 (2019). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_guide2pharma` | Harding S.D. et al. Nucleic Acids Res 46, D1091–D1106 (2018). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_hpmr` | Ben-Shlomo I. et al. Sci STKE 2003, re9 (2003). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_icellnet` | Noël F. et al. Nat Commun 12, 1089 (2021). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_italk` | Wang Y. et al. bioRxiv 507871 (2019). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_kirouac2010` | Kirouac D.C. et al. Mol Syst Biol 6, 417 (2010). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_lrdb` | Cabello-Aguilar S. et al. Nucleic Acids Res 48, e55 (2020). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_mouseconsensus` | LIANA+ mouse consensus resource (Dimitrov et al. 2022). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |
| `liana_ramilowski2015` | Ramilowski J.A. et al. Nat Commun 6, 7866 (2015). Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022); resources curated via OmniPath (Türei D. et al. Mol Syst Biol 17, e9923, 2021). | see LIANA+/OmniPath terms for this resource |

## Your own database

Any table with one ligand–receptor pair per row can be used. Complex subunits are joined with `_` (or `,`); an optional pathway column becomes the explorer's classification:

```bash
fastccc db import pairs.csv --name mydb --species mouse \
    --ligand-col ligand --receptor-col receptor --classification-col pathway \
    --source "CellChatDB v2 (mouse)" --citation "Jin et al. 2025"
fastccc run data.h5ad -c cell_type --db mydb
```

The database is written to `~/.fastccc/databases/mydb` (`FASTCCC_HOME` moves it) in the five-table layout FastCCC reads (`interaction_table.csv`, `gene_table.csv`, `protein_table.csv`, `complex_table.csv`, `complex_composition_table.csv`) plus a `database.json` with species and provenance. A directory already in that layout (for example one you exported from CellPhoneDB) is imported the same way; `fastccc db validate <dir>` checks it first. Gene symbols must match `adata.var_names` of your data: every workflow stops before scoring when fewer than 5% of the database's genes are present (wrong species or gene identifiers) and warns below 50%; the coverage is recorded in `fastccc_run.json`. The species is recorded so that human-only gene-set libraries are skipped in the explorer for other species.

In Python: `fastccc.dbtools.build_database(pairs_df, "mydb", out_dir, species="mouse")` and `fastccc.list_databases()`.
