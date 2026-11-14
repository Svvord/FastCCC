---
title: Reference panels
layout: default
parent: Reference
nav_order: 4
description: The 19 built-in healthy human tissue reference panels (cell types, cell numbers, source).
---

# Built-in healthy tissue reference panels

The panels were built from healthy human tissue data in CELLxGENE Census (Cell Ontology cell-type names). Together they hold 15,752,932 cells across 1443 cell-type entries. A panel is downloaded on first use (about 3-5 MB) and activated once (its null distributions are computed, under a minute).

Use `fastccc panels` to list them, `fastccc panels --show <tissue>` to print a panel's cell types and cell numbers, and `fastccc map` to match your labels to them.

| Tissue | Cell types | Cells | Largest cell types | Database |
|---|---:|---:|---|---|
| `adipose` | 21 | 206,363 | preadipocyte, mesothelial cell, adipocyte | CPDBv5.0.0 |
| `blood` | 61 | 494,519 | classical monocyte, natural killer cell, naive thymus-derived CD4-positive, alpha-beta T cell | CPDBv5.0.0 |
| `bone_marrow` | 79 | 355,797 | monocyte, erythrocyte, promonocyte | CPDBv5.0.0 |
| `breast` | 50 | 3,016,499 | fibroblast of mammary gland, luminal adaptive secretory precursor cell of mammary gland, luminal hormone-sensing cell of mammary gland | CPDBv5.0.0 |
| `cerebral_cortex` | 12 | 2,714,612 | neuron, oligodendrocyte, astrocyte | CPDBv5.0.0 |
| `colon` | 98 | 174,817 | enterocyte, fibroblast, plasma cell | CPDBv5.0.0 |
| `eye` | 82 | 1,972,312 | retinal rod cell, GABAergic amacrine cell, amacrine cell | CPDBv5.0.0 |
| `heart` | 50 | 1,492,187 | cardiac muscle cell, regular ventricular cardiac myocyte, mural cell | CPDBv5.0.0 |
| `kidney` | 114 | 427,896 | epithelial cell of proximal tubule, kidney loop of Henle thick ascending limb epithelial cell, endothelial cell | CPDBv5.0.0 |
| `liver` | 113 | 810,458 | erythrocyte, macrophage, natural killer cell | CPDBv5.0.0 |
| `lung` | 170 | 1,673,947 | type II pneumocyte, alveolar macrophage, macrophage | CPDBv5.0.0 |
| `lymph_node` | 85 | 214,568 | B cell, naive thymus-derived CD4-positive, alpha-beta T cell, memory B cell | CPDBv5.0.0 |
| `midbrain` | 13 | 498,610 | neuron, oligodendrocyte, central nervous system macrophage | CPDBv5.0.0 |
| `prostate_gland` | 41 | 150,605 | basal cell of prostate epithelium, epithelial cell of urethra, epithelial cell | CPDBv5.0.0 |
| `skin` | 89 | 307,335 | fibroblast, macrophage, natural killer cell | CPDBv5.0.0 |
| `small_intestine` | 130 | 241,338 | mesodermal cell, enterocyte, stromal cell | CPDBv5.0.0 |
| `spleen` | 90 | 380,855 | natural killer cell, macrophage, fibroblast | CPDBv5.0.0 |
| `stomach` | 44 | 55,474 | foveolar cell of stomach, activated CD8-positive, alpha-beta T cell, fibroblast | CPDBv5.0.0 |
| `thymus` | 101 | 564,740 | double-positive, alpha-beta thymocyte, double negative thymocyte, fibroblast | CPDBv5.0.0 |

Cell-type names follow the Cell Ontology (for example `CD4-positive, alpha-beta T cell`); `fastccc map` translates common labels (CD4 T cells, Tregs, pDC, FCGR3A+ monocytes, ...) automatically and reports what it could not match.
