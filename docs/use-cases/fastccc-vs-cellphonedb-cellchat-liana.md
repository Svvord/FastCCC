---
title: FastCCC vs CellPhoneDB, CellChat and LIANA+
layout: default
parent: Use cases
nav_order: 3
description: FastCCC compared with CellPhoneDB, CellChat, LIANA+ and MultiNicheNet for cell-cell communication analysis of single-cell RNA-seq - null distribution, scale, scoring, healthy tissue reference panels, databases, and when to choose FastCCC.
---

# FastCCC compared with CellPhoneDB, CellChat, LIANA+ and MultiNicheNet

**Short answer.** FastCCC answers what the other tools answer, which ligand–receptor pairs connect which cell types, and adds three things: it tells you what differs from healthy tissue without controls of your own, it has no permutation test and so runs on millions of cells, and it runs from one command or from a coding agent that does the work for you.

## What each tool is, in its authors' words

- **[CellPhoneDB](https://github.com/ventolab/CellphoneDB)** (Python): "a publicly available repository of HUMAN curated receptors, ligands and their interactions paired with a tool to interrogate your own single-cell transcriptomics data".
- **[CellChat](https://github.com/jinworks/CellChat)** (R): can "quantitatively characterize and compare the inferred cell-cell communication networks" by social network analysis, pattern recognition and manifold learning.
- **[LIANA+](https://github.com/saezlab/liana-py)** (Python, scverse): "a scalable framework that adapts and extends existing methods and knowledge to study cell-cell communication".
- **[MultiNicheNet](https://github.com/saeyslab/multinichenetr)** (R): "differential cell-cell communication analysis from single-cell transcriptomics data with complex multi-sample, multi-condition designs".
- **[FastCCC](https://github.com/Svvord/FastCCC)** (Python, command line): permutation-free ligand–receptor analysis with analytic p-values, sixteen scoring modules, and comparison against healthy tissue reference panels.

## Side by side

| | FastCCC | CellPhoneDB | CellChat | LIANA+ |
|---|---|---|---|---|
| Null distribution | analytic (FFT convolution), no permutations | permutation | permutation | per method |
| Datasets above a million cells | completed the 1M-cell AIDA and 2M-cell HLCA atlases in the paper's benchmark | not completed in that benchmark | not completed in that benchmark | depends on the wrapped method |
| Scoring | 16 modules combined with the Cauchy test; their pattern tells broad from sub-population signals | one score | one score | consensus of methods |
| Compare with healthy tissue | built-in reference panels for 19 tissues; every interaction Up, Down or unchanged | – | – | – |
| Two conditions that are both in your data | control group as reference (`-r condition:control:disease`) | DEG-based method with differential genes you supply | comparative analysis across conditions | targeted differential analysis, MOFA+, Tensor-cell2cell |
| Robustness | bootstrap support over cell resamples; replication across samples; effect sizes | – | – | – |
| Ligand–receptor databases | 22 packaged, including CellPhoneDB v4 and v5, CellChatDB, NicheNet and the LIANA+ consensus; import your own | CellPhoneDB | CellChatDB | OmniPath resources |
| Output | interactive self-contained HTML, TSV tables, JSON manifest | tables | R objects and plots | tables |
| Language and interface | Python, command line, agent skill | Python | R | Python |

The benchmark is in [Hou, Ma & Zhou, *Nat Commun* 16, 11428 (2025)](https://doi.org/10.1038/s41467-025-66272-z), Fig. 3: it compared FastCCC with CellPhoneDB v4 and v5, CellChat, SingleCellSignalR, ICELLNET, NicheNet and CellCall on ten datasets.

## When to choose FastCCC

**"Which communications changed compared with normal tissue, and I have no healthy controls."** `fastccc ref data.h5ad -c cell_type -r liver --auto-map` labels every interaction Up or Down against a healthy tissue panel. See [Compare with healthy tissue]({{ site.baseurl }}/use-cases/compare-cell-cell-communication-with-healthy-tissue.html).

**"My dataset has hundreds of thousands or millions of cells."** In the paper's benchmark FastCCC was the only method that completed the two atlases above one million cells. See [A million cells or more]({{ site.baseurl }}/use-cases/cell-cell-communication-million-cells.html).

**"I want CellPhoneDB's score and database, faster."** `--summary Mean --aggregation Minimum --combination Arithmetic` computes that score with an analytic p-value; in the paper the p-values correlate at 0.996 to 0.999 with CellPhoneDB's, and FastCCC is at least 40 times faster above 100,000 cells.

**"I work in R and Seurat, not Python."** You do not have to write Python. FastCCC ships an [agent skill]({{ site.baseurl }}/guides/agent.html): tell Claude Code, Codex or another coding agent what you want in plain language, and it checks the input, runs the analysis and summarises the findings. FastCCC reads `.h5ad`, so ask the agent to convert your Seurat object first.

**"I want signalling pathways, sender and receiver roles and communication patterns."** Every run writes one interactive page with pathway families per cell type, pathway enrichment, signalling roles (outgoing versus incoming), sender → ligand–receptor → receiver flows and communication patterns, each filterable and exportable as SVG. See the [result explorer]({{ site.baseurl }}/guides/explorer.html).

**"I have several samples or patients per condition."** `fastccc run cohort.h5ad -c cell_type --split-by sample` runs every sample and shows which interactions replicate across samples; `fastccc ref cohort.h5ad -c cell_type -r condition:control:disease` labels each interaction Up or Down against the control group. See [Conditions and samples]({{ site.baseurl }}/guides/conditions.html).

**"I have spatial transcriptomics data."** Run FastCCC on the expression and the cell-type labels, then keep the interactions between cell types that are neighbours in the tissue with Squidpy's neighbourhood enrichment test. The few lines of code are in the [FAQ]({{ site.baseurl }}/reference/faq.html).

**"My data is mouse."** `--db liana_mouseconsensus` uses the LIANA+ mouse consensus with MGI gene symbols; a control group or a panel built from a mouse atlas serves as the reference.

**"I do not want to depend on one score or one database."** FastCCC combines sixteen score definitions with the Cauchy test, and ships 22 databases, including CellPhoneDB, CellChatDB, NicheNet and the LIANA+ consensus, so the same data can be run on each with `--db`.

**"I want the analysis to run from one command."** Every command has `--json` output and writes a run manifest, so a result can be reproduced, scripted or handed to an agent.
