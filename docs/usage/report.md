---
title: Automated Report
layout: default
parent: Basic usage
nav_order: 5
---

# Automated HTML Report

After running FastCCC, a single function call generates a self-contained HTML report you can open in any browser. The report includes static figures plus interactive cell-type explorers organized into thematic sections. If you have two groups of cells (e.g., treated vs. control), the report adds a fourth condition-comparison tab automatically.

---

## Step 1 — Run FastCCC and get the task ID

Run FastCCC as you normally would. FastCCC saves results as files named `<task_id>_significant_results.tsv` — the task ID is the 6-character prefix in the filename.

```python
import glob, os
import scanpy as sc
import fastccc

# Load and preprocess your data
adata = sc.read_h5ad('./data/my_dataset.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Run FastCCC
fastccc.Cauchy_combination_of_statistical_analysis_methods(
    database_file_path = './db/CPDBv5.0.0',
    celltype_file_path = None,
    counts_file_path   = adata,
    convert_type       = 'hgnc_symbol',
    meta_key           = 'cell_type',   # column in adata.obs with cell type labels
    save_path          = './results/',
)

# Find the task ID from the output filename automatically
hits = glob.glob('./results/*_significant_results.tsv')
task_id = os.path.basename(hits[0]).replace('_significant_results.tsv', '')
print('task_id:', task_id)   # e.g. ab12cd
```

---

## Step 2 — Generate the report

```python
from fastccc.report import generate_report

report_path = generate_report(
    result_dir    = './results/',
    task_id       = task_id,            # from Step 1
    database_path = './db/CPDBv5.0.0',
    output_dir    = './report/',
    sample_name   = 'My Sample',        # display name shown in the report
)

print('Report saved to:', report_path)  # open this file in your browser
```

The output `report.html` is fully self-contained — all figures are embedded in the file, so you can share it as a single file with no dependencies.

---

## Two-condition comparison (optional)

If you want to compare two groups, run FastCCC separately on each group and pass all three result directories to `generate_report`. FastCCC will add a condition-comparison tab to the report.

```python
import glob, os
import scanpy as sc
import fastccc
from fastccc.report import generate_report

adata = sc.read_h5ad('./data/my_dataset.h5ad')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ── Run FastCCC on all three subsets ──────────────────────────────────────────
def run_and_get_task_id(adata_subset, save_path):
    os.makedirs(save_path, exist_ok=True)
    fastccc.Cauchy_combination_of_statistical_analysis_methods(
        database_file_path = './db/CPDBv5.0.0',
        celltype_file_path = None,
        counts_file_path   = adata_subset,
        convert_type       = 'hgnc_symbol',
        meta_key           = 'cell_type',
        save_path          = save_path,
    )
    hits = glob.glob(f'{save_path}/*_significant_results.tsv')
    return os.path.basename(hits[0]).replace('_significant_results.tsv', '')

# Full dataset
task_all = run_and_get_task_id(adata, './results/all/')

# Group A (e.g. treated)
adata_a = adata[adata.obs['condition'] == 'treated'].copy()
task_a   = run_and_get_task_id(adata_a, './results/treated/')

# Group B (e.g. control)
adata_b = adata[adata.obs['condition'] == 'control'].copy()
task_b   = run_and_get_task_id(adata_b, './results/control/')

# ── Generate report with condition-comparison tab ────────────────────────────
report_path = generate_report(
    result_dir    = './results/all/',
    task_id       = task_all,
    database_path = './db/CPDBv5.0.0',
    output_dir    = './report/',
    sample_name   = 'All Cells',

    cond_a_result_dir = './results/treated/',
    cond_a_task_id    = task_a,
    cond_a_name       = 'Treated',

    cond_b_result_dir = './results/control/',
    cond_b_task_id    = task_b,
    cond_b_name       = 'Control',
)
print('Report saved to:', report_path)
```

The report will have four tabs: **All Cells**, **Treated**, **Control**, and **Treated vs Control**.

The two-condition page is descriptive. It compares independently analysed FastCCC
result sets using significant interaction counts, communication-score changes, and
within-condition FastCCC evidence. Replicate-aware differential claims require a
between-condition statistical model matched to the study design.

---

## What the report contains

Each sample tab contains seven sections. Several global summaries are paired with interactive cell-type explorers so users can inspect the same evidence within an individual sender, receiver, or either-role context.

| Section | Figures / panels | What it shows |
|---------|------------------|---------------|
| Global CCC Overview | 1–4 | Interaction count heatmap, communication strength heatmap, chord diagram, top sender/receiver bar chart |
| L-R Pair Analysis | 5, 5B, 6, 6B, 7 | Top L-R dot plot, cell-type L-R evidence explorer, pathway classification distribution, cell-type pathway evidence explorer, pathway × cell-type heatmap |
| Pathway Enrichment | 8, 8B, 9, 9B, 10 | Ligand ORA, receptor ORA, cell-type ligand/receptor gene evidence explorers, TF regulon heatmap |
| Cell-Type Profiles | 11–13, 25–26 | I/O scatter, sender/receiver pathway specificity heatmaps, outgoing/incoming interaction flow |
| Network Analyses | 14–17, 27–28 | Network centrality, Sankey flow, autocrine/paracrine, bipartite L-R graph, asymmetry and community summaries |
| Advanced L-R | 18, 18B, 19–21, 29–30 | Pathway information flow, cell-type pathway CS explorer, per-pathway L-R multiples, L-R specificity heatmap, CS violin, L-R co-occurrence, pathway crosstalk |

The condition-comparison tab (when two conditions are provided) adds three more figures:

| Figure | What it shows |
|--------|---------------|
| 22 | Significant interaction count comparison (Condition B minus Condition A) |
| 23 | L-R comparison evidence plot with CS fold-change and within-condition evidence |
| 24 | Pathway interaction count comparison |

All sample tabs also include a figure-generation audit. If a panel cannot be
generated because a dependency, annotation field, network-backed ORA request, or
required value is unavailable, the report records the skipped/failed panel and
keeps the rest of the report usable.

---

## All parameters

```python
generate_report(
    # ── Required ──────────────────────────────────────────────────────────────
    result_dir    = './results/',      # directory containing FastCCC output files
    task_id       = 'ab12cd',          # 6-character prefix from the output filename
    database_path = './db/CPDBv5.0.0', # LRI database directory

    # ── Output ────────────────────────────────────────────────────────────────
    output_dir    = './report/',       # where to save report.html (default: results/report_<task_id>/)
    sample_name   = 'My Sample',       # display name shown in the report header

    # ── Analysis settings ─────────────────────────────────────────────────────
    pval_threshold          = 0.05,    # significance cutoff
    top_n_lr                = 30,      # number of L-R pairs shown in dot plot
    top_n_celltypes         = 20,      # number of cell types shown in bar charts
    max_chords              = 180,     # maximum directed connections shown in chord plots
    gene_sets               = ['KEGG_2021_Human', 'GO_Biological_Process_2023'],
    dpi                     = 300,     # figure resolution
    save_individual_figures = True,    # also save each figure as a separate PNG

    # ── Condition A (optional) ────────────────────────────────────────────────
    cond_a_result_dir = None,
    cond_a_task_id    = None,
    cond_a_name       = None,

    # ── Condition B (optional, requires Condition A) ──────────────────────────
    cond_b_result_dir = None,
    cond_b_task_id    = None,
    cond_b_name       = None,
)
```

---

## Choosing enrichment gene sets

The default gene sets (`KEGG_2021_Human` and `GO_Biological_Process_2023`) work for most human datasets. To see all available options:

```python
import gseapy
print(gseapy.get_library_name(organism='Human'))
```

Common choices:

| Gene set | Description |
|----------|-------------|
| `KEGG_2021_Human` | KEGG metabolic and signalling pathways |
| `GO_Biological_Process_2023` | Gene Ontology biological process terms |
| `Reactome_2022` | Reactome pathway database |
| `WikiPathway_2023_Human` | WikiPathways |
