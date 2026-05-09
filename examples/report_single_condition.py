"""
Example: generate a FastCCC HTML report for a single dataset.

Steps
-----
1. Load and preprocess your AnnData (normalize → log1p).
2. Run FastCCC to get significant interactions.
3. Call generate_report() to produce a self-contained HTML report.

Usage
-----
    python report_single_condition.py
"""

import os
import glob
import scanpy as sc
import fastccc
from fastccc.report import generate_report

# ── User-defined paths ────────────────────────────────────────────────────────
H5AD_FILE   = './data/my_dataset.h5ad'   # normalized log1p h5ad, or raw (see below)
DB_PATH     = './db/CPDBv5.0.0'          # LRI database directory
RESULT_DIR  = './results/my_sample'
REPORT_DIR  = './report/my_sample'
SAMPLE_NAME = 'My Sample'

# Cell type column in adata.obs
CELLTYPE_KEY = 'cell_type'
# ──────────────────────────────────────────────────────────────────────────────

os.makedirs(RESULT_DIR, exist_ok=True)

# ── Step 1: Load and preprocess ───────────────────────────────────────────────
adata = sc.read_h5ad(H5AD_FILE)

# FastCCC requires normalized log1p-transformed counts.
# Skip these two lines if your h5ad is already preprocessed.
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ── Step 2: Run FastCCC ───────────────────────────────────────────────────────
fastccc.Cauchy_combination_of_statistical_analysis_methods(
    database_file_path = DB_PATH,
    celltype_file_path = None,
    counts_file_path   = adata,
    convert_type       = 'hgnc_symbol',
    meta_key           = CELLTYPE_KEY,
    save_path          = RESULT_DIR,
)

# ── Step 3: Find the task ID and generate the report ─────────────────────────
hits = glob.glob(f'{RESULT_DIR}/*_significant_results.tsv')
task_id = os.path.basename(hits[0]).replace('_significant_results.tsv', '')

report_path = generate_report(
    result_dir    = RESULT_DIR,
    task_id       = task_id,
    database_path = DB_PATH,
    output_dir    = REPORT_DIR,
    sample_name   = SAMPLE_NAME,
    gene_sets     = ['KEGG_2021_Human', 'GO_Biological_Process_2023'],
)

print(f'Report saved to: {report_path}')
