"""
Example: generate a FastCCC HTML report comparing two conditions.

This script runs FastCCC in parallel on:
  - The full combined dataset  (shown as the default tab in the report)
  - Condition A                (e.g. treated, disease subtype A, time point 1)
  - Condition B                (e.g. control, disease subtype B, time point 2)

Results are cached: if a result directory already contains a
*_significant_results.tsv file, FastCCC is skipped for that condition.

The final report contains four interactive tabs:
  All Cells | Condition A | Condition B | Condition A vs Condition B (Δ)

Usage
-----
    python report_multi_condition.py

Requirements
------------
    pip install fastccc   # includes matplotlib, seaborn, jinja2, networkx, adjusttext
"""

import os
import glob
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

# ── User-defined configuration ────────────────────────────────────────────────
H5AD_FILE = './data/my_dataset.h5ad'    # normalized log1p h5ad, or raw (see below)
DB_PATH   = './db/CPDBv5.0.0'           # LRI database directory

# Column in adata.obs that identifies the condition for subsetting
CONDITION_KEY = 'condition'

# Map short labels to the exact values in adata.obs[CONDITION_KEY]
CONDITIONS = {
    'Condition A': 'group_a',   # e.g. 'treated', 'PSC', 'Day7'
    'Condition B': 'group_b',   # e.g. 'control', 'PBC', 'Day0'
}

# Cell type column in adata.obs
CELLTYPE_KEY = 'cell_type'

RESULT_DIR_ALL = './results/all'
RESULT_DIR_A   = './results/condition_a'
RESULT_DIR_B   = './results/condition_b'
REPORT_DIR     = './report'
# ──────────────────────────────────────────────────────────────────────────────


def _find_task_id(result_dir: str):
    """Return the FastCCC task ID if results already exist, else None."""
    hits = glob.glob(f'{result_dir}/*_significant_results.tsv')
    return os.path.basename(hits[0]).replace('_significant_results.tsv', '') if hits else None


def _run_fastccc(label: str, result_dir: str) -> str:
    """
    Load the h5ad, optionally subset to one condition, and run FastCCC.
    Returns the task ID.
    """
    import scanpy as sc
    import fastccc

    adata = sc.read_h5ad(H5AD_FILE)

    # FastCCC requires normalized log1p-transformed counts.
    # Skip these two lines if your h5ad is already preprocessed.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    if label in CONDITIONS:
        adata = adata[adata.obs[CONDITION_KEY] == CONDITIONS[label]].copy()

    print(f'[{label}] {adata.n_obs:,} cells → FastCCC…', flush=True)
    t0 = time.perf_counter()
    fastccc.Cauchy_combination_of_statistical_analysis_methods(
        database_file_path = DB_PATH,
        celltype_file_path = None,
        counts_file_path   = adata,
        convert_type       = 'hgnc_symbol',
        meta_key           = CELLTYPE_KEY,
        save_path          = result_dir,
    )
    task_id = _find_task_id(result_dir)
    print(f'[{label}] done in {time.perf_counter() - t0:.1f}s  task_id={task_id}', flush=True)
    return task_id


if __name__ == '__main__':
    wall_start = time.perf_counter()

    runs = [
        ('All Cells',   RESULT_DIR_ALL),
        ('Condition A', RESULT_DIR_A),
        ('Condition B', RESULT_DIR_B),
    ]
    for _, rdir in runs:
        os.makedirs(rdir, exist_ok=True)

    # ── Step 1: Run FastCCC (skip cached conditions) ──────────────────────────
    pending = [(lbl, rdir) for lbl, rdir in runs if not _find_task_id(rdir)]
    cached  = [(lbl, rdir) for lbl, rdir in runs if _find_task_id(rdir)]

    for lbl, rdir in cached:
        print(f'[{lbl}] cached — task_id: {_find_task_id(rdir)}')

    if pending:
        print(f'Running {len(pending)} condition(s) in parallel…')
        t0 = time.perf_counter()
        with ProcessPoolExecutor(max_workers=len(pending)) as pool:
            futures = {pool.submit(_run_fastccc, lbl, rdir): lbl for lbl, rdir in pending}
            for fut in as_completed(futures):
                try:
                    fut.result()
                except Exception as e:
                    print(f'ERROR [{futures[fut]}]: {e}')
                    raise
        print(f'FastCCC wall time: {time.perf_counter() - t0:.1f}s')

    task_all = _find_task_id(RESULT_DIR_ALL)
    task_a   = _find_task_id(RESULT_DIR_A)
    task_b   = _find_task_id(RESULT_DIR_B)
    assert all([task_all, task_a, task_b]), 'Missing results for one or more conditions.'

    # ── Step 2: Generate the report ───────────────────────────────────────────
    from fastccc.report import generate_report

    name_a, name_b = list(CONDITIONS.keys())

    print('\nGenerating report…')
    t0 = time.perf_counter()
    report_path = generate_report(
        result_dir    = RESULT_DIR_ALL,
        task_id       = task_all,
        database_path = DB_PATH,
        output_dir    = REPORT_DIR,
        sample_name   = 'All Cells',
        gene_sets     = ['KEGG_2021_Human', 'GO_Biological_Process_2023'],
        cond_a_result_dir = RESULT_DIR_A,
        cond_a_task_id    = task_a,
        cond_a_name       = name_a,
        cond_b_result_dir = RESULT_DIR_B,
        cond_b_task_id    = task_b,
        cond_b_name       = name_b,
    )
    print(f'Report done in {time.perf_counter() - t0:.1f}s')
    print(f'Total wall time: {time.perf_counter() - wall_start:.1f}s')
    print(f'Report: {report_path}')
