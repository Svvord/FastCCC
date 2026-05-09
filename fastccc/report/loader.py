"""Load and integrate FastCCC outputs with database annotations."""

import glob
import os
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from loguru import logger


@dataclass
class CCCData:
    significant: pd.DataFrame       # significant LR pairs with pathway annotations
    pvals: pd.DataFrame             # full p-value matrix (celltype_pair × LRI)
    strength: pd.DataFrame          # full strength matrix
    counts_matrix: pd.DataFrame     # celltype × celltype count matrix
    strength_matrix: pd.DataFrame   # celltype × celltype mean strength matrix
    celltypes: list
    task_id: str
    sample_name: str
    database_path: str
    receptor_tf: Optional[pd.DataFrame] = None


def _find_file(result_dir: str, task_id: str, pattern: str, fallback_glob: str) -> str:
    primary = os.path.join(result_dir, f"{task_id}_{pattern}")
    if os.path.exists(primary):
        return primary
    matches = sorted(glob.glob(os.path.join(result_dir, f"{task_id}{fallback_glob}")))
    if matches:
        logger.warning(f"Using fallback file: {matches[0]}")
        return matches[0]
    raise FileNotFoundError(f"Cannot find {pattern} for task {task_id} in {result_dir}")


def load_results(
    result_dir: str,
    task_id: str,
    database_path: str,
    sample_name: str = "Sample",
    pval_threshold: float = 0.05,
) -> CCCData:

    # --- Load core files ---
    sig_path = _find_file(result_dir, task_id, "significant_results.tsv", "*significant_results.tsv")
    pval_path = _find_file(result_dir, task_id, "Cauchy_pvals.tsv", "*pvals.tsv")
    str_path  = _find_file(result_dir, task_id, "average_interactions_strength.tsv", "*interactions_strength.tsv")

    significant = pd.read_csv(sig_path, sep='\t')
    pvals       = pd.read_csv(pval_path, sep='\t', index_col=0)
    strength    = pd.read_csv(str_path,  sep='\t', index_col=0)

    # Align columns (strength may have more columns than pvals after Cauchy merge)
    shared_cols = pvals.columns.intersection(strength.columns)
    pvals    = pvals[shared_cols]
    strength = strength[shared_cols]

    # --- Annotate with pathway classification from DB ---
    itbl = pd.read_csv(os.path.join(database_path, 'interaction_table.csv'))
    itbl = itbl.set_index('id_cp_interaction')[['classification', 'directionality', 'is_ppi']]
    significant = significant.merge(itbl, left_on='LRI_ID', right_index=True, how='left')
    significant['classification'] = significant['classification'].fillna('Unknown')
    significant['directionality'] = significant['directionality'].fillna('Unknown')

    # --- Build celltype × celltype matrices ---
    celltypes = sorted(set(
        [idx.split('|')[0] for idx in pvals.index] +
        [idx.split('|')[1] for idx in pvals.index]
    ))

    sig_mask = (pvals < pval_threshold).astype(float)

    counts_matrix   = pd.DataFrame(0,   index=celltypes, columns=celltypes, dtype=int)
    strength_matrix = pd.DataFrame(0.0, index=celltypes, columns=celltypes)

    for row_idx in pvals.index:
        ct1, ct2 = row_idx.split('|')
        n_sig = int(sig_mask.loc[row_idx].sum())
        mean_s = float((strength.loc[row_idx] * sig_mask.loc[row_idx]).sum())
        counts_matrix.loc[ct1, ct2]   = n_sig
        strength_matrix.loc[ct1, ct2] = mean_s

    # --- Load receptor → TF table (CPDBv5 only) ---
    tf_path = os.path.join(database_path, 'receptor_to_transcription_factor.csv')
    receptor_tf = None
    if os.path.exists(tf_path):
        receptor_tf = pd.read_csv(tf_path)

    logger.info(
        f"Loaded {len(significant)} significant interactions across "
        f"{len(celltypes)} cell types for sample '{sample_name}'."
    )

    return CCCData(
        significant=significant,
        pvals=pvals,
        strength=strength,
        counts_matrix=counts_matrix,
        strength_matrix=strength_matrix,
        celltypes=celltypes,
        task_id=task_id,
        sample_name=sample_name,
        database_path=database_path,
        receptor_tf=receptor_tf,
    )
