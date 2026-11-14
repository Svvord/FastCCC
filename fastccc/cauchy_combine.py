import sys
import os, glob
import math
import getopt
import numpy as np
import pandas as pd
import pickle as pkl
from loguru import logger

SMALL_P_METHODS = ('acat', 'legacy')


def combine_pvalues(pval_mat, weights=None, small_p_method='acat'):
    """Cauchy combination of module p-values.

    ``pval_mat`` has the modules on axis 1 (shape ``(pairs, modules, ...)``). Returns the
    combined statistic and p-value with axis 1 removed. P-values are first clamped to
    [0, 1]; NaN entries propagate.

    The statistic averages ``tan(pi * (0.5 - p))`` over the modules, whose extremes the
    formula does not define (p = 0 gives +inf, p = 1 gives -inf). ``small_p_method``
    selects how they are handled:

    - ``'acat'`` (default): a p-value of exactly 1 (a summary of 0, e.g. the median of a
      gene expressed in fewer than half of the cells) is truncated to ``1 - 1/d`` so that
      it counts as the weakest possible evidence rather than as an infinite veto; p-values
      below 1e-15 enter as ``1 / (p * pi)``, the tail approximation of the ACAT
      implementation, so an exact 0 gives a combined p-value of 0.
      With a single module (``d = 1``) nothing is truncated and the combined p-value is
      that module's p-value.
    - ``'legacy'``: the rule of the published method. P-values at or below 1e-16 are
      treated as 0 and enter through ``tan``, which caps their term at ``tan(pi / 2)`` in
      floating point; p = 1 enters as ``tan(-pi / 2)``.
    """
    if small_p_method not in SMALL_P_METHODS:
        raise ValueError(f"small_p_method must be one of {SMALL_P_METHODS}, got {small_p_method!r}")
    P = np.clip(np.asarray(pval_mat, dtype=float), 0.0, 1.0)
    d = P.shape[1]
    w = np.ones(d) / d if weights is None else np.asarray(weights, dtype=float) / np.sum(weights)
    T = np.empty_like(P)
    if small_p_method == 'acat':
        if d > 1:   # a single module has nothing to combine with: its p-value is returned as is (1 - 1/1 would turn p = 1 into 0)
            P = np.where(P >= 1.0, 1.0 - 1.0 / d, P)
        small = P < 1e-15
        T[~small] = np.tan(np.pi * (0.5 - P[~small]))
        with np.errstate(divide='ignore'):
            T[small] = 1.0 / (P[small] * np.pi)
    else:
        Ps = np.where(P <= 1e-16, 0.0, P)
        T[...] = np.tan(np.pi * (0.5 - Ps))
    S = np.tensordot(w, np.moveaxis(T, 1, 0), axes=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        # the upper tail of the standard Cauchy; 1 / (pi * S) keeps precision for very large statistics
        Pc = np.where(S > 1e15, 1.0 / (np.pi * S), 0.5 - np.arctan(S) / np.pi)
    return S, Pc


def cauthy_combine(fastCCC_dir, task_id=None, small_p_method='acat'):
    if task_id is None:
        logger.warning("No task_id is provided, all pvals files will be combined.")
        pval_paths = glob.glob(fastCCC_dir+os.sep+'*pvals.tsv')
    else:
        logger.info(f"Task ID for combining is :{task_id}")
        pval_paths = glob.glob(fastCCC_dir+os.sep+f'{task_id}*pvals.tsv')
    logger.info(f"There are {len(pval_paths)} pval files.")
    joined_path = '\n'.join(pval_paths)
    logger.debug(f"\n{joined_path}")

    ct_pairs, cpis = None, None
    comb_dict = dict()
    for pval_path in pval_paths:
        comb = os.path.basename(pval_path)
        comb = comb.replace('pvals.tsv', '')
        pval_df = pd.read_csv(pval_path, header=0, index_col=0, sep='\t')
        if ct_pairs is None:
            ct_pairs = pval_df.index.tolist()
        if cpis is None:
            cpis = pval_df.columns.tolist()
        comb_dict[comb] = pval_df.values

    pval_mat = []
    for comb, values in comb_dict.items():
        pval_mat.append(np.expand_dims(values, axis=1))
    pval_mat = np.concatenate(pval_mat, axis=1)
    T, P = combine_pvalues(pval_mat, small_p_method=small_p_method)

    T_df = pd.DataFrame(T, index=ct_pairs, columns=cpis)
    P_df = pd.DataFrame(P, index=ct_pairs, columns=cpis)

    if task_id is None:
        T_df.to_csv(fastCCC_dir+os.sep+'Cauchy_stats.tsv', sep='\t')
        P_df.to_csv(fastCCC_dir+os.sep+'Cauchy_pvals.tsv', sep='\t')
    else:
        T_df.to_csv(fastCCC_dir+os.sep+f'{task_id}_Cauchy_stats.tsv', sep='\t')
        P_df.to_csv(fastCCC_dir+os.sep+f'{task_id}_Cauchy_pvals.tsv', sep='\t')

# if __name__ == "__main__":
#     cauthy_combine(fastCCC_dir)