import sys
import os, glob
import math
import getopt
import numpy as np
import pandas as pd
import pickle as pkl
from loguru import logger

def cauthy_combine(fastCCC_dir, task_id=None):
    if task_id is None:
        logger.warning("No task_id is provided, all pvals files will be combined.")
        pval_paths = glob.glob(fastCCC_dir+os.sep+'*pvals.csv')
    else:
        logger.info(f"Task ID for combining is :{task_id}")
        pval_paths = glob.glob(fastCCC_dir+os.sep+f'{task_id}*pvals.csv')
    logger.info(f"There are {len(pval_paths)} pval files.")
    joined_path = '\n'.join(pval_paths)
    logger.debug(f"\n{joined_path}")

    ct_pairs, cpis = None, None
    comb_dict = dict()
    for pval_path in pval_paths:
        comb = os.path.basename(pval_path)
        comb = comb.replace('pvals.csv', '')
        pval_df = pd.read_csv(pval_path, header=0, index_col=0)
        if ct_pairs is None:
            ct_pairs = pval_df.index.tolist()
        if cpis is None:
            cpis = pval_df.columns.tolist()
        comb_dict[comb] = pval_df.values

    pval_mat = []
    for comb, values in comb_dict.items():
        pval_mat.append(np.expand_dims(values, axis=1))
    pval_mat = np.concatenate(pval_mat, axis=1)
    weight = np.ones(len(comb_dict)) / len(comb_dict)
    T = pval_mat.copy()
    T[np.where(pval_mat == 1)] = np.tan(-np.pi/2)
    T[foo] = np.tan(np.pi*(0.5 - T[(foo:=np.where(pval_mat != 1))]))
    T =  weight @ T
    P = 0.5 - np.arctan(T) / np.pi

    T_df = pd.DataFrame(T, index=ct_pairs, columns=cpis)
    P_df = pd.DataFrame(P, index=ct_pairs, columns=cpis)

    if task_id is None:
        T_df.to_csv(fastCCC_dir+os.sep+'Cauchystats.csv')
        P_df.to_csv(fastCCC_dir+os.sep+'Cauchypvals.csv')
    else:
        T_df.to_csv(fastCCC_dir+os.sep+f'{task_id}_Cauchystats.csv')
        P_df.to_csv(fastCCC_dir+os.sep+f'{task_id}_Cauchypvals.csv')

# if __name__ == "__main__":
#     cauthy_combine(fastCCC_dir)