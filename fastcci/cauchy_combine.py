import sys
import os, glob
import math
import getopt
import numpy as np
import pandas as pd
import pickle as pkl

usage='Usage:'+sys.argv[0]
usage+='''<Required>[Options]
    <Required>
    -d --fastCCI_dir directory that contains fastCCI results
'''

if len(sys.argv) < 2 or not sys.argv[1].startswith('-'):
    sys.exit(usage)


optlist,alist=getopt.getopt(sys.argv[1:],
                            'hd:',
                            [
                                'help=',
                                'fastCCI_dir='
                            ])
for opt, arg in optlist:
    if opt in ['-h', '--help']:
        sys.exit(usage)
    elif opt in ['-d', '--fastCCI_dir']:
        fastCCI_dir = arg


# TODO: for debug purpose
#fastCCI_dir = "/net/mulan/home/wenjinma/projects/FastCCI-benchmark/simdata/PDAC/simulated_dataset1/tools/fastCCI/output"

def cauthy_combine(fastCCI_dir):
    pval_paths = glob.glob(fastCCI_dir+os.sep+'*pvals.csv')

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

    T_df.to_csv(fastCCI_dir+os.sep+'Cauchystats.csv')
    P_df.to_csv(fastCCI_dir+os.sep+'Cauchypvals.csv')

if __name__ == "__main__":
    cauthy_combine(fastCCI_dir)