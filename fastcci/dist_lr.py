import pandas as pd
import numpy as np
import itertools
import timeit
from .distrib import get_pvalue_from_pmf

def calculate_key_interactions_pvalue(mean_pmf, interactions, interactions_strength, percent_analysis, method='Arithmetic'):
    p1_index = []
    p2_index = []
    all_index = []
    for i in itertools.product(sorted(mean_pmf.index), sorted(mean_pmf.index)):
        p1_index.append(i[0])
        p2_index.append(i[1])
        all_index.append('|'.join(i))
        
    p1 = mean_pmf.loc[p1_index, interactions['multidata_1_id']]
    p2 = mean_pmf.loc[p2_index, interactions['multidata_2_id']]
    p1.columns = interactions.index
    p2.columns = interactions.index
    p1.index = all_index
    p2.index = all_index
    
    p1_items = p1.values[np.where(percent_analysis)]
    p2_items = p2.values[np.where(percent_analysis)]
    
    
    start = timeit.default_timer()
    if method == 'Arithmetic':
        pval_pmfs = (p1_items & p2_items) / 2
    elif method == 'Geometric':
        pval_pmfs = p1_items * p2_items
    stop = timeit.default_timer()
    
    mean_gt = interactions_strength.values[np.where(percent_analysis)]
    est = []
    for i, value in enumerate(mean_gt):
        pval_est = get_pvalue_from_pmf(value, pval_pmfs[i])
        est.append(pval_est)
        
    pvalues = np.ones_like(interactions_strength)
    pvalues[np.where(percent_analysis)] = est
    pvalues = pd.DataFrame(pvalues, index=interactions_strength.index, columns=interactions_strength.columns)
    return pvalues