import copy
import numpy as np
import pandas as pd
from .distrib import Distribution, get_norm_minimum_cdf_func, pmf_bins

def combine_complex_distribution_df(set_distribution_df, complex_table, complex_process_function):
    def sub_func(x):
        return complex_process_function(*x)

    def func(x):
        x = [sub_x for sub_x in x if sub_x in set_distribution_df.columns]
        if len(x) == 0:
            # 如果出现这种情况，其实应该报错，前面数据处理应该排除这种情况。
            return pd.Series(index=mean_pmf.index)
        return set_distribution_df.loc[:,x].apply(sub_func, axis=1)

    if not complex_table.empty:
        complex_pmf = complex_table.groupby('complex_multidata_id').apply(lambda x: x['protein_multidata_id'].values, include_groups=False).apply(lambda x:func(x)).T
        complex_pmf = complex_pmf.dropna(axis=1)
        set_distribution_df = pd.concat((set_distribution_df, complex_pmf), axis=1)
        
    return set_distribution_df


def get_minimum_distribution(*pmf_list):
    '''
    distribution_list = k distribution 
    F = k * len(distribuion) matrix (pmf -> cdf use np.cumsum)
    minimum F(x) = 1 - \\product (1 - F_i(x))
    so f(x) = np.diff(F(x))
    '''
    for pmf in pmf_list:
        assert pmf.is_align, "PMF is not aligned"
    if len(pmf_list) == 1:
        return pmf_list[0]
    
    ### update 0203 ###
    loc_list, scale_list = [], []
    flag = True
    for pmf in pmf_list:
        if pmf.is_normal_type:
            loc_list.append(pmf.loc)
            scale_list.append(pmf.scale)
        else:
            flag = False
            break
    if flag:
        
        ### update 0204 ###
        arg_loc = np.argsort(loc_list)
        loc_min = loc_list[arg_loc[0]]
        loc_other = np.array(loc_list)[arg_loc[1:]]
        scale_min = scale_list[arg_loc[0]]
        scale_other = np.array(scale_list)[arg_loc[1:]]
        if loc_min + 5* scale_min <= np.min(loc_other - 5*scale_other):
            return Distribution(dtype='normal', loc=loc_min, scale=scale_min)
        
        ### ### ### ### ###
        
        analytic_func = get_norm_minimum_cdf_func(locs = loc_list, scales=scale_list)
        f = np.diff(analytic_func(pmf_bins), prepend=0)  # 没矫正和为1
        pmf = Distribution(dtype='other', pmf_array=f, is_align=True)
        pmf.is_analytic = True
        pmf.cdf_analytic_func = analytic_func
        return pmf
    ### ### ### ### ###
    
    distribution_list = [pmf.get_pmf_array() for pmf in pmf_list]
    max_length = max(arr.shape[0] for arr in distribution_list)
    padded_arrays = [np.pad(arr, (0, max_length - arr.shape[0]), 'constant') for arr in distribution_list]
    F = np.vstack(padded_arrays)
    # F = np.stack(distribution_list)
    F = np.cumsum(F, axis=1)
    minF = 1 - np.prod(1-F, axis=0)
    f = np.diff(minF, prepend=0)
    return Distribution(dtype='other', pmf_array=f, is_align=True)


def get_average_distribution(*distribution_list):

    assert len(distribution_list) > 0, "Distribution list is empty."
    for distribution in distribution_list:
        assert distribution.is_align, "PMF is not aligned"
    if len(distribution_list) == 1:
        return distribution_list[0]
    
    new_distribution = distribution_list[0]
    for i in range(1, len(distribution_list)):
        new_distribution += distribution_list[i]
    new_distribution /= len(distribution_list)
    return new_distribution