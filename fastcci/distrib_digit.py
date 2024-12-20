import os
import numpy as np
import pandas as pd
from scipy.special import loggamma, gamma
from scipy.optimize import minimize

from scipy.signal import fftconvolve
from scipy.stats import norm

from typing import Literal, Union
from functools import lru_cache

from tqdm import tqdm
from loguru import logger

n_bins = 50
precision = 0.01
n_fft = 30
pmf_bins = np.arange(0, n_bins+precision - 1e-10, precision)
pmf_int_bins = np.arange(len(pmf_bins))


def get_norm_cdf_func(loc, scale):
    def analytic_func(x):
        return norm.cdf(x, loc, scale)
    return analytic_func

def get_norm_minimum_cdf_func(locs, scales):
    def analytic_func(x, locs=locs, scales=scales):
        F = []
        for loc, scale in zip(locs, scales):
            F.append(norm.cdf(x, loc, scale))
        minF = 1 - np.product(1-np.array(F), axis=0)
        return minF
    return analytic_func


class Distribution_digit:
    def __init__(
        self, 
        dtype: Literal['normal', 'gaussian', 'other'] = 'other', 
        loc: Union[float, None] = None, 
        scale: Union[float, None] = None, 
        samples: Union[np.ndarray, None] = None,
        pmf_array: Union[np.ndarray, None] = None, 
        is_align: Union[bool, None] = None,
        eps: float = 1e-4
    ):

        if dtype in ['normal', 'gaussian']:
            assert loc is not None and scale is not None
            is_align = True
            dtype = 'normal'
        else:
            if samples is not None:
                pmf_array = get_pmf_array_from_samples(samples)
                is_align = True
            assert pmf_array is not None and is_align is not None
            pmf_array = np.trim_zeros(pmf_array, trim='b')
            assert np.abs(1 - np.sum(pmf_array)) < eps, f"The sum of PMF is {np.sum(pmf_array)} not 1."
            # if np.abs(1 - np.sum(pmf_array)) >= 1e-7:
            #     print(f"The sum of PMF is {np.sum(pmf_array)} not 1.")
            
        self.loc = loc
        self.scale = scale
        self.dtype = dtype
        self.is_align = is_align
        self.pmf_array = pmf_array
        self.is_normal_type = self.dtype == 'normal'
        self.is_analytic = self.is_normal_type # normal or combine of normal
        self.is_complex_analytic = False
        if self.is_normal_type:
            self.cdf_analytic_func = get_norm_cdf_func(loc, scale)
            self.min_cdf_non_zero = self.loc - 5 * self.scale
            self.min_cdf_one = self.loc + 5 * self.scale
            self.support_length = self.min_cdf_one - self.min_cdf_non_zero
            
        # supplementary property
        self.cdf_array = None
    
    @property
    def mean(self):
        return self.get_mean()
    
    @property
    def std(self):
        return self.get_std()
        
    @property    
    def var(self):
        if self.is_normal_type:
            return np.square(self.scale)
        return np.sum(np.square(np.arange(len(self.pmf_array)) * precision - self.mean) 
                      * self.pmf_array )
    @property
    def pmf(self):
        return self.get_pmf_array()
    
    @property
    def cdf(self):
        return self.get_cdf_array()
        
    def get_mean(self):
        if self.loc is None:
            self.loc = np.sum(self.pmf_array * np.arange(len(self.pmf_array))) * precision
        return self.loc
    
    def get_std(self):
        if self.scale is None:
            self.scale = np.sqrt(self.var)
        return self.scale
            
    def get_pmf_array(self):
        if self.pmf_array is None: # which is a normal distribution
            # >>>>>>>>>>> update precision >>>>>>>>>> #
            # old version:
            # pdf = norm.pdf(pmf_bins, loc=self.loc, scale=self.scale)
            # self.pmf_array = pdf / np.sum(pdf)
            # new version: 
            cdf = norm.cdf(pmf_bins, loc=self.loc, scale=self.scale)
            self.pmf_array = np.diff(cdf, prepend=0)
            # <<<<<<<<<<< update precision <<<<<<<<<< #
            self.pmf_array = np.trim_zeros(self.pmf_array, trim='b')
        return self.pmf_array
        
    def get_cdf_array(self):
        if self.cdf_array is None:
            self.cdf_array = np.cumsum(self.get_pmf_array())
        return self.cdf_array
        
    def __add__(self, pmf2):
        if self.is_normal_type and pmf2.is_normal_type:
            new_loc = self.loc + pmf2.loc
            new_scale = np.sqrt(np.square(self.scale) + np.square(pmf2.scale))
            return Distribution_digit(dtype='normal', loc=new_loc, scale=new_scale)
        
        new_pmf_array = fftconvolve(self.get_pmf_array(), pmf2.get_pmf_array())
        new_pmf = Distribution_digit(dtype='other', pmf_array=new_pmf_array, is_align=False)
        
        # 为了精确 打得补丁
        if self.is_analytic and pmf2.is_analytic:
            self.is_complex_analytic = True
            new_pmf.ligand = self
            new_pmf.receptor = pmf2
        ###
 
        return new_pmf
        
    def __pow__(self, n):
        assert type(n) == int
        assert n >= 1
        if n == 1:
            return self
        
        if n >= n_fft or self.is_normal_type:
            new_loc = self.loc * n
            new_scale = np.sqrt(n) * self.scale
            return Distribution_digit(dtype='normal', loc=new_loc, scale=new_scale)
        
        else:
            new_pmf_array = self.get_pmf_array()
            for i in range(1, n):
                new_pmf_array = fftconvolve(new_pmf_array, self.get_pmf_array())
                # new_pmf_array = np.trim_zeros(new_pmf_array, trim='b') 应该是不会进入这个优化，因为前面都是有非0，卷积仍然非0，除非太小？
            return Distribution_digit(dtype='other', pmf_array=new_pmf_array, is_align=False)

    def __truediv__(self, n):
        assert type(n) == int
        assert n >= 1
        if n == 1:
            return self
        
        if self.is_normal_type:
            new_loc = self.loc / n
            new_scale = self.scale / n
            return Distribution_digit(dtype='normal', loc=new_loc, scale=new_scale)

        cdf = self.get_cdf_array()
        bins = pmf_int_bins * n
        n_cut = np.sum(bins < len(cdf))
        # 增加 append 1.0，修正可能的bug
        average_pmf_array = np.diff(np.append(cdf[bins[:n_cut]], 1.0), prepend=0)
        average_pmf =  Distribution_digit(dtype='other', pmf_array=average_pmf_array, is_align=True)
        
        # 为了精确 打得补丁
        if self.is_complex_analytic:
            average_pmf.is_complex_analytic = True
            average_pmf.ligand = self.ligand
            average_pmf.receptor = self.receptor
        ####
            
        return average_pmf

    def __and__(self, pmf2):
        if self.is_normal_type and pmf2.is_normal_type:
            return self.__add__(pmf2)
        p1 = self.get_pmf_array()
        p2 = pmf2.get_pmf_array()
        pr0 = p1[0] + p2[0] - p1[0]*p2[0]
        new_pmf_array = fftconvolve(p1[1:], p2[1:])
        new_pmf_array = np.insert(new_pmf_array, 0, [pr0, 0])
        return Distribution_digit(dtype='other', pmf_array=new_pmf_array, is_align=False)



def get_pmf_array_from_samples_for_digitized_bins(samples, n_bins=50):
    assert type(n_bins) == int and n_bins > 0, "n_bins should be positive integer"
    assert np.max(samples) <= n_bins, "error: n_bins < max(samples)"
    pmf_array, _ = np.histogram(samples, np.arange(n_bins+2))
    pmf_array = pmf_array / np.sum(pmf_array)
    return pmf_array


def get_minimum_distribution_for_digit(*pmf_list):
    '''
    distribution_list = k distribution 
    F = k * len(distribuion) matrix (pmf -> cdf use np.cumsum)
    minimum F(x) = 1 - \product (1 - F_i(x))
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
            return Distribution_digit(dtype='normal', loc=loc_min, scale=scale_min)
        
        ### ### ### ### ###
        
        analytic_func = get_norm_minimum_cdf_func(locs = loc_list, scales=scale_list)
        f = np.diff(analytic_func(pmf_bins), prepend=0)  # 没矫正和为1
        pmf = Distribution_digit(dtype='other', pmf_array=f, is_align=True)
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
    minF = 1 - np.product(1-F, axis=0)
    f = np.diff(minF, prepend=0)
    return Distribution_digit(dtype='other', pmf_array=f, is_align=True)