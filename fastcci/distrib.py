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

precision = 0.01 
max_value_4_log1p = 14
pmf_bin_edge = np.arange(0, max_value_4_log1p + 1e-10, precision)
pmf_bins = np.arange(0, max_value_4_log1p - 1e-10, precision)
pmf_int_bins = np.arange(len(pmf_bins))

mul_index_file_name = os.path.dirname(os.path.abspath(__file__))+os.sep+str(precision)+'index_mat.npy'
if os.path.exists(mul_index_file_name):
    print("Loading cached files.")
    index_mat = np.load(mul_index_file_name)
else:
    ############################################
    product_mat = np.outer(pmf_bins, pmf_bins)
    index_mat = np.zeros((len(pmf_bins), len(pmf_bins)), dtype=np.int16)
    last_value = -1
    for i, item in enumerate(pmf_bins):
        value = item ** 2
        index =np.where(np.logical_and(product_mat > last_value, product_mat <= value))
        index_mat[index] = i
        last_value = value
    np.save(mul_index_file_name, index_mat) 
    ############################################

def update_precision(new_precision):
    global precision, pmf_bin_edge, pmf_bins, pmf_int_bins
    precision = new_precision
    pmf_bin_edge = np.arange(0, max_value_4_log1p + 1e-10, precision)
    pmf_bins = np.arange(0, max_value_4_log1p - 1e-10, precision)
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


class Distribution:
    def __init__(
        self, 
        dtype: Literal['normal', 'gaussian', 'other'] = 'other', 
        loc: Union[float, None] = None, 
        scale: Union[float, None] = None, 
        samples: Union[np.ndarray, None] = None,
        pmf_array: Union[np.ndarray, None] = None, 
        is_align: Union[bool, None] = None,
        eps: float = 1e-7
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
            if np.abs(1 - np.sum(pmf_array)) >= eps:
                print(f"The sum of PMF is {np.sum(pmf_array)} not 1.")
            
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
            pdf = norm.pdf(pmf_bins, loc=self.loc, scale=self.scale)
            self.pmf_array = pdf / np.sum(pdf)
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
            return Distribution(dtype='normal', loc=new_loc, scale=new_scale)
        
        new_pmf_array = fftconvolve(self.get_pmf_array(), pmf2.get_pmf_array())
        new_pmf = Distribution(dtype='other', pmf_array=new_pmf_array, is_align=False)
        
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
        
        if n >= 30 or self.is_normal_type:
            new_loc = self.loc * n
            new_scale = np.sqrt(n) * self.scale
            return Distribution(dtype='normal', loc=new_loc, scale=new_scale)
        
        else:
            new_pmf_array = self.get_pmf_array()
            for i in range(1, n):
                new_pmf_array = fftconvolve(new_pmf_array, self.get_pmf_array())
                # new_pmf_array = np.trim_zeros(new_pmf_array, trim='b') 应该是不会进入这个优化，因为前面都是有非0，卷积仍然非0，除非太小？
            return Distribution(dtype='other', pmf_array=new_pmf_array, is_align=False)

            
    def __truediv__(self, n):
        assert type(n) == int
        assert n >= 1
        if n == 1:
            return self
        
        if self.is_normal_type:
            new_loc = self.loc / n
            new_scale = self.scale / n
            return Distribution(dtype='normal', loc=new_loc, scale=new_scale)

        cdf = self.get_cdf_array()
        bins = pmf_int_bins * n
        n_cut = np.sum(bins < len(cdf))
        average_pmf_array = np.diff(cdf[bins[:n_cut]], prepend=0)
        average_pmf =  Distribution(dtype='other', pmf_array=average_pmf_array, is_align=True)
        
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
        return Distribution(dtype='other', pmf_array=new_pmf_array, is_align=False)
    
    
    #######################################
    # === @marvinquiet: add more transformation of PMF === 
    '''locate_bin() and calculate_probs() can be shared between log & expotential transformation'''
    def _locate_bin(self, left_ptr, right_ptr, bin_intervals):
        '''Locate the corresponding log-transformed bins in pmf array'''
        interval_list = []
        for idx, (interval_left, interval_right) in enumerate(bin_intervals):
            # left pointer locate in the interval
            if left_ptr >= interval_left and left_ptr <= interval_right:
                interval_list.append(idx)
            # right pointer locate in the interval
            if right_ptr >= interval_left and right_ptr <= interval_right:
                interval_list.append(idx)
        interval_list = list(set(interval_list)) # in case of duplicated index
        interval_list.sort()
        return interval_list

    def _calculate_probs(self, left_ptr, right_ptr, interval_list, bin_intervals, prob_array, stepsize):
        res = 0.0
        if len(interval_list) == 0:
            res = 0.0
        elif len(interval_list) == 1: # both left/right locate in one interval
            prob = prob_array[interval_list[0]]
            bin_len = bin_intervals[interval_list[0]][1] - bin_intervals[interval_list[0]][0]
            # evenly distributed probablities
            res = prob/bin_len*stepsize
        else: # left locates in one interval and right locates in another interval
            # left prob remained in interval
            left_bin = bin_intervals[interval_list[0]]
            left_prob = prob_array[interval_list[0]]
            left_bin_len = bin_intervals[interval_list[0]][1] - bin_intervals[interval_list[0]][0]
            left_bin_len_remained = left_bin[1]-left_ptr # remained in left length
            left_remained_prob = left_prob/left_bin_len*left_bin_len_remained
            right_bin = bin_intervals[interval_list[-1]]
            right_prob = prob_array[interval_list[-1]]
            right_bin_len = bin_intervals[interval_list[-1]][1] - bin_intervals[interval_list[-1]][0]
            right_bin_len_remained = right_ptr-right_bin[0]
            right_remained_prob = right_prob/right_bin_len*right_bin_len_remained
            res = left_remained_prob+right_remained_prob
            if interval_list[-1]-interval_list[0] > 1:
                res += prob_array[(interval_list[0]+1):interval_list[-1]].sum()
        return res
    def __log__(self, log_bin_size=100):
        '''Perform log transformation on PMF. Coordinates and bin size adjusted.
        1. Set aside the first element as it involves -Inf
        2. Log-transform the rest of the PMF array
        3. Assign probabilities according to the length of the bin
        '''
        # ignore the first element as it involves -Inf
        pmf_bin_coords = np.log10(pmf_bins[1:])
        # new log-transformed coordinates
        log_pmf_bin_coords, log_pmf_bin_step = np.linspace(pmf_bin_coords[0], pmf_bin_coords[-1], num=log_bin_size, retstep=True) 
        # generate new array
        pmf_array = self.pmf_array[1:] 
        pmf_bin_coords = pmf_bin_coords[:len(pmf_array)+1]
        pmf_bin_intervals = [(pmf_bin_coords[_], pmf_bin_coords[_ + 1]) for _ in range(len(pmf_bin_coords) - 1)]
        if len(pmf_bin_intervals) < len(pmf_bins) - 1:
            remained_intervals = (pmf_bin_intervals[-1][1], np.log10(pmf_bins[-1]))
            pmf_bin_intervals.append(remained_intervals)
            pmf_array = np.append(pmf_array, 0.0)
        log_pmf_bins = np.zeros_like(log_pmf_bin_coords)
        log_pmf_bins[0] = self.pmf_array[0]
        #TODO: can be accelerated by two pointers algorithm
        for log_pmf_ptr in range(1, len(log_pmf_bins)):
            log_pmf_left, log_pmf_right = log_pmf_bin_coords[log_pmf_ptr-1], log_pmf_bin_coords[log_pmf_ptr]
            pmf_ptr_list = self._locate_bin(log_pmf_left, log_pmf_right, pmf_bin_intervals)
            log_pmf_bins[log_pmf_ptr] = self._calculate_probs(log_pmf_left, log_pmf_right, pmf_ptr_list, pmf_bin_intervals, pmf_array, log_pmf_bin_step)
        return Distribution(dtype='other', pmf_array=log_pmf_bins, is_align=False)

    def __logtruediv__(self, n, log_bin_size=100):
        assert type(n) == int
        assert n >= 1
        if n == 1:
            return self
        if self.is_normal_type:
            new_loc = self.loc / n
            new_scale = self.scale / n
            return Distribution(dtype='normal', loc=new_loc, scale=new_scale)
        cdf = self.get_cdf_array()
        log_pmf_int_bins = np.arange(log_bin_size)
        bins = log_pmf_int_bins * n
        n_cut = np.sum(bins < len(cdf))
        average_pmf_array = np.diff(cdf[bins[:n_cut]], prepend=0)
        average_pmf =  Distribution(dtype='other', pmf_array=average_pmf_array, is_align=True)
        # 为了精确 打得补丁
        if self.is_complex_analytic:
            average_pmf.is_complex_analytic = True
            average_pmf.ligand = self.ligand
            average_pmf.receptor = self.receptor
        ####
        return average_pmf

    def __exp__(self, log_bin_size=100):
        '''Perform exponential transform on log-transformed PMF. Coordinates and bin size adjusted.
        1. Set aside the first element as it involves -Inf
        2. Expotential the rest of log-transformed PMF
        3. Assign probabilities according to the length of the bin
        '''
        # ignore the first element as it involves -Inf
        log_pmf_bin_coords, log_pmf_bin_step = np.linspace(np.log10(pmf_bins[1]), np.log10(pmf_bins[-1]), num=log_bin_size, retstep=True)
        log_pmf_bin_coords = np.power(10, log_pmf_bin_coords) # transform back to original scale
        log_pmf_array = self.pmf_array[1:]
        log_pmf_bin_coords = log_pmf_bin_coords[:len(log_pmf_array)+1]
        log_pmf_bin_intervals = [(log_pmf_bin_coords[_], log_pmf_bin_coords[_ + 1]) for _ in range(len(log_pmf_bin_coords) - 1)]
        if len(log_pmf_bin_intervals) < log_bin_size-1:
            remained_intervals = (log_pmf_bin_intervals[-1][1], pmf_bins[-1])
            log_pmf_bin_intervals.append(remained_intervals)
            log_pmf_array = np.append(log_pmf_array, 0.0)
        # new expotential-transformed coordinates
        pmf_bin_coords = pmf_bins[1:]
        new_pmf_bins = np.zeros_like(pmf_bin_coords)
        new_pmf_bins[0] = self.pmf_array[0]
        for pmf_ptr in range(1, len(new_pmf_bins)):
            pmf_left, pmf_right = pmf_bin_coords[pmf_ptr-1], pmf_bin_coords[pmf_ptr]
            log_pmf_ptr_list = self._locate_bin(pmf_left, pmf_right, log_pmf_bin_intervals)
            new_pmf_bins[pmf_ptr] = self._calculate_probs(pmf_left, pmf_right, log_pmf_ptr_list, log_pmf_bin_intervals, log_pmf_array, precision)
        return Distribution(dtype='other', pmf_array=new_pmf_bins, is_align=False)

    def __wenjing_mul__(self, pmf2):
        '''log-transform and the exponential back
        '''
        # if self.is_normal_type and pmf2.is_normal_type:
        #     # TODO: return results with formula?
        #     return self.__add__(pmf2)
        log_p1 = self.__log__()
        log_p2 = pmf2.__log__()
        log_p1p2 = log_p1 & log_p2 # log x + log y
        avg_pmf = log_p1p2.__logtruediv__(n=2) # log (xy) / 2
        return avg_pmf.__exp__() # exp(log(xy)/2) = sqrt(xy)

    # === @marvinquiet: end ===
    
    
    def __mul__(self, distribution2):
        '''
        Output is the distribution of sqrt(XY) not XY
        '''
        foo = np.outer(self.pmf, distribution2.pmf)
        h, w = foo.shape
        local_index_mat = index_mat[:h, :w]
        data_flat = foo.ravel()
        labels_flat = local_index_mat.ravel()
        pmf_array = np.bincount(labels_flat, weights=data_flat)
        return Distribution(dtype='other', pmf_array=pmf_array, is_align=True)
    
    #######################################
    
    
        


def get_pmf_array_from_samples(samples, is_log1p=True):
    assert np.min(samples) >= 0, "Expression value should be non-negative"
    if is_log1p:
        assert np.max(samples) < max_value_4_log1p, "Support domain is not valid."
        pmf_array, _ = np.histogram(samples, pmf_bin_edge)
        pmf_array = pmf_array / np.sum(pmf_array)
        return pmf_array
    else:
        # count data
        return get_pmf_array_from_samples(np.log1p(samples))
    
def get_distribution_from_samples(samples, is_log1p=True):
    pmf_array = get_pmf_array_from_samples(samples, is_log1p)
    loc = np.mean(samples)
    scale = np.std(samples)
    return Distribution('other', pmf_array=pmf_array, loc=loc, scale=scale, is_align=True)



def get_pvalue_from_pmf(value, pmf):
    y = pmf.get_pmf_array()
    if pmf.is_analytic:
        pvalue =  1 - pmf.cdf_analytic_func(value)
    elif pmf.is_complex_analytic:
        pvalue = get_pvalue_from_complex_pmf(value, pmf)
    else:
        pvalue = 1 - np.sum(y[:int(np.ceil(value / precision))])
    return pvalue


def get_precise_pmf_array(pmf, n_segment=100):
    '''
    只在支持域上划分n_segment个小段，然后算cdf准确数值，然后算差值作为精准pmf
    支持域从真实pmf的第一个非0起到最后一个非0，对应cdf的第一个非0和第一个1
    '''
    assert pmf.is_analytic    
    assert pmf.min_cdf_non_zero < pmf.min_cdf_one
    support_domain = np.linspace(pmf.min_cdf_non_zero, pmf.min_cdf_one, n_segment+1)
    precise_cdf_array = pmf.cdf_analytic_func(support_domain)
    precise_pmf_array = np.diff(precise_cdf_array, prepend=0)
    return precise_pmf_array, support_domain
    

def get_pvalue_from_complex_pmf(value, pmf):
    # assert pmf.is_from_complex, "This pmf is not calculated from complex proteins."
    assert pmf.is_complex_analytic, "This pmf is not \"complex_analytic\"."
    assert pmf.ligand.is_analytic
    assert pmf.receptor.is_analytic
    
    value *= 2 # we calculate the sum of ligand and receptor instead
    pmf.min_cdf_non_zero = (pmf.ligand.min_cdf_non_zero + pmf.receptor.min_cdf_non_zero) / 2
    pmf.min_cdf_one = (pmf.ligand.min_cdf_one + pmf.receptor.min_cdf_one) / 2
    print(pmf.min_cdf_non_zero, pmf.min_cdf_one)
    if value <= pmf.min_cdf_non_zero * 2:
        return 1.0
    elif value >= pmf.min_cdf_one * 2:
        return 0.0
    else:
        p_value = 0.0
        if pmf.ligand.support_length > pmf.receptor.support_length:
            accurate_pmf = pmf.ligand
            discrete_pmf = pmf.receptor
        else:
            accurate_pmf = pmf.receptor
            discrete_pmf = pmf.ligand

        precise_pmf_array, support_domain = get_precise_pmf_array(discrete_pmf)
        print(precise_pmf_array[49:54])
        for x, p in zip(support_domain, precise_pmf_array):
            p_value += p * (1 - accurate_pmf.cdf_analytic_func(value - x))
        return p_value
    
    
# def get_quantile_pmf_for_n_iid_distribution(distribution, n, quantile=0.9):
#     '''
#     $
#     f_{X_i} = \frac{n!}{(n-i)!(i-1)!} 
#               * [F(x_i)]^(i-1)
#               * f(x_i)
#               * [1 - F(x_i)]^(n-i)
#     $
#     '''
#     assert 0 < quantile and quantile < 1, "Quantile should be 0-1."
#     fx = distribution.pmf
#     Fx = distribution.cdf.copy()
#     Qx = 1 - Fx
#     # 0606 update
#     # Qx = np.insert(Qx, 0, 1.0)[:-1]
#     Fx = np.insert(Fx, 0, np.PZERO)[:-1]
    
    
#     i = max(int(quantile*n), 1)
#     with np.errstate(divide='ignore', invalid='ignore'):
#         logfx = np.log(fx)
#         logFx = (i-1) * np.log(Fx) if i > 1 else 0.0
#         logQx = (n-i) * np.log(Qx) if i < n else 0.0
#         logp = logfx + logFx + logQx
#         logp = logp - np.max(logp)
#         fxi = np.exp(logp)
#     fxi /= np.sum(fxi)
#     if np.abs(1 - np.sum(fxi)) >= 1e-7:
#         print(np.sum(fxi))
#     return fxi

#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
m = 100000
with np.errstate(divide='ignore', invalid='ignore'):
    logFx = np.log(np.arange(1, m+1))
    logQx = np.log(m+1-np.arange(1, m+1))
    m -= 1
    
@lru_cache(maxsize=1)
def get_fxi_cdf(n, i):
    with np.errstate(divide='ignore', invalid='ignore'):
        logp = (i-1) * logFx + (n-i) * logQx
        logp = logp - np.max(logp)
        fxi = np.exp(logp)
        fxi /= np.sum(fxi)
    return np.cumsum(fxi)

            
def get_quantile_pmf_for_n_iid_distribution(distribution, n, quantile=0.9):
    '''
    $
    f_{X_i} = \frac{n!}{(n-i)!(i-1)!} 
              * [F(x_i)]^(i-1)
              * f(x_i)
              * [1 - F(x_i)]^(n-i)
    $
    '''
    assert 0 < quantile and quantile < 1, "Quantile should be 0-1."
    i = max(int(quantile*(n-1) + 1), 1)
    fxi_cdf = get_fxi_cdf(n, i)

    section = np.round((m+1) * distribution.pmf)
    section[0] = m + 1 - np.sum(section[1:])
    section = np.int32(np.cumsum(section))
    cdf = fxi_cdf[section-1]
    pmf = np.diff(cdf, prepend=0)
    pmf = np.clip(pmf, a_min=0, a_max=None)
    return pmf
#=================================================================
