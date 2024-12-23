import scanpy as sc
import numpy as np
import pandas as pd
from .build_reference import rank_preprocess, get_fastcci_input
from .build_reference import calculate_L_R_and_IS_score, calculate_L_R_and_IS_percents
from .build_reference import calculate_mean_pmfs, record_hk_genes
from . import score
from .core import calculate_cluster_percents
from loguru import logger
import pickle
import itertools
from scipy.stats import norm
import os

precision = 0.01

def get_norm_minimum_ppf_func(locs, scales):
    def analytic_func(x, locs=locs, scales=scales):
        F = []
        for loc, scale in zip(locs, scales):
            F.append(norm.cdf(x, loc, scale))
        minF = 1 - np.prod(1-np.array(F), axis=0)
        return minF
    def inverse_f(y, low=0, high=50, tol=1e-6): # high=n_bins
        assert y > tol and y < 1 - tol
        while high - low > tol:
            mid = (low + high) / 2
            if analytic_func(mid) < y:
                low = mid
            else:
                high = mid
        return (low + high) / 2
    return inverse_f

def get_norm_ppf_func(loc, scale):
    def analytic_func(x):
        return norm.ppf(x, loc, scale)
    return analytic_func
    
def get_threshold_from_pmf(pmf, alpha=0.025):
    '''
    alpha: 概率分布右侧面积
    '''
    if pmf.is_analytic:
        threshold = get_norm_ppf_func(pmf.loc, pmf.scale)(1-alpha)
    elif pmf.is_complex_analytic:
        raise RuntimeError
    else:
        y = pmf.get_cdf_array()
        threshold = np.argmax((1 - y) <= alpha) * precision
    return threshold

def compare_strategy():
    pass

def get_valid_pval_pmfs(mean_pmfs, interactions_strength, interactions, percents_analysis):
    p1_index = []
    p2_index = []
    all_index = []
    for i in interactions_strength.index:
        p1_index.append(i.split('|')[0])
        p2_index.append(i.split('|')[1])
        all_index.append(i)
        
    p1 = mean_pmfs.loc[p1_index, interactions['multidata_1_id']]
    p2 = mean_pmfs.loc[p2_index, interactions['multidata_2_id']]
    p1.columns = interactions.index
    p2.columns = interactions.index
    p1.index = all_index
    p2.index = all_index

    p1_items = p1.values[np.where(percents_analysis)]
    p2_items = p2.values[np.where(percents_analysis)]

    pval_pmfs = (p1_items & p2_items) / 2

    return pval_pmfs

def compare_with_reference(counts_df, labels_df, complex_table, interactions, reference_path, save_path, k=2.59, debug_mode=False):
    assert os.path.exists(reference_path), "Reference dir doesn't exist."
    logger.info("Loading reference data.")
    with open(f'{reference_path}/complex_table.pkl', 'rb') as f:
        ref_complex_table = pickle.load(f)
    with open(f'{reference_path}/interactions.pkl', 'rb') as f:
        ref_interactions = pickle.load(f)
    with open(f'{reference_path}/ref_gene_pmf_dict.pkl', 'rb') as f:
        ref_gene_pmf_dict = pickle.load(f)
    with open(f'{reference_path}/ref_min_percentile.txt', 'rt') as f:
        min_percentile = float(f.readline().strip())
        logger.info(f"Ref min_percentile is loaded as {min_percentile}")
    ref_pvals = pd.read_csv(f'{reference_path}/ref_pvals.txt', sep='\t', index_col=0)
    ref_p1 = pd.read_csv(f'{reference_path}/ref_interactions_strength_L.txt', sep='\t', index_col=0)
    ref_p2 = pd.read_csv(f'{reference_path}/ref_interactions_strength_R.txt', sep='\t', index_col=0)
    ref_percents_analysis = pd.read_csv(f'{reference_path}/ref_percents_analysis.txt', sep='\t', index_col=0)
    ref_L_perc = pd.read_csv(f'{reference_path}/ref_percents_L.txt', sep='\t', index_col=0)
    ref_R_perc = pd.read_csv(f'{reference_path}/ref_percents_R.txt', sep='\t', index_col=0)
    logger.success("Reference data is loaded.")

    logger.info("Calculating IS score for query data.")
    mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
    null_mean_counts = mean_counts.copy()
    null_mean_counts.values[:,:] = np.repeat(np.array([counts_df.mean(axis=0)]), len(mean_counts), axis=0).reshape(*mean_counts.shape)
    complex_func = score.calculate_complex_min_func
    mean_counts = score.combine_complex_distribution_df(mean_counts, complex_table, complex_func)
    null_mean_counts = score.combine_complex_distribution_df(null_mean_counts, complex_table, complex_func)
    null_interactions_strength = score.calculate_interactions_strength(null_mean_counts, interactions, method='Arithmetic')
    p1, p2, interactions_strength = calculate_L_R_and_IS_score(mean_counts, interactions)
    percents = calculate_cluster_percents(counts_df, labels_df, complex_table)

    logger.info("Filtering reference data.")
    common_ind = sorted(set(ref_pvals.index) & set(interactions_strength.index))
    common_col = sorted(set(ref_pvals.columns) & set(interactions_strength.columns))
    ref_pvals = ref_pvals.loc[common_ind, common_col]
    ref_percents_analysis= ref_percents_analysis.loc[common_ind, common_col]
    ref_p1= ref_p1.loc[common_ind, common_col]
    ref_p2= ref_p2.loc[common_ind, common_col]
    ref_L_perc = ref_L_perc.loc[common_ind, common_col]
    ref_R_perc = ref_R_perc.loc[common_ind, common_col]

    logger.info("Filtering by using reference.")
    # filtering use ref
    complex_table = complex_table.loc[[item for item in complex_table.index if item in ref_complex_table.index]]
    interactions = interactions.loc[[item for item in interactions.index if item in ref_interactions.index]]
    L_perc, R_perc, percents_analysis = calculate_L_R_and_IS_percents(percents, interactions, threshold=min_percentile)
    interactions_strength = interactions_strength.loc[percents_analysis.index, percents_analysis.columns]
    null_interactions_strength = null_interactions_strength.loc[percents_analysis.index, percents_analysis.columns]
    p1 = p1.loc[percents_analysis.index, percents_analysis.columns]
    p2 = p2.loc[percents_analysis.index, percents_analysis.columns]

    p1_list = p1.values[np.where(percents_analysis)]
    p2_list = p2.values[np.where(percents_analysis)]

    logger.info("Inferring sig. boundaries.")
    mean_pmfs = calculate_mean_pmfs(counts_df, labels_df, complex_table, ref_gene_pmf_dict)
    pval_pmfs = get_valid_pval_pmfs(mean_pmfs, interactions_strength, interactions, percents_analysis)

    up_bound_list = []
    low_bound_list = []
    for i in range(len(pval_pmfs)):
        threshold = get_threshold_from_pmf(pval_pmfs[i], 0.05)
        up_bound = threshold + threshold / k * 1.96 
        low_bound = threshold - threshold / k * 1.96 
        up_bound_list.append(up_bound)
        low_bound_list.append(low_bound)

    valid_IS_list = interactions_strength.values[np.where(percents_analysis)]
    null_IS_list = null_interactions_strength.values[np.where(percents_analysis)]
    L_perc_list = L_perc.values[np.where(percents_analysis)]
    R_perc_list = R_perc.values[np.where(percents_analysis)]

    # strategy
    prediction = []
    est_by_ref = []
    results = []
    for i, (index, col) in enumerate(zip(*np.where(percents_analysis))):
        index = percents_analysis.index[index]
        col = percents_analysis.columns[col]
        IS = interactions_strength.loc[index, col]
        assert valid_IS_list[i] == IS
        low_IS = low_bound_list[i]
        up_IS = up_bound_list[i]
        null_IS = null_IS_list[i]

        ligand_IS = p1_list[i]
        receptor_IS = p2_list[i]
        ligand_perc = L_perc_list[i]
        receptor_perc = R_perc_list[i]

        ref_flag = index in ref_pvals.index
        is_in_ref = True if ref_flag else False
        ref_perc = ref_percents_analysis.loc[index, col] if ref_flag else np.nan
        ref_ligand_perc = ref_L_perc.loc[index, col] if ref_flag else np.nan
        ref_receptor_perc = ref_R_perc.loc[index, col] if ref_flag else np.nan
        ligand_low = np.max(ref_p1.loc[index, col] * (1-1/k*1.96), 0) if ref_flag else np.nan
        ligand_high = ref_p1.loc[index, col] * (1+1/k*1.96) if ref_flag else np.nan
        receptor_low = np.max(ref_p2.loc[index, col], 0) * (1-1/k*1.96) if ref_flag else np.nan
        receptor_high = ref_p2.loc[index, col] * (1+1/k*1.96) if ref_flag else np.nan
        ref_sig = (ref_pvals.loc[index, col] < 0.05) if ref_flag else np.nan # 参考是否显著
        ligand_range = f"{ligand_low}-{ligand_high}" if ref_flag else np.nan
        receptor_range = f"{receptor_low}-{receptor_high}" if ref_flag else np.nan


        if IS > up_IS:
            prediction.append(1)
            est_by_ref.append(True)
        elif IS < low_IS:
            prediction.append(0)
            est_by_ref.append(True)
        else:
            if ref_flag:
                if ref_sig:
                    if ligand_IS >= ligand_low and receptor_IS >= receptor_low:
                        prediction.append(1)
                        est_by_ref.append(True)
                    else:
                        prediction.append(1 if IS > null_IS else 0)
                        est_by_ref.append(False)
                else:
                    if ligand_IS <= ligand_high and receptor_IS <= receptor_high:
                        prediction.append(0)
                        est_by_ref.append(True)
                    else:
                        prediction.append(1 if IS > null_IS else 0)
                        est_by_ref.append(False)
            else:
                prediction.append(1 if IS > null_IS else 0)
                est_by_ref.append(False)

        if ref_flag:
            if prediction[-1] and ref_sig:
                comparison = "Both Sig"
            elif prediction[-1] and not ref_sig:
                comparison = "Up"
            elif not prediction[-1] and ref_sig:
                comparison = "Down"
            else:
                comparison = "Both NS"
        else:
            comparison = np.nan

        results.append((
            index, is_in_ref, col, f"{IS}", null_IS, f"{low_IS}-{up_IS}", 
            True, ligand_perc, receptor_perc, 
            ref_perc, ref_ligand_perc, ref_receptor_perc, 
            ligand_IS, ligand_range, 
            receptor_IS, receptor_range,
            bool(prediction[-1]), ref_sig, comparison
        ))
    
    for index, col in zip(*np.where(np.logical_and(ref_pvals < 0.05, percents_analysis.loc[common_ind, common_col] == False))):
        index = ref_pvals.index[index]
        col = ref_pvals.columns[col]
        IS = interactions_strength.loc[index, col]
        null_IS = null_interactions_strength.loc[index, col]
        ligand_perc = L_perc.loc[index, col]
        receptor_perc = R_perc.loc[index, col]
        ref_ligand_perc = ref_L_perc.loc[index, col]
        ref_receptor_perc = ref_R_perc.loc[index, col]
        ligand_IS = p1.loc[index, col]
        receptor_IS = p2.loc[index, col]

        results.append((
            index, True, col, f"{IS}", null_IS, np.nan, 
            False, ligand_perc, receptor_perc, 
            True, ref_ligand_perc, ref_receptor_perc, 
            ligand_IS, np.nan, 
            receptor_IS, np.nan,
            False, True, "Down"
        ))


    results_df = pd.DataFrame(
        results, 
        columns=[
            'sender|receiver', 'is_in_ref', 'LRI_id', 'IS_score', 'null_IS',
            'Sig_thres_interval_by_ref',
            '>min_percentile', 'L_perc', 'R_perc',
            '>min_percentile_ref', 'L_perc_ref', 'R_perc_ref',
            'IS_L_part','IS_L_interval_by_ref',
            'IS_R_part', 'IS_R_interval_by_ref',
            'is_sig', 'is_sig_ref', 'comparison'
        ]
    )

    if debug_mode:
        logger.debug("Entering debug process")
        real_pvals = pd.read_csv(f'{save_path}/debug_pvals.txt', sep='\t', index_col=0)
        real_pvals = real_pvals.loc[interactions_strength.index, interactions_strength.columns]
        real_pvals_values = real_pvals.values[np.where(percents_analysis)]
        intersection = np.sum(np.logical_and(prediction, real_pvals_values < 0.05))
        union = np.sum(np.logical_or(prediction, real_pvals_values < 0.05))
        print(f"Intersection:{intersection}, Union:{union}")
        print(f"#Pred.Sig:{np.sum(prediction)}, #Real.Sig:{np.sum(real_pvals_values < 0.05)}, #valid:{len(real_pvals_values)}")
        print(f"Precision:{intersection/np.sum(prediction)}, IoU:{intersection/union}, Recall:{intersection/np.sum(real_pvals_values < 0.05)}")

        est_by_ref = np.array(est_by_ref)
        print(f"#by_null: {np.sum(1 - est_by_ref)}, #by_ref: {np.sum(est_by_ref)}")
        real_pvals_values = real_pvals_values[est_by_ref]
        prediction = np.array(prediction)[est_by_ref]
        intersection = np.sum(np.logical_and(prediction, real_pvals_values < 0.05))
        union = np.sum(np.logical_or(prediction, real_pvals_values < 0.05))
        print(f"Real Prec:{intersection/np.sum(prediction)}, Real IoU:{intersection/union}, Real Recall:{intersection/np.sum(real_pvals_values < 0.05)}")
        logger.success("Debug ends.")


    ####### save reference results #######
    logger.info("Saving inference results.")
    percents_analysis.to_csv(f'{save_path}/query_percents_analysis.txt', sep='\t')
    interactions_strength.to_csv(f'{save_path}/query_interactions_strength.txt', sep='\t')
    results_df.to_csv(f'{save_path}/query_infer_results.txt', sep='\t')


def calculate_adjust_factor(query, reference_path, save_path, debug_mode=False):
    mean_hk_rnk, gene_index = record_hk_genes(query)
    query_hk = pd.DataFrame(np.array(mean_hk_rnk).flatten(), index=gene_index, columns=['query_hk'])
    ref_hk = pd.read_csv(f'{reference_path}/ref_hk.txt', sep='\t', index_col=0)
    hk_df = pd.merge(query_hk, ref_hk, left_index=True, right_index=True)
    k = np.sqrt(1 / np.mean(np.square((hk_df.ref_hk - hk_df.query_hk) / hk_df.ref_hk)))

    if debug_mode:
        import matplotlib.pyplot as plt 
        plt.figure(figsize=(6,3))
        plt.subplot(1,2,1)
        plt.scatter(hk_df.query_hk, hk_df.ref_hk, s=3, alpha=0.1)
        plt.subplot(1,2,2)
        plt.scatter(np.log(hk_df.query_hk), np.log(hk_df.ref_hk), s=3, alpha=0.1)
        plt.savefig(f'{save_path}/hk_scatter.jpg')

        mean_list = []
        std_list = []
        mean_by_std_list = []
        for i in np.arange(0,30):
            foo = np.where(np.logical_and(hk_df.ref_hk >=i, hk_df.ref_hk < i+5))
            mean_ = np.mean(np.array(hk_df.ref_hk)[foo])
            std_ = np.sqrt(np.mean(np.square(np.array(hk_df.query_hk)[foo] - np.array(hk_df.ref_hk)[foo])))
            mean_list.append(mean_)
            std_list.append(std_)
            mean_by_std_list.append(mean_/std_)
        plt.figure(figsize=(3,3))
        plt.scatter(mean_list, std_list)
        plt.title(np.nanmean(mean_by_std_list))
        print(np.nanmean(mean_by_std_list))
        plt.savefig(f'{save_path}/hk_mean_by_std.jpg')

        return np.nanmean(mean_by_std_list)

    return k
    



def infer_query_workflow(database_file_path, reference_path, query_counts_file_path, celltype_file_path, save_path, meta_key=None, debug_mode=False):
    logger.info(f"Start inferring by using CCI reference: {reference_path}")

    if not os.path.exists(save_path):
        os.makedirs(save_path)
        logger.success(f"Query save dir is created.")
    else:
        logger.warning(f"{save_path} already exists, all files will be overwritten")


    query = sc.read_h5ad(query_counts_file_path)
    sc.pp.filter_cells(query, min_genes=50) # basic QC
    logger.info(f"Reading query adata, {query.shape[0]} cells x {query.shape[1]} genes")
    if debug_mode:
        query.write_h5ad(f"{save_path}/debug_digit.h5ad")

    if meta_key is not None:
        labels_df = pd.DataFrame(query.obs[meta_key])
        labels_df.columns = ['cell_type']
        labels_df.index.name = 'barcode_sample'
    else:
        labels_df = pd.read_csv(meta_file_path, sep='\t', index_col=0)
        for barcode in query.obs_names:
            assert barcode in labels_df.index, "The index of query data doesn't match the index of labels"
        labels_df = labels_df.loc[query.obs_names, :]

    query = rank_preprocess(query)
    k = calculate_adjust_factor(query, reference_path, save_path, debug_mode)
    logger.info(f"k={k}")
    if k < 2.6:
        k = 2.6
    counts_df, complex_table, interactions = get_fastcci_input(query, database_file_path)
    if debug_mode:
        logger.debug("Entering debug process")
        from .build_reference import fastcci_for_reference
        fastcci_for_reference('', save_path, counts_df, labels_df, complex_table, interactions, min_percentile = 0.1, debug_mode=True)
        logger.debug("Debug ends.")
    compare_with_reference(counts_df, labels_df, complex_table, interactions, reference_path, save_path, k=k, debug_mode=debug_mode)
    logger.success("Inference workflow done.")