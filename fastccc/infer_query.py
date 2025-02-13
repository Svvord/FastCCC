import scanpy as sc
import numpy as np
import pandas as pd
from .build_reference import rank_preprocess, get_fastccc_input
from .build_reference import calculate_L_R_and_IS_score, calculate_L_R_and_IS_percents
from .build_reference import calculate_mean_pmfs, record_hk_genes
from . import score
from .core import calculate_cluster_percents
from loguru import logger
import pickle
import itertools
from scipy.stats import norm
import os
import tomllib
from collections import Counter
from .distrib_digit import get_minimum_distribution_for_digit
from . import dist_complex
from . import dist_lr
import json

# logger.remove()  # 移除默认的日志处理器
# logger.add(
#     sink=lambda msg: print(msg.strip()),
#     format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | <cyan>{message}</cyan>",
#     colorize=True
# )



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

def get_ligand_receptor_summary(mean_counts, interactions):
    p1_index = []
    p2_index = []
    all_index = []
    for i in itertools.product(sorted(mean_counts.index), sorted(mean_counts.index)):
        p1_index.append(i[0])
        p2_index.append(i[1])
        all_index.append('|'.join(i))
    p1 = mean_counts.loc[p1_index, interactions['multidata_1_id']]
    p2 = mean_counts.loc[p2_index, interactions['multidata_2_id']]
    p1.columns = interactions.index
    p2.columns = interactions.index
    p1.index = all_index
    p2.index = all_index
    interactions_strength = (p1 + p2)/2 * (p1 > 0) * (p2>0)
    return p1, p2, interactions_strength

def get_celltype_mapping_dict(celltype_mapping_dict):
    match celltype_mapping_dict:
        case None:
            logger.info("Reference cell types label will be used directly.")
            return None
        case str():
            logger.info("Parsing celltype_mapping_dict file.")
            try:
                with open(celltype_mapping_dict, 'r', encoding='utf-8') as file:
                    celltype_mapping_dict = json.load(file)
                return celltype_mapping_dict
            except FileNotFoundError:
                logger.error(f"The file '{celltype_mapping_dict}' was not found.")
                logger.info("Reference cell types label will be used directly.")
                return None
            except json.JSONDecodeError:
                logger.error(f"The file '{celltype_mapping_dict}' is not a valid JSON file.")
                logger.info("Reference cell types label will be used directly.")
                return None
            return "The variable is a string."
        case dict():
            logger.info("The variable celltype_mapping_dict is used.")
            return celltype_mapping_dict
        case _:
            logger.error(f"The variable celltype_mapping_dict is of an unknown type.")
            logger.info("Reference cell types label will be used directly.")
            return None

def aggregate_by_weight(df, merge_dict, weight_dict):
    results = []
    index = []
    for key in sorted(merge_dict):
        index.append(key)
        row = []
        sum_ = 0
        if len(merge_dict[key]) == 1:
            results.append(df.loc[merge_dict[key][0]].values)
        else:
            for ct in merge_dict[key]:
                assert weight_dict[ct] > 0
                sum_ += weight_dict[ct]
                row.append(df.loc[ct].values * weight_dict[ct])
            results.append(np.array(row).sum(axis=0) / sum_)
    return pd.DataFrame(results, index=index, columns=df.columns)


def load_LRI_id_to_gene_symbol_dict(LRI_db_path):
    from . import preprocess
    interactions = preprocess.get_interactions(LRI_db_path)
    id2symbol_dict = (pd.read_csv(f'{LRI_db_path}/gene_table.csv')[['protein_id', 'hgnc_symbol']]\
    .set_index('protein_id')['hgnc_symbol']).to_dict()
    ##### complex_table ######
    complex_composition = pd.read_csv(os.path.join(LRI_db_path, 'complex_composition_table.csv'))
    complex_table = pd.read_csv(os.path.join(LRI_db_path, 'complex_table.csv'))
    complex_table = complex_table.merge(complex_composition, left_on='complex_multidata_id', right_on='complex_multidata_id')
    foo_dict = complex_table.groupby('complex_multidata_id').apply(lambda x: list(x['protein_multidata_id'].values)).to_dict()
    for key in foo_dict:
        value = ','.join([id2symbol_dict[item] for item in foo_dict[key]])
        id2symbol_dict[key] = value
    return id2symbol_dict


def get_null_ref_mean_counts(ref_gene_pmf_dict, ref_complex_table):
    # 使用 reference 数据，在没有任何条件下，gene 的均值应该是多少，complex 也同样被考虑。
    # 返回一个 dataframe
    null_ref_mean_counts = []
    columns = []
    for multi_id in ref_gene_pmf_dict:
        columns.append(multi_id)
        null_ref_mean_counts.append(ref_gene_pmf_dict[multi_id][1].mean)
    null_ref_mean_counts = pd.DataFrame([null_ref_mean_counts], index=['ref_null'], columns=columns)
    complex_func = score.calculate_complex_min_func
    null_ref_mean_counts = score.combine_complex_distribution_df(null_ref_mean_counts, ref_complex_table, complex_func)
    return null_ref_mean_counts


def compare_with_reference(counts_df, labels_df, complex_table, interactions, reference_path, save_path, config, celltype_mapping_dict, database_file_path, k=2.59, debug_mode=False):
    assert os.path.exists(reference_path), "Reference dir doesn't exist."

    logger.info("Loading reference data.")
    with open(f'{reference_path}/complex_table.pkl', 'rb') as f:
        ref_complex_table = pickle.load(f)
    with open(f'{reference_path}/interactions.pkl', 'rb') as f:
        ref_interactions = pickle.load(f)
    with open(f'{reference_path}/ref_gene_pmf_dict.pkl', 'rb') as f:
        ref_gene_pmf_dict = pickle.load(f)

    with open(f'{reference_path}/ref_percents.pkl', 'rb') as f:
        ref_percents = pickle.load(f)
    with open(f'{reference_path}/ref_mean_counts.pkl', 'rb') as f:
        ref_mean_counts = pickle.load(f)

    gene_list = [gene for gene in ref_mean_counts.columns if gene in counts_df.columns]
    ref_complex_table = ref_complex_table.loc[[item for item in ref_complex_table.index if item in complex_table.index]]
    ref_interactions = ref_interactions.loc[[item for item in ref_interactions.index if item in interactions.index]]

    null_ref_mean_counts = get_null_ref_mean_counts(ref_gene_pmf_dict, ref_complex_table)

    ref_label_counter = config['celltype']
    label_counter = Counter(labels_df['cell_type'])

    # 25-01-08 add
    # Combine cell types
    celltype_mapping_dict = get_celltype_mapping_dict(celltype_mapping_dict)
    if celltype_mapping_dict is None:

        ref_meta_dict = {}
        for item in ref_label_counter:
            if item not in label_counter:
                continue
            ref_meta_dict[item] = ref_label_counter[item]
        logger.debug(f"Valid reference cell type data:\n{str(ref_meta_dict)}")

    else:
        for key, value in celltype_mapping_dict.items():
            assert key in ref_label_counter, 'Please check the spelling, capitalization, spacing, and format of “your cell type”, as well as whether it is included in the reference.'
            assert value in label_counter, 'Please check the spelling, capitalization, spacing, and format of “your cell type”, as well as whether it is included in the query.'
        
        query_to_reference_dict = {}
        for key, value in celltype_mapping_dict.items():
            if value not in query_to_reference_dict:
                query_to_reference_dict[value] = [key]
            else:
                query_to_reference_dict[value].append(key)

        for item in ref_label_counter:
            if item not in celltype_mapping_dict and item in label_counter:
                query_to_reference_dict[item] = [item]

        ref_mean_counts = aggregate_by_weight(ref_mean_counts, query_to_reference_dict, ref_label_counter)
        ref_percents = aggregate_by_weight(ref_percents, query_to_reference_dict, ref_label_counter)

        ref_meta_dict = {}
        for key, values in query_to_reference_dict.items():
            sum_ = 0
            for item in values:
                sum_ += ref_label_counter[item]
            ref_meta_dict[key] = sum_
        logger.debug(f"Valid reference cell type data:\n{str(ref_meta_dict)}")

    unique_meta_dict = {}
    for item in label_counter:
        if item not in ref_meta_dict:
            unique_meta_dict[item] = label_counter[item]
    logger.debug(f"Unique query cell type data:\n{str(unique_meta_dict)}")
    

    # 25-01-08 end


    
    # ref_meta_dict = {}
    # for item in ref_label_counter:
    #     if item not in label_counter:
    #         continue
    #     ref_meta_dict[item] = ref_label_counter[item]

    ref_mean_counts = ref_mean_counts.loc[[item for item in ref_mean_counts.index if item in ref_meta_dict]]
    ref_percents = ref_percents.loc[[item for item in ref_mean_counts.index if item in ref_meta_dict]]

    ref_p1, ref_p2, ref_CS = get_ligand_receptor_summary(ref_mean_counts, ref_interactions)
    ref_L_perc, ref_R_perc, ref_percents_analysis = calculate_L_R_and_IS_percents(ref_percents, ref_interactions, threshold=config['min_percentile'])

    ref_clusters_mean_dict = {}
    for celltype in sorted(ref_meta_dict):
        ref_clusters_mean_dict[celltype]  = {}
        n_sum = ref_meta_dict[celltype]
        if n_sum < 100:
            for gene in gene_list:
                ref_clusters_mean_dict[celltype][gene] = ref_gene_pmf_dict[gene][n_sum]
        else:
            for gene in gene_list:
                ref_clusters_mean_dict[celltype][gene] = ref_gene_pmf_dict[gene][1] ** n_sum / n_sum
    ref_mean_pmfs = pd.DataFrame(ref_clusters_mean_dict).T
    complex_func = get_minimum_distribution_for_digit
    if len(ref_mean_pmfs):
        ref_mean_pmfs = dist_complex.combine_complex_distribution_df(ref_mean_pmfs, ref_complex_table, complex_func)
    logger.success("Reference data is loaded.")

    logger.info("Calculating CS score for query data.")
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
    common_ind = sorted(set(ref_percents_analysis.index) & set(interactions_strength.index))
    common_col = sorted(set(ref_percents_analysis.columns) & set(interactions_strength.columns))
    if len(ref_mean_pmfs):
        ref_pvals = dist_lr.calculate_key_interactions_pvalue(
            ref_mean_pmfs, ref_interactions, ref_CS, ref_percents_analysis, method='Arithmetic'
        )
        ref_pvals = ref_pvals.loc[common_ind, common_col]
    else:
        ref_pvals = pd.DataFrame([], index=common_ind, columns=common_col)
    # real_ref_pvals = pd.read_csv(f'{reference_path}/ref_pvals.txt', sep='\t', index_col=0)
    # real_ref_pvals = real_ref_pvals.loc[common_ind, common_col]
    # print(f"Error: {np.sum(np.abs(real_ref_pvals.values - ref_pvals.values))}")
    # print(f"Equal: {np.array_equal(real_ref_pvals.values, ref_pvals.values)}")
    # tmp = np.where(real_ref_pvals.values!= ref_pvals.values)
    # print('real')
    # print(real_ref_pvals.values[tmp])
    # print('ref')
    # print(ref_pvals.values[tmp])
    # print(f'ref_pvals shape: {ref_pvals.shape}')
    ref_percents_analysis= ref_percents_analysis.loc[common_ind, common_col]
    ref_p1= ref_p1.loc[common_ind, common_col]
    ref_p2= ref_p2.loc[common_ind, common_col]
    ref_L_perc = ref_L_perc.loc[common_ind, common_col]
    ref_R_perc = ref_R_perc.loc[common_ind, common_col]

    

    logger.info("Filtering by using reference.")
    # filtering use ref
    complex_table = complex_table.loc[[item for item in complex_table.index if item in ref_complex_table.index]]
    interactions = interactions.loc[[item for item in interactions.index if item in ref_interactions.index]]
    L_perc, R_perc, percents_analysis = calculate_L_R_and_IS_percents(percents, interactions, threshold=config['min_percentile'])
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

    id2symbol_dict = load_LRI_id_to_gene_symbol_dict(database_file_path)

    # strategy
    prediction = []
    est_by_ref = []
    results = []
    for i, (index, col) in enumerate(zip(*np.where(percents_analysis))):
        index = percents_analysis.index[index]
        col = percents_analysis.columns[col]

        ligand_gene_name = interactions.loc[col].multidata_1_id
        receptor_gene_name = interactions.loc[col].multidata_2_id
        ligand_gene_name = id2symbol_dict[ligand_gene_name]
        receptor_gene_name =  id2symbol_dict[receptor_gene_name]

        #01-11 add
        if debug_mode == True:
            mid1 = interactions.loc[col].multidata_1_id
            mid2 = interactions.loc[col].multidata_2_id
            # L_null_ref = ref_mean_pmfs.iloc[0][mid1].mean
            # R_null_ref = ref_mean_pmfs.iloc[0][mid2].mean
            L_null_ref = null_ref_mean_counts[mid1].item()
            R_null_ref = null_ref_mean_counts[mid2].item()

        #01-11 end

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

        if debug_mode:
            results.append((
                index, is_in_ref, col,  ligand_gene_name, receptor_gene_name,
                f"{IS}", null_IS, f"{low_IS}-{up_IS}", L_null_ref, R_null_ref,
                True, ligand_perc, receptor_perc, 
                ref_perc, ref_ligand_perc, ref_receptor_perc, 
                ligand_IS, ligand_range, 
                receptor_IS, receptor_range,
                bool(prediction[-1]), ref_sig, comparison
            ))
        else:   
            results.append((
                index, is_in_ref, col,  ligand_gene_name, receptor_gene_name,
                f"{IS}", null_IS, f"{low_IS}-{up_IS}", 
                True, ligand_perc, receptor_perc, 
                ref_perc, ref_ligand_perc, ref_receptor_perc, 
                ligand_IS, ligand_range, 
                receptor_IS, receptor_range,
                bool(prediction[-1]), ref_sig, comparison
            ))
    
    for index, col in zip(*np.where(np.logical_and(ref_pvals < 0.05, percents_analysis.loc[common_ind, common_col] == False))):
        index = ref_pvals.index[index]
        col = ref_pvals.columns[col]

        #01-11 add
        if debug_mode == True:
            mid1 = interactions.loc[col].multidata_1_id
            mid2 = interactions.loc[col].multidata_2_id
            # L_null_ref = ref_mean_pmfs.iloc[0][mid1].mean
            # R_null_ref = ref_mean_pmfs.iloc[0][mid2].mean
            L_null_ref = null_ref_mean_counts[mid1].item()
            R_null_ref = null_ref_mean_counts[mid2].item()

        #01-11 end

        # 01-10 add
        sender, receiver = index.split('|')
        mid_1 = interactions.loc[col, 'multidata_1_id']
        mid_2 = interactions.loc[col, 'multidata_2_id']
        pmf_1 = mean_pmfs.loc[sender, mid_1]
        pmf_2 = mean_pmfs.loc[receiver, mid_2]
        pmf = (pmf_1 & pmf_2) / 2
        threshold = get_threshold_from_pmf(pmf, 0.05)
        up_IS = threshold + threshold / k * 1.96 
        low_IS = threshold - threshold / k * 1.96 
        # 01-10 end

        ligand_gene_name = interactions.loc[col].multidata_1_id
        receptor_gene_name = interactions.loc[col].multidata_2_id
        ligand_gene_name = id2symbol_dict[ligand_gene_name]
        receptor_gene_name =  id2symbol_dict[receptor_gene_name]

        IS = interactions_strength.loc[index, col]
        null_IS = null_interactions_strength.loc[index, col]
        ligand_perc = L_perc.loc[index, col]
        receptor_perc = R_perc.loc[index, col]
        ref_ligand_perc = ref_L_perc.loc[index, col]
        ref_receptor_perc = ref_R_perc.loc[index, col]
        ligand_IS = p1.loc[index, col]
        receptor_IS = p2.loc[index, col]

        ligand_low = np.max(ref_p1.loc[index, col] * (1-1/k*1.96), 0)
        ligand_high = ref_p1.loc[index, col] * (1+1/k*1.96)
        receptor_low = np.max(ref_p2.loc[index, col], 0) * (1-1/k*1.96)
        receptor_high = ref_p2.loc[index, col] * (1+1/k*1.96)
        ligand_range = f"{ligand_low}-{ligand_high}"
        receptor_range = f"{receptor_low}-{receptor_high}"

        if debug_mode:
            results.append((
                index, True, col, ligand_gene_name, 
                receptor_gene_name,
                f"{IS}", null_IS, f"{low_IS}-{up_IS}", L_null_ref, R_null_ref,
                False, ligand_perc, receptor_perc, 
                True, ref_ligand_perc, ref_receptor_perc, 
                ligand_IS, ligand_range, 
                receptor_IS, receptor_range,
                False, True, "Down"
            ))
        else:
            results.append((
                index, True, col, ligand_gene_name, 
                receptor_gene_name,
                f"{IS}", null_IS, f"{low_IS}-{up_IS}", 
                False, ligand_perc, receptor_perc, 
                True, ref_ligand_perc, ref_receptor_perc, 
                ligand_IS, ligand_range, 
                receptor_IS, receptor_range,
                False, True, "Down"
            ))

    if debug_mode:
        results_df = pd.DataFrame(
            results, 
            columns=[
                'sender|receiver', 'is_in_ref', 'LRI_id', 
                'ligand', 'receptor', 'IS_score', 'null_IS',
                'Sig_thres_interval_by_ref', 'L_null_ref', 'R_null_ref',
                '>min_percentile', 'L_perc', 'R_perc',
                '>min_percentile_ref', 'L_perc_ref', 'R_perc_ref',
                'IS_L_part','IS_L_interval_by_ref',
                'IS_R_part', 'IS_R_interval_by_ref',
                'is_sig', 'is_sig_ref', 'comparison'
            ]
        )
    else:
        results_df = pd.DataFrame(
            results, 
            columns=[
                'sender|receiver', 'is_in_ref', 'LRI_id', 
                'ligand', 'receptor', 'IS_score', 'null_IS',
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
        with open(f'{save_path}/debug_results.txt', 'wt') as f:
            f.write(f"Intersection:{intersection}, Union:{union}\n")
            f.write(f"#Pred.Sig:{np.sum(prediction)}, #Real.Sig:{np.sum(real_pvals_values < 0.05)}, #valid:{len(real_pvals_values)}\n")
            f.write(f"Precision:{intersection/np.sum(prediction)}, IoU:{intersection/union}, Recall:{intersection/np.sum(real_pvals_values < 0.05)}\n")

        est_by_ref = np.array(est_by_ref)
        print(f"#by_null: {np.sum(1 - est_by_ref)}, #by_ref: {np.sum(est_by_ref)}")
        real_pvals_values = real_pvals_values[est_by_ref]
        prediction = np.array(prediction)[est_by_ref]
        intersection = np.sum(np.logical_and(prediction, real_pvals_values < 0.05))
        union = np.sum(np.logical_or(prediction, real_pvals_values < 0.05))
        print(f"Real Prec:{intersection/np.sum(prediction)}, Real IoU:{intersection/union}, Real Recall:{intersection/np.sum(real_pvals_values < 0.05)}")
        with open(f'{save_path}/debug_results.txt', 'at') as f:
            f.write(f"#by_null: {np.sum(1 - est_by_ref)}, #by_ref: {np.sum(est_by_ref)}\n")
            f.write(f"Real Prec:{intersection/np.sum(prediction)}, Real IoU:{intersection/union}, Real Recall:{intersection/np.sum(real_pvals_values < 0.05)}")
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
    k = np.nanmean(mean_by_std_list)
    logger.debug(f"k={k}")
    k = max(3, k)
    
    if debug_mode:
        import matplotlib.pyplot as plt 
        plt.figure(figsize=(6,3))
        plt.subplot(1,2,1)
        plt.scatter(hk_df.query_hk, hk_df.ref_hk, s=3, alpha=0.1)
        plt.subplot(1,2,2)
        plt.scatter(np.log(hk_df.query_hk), np.log(hk_df.ref_hk), s=3, alpha=0.1)
        plt.savefig(f'{save_path}/hk_scatter.jpg')

        plt.figure(figsize=(3,3))
        plt.scatter(mean_list, std_list)
        plt.title(np.nanmean(mean_by_std_list))
        plt.savefig(f'{save_path}/hk_mean_by_std.jpg')
        
    return k
    

def update_reference_for_first_time_activation(reference_path):
    with open(f'{reference_path}/basic_info_dict.pkl', 'rb') as f:
        basic_info_dict = pickle.load(f)
    
    from scipy.signal import fftconvolve
    from .distrib_digit import Distribution_digit
    
    n_bins = 50
    precision_digit = 0.01
    pmf_bins_digit = np.arange(0, n_bins+precision_digit - 1e-10, precision_digit)

    n_fft = 100
    gene_sum_pmf_dict = {}
    for gene in basic_info_dict:
        loc = basic_info_dict[gene]['loc'] 
        scale =  basic_info_dict[gene]['scale']
        gene_sum_pmf_dict[gene] = {1: basic_info_dict[gene]['expr_dist']}
        for item in range(2,n_fft):
            gene_sum_pmf_dict[gene][item] = fftconvolve(gene_sum_pmf_dict[gene][item-1], gene_sum_pmf_dict[gene][1])

    gene_pmf_dict = {}
    for gene in basic_info_dict:
        gene_pmf_dict[gene] = {}
        loc = basic_info_dict[gene]['loc']
        scale = basic_info_dict[gene]['scale']
        for item in range(1,n_fft):
            pmf = gene_sum_pmf_dict[gene][item]
            cdf = np.cumsum(pmf)
            pmf_array = np.diff(cdf[np.int64(pmf_bins_digit * item)],prepend=0)
            if item == 1:
                gene_pmf_dict[gene][item] = Distribution_digit('other', pmf_array=pmf_array, loc=loc, scale=scale, is_align=True)
            else:
                gene_pmf_dict[gene][item] = Distribution_digit('other', pmf_array=pmf_array, is_align=True)

    with open(f'{reference_path}/ref_gene_pmf_dict.pkl', 'wb') as f:
        pickle.dump(gene_pmf_dict, f)

    logger.success(f"Reference panel data updated.")


def infer_query_workflow(database_file_path, reference_path, query_counts_file_path, celltype_file_path, save_path, celltype_mapping_dict=None, meta_key=None, debug_mode=False):
    
    with open(f'{reference_path}/config.toml', 'rb') as f:
        config = tomllib.load(f)
        logger.info(f"Start inferring by using CCC reference: {config['reference_name']}")
        logger.info(f"Reference min_percentile = {config['min_percentile']}")
        logger.info(f"Reference LRI DB = {config['LRI_database']}")

    if not os.path.exists(f'{reference_path}/ref_gene_pmf_dict.pkl'):
        logger.info(f"Updating reference panel data for first time activation.")
        update_reference_for_first_time_activation(reference_path)

    if not os.path.exists(save_path):
        os.makedirs(save_path)
        logger.success(f"Query save dir is created.")
    else:
        logger.warning(f"{save_path} already exists, all files will be overwritten")


    query = sc.read_h5ad(query_counts_file_path)
    sc.pp.filter_cells(query, min_genes=50) # basic QC
    logger.info(f"Reading query adata, {query.shape[0]} cells x {query.shape[1]} genes")
    
    if meta_key is not None:
        labels_df = pd.DataFrame(query.obs[meta_key])
        labels_df.columns = ['cell_type']
        labels_df.index.name = 'barcode_sample'
    else:
        labels_df = pd.read_csv(celltype_file_path, sep='\t', index_col=0)
        for barcode in query.obs_names:
            assert barcode in labels_df.index, "The index of query data doesn't match the index of labels"
        labels_df = labels_df.loc[query.obs_names, :]

    query = rank_preprocess(query)
    if debug_mode:
        query.write_h5ad(f"{save_path}/debug_digit.h5ad")

    k = calculate_adjust_factor(query, reference_path, save_path, debug_mode)

    counts_df, complex_table, interactions = get_fastccc_input(query, database_file_path)
    
    if debug_mode:
        logger.debug("Entering debug process")
        from .build_reference import fastccc_for_reference
        fastccc_for_reference('', save_path, counts_df, labels_df, complex_table, interactions, min_percentile = config['min_percentile'], query_debug_mode=True)
        logger.debug("Debug ends.")

    compare_with_reference(
        counts_df, labels_df, complex_table, interactions, reference_path, save_path, 
        config=config, celltype_mapping_dict = celltype_mapping_dict, 
        database_file_path = database_file_path,
        k=k, debug_mode=debug_mode
    )
    logger.success("Inference workflow done.")