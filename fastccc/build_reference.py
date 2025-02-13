import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from tqdm import tqdm
from loguru import logger
from .preprocess import get_interactions
from . import preproc_utils
from .core import calculate_cluster_percents, analyze_interactions_percents
from .distrib_digit import Distribution_digit, get_pmf_array_from_samples_for_digitized_bins,  get_minimum_distribution_for_digit
from . import dist_complex
from . import dist_lr
from . import score
import itertools
from scipy.signal import fftconvolve
from collections import Counter
import os
import pickle

def digitize_transform(x, n_bins=50):
    def _digitize(x: np.ndarray, bins: np.ndarray, side="both") -> np.ndarray:
        assert x.ndim == 1 and bins.ndim == 1
        left_digits = np.digitize(x, bins)
        if side == "one":
            return left_digits
        right_digits = np.digitize(x, bins, right=True)
        rands = np.random.rand(len(x))  # uniform random numbers
        digits = rands * (right_digits - left_digits) + left_digits
        digits = np.ceil(digits).astype(np.int64)
        return digits
        
    # non_zero_ids = x.nonzero()
    # non_zero_row = x[non_zero_ids]
    '''
    input x 就是 非0的，直接针对csr，coo数据的
    '''
    bins = np.quantile(x, np.linspace(0, 1, n_bins - 1))
    non_zero_digits = _digitize(x, bins)
    return non_zero_digits
    

def calculate_L_R_and_IS_score(mean_counts, interactions):
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

def calculate_L_R_and_IS_percents(cluster_percents, interactions, threshold=0.1, sep='|'):
    
    p1_index = []
    p2_index = []
    all_index = []
    for i in itertools.product(sorted(cluster_percents.index), sorted(cluster_percents.index)):
        p1_index.append(i[0])
        p2_index.append(i[1])
        all_index.append(sep.join(i))
        
    p1 = cluster_percents.loc[p1_index, interactions['multidata_1_id']]
    p2 = cluster_percents.loc[p2_index, interactions['multidata_2_id']]
    p1.columns = interactions.index
    p2.columns = interactions.index
    p1.index = all_index
    p2.index = all_index
    
    interactions_strength = (p1>threshold) * (p2>threshold)
    # print((p1>threshold) * (p2>threshold))
    
    return p1, p2, interactions_strength

def calculate_mean_pmfs(counts_df, labels_df, complex_table, gene_pmf_dict, n_fft=100):
    meta_dict = Counter(labels_df.cell_type)
    ####### clusters_mean #######
    clusters_mean_dict = {}
    for celltype in sorted(meta_dict):
        clusters_mean_dict[celltype]  = {}
        n_sum = meta_dict[celltype]
        if n_sum < n_fft:
            for gene in counts_df.columns:
                if gene not in gene_pmf_dict:
                    continue
                else:
                    clusters_mean_dict[celltype][gene] = gene_pmf_dict[gene][n_sum]
        else:
            for gene in counts_df.columns:
                if gene not in gene_pmf_dict:
                    continue
                else:
                    clusters_mean_dict[celltype][gene] = gene_pmf_dict[gene][1] ** n_sum / n_sum
    mean_pmfs = pd.DataFrame(clusters_mean_dict).T
    complex_func = get_minimum_distribution_for_digit
    mean_pmfs = dist_complex.combine_complex_distribution_df(mean_pmfs, complex_table, complex_func)
    return mean_pmfs


def rank_preprocess(adata):
    np.random.seed(42) # add seed to ensure reproduiablity
    assert issparse(adata.X), "Anndata.X should be a sparse matrix format."
    if adata.shape[1] < 5000:
        logger.warning("Do you use whole transcriptomes? Raw data w\o filtering genes should work better.")

    for i in tqdm(range(adata.shape[0]), desc="Ranking genes for cells", unit="cell", 
              bar_format="{l_bar}{bar} | {n_fmt}/{total_fmt} cells completed", leave=False):
        indices = slice(adata.X.indptr[i], adata.X.indptr[i+1])
        x = adata.X.data[indices]
        adata.X.data[indices] = digitize_transform(x)
    logger.success("Rank preprocess done.")
    return adata

def get_fastccc_input(adata, lrdb_file_path, convert_type = 'hgnc_symbol'):
    logger.info("Loading LRIs database. hgnc_symbol as gene name is requested.")
    interactions = get_interactions(lrdb_file_path)
    ##### gene_table ########
    gene_table = pd.read_csv(os.path.join(lrdb_file_path, 'gene_table.csv'))
    protein_table = pd.read_csv(os.path.join(lrdb_file_path, 'protein_table.csv'))
    gene_table = gene_table.merge(protein_table, left_on='protein_id', right_on='id_protein')
    #########################


    ##### complex_table ######
    complex_composition = pd.read_csv(os.path.join(lrdb_file_path, 'complex_composition_table.csv'))
    complex_table = pd.read_csv(os.path.join(lrdb_file_path, 'complex_table.csv'))
    complex_table = complex_table.merge(complex_composition, left_on='complex_multidata_id', right_on='complex_multidata_id')

    # 让我们只关注 'complex_multidata_id'，'protein_multidata_id'
    complex_table = complex_table[['complex_multidata_id','protein_multidata_id']]
    '''
    complex_table(pandas.DataFrame):
    =======================================================
            | complex_multidata_id   |  protein_multidata_id
    -------------------------------------------------------
    0      |             1355       |          1134
    1      |             1356       |          1175
    2      |             1357       |          1167
    =======================================================
    '''
    ##########################

    ##### feature to id conversion  ######
    # 不在 标准列表 里的 gene 就不要了
    tmp = gene_table[[convert_type, 'protein_multidata_id']]
    tmp = tmp.drop_duplicates()
    tmp.set_index('protein_multidata_id', inplace=True)

    select_columns = []
    columns_names = []
    for foo, boo in zip(tmp.index, tmp[convert_type]):
        if boo in adata.var_names:#counts.columns:
            select_columns.append(boo)
            columns_names.append(foo)

    reduced_counts = adata[:, select_columns].to_df()
    reduced_counts.columns = columns_names
    reduced_counts = reduced_counts.T.groupby(reduced_counts.columns).mean().T
    # FutureWarning: DataFrame.groupby with axis=1 is deprecated. Do `frame.T.groupby(...)` without axis instead.
    # reduced_counts = reduced_counts.groupby(reduced_counts.columns, axis=1).mean()
    ######################################

    ########## filter genes  ############
    # gene 在 所有 cell 上为 0 不要
    reduced_counts = preproc_utils.filter_empty_genes(reduced_counts)
    ######################################
    
    ######################################################################
    #                      3.Other DF filtered                       #
    ######################################################################
    # 一个 interaction 可能只有 partA 存在，但是 partB 不存在
    # 只有 如果 任意一部分不存在， 另一部分没必要参与后续计算

    ##### delete item not involved interactions ####

    foo_dict = complex_table.groupby('complex_multidata_id').apply(lambda x: list(x['protein_multidata_id'].values)).to_dict()
    '''
    dictionary complex_id: [protein_id_1, pid2, pid3, ...]
    foo_dict = {
        1355: [1134],
        1356: [1175],
        xxxx: [AAAA, BBBB, CCCC],
    }
    '''

    def __content__(key):
        if key not in foo_dict:
            return [key]
        else:
            return foo_dict[key]

    def __exist__(key, df):
        # 目前的 complex 策略就是 全部都要有
        # 经测试，这是cpdb用的策略
        for item in __content__(key):
            if item not in df.columns:
                return False
        return True

    temp_list = []
    temp_dict = {}
    for item in reduced_counts.columns:
        temp_dict[item] = False

    # 注释是为了验证 interactions 的过滤策略， 完全一致
    # print(interactions)
    select_index = []
    for partA, partB in zip(interactions.multidata_1_id, interactions.multidata_2_id):
        if __exist__(partA, reduced_counts) and __exist__(partB, reduced_counts):
            temp_list.extend([partA, partB])
            select_index.append(True)
        else:
            select_index.append(False)
    interactions_filtered = interactions[select_index]

    for item in temp_list:
        for subitem in __content__(item):
            if subitem in temp_dict:
                temp_dict[subitem] = True
    select_index = [key for key in temp_dict if temp_dict[key]]
    reduced_counts = reduced_counts[select_index]

    counts_df = reduced_counts
    temp_list = set(temp_list)
    select_index = [True if item in temp_list else False for item in complex_table.complex_multidata_id]
    complex_table = complex_table[select_index]
    interactions = interactions_filtered
    logger.success("Requested data for fastccc is prepared.")
    return counts_df, complex_table, interactions


def fastccc_for_reference(reference_name, save_path, counts_df, labels_df, complex_table, interactions, min_percentile = 0.1, ref_debug_mode=False, query_debug_mode=False, for_uploading=False):
    logger.info("Running FastCCC.")
    mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
    complex_func = score.calculate_complex_min_func
    mean_counts = score.combine_complex_distribution_df(mean_counts, complex_table, complex_func)
    percents = calculate_cluster_percents(counts_df, labels_df, complex_table)

    n_bins = 50
    precision_digit = 0.01
    pmf_bins_digit = np.arange(0, n_bins+precision_digit - 1e-10, precision_digit)

    ####### 
    logger.info("Calculating null distributions.")
    n_fft = 100
    gene_sum_pmf_dict = {}
    basic_info_dict = {}
    for gene in counts_df.columns:
        samples = counts_df[gene].values

        loc = np.mean(samples)
        scale = np.std(samples)
        basic_info_dict[gene] = {'loc':loc, 'scale':scale}

        gene_sum_pmf_dict[gene] = {1: get_pmf_array_from_samples_for_digitized_bins(samples)}
        basic_info_dict[gene]['expr_dist'] = gene_sum_pmf_dict[gene][1]

        for item in range(2,n_fft):
            gene_sum_pmf_dict[gene][item] = fftconvolve(gene_sum_pmf_dict[gene][item-1], gene_sum_pmf_dict[gene][1])

    gene_pmf_dict = {}
    for gene in counts_df.columns:
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

    if ref_debug_mode or query_debug_mode:
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
        L_perc, R_perc, percents_analysis = calculate_L_R_and_IS_percents(percents, interactions, threshold=min_percentile)
        
        meta_dict = Counter(labels_df.cell_type)
        ####### clusters_mean #######
        clusters_mean_dict = {}
        for celltype in sorted(meta_dict):
            clusters_mean_dict[celltype]  = {}
            n_sum = meta_dict[celltype]
            if n_sum < n_fft:
                for gene in counts_df.columns:
                    clusters_mean_dict[celltype][gene] = gene_pmf_dict[gene][n_sum]
            else:
                for gene in counts_df.columns:
                    clusters_mean_dict[celltype][gene] = gene_pmf_dict[gene][1] ** n_sum / n_sum
        mean_pmfs = pd.DataFrame(clusters_mean_dict).T
        complex_func = get_minimum_distribution_for_digit
        mean_pmfs = dist_complex.combine_complex_distribution_df(mean_pmfs, complex_table, complex_func)

        logger.info("Calculating sig. LRIs.")
        pvals = dist_lr.calculate_key_interactions_pvalue(
            mean_pmfs, interactions, interactions_strength, percents_analysis, method='Arithmetic'
        )

        if query_debug_mode:
            pvals.to_csv(f'{save_path}/debug_pvals.txt', sep='\t')
            return

        if ref_debug_mode:
            pvals.to_csv(f'{save_path}/ref_pvals.txt', sep='\t')
            percents_analysis.to_csv(f'{save_path}/ref_percents_analysis.txt', sep='\t')
            L_perc.to_csv(f'{save_path}/ref_percents_L.txt', sep='\t')
            R_perc.to_csv(f'{save_path}/ref_percents_R.txt', sep='\t')
            # interactions_strength.to_csv(f'{save_path}/ref_interactions_strength.csv')
            p1.to_csv(f'{save_path}/ref_interactions_strength_L.txt', sep='\t')
            p2.to_csv(f'{save_path}/ref_interactions_strength_R.txt', sep='\t')

    ####### save reference results #######
    logger.info("Saving reference.")
    if for_uploading:
        with open(f'{save_path}/basic_info_dict.pkl', 'wb') as f:
            pickle.dump(basic_info_dict, f)
    else:
        with open(f'{save_path}/ref_gene_pmf_dict.pkl', 'wb') as f:
            pickle.dump(gene_pmf_dict, f)
    with open(f'{save_path}/ref_percents.pkl', 'wb') as f:
        pickle.dump(percents, f)
    with open(f'{save_path}/ref_mean_counts.pkl', 'wb') as f:
        pickle.dump(mean_counts, f)
    with open(f'{save_path}/complex_table.pkl', 'wb') as f:
        pickle.dump(complex_table, f)
    with open(f'{save_path}/interactions.pkl', 'wb') as f:
        pickle.dump(interactions, f)


def record_hk_genes(adata):
    from .hk_genes import housekeeping_genes
    select_index = [item for item in adata.var_names if item in housekeeping_genes]
    hk_adata = adata[:, select_index]
    mean_hk_rnk = hk_adata.X.mean(axis=0)
    return mean_hk_rnk, select_index

def record_adjustment_info(adata, save_path):
    mean_hk_rnk, gene_index = record_hk_genes(adata)
    ref_hk = pd.DataFrame(np.array(mean_hk_rnk).flatten(), index=gene_index, columns=['ref_hk'])
    ref_hk.to_csv(f'{save_path}/ref_hk.txt', sep='\t')


reference_config = {}

from datetime import date

def dumps(toml_dict, table=""):
    document = []
    for key, value in toml_dict.items():
        match value:
            case dict():
                table_key = f"{table}.{key}" if table else key
                document.append(
                    f"\n[{table_key}]\n{_dumps_dict(value)}"
                )
            case _:
                document.append(f"{key} = {_dumps_value(value)}")
    return "\n".join(document)

def _dumps_dict(toml_dict):
    document = []
    for key, value in toml_dict.items():
        key = f'"{key}"'
        document.append(f"{key} = {_dumps_value(value)}")
    return "\n".join(document)

def _dumps_value(value):
    match value:
        case bool():
            return "true" if value else "false"
        case float() | int():
            return str(value)
        case str():
            return f'"{value}"'
        case date():
            return value.isoformat()
        case list():
            return f"[{', '.join(_dumps_value(v) for v in value)}]"
        case _:
            raise TypeError(
                f"{type(value).__name__} {value!r} is not supported"
            )

def save_config(save_path):
    logger.info("Saving reference config.")
    save_content = dumps(reference_config)
    with open(f'{save_path}/config.toml', 'w') as f:
        f.write(save_content) 


def build_reference_workflow(database_file_path, reference_counts_file_path, celltype_file_path, reference_name, save_path, meta_key=None, min_percentile = 0.1, debug_mode=False, for_uploading=False):
    logger.info(f"Start building CCC reference.")

    reference_config['reference_name'] = reference_name
    reference_config['min_percentile'] = min_percentile
    if database_file_path.endswith('/'):
        reference_config['LRI_database'] = database_file_path[:-1].split('/')[-1]
    else:
        reference_config['LRI_database'] = database_file_path.split('/')[-1]

    logger.info(f"Reference_name = {reference_config['reference_name']}")
    logger.info(f"min_percentile = {reference_config['min_percentile']}")
    logger.info(f"LRI database = {reference_config['LRI_database']}")

    save_path = os.path.join(save_path, reference_name)
    if not os.path.exists(save_path):
        os.makedirs(save_path)
        logger.success(f"Reference save dir {save_path} is created.")
    else:
        logger.warning(f"{save_path} already exists, all files will be overwritten")

    reference = sc.read_h5ad(reference_counts_file_path)
    sc.pp.filter_cells(reference, min_genes=50)
    logger.info(f"Reading reference adata, {reference.shape[0]} cells x {reference.shape[1]} genes.")

    if meta_key is not None:
        labels_df = pd.DataFrame(reference.obs[meta_key])
        labels_df.columns = ['cell_type']
        labels_df.index.name = 'barcode_sample'
    else:
        labels_df = pd.read_csv(celltype_file_path, sep='\t', index_col=0)
        for barcode in reference.obs_names:
            assert barcode in labels_df.index, "The index of query data doesn't match the index of labels"
        labels_df = labels_df.loc[reference.obs_names, :]
    
    ct_counter = Counter(labels_df['cell_type'])
    reference_config['celltype'] = ct_counter
    
    reference = rank_preprocess(reference)
    record_adjustment_info(reference, save_path)
    counts_df, complex_table, interactions = get_fastccc_input(reference, database_file_path)
    fastccc_for_reference(reference_name, save_path, counts_df, labels_df, complex_table, interactions, min_percentile, ref_debug_mode=debug_mode, for_uploading=for_uploading)
    save_config(save_path)
    logger.success(f"Reference '{reference_name}' is built.")

