import pandas as pd
import numpy as np
from .distrib import Distribution, get_pvalue_from_pmf
from collections import Counter
import itertools
import timeit
from .ccc_utils import get_current_memory
from . import preprocess
from . import score
from . import dist_iid_set
from . import dist_complex 
from . import dist_lr
import warnings




def statistical_analysis_method(
    database_file_path,
    celltype_file_path,
    counts_file_path,
    convert_type = 'hgnc_symbol',
    cluster_distrib_method = 'Mean',
    complex_distrib_method = 'Minimum',
    LR_distrib_method = 'Arithmetic',
    quantile = 0.9,
    min_percentile = 0.1,
    style = None
    
):
    method_qt_dict = {
        "Median": 0.5,
        "Q2": 0.5,
        "Q3": 0.75,
        "Quantile": quantile
    }
    if cluster_distrib_method in ["Median", "Q2", "Q3", "Quantile"]:
        quantile = method_qt_dict[cluster_distrib_method]
        cluster_distrib_method = 'Quantile'
        assert quantile < 1 and quantile > 0, "The quantile should be within the range (0, 1)."
        if quantile < 0.5:
            warnings.warn(f"Invalid quantile: {quantile}. The quantile is too low to be effective.", UserWarning)
        min_percentile = max(min_percentile, 1-quantile)
    
    if style not in {None, "cpdb"}:
        raise ValueError(f"Invalid style: {style}. Must be one of [None, 'cpdb']")
    if style == 'cpdb':
        cluster_distrib_method = 'Mean'
        complex_distrib_method = 'Minimum'
        LR_distrib_method = 'Arithmetic'
        
    counts_df, labels_df, complex_table, interactions = preprocess.get_input_data(
        database_file_path, 
        celltype_file_path,
        counts_file_path,
        convert_type
    )
    
    # Stage I : calculate L-R expression score:
    ## Scores for clusters
    if cluster_distrib_method == 'Mean':
        mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
    elif cluster_distrib_method == 'Quantile':
        mean_counts = score.calculate_cluster_quantile(counts_df, labels_df, quantile)
    ## Scores for complex
    if complex_distrib_method == 'Minimum':
        complex_func = score.calculate_complex_min_func
    elif complex_distrib_method == 'Average':
        complex_func = score.calculate_complex_mean_func
    mean_counts = score.combine_complex_distribution_df(mean_counts, complex_table, score.calculate_complex_min_func)
    ## Scores for L-R expression
    interactions_strength = score.calculate_interactions_strength(mean_counts, interactions, method=LR_distrib_method)
    
    # Stage II : Percentages
    percents = calculate_cluster_percents(counts_df, labels_df, complex_table)
    percents_analysis = analyze_interactions_percents(percents, interactions, threshold=min_percentile)
    
    # Stage III : Null Distribution
    if cluster_distrib_method == 'Mean':
        mean_pmfs = dist_iid_set.calculate_cluster_mean_distribution(counts_df, labels_df)
    elif cluster_distrib_method == 'Quantile':
        mean_pmfs = dist_iid_set.calculate_cluster_quantile_distribution(counts_df, labels_df, quantile)
    if complex_distrib_method == 'Minimum':
        complex_func = dist_complex.get_minimum_distribution
    elif complex_distrib_method == 'Average':
        complex_func = dist_complex.get_average_distribution
    mean_pmfs = dist_complex.combine_complex_distribution_df(mean_pmfs, complex_table, complex_func)
    
    # Stage IV : P-values
    pvals = dist_lr.calculate_key_interactions_pvalue(
        mean_pmfs, interactions, interactions_strength, percents_analysis, method=LR_distrib_method
    )
    return interactions_strength, pvals, percents_analysis
    
        
    



def calculate_cluster_mean(counts_df, labels_df, complex_table):
    '''
    input:
        counts_df(pandas.DataFrame): 
            sample * feature matrix. 
            The expression level of every gene in every cell.
            i.e.:
                =======================================
                      | gene1 | gene2 | gene3 | gene4 
                =======================================
                cell1 |  0.4  |  0.1  |  0.8  |  0.0
                cell2 |  0.0  |  0.6  |  0.0  |  0.4  
                cell3 |  0.0  |  0.1  |  0.0  |  0.0  
                =======================================

        labels_df(pandas.DataFrame): 
            sample * 2 matrix. 
            col_names must be (index, cell_type)
            Col1 is sample name the same as counts_df. 
            Col2 is celltypes or labels
            i.e.:
                =================
                sample | celltype
                =================
                cell1  |  B cell 
                cell2  |  T cell 
                cell3  |  B cell
                =================



        complex_table(pandas.DataFrame): 
            n * 2 matrix. 
            Col1 is complex compound ID.
            Col2 protein ID.
            i.e.:
                C1 is composed of gene1, gene2, gene3. The data is:
                ===========
                cid |  pid
                ===========
                C1  | gene1
                C1  | gene2
                C1  | gene3
                ===========


    output:
        mean_counts(pandas.DataFrame):
            celltype * feature(with complex feature) matrix. 
            item = average level of expression for each feature within each celltype. 
    '''

    ###################### Check Data ##########################



    ############################################################

    counts_df_with_labels = counts_df.merge(labels_df, left_index=True, right_index=True)
    mean_counts = counts_df_with_labels.groupby('cell_type').mean()

    # complex count dataframe 
    def create_complex_count_func(x):
        x = [sub_x for sub_x in x if sub_x in mean_counts.columns]
        if len(x) == 0:
            return pd.Series(index=mean_counts.index)
        return mean_counts.loc[:,x].T.min()
    
    if not complex_table.empty:
        complex_counts = complex_table.groupby('complex_multidata_id').apply(
            lambda x: x['protein_multidata_id'].values).apply(
            lambda x:create_complex_count_func(x)).T

        # 合成 最终 mean_counts dataframe
        mean_counts = pd.concat((mean_counts, complex_counts), axis=1)
    
    ##############################################################
    #                         Return                             #
    ##############################################################
    mean_counts = mean_counts.dropna(axis=1)
    return mean_counts


def calculate_cluster_quantile(counts_df, labels_df, complex_table, qt=0.9):
    '''
    input:
        counts_df(pandas.DataFrame): 
            sample * feature matrix. 
            The expression level of every gene in every cell.
            i.e.:
                =======================================
                      | gene1 | gene2 | gene3 | gene4 
                =======================================
                cell1 |  0.4  |  0.1  |  0.8  |  0.0
                cell2 |  0.0  |  0.6  |  0.0  |  0.4  
                cell3 |  0.0  |  0.1  |  0.0  |  0.0  
                =======================================

        labels_df(pandas.DataFrame): 
            sample * 2 matrix. 
            col_names must be (index, cell_type)
            Col1 is sample name the same as counts_df. 
            Col2 is celltypes or labels
            i.e.:
                =================
                sample | celltype
                =================
                cell1  |  B cell 
                cell2  |  T cell 
                cell3  |  B cell
                =================



        complex_table(pandas.DataFrame): 
            n * 2 matrix. 
            Col1 is complex compound ID.
            Col2 protein ID.
            i.e.:
                C1 is composed of gene1, gene2, gene3. The data is:
                ===========
                cid |  pid
                ===========
                C1  | gene1
                C1  | gene2
                C1  | gene3
                ===========


    output:
        mean_counts(pandas.DataFrame):
            celltype * feature(with complex feature) matrix. 
            item = average level of expression for each feature within each celltype. 
    '''

    ###################### Check Data ##########################



    ############################################################

    counts_df_with_labels = counts_df.merge(labels_df, left_index=True, right_index=True)
    # mean_counts = counts_df_with_labels.groupby('cell_type').quantile(qt)
    mean_counts = counts_df_with_labels.groupby('cell_type').apply(lambda x: pd.Series(np.quantile(x,qt, axis=0, method='lower'), index=x.columns))

    # complex count dataframe 
    def create_complex_count_func(x):
        x = [sub_x for sub_x in x if sub_x in mean_counts.columns]
        if len(x) == 0:
            return pd.Series(index=mean_counts.index)
        return mean_counts.loc[:,x].T.min()
    
    if not complex_table.empty:
        complex_counts = complex_table.groupby('complex_multidata_id').apply(
            lambda x: x['protein_multidata_id'].values).apply(
            lambda x:create_complex_count_func(x)).T

        # 合成 最终 mean_counts dataframe
        mean_counts = pd.concat((mean_counts, complex_counts), axis=1)
    
    ##############################################################
    #                         Return                             #
    ##############################################################
    mean_counts = mean_counts.dropna(axis=1)
    return mean_counts




def calculate_cluster_percents(counts_df, labels_df, complex_table):
    '''
    input:
        counts_df(pandas.DataFrame): 
            sample * feature matrix. 
            The expression level of every gene in every cell.
            i.e.:
                =======================================
                      | gene1 | gene2 | gene3 | gene4 
                =======================================
                cell1 |  0.4  |  0.1  |  0.8  |  0.0
                cell2 |  0.0  |  0.6  |  0.0  |  0.4  
                cell3 |  0.0  |  0.1  |  0.0  |  0.0  
                =======================================

        labels_df(pandas.DataFrame): 
            sample * 2 matrix. 
            col_names must be (index, cell_type)
            Col1 is sample name the same as counts_df. 
            Col2 is celltypes or labels
            i.e.:
                =================
                sample | celltype
                =================
                cell1  |  B cell 
                cell2  |  T cell 
                cell3  |  B cell
                =================



        complex_table(pandas.DataFrame): 
            n * 2 matrix. 
            Col1 is complex compound ID.
            Col2 protein ID.
            i.e.:
                C1 is composed of gene1, gene2, gene3. The data is:
                ===========
                cid |  pid
                ===========
                C1  | gene1
                C1  | gene2
                C1  | gene3
                ===========


    output:
        mean_counts(pandas.DataFrame):
            celltype * feature(with complex feature) matrix. 
            item = average level of expression for each feature within each celltype. 
    '''

    ###################### Check Data ##########################



    ############################################################
    
    counts_df = counts_df > 0
    counts_df_with_labels = counts_df.merge(labels_df, left_index=True, right_index=True)
    mean_counts = counts_df_with_labels.groupby('cell_type').mean()

    # complex count dataframe 
    def create_complex_count_func(x):
        x = [sub_x for sub_x in x if sub_x in mean_counts.columns]
        if len(x) == 0:
            return pd.Series(index=mean_counts.index)
        return mean_counts.loc[:,x].T.min()
    if not complex_table.empty:
        complex_counts = complex_table.groupby('complex_multidata_id').apply(
            lambda x: x['protein_multidata_id'].values).apply(
            lambda x:create_complex_count_func(x)).T

        # 合成 最终 mean_counts dataframe
        mean_counts = pd.concat((mean_counts, complex_counts), axis=1)
    
    ##############################################################
    #                         Return                             #
    ##############################################################
    mean_counts = mean_counts.dropna(axis=1)
    return mean_counts




# def calculate_cluster_mean_distribution(counts_df, labels_df, complex_table):
#     meta_dict = Counter(labels_df.cell_type)
    
#     ####### 
#     gene_sum_pmf_dict = {}
#     for gene in counts_df.columns:
#         samples = counts_df[gene].values
#         gene_sum_pmf_dict[gene] = {1: get_distribution_from_samples(samples)}
#         for item in range(2,30):
#             gene_sum_pmf_dict[gene][item] = gene_sum_pmf_dict[gene][item-1] + gene_sum_pmf_dict[gene][1]
            
#     for gene in counts_df.columns:
#         for item in range(2,30):
#             gene_sum_pmf_dict[gene][item] = gene_sum_pmf_dict[gene][item] / item
            
#     ####### clusters_mean #######
#     clusters_mean_dict = {}
#     for celltype in meta_dict:
#         clusters_mean_dict[celltype]  = {}
#         n_sum = meta_dict[celltype]
#         if n_sum < 30:
#             for gene in counts_df.columns:
#                 clusters_mean_dict[celltype][gene] = gene_sum_pmf_dict[gene][n_sum]
#         else:
#             for gene in counts_df.columns:
#                 clusters_mean_dict[celltype][gene] = gene_sum_pmf_dict[gene][1] ** n_sum / n_sum
#     mean_pmf = pd.DataFrame(clusters_mean_dict).T
    
#     def sub_func(x):
#         return get_minimum_distribution(*x)

#     def func(x):
#         x = [sub_x for sub_x in x if sub_x in mean_pmf.columns]
#         if len(x) == 0:
#             return pd.Series(index=mean_pmf.index)
#         return mean_pmf.loc[:,x].apply(sub_func, axis=1)
    
#     if not complex_table.empty:
#         complex_pmf = complex_table.groupby('complex_multidata_id').apply(lambda x: x['protein_multidata_id'].values).apply(lambda x:func(x)).T
#         complex_pmf = complex_pmf.dropna(axis=1)
#         mean_pmf = pd.concat((mean_pmf, complex_pmf), axis=1)
        
#     return mean_pmf


def calculate_interactions_strength(mean_counts, interactions, sep='|'):
    
    # 
    p1_index = []
    p2_index = []
    all_index = []
    for i in itertools.product(sorted(mean_counts.index), sorted(mean_counts.index)):
        p1_index.append(i[0])
        p2_index.append(i[1])
        all_index.append(sep.join(i))
        
    p1 = mean_counts.loc[p1_index, interactions['multidata_1_id']]
    p2 = mean_counts.loc[p2_index, interactions['multidata_2_id']]
    p1.columns = interactions.index
    p2.columns = interactions.index
    p1.index = all_index
    p2.index = all_index
    
    interactions_strength = (p1 + p2)/2 * (p1 > 0) * (p2>0)
    
    return interactions_strength

def calculate_interactions_strength_multiply(mean_counts, interactions, sep='|'):
    
    # 
    p1_index = []
    p2_index = []
    all_index = []
    for i in itertools.product(sorted(mean_counts.index), sorted(mean_counts.index)):
        p1_index.append(i[0])
        p2_index.append(i[1])
        all_index.append(sep.join(i))
        
    p1 = mean_counts.loc[p1_index, interactions['multidata_1_id']]
    p2 = mean_counts.loc[p2_index, interactions['multidata_2_id']]
    p1.columns = interactions.index
    p2.columns = interactions.index
    p1.index = all_index
    p2.index = all_index
    
    interactions_strength = np.sqrt(p1*p2)#(p1 + p2)/2 * (p1 > 0) * (p2>0)
    
    return interactions_strength


def analyze_interactions_percents(cluster_percents, interactions, threshold=0.1, sep='|'):
    
    # 
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
    
    return interactions_strength


    
    
def calculate_key_interactions_pvalue(mean_pmf, interactions, interactions_strength, percent_analysis):
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
    pval_pmfs = (p1_items & p2_items) / 2
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


########################
def calculate_key_interactions_pvalue_multiply_version(mean_pmf, interactions, interactions_strength, percent_analysis):
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
########################



def calculate_key_interactions_pvalue_by_part(mean_counts, mean_pmf, interactions, interactions_strength, percent_analysis):
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
    
    mean_count_p1 = mean_counts.loc[p1_index, interactions['multidata_1_id']]
    mean_count_p2 = mean_counts.loc[p2_index, interactions['multidata_2_id']]
    mean_count_p1.columns = interactions.index
    mean_count_p2.columns = interactions.index
    mean_count_p1.index = all_index
    mean_count_p2.index = all_index
    mean_count_p1_item = mean_count_p1.values[np.where(percent_analysis)]
    mean_count_p2_item = mean_count_p2.values[np.where(percent_analysis)]
    
    
    est1 = []
    for i, value in enumerate(mean_count_p1_item):
        pval_est = get_pvalue_from_pmf(value, p1_items[i])
        est1.append(pval_est)
    est1 = np.array(est1) < 0.05
        
    est2 = []
    for i, value in enumerate(mean_count_p2_item):
        pval_est = get_pvalue_from_pmf(value, p2_items[i])
        est2.append(pval_est)
    est2 = np.array(est2) < 0.05
        
    significant_flag_mat = np.zeros_like(interactions_strength, dtype=bool)
    significant_flag_mat[np.where(percent_analysis)] = np.logical_and(est1, est2)
    significant_flag_mat = pd.DataFrame(significant_flag_mat, index=interactions_strength.index, columns=interactions_strength.columns)
    return significant_flag_mat