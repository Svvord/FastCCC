import pandas as pd
import numpy as np
from .distrib import Distribution, get_pvalue_from_pmf
from collections import Counter
import itertools
import timeit
from .ccc_utils import get_current_memory, create_significant_interactions_df
from . import preprocess
from . import score
from . import dist_iid_set
from . import dist_complex 
from . import dist_lr
import warnings
import datetime
import os
from .cauchy_combine import cauthy_combine
import uuid
from loguru import logger
import glob

def generate_task_id(length=6):
    unique_part = uuid.uuid4().hex[:length]
    return unique_part

def mkdir(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        logger.info(f"Directory created: {directory_path}")
    else:
        logger.info(f"Directory already exists: {directory_path}")


def __save(dataframe, save_path, task_id, timestamp, method_key, suffix, index=True):
    if method_key == '':
        save_file = os.path.join(save_path, f"{task_id}_{suffix}.tsv")
    else:
        save_file = os.path.join(save_path, f"{task_id}_{timestamp}_{method_key}_{suffix}.tsv")
    dataframe.to_csv(save_file, sep='\t', index=index)

def __save_file(
    interactions_strength, 
    pvals, 
    percents_analysis, 
    save_path, 
    task_id,
    timestamp = None,
    method_key = ''
):
    assert os.path.isdir(save_path), "{save_path} doesn't exist or not a dir"
    if timestamp is None:
        now = datetime.datetime.now()
        timestamp = now.strftime("%Y%m%d_%H%M%S")
    __save(interactions_strength, save_path, task_id, timestamp, method_key, 'interactions_strength')
    __save(pvals, save_path, task_id, timestamp, method_key, 'pvals')
    __save(percents_analysis, save_path, task_id, timestamp, method_key, 'percents_analysis')


def Cauchy_combination_of_statistical_analysis_methods(
    database_file_path,
    celltype_file_path,
    counts_file_path,
    convert_type = 'hgnc_symbol',
    single_unit_summary_list = ['Mean', 'Median', 'Q3', 'Quantile_0.9'],
    complex_aggregation_list = ['Minimum', 'Average'],
    LR_combination_list = ['Arithmetic', 'Geometric'],
    min_percentile = 0.1,
    save_path = None,
    meta_key = None, 
    select_list = [], 
    filter_ = False,
    use_DEG = False
):
    # Generate unique task ID.
    task_id = generate_task_id()
    logger.info(f"Task id is {task_id}.")

    # Set save directory.
    if save_path is None:
        save_path = './results/'
        logger.info(f"Results will be saved to \"{save_path}\".")
    mkdir(save_path)

    method_qt_dict = {
        "Median": 0.5,
        "Q2": 0.5,
        "Q3": 0.75,
    }

    counts_df, labels_df, complex_table, interactions = preprocess.get_input_data(
        database_file_path, 
        celltype_file_path,
        counts_file_path,
        convert_type,
        meta_key = meta_key,
        select_list = select_list,
        filter_ = filter_
    )
    logger.success("Data preprocessing done.")

    for cluster_distrib_method in single_unit_summary_list:
        cluster_distrib_key = cluster_distrib_method
        current_min_percentile = min_percentile
        # check input
        if cluster_distrib_method.startswith('Quantile_'):
            try:
                quantile = float(cluster_distrib_method.split('_')[-1])
                method_qt_dict["Quantile"] = quantile
                cluster_distrib_method = "Quantile"
            except ValueError:
                print('Quantile parameter in cluster_distrib_metho_list should be e.g. "Quantile_0.9"')
                raise ValueError
        if cluster_distrib_method in ["Median", "Q2", "Q3", "Quantile"]:
            quantile = method_qt_dict[cluster_distrib_method]
            cluster_distrib_method = 'Quantile'
            assert quantile < 1 and quantile > 0, "The quantile should be within the range (0, 1)."
            if quantile < 0.5:
                warnings.warn(f"Invalid quantile: {quantile}. The quantile is too low to be effective.", UserWarning)
            current_min_percentile = max(current_min_percentile, 1-quantile)
            logger.debug(f"Adjusted percentile is {current_min_percentile}")

        # Stage I : Percentages
        percents = calculate_cluster_percents(counts_df, labels_df, complex_table)
        percents_analysis = analyze_interactions_percents(percents, interactions, threshold=current_min_percentile)

        # Stage II : Counts Distribution
        ## Scores for clusters
        if cluster_distrib_method == 'Mean':
            mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
        elif cluster_distrib_method == 'Quantile':
            mean_counts = score.calculate_cluster_quantile(counts_df, labels_df, quantile)

        # Stage III : Null Distribution
        if cluster_distrib_method == 'Mean':
            mean_pmfs = dist_iid_set.calculate_cluster_mean_distribution(counts_df, labels_df)
        elif cluster_distrib_method == 'Quantile':
            mean_pmfs = dist_iid_set.calculate_cluster_quantile_distribution(counts_df, labels_df, quantile)

        for complex_distrib_method in complex_aggregation_list:
            ## Counts for complex
            if complex_distrib_method == 'Minimum':
                complex_func = score.calculate_complex_min_func
            elif complex_distrib_method == 'Average':
                complex_func = score.calculate_complex_mean_func
            mean_counts_with_complex = score.combine_complex_distribution_df(mean_counts, complex_table, complex_func)
            ## Pmfs for complex
            if complex_distrib_method == 'Minimum':
                complex_func = dist_complex.get_minimum_distribution
            elif complex_distrib_method == 'Average':
                complex_func = dist_complex.get_average_distribution
            mean_pmfs_with_complex = dist_complex.combine_complex_distribution_df(mean_pmfs, complex_table, complex_func)
        

            for LR_distrib_method in LR_combination_list:
                method_key = f'{cluster_distrib_key}_{complex_distrib_method}_{LR_distrib_method}'

                logger.info(f"Running:\n-> {cluster_distrib_key} for single-unit summary function.\n"
                    + f"-> {complex_distrib_method} for multi-unit complex aggregation.\n"
                    + f"-> {LR_distrib_method} for L-R combination to compute the CS.\n"
                    + f"-> Percentile is {current_min_percentile}.")

                # Stage IV : P-values
                ## calculate L-R expression score:
                interactions_strength = score.calculate_interactions_strength(mean_counts_with_complex, interactions, method=LR_distrib_method)
                ## P-values for L-R expression:
                pvals = dist_lr.calculate_key_interactions_pvalue(
                    mean_pmfs_with_complex, interactions, interactions_strength, percents_analysis, method=LR_distrib_method
                )
                __save_file(interactions_strength, pvals, percents_analysis, save_path, task_id=task_id, method_key=method_key)
                logger.success("CS scoring module calculation done.")      
    logger.success("All scoring modules calculation done.")  
    cauthy_combine(save_path, task_id=task_id)
    logger.success("Cauthy combination done.") 

    pvals = pd.read_csv(os.path.join(save_path, f'{task_id}_Cauchy_pvals.tsv'), index_col=0, sep='\t')
    if use_DEG:
        mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
        complex_func = score.calculate_complex_min_func
        mean_counts = score.combine_complex_distribution_df(mean_counts, complex_table, complex_func)
        mean_pmfs = dist_iid_set.calculate_cluster_mean_distribution(counts_df, labels_df)
        complex_func = dist_complex.get_minimum_distribution
        mean_pmfs = dist_complex.combine_complex_distribution_df(mean_pmfs, complex_table, complex_func)
        pvals = check_key_interactions_pvalue_by_DEG(mean_counts, mean_pmfs, interactions, pvals)
        pvals.to_csv(os.path.join(save_path, f'{task_id}_Cauchy_with_DEG_pvals.tsv'), sep='\t')
        logger.success("DEGs selection done.")
    
    significant_results = create_significant_interactions_df(pvals, database_file_path)
    __save(
        significant_results, 
        save_path = save_path,
        task_id = task_id, 
        method_key = '', 
        timestamp = '',
        suffix = 'significant_results',
        index = False
    )

    file_list = glob.glob(f'{save_path}/{task_id}*interactions_strength.tsv')
    logger.info(f"Integrating {len(file_list)} interactions_strength files.")
    average_strength_df = None
    for item in file_list:
        strength_df = pd.read_csv(item, index_col=0, sep='\t')
        if average_strength_df is None:
            average_strength_df = strength_df.copy()
        else:
            average_strength_df += strength_df.copy()
    average_strength_df = average_strength_df / len(file_list)
    average_strength_df.to_csv(os.path.join(save_path, f'{task_id}_average_interactions_strength.tsv'), sep='\t')
    logger.success("Average CS across all scoring methods calculation done.")






def statistical_analysis_method(
    database_file_path,
    celltype_file_path,
    counts_file_path,
    convert_type = 'hgnc_symbol',
    single_unit_summary = 'Mean',
    complex_aggregation = 'Minimum',
    LR_combination = 'Arithmetic',
    min_percentile = 0.1,
    style = None,
    meta_key=None, 
    select_list=[], 
    filter_= False,
    use_DEG = False,
    save_path = None
):
    # Generate unique task ID.
    task_id = generate_task_id()
    logger.info(f"Task id is {task_id}.")

    # Set save directory.
    if save_path is None:
        save_path = './results/'
        logger.info(f"Results will be saved to \"{save_path}\".")
    mkdir(save_path)

    method_qt_dict = {
        "Median": 0.5,
        "Q2": 0.5,
        "Q3": 0.75,
    }

    single_unit_key = single_unit_summary

    if single_unit_summary.startswith("Quantile_"):
        quantile = float(single_unit_summary.split("Quantile_")[1])
        single_unit_summary = "Quantile"

    if single_unit_summary in ["Median", "Q2", "Q3", "Quantile"]:
        if single_unit_summary in method_qt_dict:
            quantile = method_qt_dict[single_unit_summary]
            single_unit_summary = 'Quantile'
        assert quantile < 1 and quantile > 0, "The quantile should be within the range (0, 1)."
        if quantile < 0.5:
            warnings.warn(f"Invalid quantile: {quantile}. The quantile is too low to be effective.", UserWarning)
        min_percentile = max(min_percentile, 1-quantile)
    
    if style not in {None, "cpdb"}:
        raise ValueError(f"Invalid style: {style}. Must be one of [None, 'cpdb']")
    if style == 'cpdb':
        single_unit_summary = 'Mean'
        complex_aggregation = 'Minimum'
        LR_combination = 'Arithmetic'
        min_percentile = 0.1

    method_key = f'{single_unit_key}_{complex_aggregation}_{LR_combination}'
    
        
    counts_df, labels_df, complex_table, interactions = preprocess.get_input_data(
        database_file_path, 
        celltype_file_path,
        counts_file_path,
        convert_type,
        meta_key = meta_key,
        select_list = select_list,
        filter_ = filter_
    )
    logger.success("Data preprocessing done.")


    logger.info(f"Running:\n-> {single_unit_key} for single-unit summary function.\n"
                    + f"-> {complex_aggregation} for multi-unit complex aggregation.\n"
                    + f"-> {LR_combination} for L-R combination to compute the CS.\n"
                    + f"-> Percentile is {min_percentile}.")
    
    # Stage I : calculate L-R expression score:
    ## Scores for clusters
    if single_unit_summary == 'Mean':
        mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
    elif single_unit_summary == 'Quantile':
        mean_counts = score.calculate_cluster_quantile(counts_df, labels_df, quantile)
    ## Scores for complex
    if complex_aggregation == 'Minimum':
        complex_func = score.calculate_complex_min_func
    elif complex_aggregation == 'Average':
        complex_func = score.calculate_complex_mean_func
    mean_counts = score.combine_complex_distribution_df(mean_counts, complex_table, complex_func)
    ## Scores for L-R expression
    interactions_strength = score.calculate_interactions_strength(mean_counts, interactions, method=LR_combination)
    
    # Stage II : Percentages
    percents = calculate_cluster_percents(counts_df, labels_df, complex_table)
    percents_analysis = analyze_interactions_percents(percents, interactions, threshold=min_percentile)

    # Stage III : Null Distribution
    if single_unit_summary == 'Mean':
        mean_pmfs = dist_iid_set.calculate_cluster_mean_distribution(counts_df, labels_df)
    elif single_unit_summary == 'Quantile':
        mean_pmfs = dist_iid_set.calculate_cluster_quantile_distribution(counts_df, labels_df, quantile)
    if complex_aggregation == 'Minimum':
        complex_func = dist_complex.get_minimum_distribution
    elif complex_aggregation == 'Average':
        complex_func = dist_complex.get_average_distribution
    mean_pmfs = dist_complex.combine_complex_distribution_df(mean_pmfs, complex_table, complex_func)
    
    # Stage IV : P-values
    pvals = dist_lr.calculate_key_interactions_pvalue(
        mean_pmfs, interactions, interactions_strength, percents_analysis, method=LR_combination
    )
    logger.success("FastCCC calculation done.")

    __save_file(interactions_strength, pvals, percents_analysis, save_path, task_id=task_id, method_key=method_key)

    if use_DEG:
        pvals = check_key_interactions_pvalue_by_DEG(mean_counts, mean_pmfs, interactions, pvals)
        pvals.to_csv(os.path.join(save_path, f'{task_id}_DEG_pvals.tsv'), sep='\t')
        logger.success("DEGs selection done.")

    significant_results = create_significant_interactions_df(pvals, database_file_path)
    __save(
        significant_results, 
        save_path = save_path,
        task_id = task_id, 
        method_key = '', 
        timestamp = '',
        suffix = 'significant_results',
        index = False
    )

    return interactions_strength, pvals, percents_analysis

# --- 06/27/2024 @marvinquiet
def cluster_markers_method(counts_file_path, celltype_file_path,
                           cluster_distrib_method = 'Mean',
                           quantile=0.9):
    '''Calculate cluster markers pvalues based on mean counts and null distribution
    input:
        counts_file_path(str): log-normalized count matrix path
        celltype_file_path(str): metadata file path that contains cell labels
        cluster_distrib_method(str): Mean/Quantile, method to calculate cluster distribution
        quantile(float): quantile to calculate the cluster distribution
    output:
        mean_counts(pandas.DataFrame): celltype * gene mean counts matrix.
        mean_pmfs(pandas.DataFrame): celltype * gene distribution matrix.
        clustermarker_pvals(pandas.DataFrame): celltype * gene marker pvalues matrix.
    '''
    counts_df, labels_df = preprocess.get_count_data(counts_file_path, celltype_file_path)
    # Preprocess input parameters
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
    # Cluster Mean Counts
    if cluster_distrib_method == 'Mean':
        mean_counts = score.calculate_cluster_mean(counts_df, labels_df)
    elif cluster_distrib_method == 'Quantile':
        mean_counts = score.calculate_cluster_quantile(counts_df, labels_df, quantile)
    # Null Distribution
    if cluster_distrib_method == 'Mean':
        mean_pmfs = dist_iid_set.calculate_cluster_mean_distribution(counts_df, labels_df)
    elif cluster_distrib_method == 'Quantile':
        mean_pmfs = dist_iid_set.calculate_cluster_quantile_distribution(counts_df, labels_df, quantile)
    # get pvalues
    clustermarker_pvals = pd.DataFrame(
        list(map(get_pvalue_from_pmf, mean_counts.stack(), mean_pmfs.stack())),
        index=mean_counts.stack().index
    ).unstack()
    clustermarker_pvals = clustermarker_pvals.droplevel(0, axis=1) # drop the additional level generated by stack/unstack
    clustermarker_pvals = clustermarker_pvals.fillna(1.0) # fill NaN with 1.0
    return mean_counts, mean_pmfs, clustermarker_pvals
# --- end


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
    mean_counts = counts_df_with_labels.groupby('cell_type', as_index=True, observed=True).mean()

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


def check_key_interactions_pvalue_by_DEG(mean_counts, mean_pmf, interactions, pvals):
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
    
    p1_items = p1.values[np.where(pvals.values < 0.05)]
    p2_items = p2.values[np.where(pvals.values < 0.05)]
    
    mean_count_p1 = mean_counts.loc[p1_index, interactions['multidata_1_id']]
    mean_count_p2 = mean_counts.loc[p2_index, interactions['multidata_2_id']]
    mean_count_p1.columns = interactions.index
    mean_count_p2.columns = interactions.index
    mean_count_p1.index = all_index
    mean_count_p2.index = all_index
    mean_count_p1_item = mean_count_p1.values[np.where(pvals.values < 0.05)]
    mean_count_p2_item = mean_count_p2.values[np.where(pvals.values < 0.05)]
    
    
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
        
    mask = np.float64(np.logical_and(est1, est2))
    masked_values = mask * pvals.values[np.where(pvals.values < 0.05)] + (1 - mask)
    pvals.values[np.where(pvals.values < 0.05)] = masked_values
    return pvals
