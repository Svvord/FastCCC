import pandas as pd
import numpy as np
import itertools
from functools import partial


def calculate_cluster_mean(counts_df, labels_df):
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
    mean_counts = counts_df_with_labels.groupby('cell_type', as_index=True, observed=True).mean()
    return mean_counts


def calculate_cluster_quantile(counts_df, labels_df, qt=0.9):
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
    mean_counts = counts_df_with_labels.groupby('cell_type', observed=True).apply(lambda x: pd.Series(np.quantile(x,qt, axis=0, method='lower'), index=x.columns), include_groups=False)
    if 'cell_type' in mean_counts:
        mean_counts.drop(columns=['cell_type'], inplace=True)
    return mean_counts


def combine_complex_distribution_df(mean_counts, complex_table, complex_count_func):
    func = partial(complex_count_func, mean_counts=mean_counts)
    # complex count dataframe 
    if not complex_table.empty:
        complex_counts = complex_table.groupby('complex_multidata_id').apply(
            lambda x: x['protein_multidata_id'].values, include_groups=False).apply(
            lambda x:func(x)).T
        # 合成 最终 mean_counts dataframe
        mean_counts = pd.concat((mean_counts, complex_counts), axis=1)
    ##############################################################
    #                         Return                             #
    ##############################################################
    mean_counts = mean_counts.dropna(axis=1)
    return mean_counts

def calculate_complex_min_func(x, mean_counts):
    x = [sub_x for sub_x in x if sub_x in mean_counts.columns]
    if len(x) == 0:
        return pd.Series(index=mean_counts.index)
    return mean_counts.loc[:,x].T.min()

def calculate_complex_mean_func(x, mean_counts):
    x = [sub_x for sub_x in x if sub_x in mean_counts.columns]
    if len(x) == 0:
        return pd.Series(index=mean_counts.index)
    return mean_counts.loc[:,x].T.mean()


def calculate_interactions_strength(mean_counts, interactions, method='Arithmetic', sep='|'):
    assert method in ['Arithmetic', 'Geometric'], "Only support 'Arithmetic' or 'Geometric' mean."
    
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
    
    if method == 'Arithmetic':
        interactions_strength = (p1 + p2)/2 * (p1 > 0) * (p2>0)
    elif method == 'Geometric':
        interactions_strength = np.sqrt(p1*p2)
        
    return interactions_strength