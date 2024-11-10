cpdb_file_path = './db/mydb'
meta_file_path = './data/data_tutorial/metadata.tsv'
counts_file_path = './data/data_tutorial/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

microenvs_file_path = 'data/data_tutorial/microenvironment.tsv'
out_path = 'results/method1'

import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import os
from . import preproc_utils
from .ccc_utils import get_current_memory
import timeit


def get_interactions(cpdb_file_path, select_list=[]):
    interactions = pd.read_csv(os.path.join(cpdb_file_path, 'interaction_table.csv'))[['id_cp_interaction', 'multidata_1_id', 'multidata_2_id']]
    interactions.set_index('id_cp_interaction', inplace=True)
    interactions = interactions[['multidata_1_id', 'multidata_2_id']].drop_duplicates()
    if select_list:
        interactions = interactions.loc[select_list]
    return interactions


def get_input_data(cpdb_file_path, meta_file_path, counts_file_path, convert_type, meta_key=None, select_list=[], filter_=False):

    '''
    input:
        1. cpdb_file_path(str):
            a path of dir which contains all interaction data
            details see temp
        2. coutns_file_path(str):


    temp:
        I1 -> cpdb_file_path/interaction_table.csv -> interactions(pandas.DataFrame):
        =======================================================
        id_cp_interaction  |  multidata_1_id  |  multidata_2_id
        =======================================================
        CPI-CS0A5B6BD7A    |       1364       |       797
        CPI-CS047D9C0D7    |       1532       |       797
        CPI-CS04A56D5BE    |       1364       |       1004
        CPI-CS0F5B070C5    |       1532       |       1004
        =======================================================

        I2 -> counts(adata):
        main_frame is:
        ===============================
              | gene1 | gene2 | gene3
        ===============================
        cell1 |  0.5  |  0.6  |  0.0
        cell2 |  0.0  |  0.2  |  0.0   
        cell3 |  0.5  |  0.0  |  0.4   
        cell4 |  0.0  |  0.0  |  0.1  
        ===============================


    output:

    '''
    


    ######################################################################
    #                          1.Data Input                              #
    ######################################################################

    ### interactions, counts, labels  ##
    interactions = get_interactions(cpdb_file_path, select_list)
    start = timeit.default_timer()
    counts = anndata.read_h5ad(counts_file_path)#.to_df()
    counts.var_names_make_unique()
    stop = timeit.default_timer()
    print('Read Time: ', stop - start) 
    
    #MMMMM
    current_memory = get_current_memory()
    print("reading_count\t{:.2f}MB".format(current_memory))
    #MMMMM
    # 0815 add:
    if filter_:
        sc.pp.filter_cells(counts, min_counts=1)
        sc.pp.filter_genes(counts, min_cells=1)

    
    if meta_key is not None:
        labels_df = pd.DataFrame(counts.obs[meta_key])
        labels_df.columns = ['cell_type']
        labels_df.index.name = 'barcode_sample'
    else:
        labels_df = pd.read_csv(meta_file_path, sep='\t', index_col=0)
    
    '''
    interactions(pandas.DataFrame):
    =======================================================
                      |  multidata_1_id  |  multidata_2_id
    id_cp_interaction |                  |            
    =======================================================
    CPI-CS0A5B6BD7A   |        1364      |       797
    CPI-CS047D9C0D7   |        1532      |       797
    CPI-CS04A56D5BE   |        1364      |      1004
    =======================================================
    '''
    ####################################



    ##### gene_table ########
    gene_table = pd.read_csv(os.path.join(cpdb_file_path, 'gene_table.csv'))
    protein_table = pd.read_csv(os.path.join(cpdb_file_path, 'protein_table.csv'))
    gene_table = gene_table.merge(protein_table, left_on='protein_id', right_on='id_protein')
    #########################


    ##### complex_table ######
    complex_composition = pd.read_csv(os.path.join(cpdb_file_path, 'complex_composition_table.csv'))
    complex_table = pd.read_csv(os.path.join(cpdb_file_path, 'complex_table.csv'))
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


   


    ######################################################################
    #                      2.Counts DF preprocess                        #
    ######################################################################
    
    ##### feature to id conversion  ######
    # 不在 标准列表 里的 gene 就不要了
    tmp = gene_table[[convert_type, 'protein_multidata_id']]
    tmp = tmp.drop_duplicates()
    tmp.set_index('protein_multidata_id', inplace=True)

    select_columns = []
    columns_names = []
    for foo, boo in zip(tmp.index, tmp[convert_type]):
        if boo in counts.var_names:#counts.columns:
            select_columns.append(boo)
            columns_names.append(foo)

    reduced_counts = counts[:, select_columns].to_df()
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

    # def __exist__(key, df):
    #     # 目前的 complex 策略就是 有一个 就可以了
    #     for item in __content__(key):
    #         if item in df.columns:
    #             return True
    #     return False
    
    def __exist__(key, df):
        # 目前的 complex 策略就是 全部都要有
        # 经测试，这是cpdb用的策略
        for item in __content__(key):
            if item not in df.columns:
                return False
        return True
    
    # ###### 检验 __exist__ ####
    # select_index = []
    # for item in complex_table['complex_multidata_id']:
    #     if __exist__(item, reduced_counts):
    #         select_index.append(True)
    #     else:
    #         select_index.append(False)
    # complex_composition_filtered = complex_table[select_index]
    # #########################

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
    
    ###### cpdb verification #######
#     select_cols = []
#     for item in set(temp_list):
#         if item not in foo_dict:
#             select_cols.append(item)
#     counts_simple = reduced_counts[select_cols]

#     select_cols = []
#     for item in set(temp_list):
#         if item in foo_dict:
#             select_cols.extend(__content__(item))
#     # 从这里开始出现了不一致， cpdb先看 complex是否存在 （所有protein组成部分都在），
#     # 然后生成了 counts_complex 部分，然后 和 simplex 拼接，
#     # 必然 出现多了一些没用的部分 
#     counts_complex = reduced_counts[list(set(select_cols))]
    ################################
    
            
    for item in temp_list:
        for subitem in __content__(item):
            if subitem in temp_dict:
                temp_dict[subitem] = True
    select_index = [key for key in temp_dict if temp_dict[key]]
    reduced_counts = reduced_counts[select_index]
    ################################################
    
    #MMMMM
    current_memory = get_current_memory()
    print("get_input_data_peak\t{:.2f}MB".format(current_memory))
    #MMMMM
    
    ######################################################################
    #                               Return                               #
    ######################################################################
    counts_df = reduced_counts
    temp_list = set(temp_list)
    select_index = [True if item in temp_list else False for item in complex_table.complex_multidata_id]
    complex_table = complex_table[select_index]
    interactions = interactions_filtered
    return counts_df, labels_df, complex_table, interactions
    
    
    
    