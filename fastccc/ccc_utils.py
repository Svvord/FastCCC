import numpy as np
import pandas as pd
import os
import psutil


def create_significant_interactions_df(
    pvals, 
    LRI_db_path, 
    save = False, 
    save_path = './temp/', 
    save_file = 'significant_interaction_list', 
    seperator = '|'
):

    from . import preprocess
    interactions = preprocess.get_interactions(LRI_db_path)
    id2symbol_dict = (pd.read_csv(f'{LRI_db_path}/gene_table.csv')[['protein_id', 'hgnc_symbol']]\
    .set_index('protein_id')['hgnc_symbol']).to_dict()
    ##### complex_table ######
    complex_composition = pd.read_csv(os.path.join(LRI_db_path, 'complex_composition_table.csv'))
    complex_table = pd.read_csv(os.path.join(LRI_db_path, 'complex_table.csv'))
    complex_table = complex_table.merge(complex_composition, left_on='complex_multidata_id', right_on='complex_multidata_id')
    foo_dict = complex_table.groupby('complex_multidata_id').apply(lambda x: list(x['protein_multidata_id'].values), include_groups=False).to_dict()
    for key in foo_dict:
        value = ','.join([id2symbol_dict[item] for item in foo_dict[key]])
        id2symbol_dict[key] = value

    x, y = np.where(pvals.values<0.05)
    lines = []
    for subx, suby in zip(x,y):
        celltype_A, celltype_B = pvals.index[subx].split(seperator)
        interaction_ID = pvals.columns[suby]
        lines.append((celltype_A, celltype_B, interaction_ID, pvals.iloc[subx, suby]))

    output_df = pd.DataFrame(lines, columns=['sender_celltype', 'receiver_celltype', 'LRI_ID', 'p-value'])
    output_df = output_df.merge(interactions, left_on='LRI_ID', right_index=True, how='left')

    output_df.multidata_1_id = [id2symbol_dict[ligand_gene_name] for ligand_gene_name in output_df.multidata_1_id.tolist()]
    output_df.multidata_2_id = [id2symbol_dict[receptor_gene_name] for receptor_gene_name in output_df.multidata_2_id.tolist()]

    output_df = output_df[['sender_celltype', 'receiver_celltype', 'LRI_ID', 'multidata_1_id', 'multidata_2_id', 'p-value']]
    output_df = output_df.rename(columns={'multidata_1_id': 'ligand', 'multidata_2_id': 'receptor'})

    save_file = os.path.join(save_path, save_file)
    if save:
        output_df.to_excel(f'{save_file}.xlsx')
    return output_df


def create_significant_interactions_with_flag_df(pvals, significant_flag, interactions, save_path='./temp/', save_file='significant_interaction_list'):
    seperator = '|'
    x, y = np.where(pvals.values<0.05)
    lines = []
    for subx, suby in zip(x,y):
        if not significant_flag.iloc[subx, suby]:
            continue
        celltype_A, celltype_B = pvals.index[subx].split(seperator)
        interaction_ID = pvals.columns[suby]
        lines.append((celltype_A, celltype_B, interaction_ID, pvals.iloc[subx, suby]))
    output_df = pd.DataFrame(lines, columns=['Ligand_celltype', 'Receptor_celltype', 'Interaction_ID', 'P-val'])
    output_df = output_df.merge(interactions, left_on='Interaction_ID', right_index=True, how='left')
    save_file = os.path.join(save_path, save_file)
    output_df.to_excel(f'{save_file}.xlsx')
    return output_df

def get_current_memory():
    """
    获取当前内存占用
    usage:
    current_memory = get_current_memory()
    print("当前内存占用: {:.2f} MB".format(current_memory))
    """
    process = psutil.Process()
    return process.memory_info().rss / (1024 * 1024)  # 转换为MB
