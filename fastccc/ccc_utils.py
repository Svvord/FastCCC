import numpy as np
import pandas as pd
import os
import psutil

def create_significant_interactions_df(pvals, interactions, save=False, save_path='./temp/', save_file='significant_interaction_list'):
    seperator = '|'
    x, y = np.where(pvals.values<0.05)
    lines = []
    for subx, suby in zip(x,y):
        celltype_A, celltype_B = pvals.index[subx].split(seperator)
        interaction_ID = pvals.columns[suby]
        lines.append((celltype_A, celltype_B, interaction_ID, pvals.iloc[subx, suby]))
    output_df = pd.DataFrame(lines, columns=['Ligand_celltype', 'Receptor_celltype', 'Interaction_ID', 'P-val'])
    output_df = output_df.merge(interactions, left_on='Interaction_ID', right_index=True, how='left')
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
