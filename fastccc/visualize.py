import pandas as pd
from . import ccc_utils
from . import preprocess
from collections import Counter
import numpy as np

def extract_unique_identifier(significant_interactions_df):
    unis = []
    for c1, c2, iid in significant_interactions_df[['Ligand_celltype','Receptor_celltype','Interaction_ID']].values:
        unis.append('|'.join([c1,c2,iid]))
    return set(unis)


# def get_net_for_number_of_interactions_old_version(pval_df_path, interaction_database):
#     pval_df = pd.read_csv(pval_df_path, index_col=0)
#     interactions = preprocess.get_interactions(interaction_database)
#     significant_interactions = ccc_utils.create_significant_interactions_df(pval_df, interactions, save=False)
#     cci_set = extract_unique_identifier(significant_interactions)

#     counter = Counter(['|'.join(item.split("|")[:2]) for item in sorted(cci_set)])
#     result = []
#     for key in counter:
#         result.append(key.split('|') + [counter[key]])
#     edge_df = pd.DataFrame(result, columns=['from', 'to', 'weight'])

#     node_size_dict = {}
#     for key in counter:
#         c1, c2 = key.split('|')
#         score = counter[key]
#         if c1 not in node_size_dict:
#             node_size_dict[c1] = 0
#         if c2 not in node_size_dict:
#             node_size_dict[c2] = 0
#         node_size_dict[c1] += score
#         node_size_dict[c2] += score
#     result = []
#     for key in sorted(node_size_dict):
#         result.append((key, node_size_dict[key]))
#     node_df = pd.DataFrame(result, columns=['id', 'all.size'])
#     return node_df, edge_df


def get_net_for_number_of_interactions(pval_df_path):
    pval_df = pd.read_csv(pval_df_path, index_col=0)
    weight = list((pval_df.values < 0.05).sum(axis=1).flatten())
    from_list, to_list = [], []
    for item in pval_df.index:
        item = item.split('|')
        from_list.append(item[0])
        to_list.append(item[1])
    result = np.c_[from_list, to_list, weight]
    edge_df = pd.DataFrame(result, columns=['from', 'to', 'weight'])
    edge_df['weight'] = np.int64(edge_df['weight'])
    edge_df = edge_df[edge_df.weight > 0]

    node_size_dict = {}
    for i in range(edge_df.shape[0]):
        c1 = edge_df.iloc[i,0]
        c2 = edge_df.iloc[i,1]
        score = edge_df.iloc[i,2]
        if c1 not in node_size_dict:
            node_size_dict[c1] = 0
        if c2 not in node_size_dict:
            node_size_dict[c2] = 0
        node_size_dict[c1] += score
        if c1 != c2:
            node_size_dict[c2] += score
    result = []
    for key in sorted(node_size_dict):
        result.append((key, node_size_dict[key]))
    node_df = pd.DataFrame(result, columns=['id', 'all.size'])
    return node_df, edge_df

def get_net_for_interaction_strength(pval_df_path, interaction_strength_path):
    interaction_strength = pd.read_csv(interaction_strength_path, index_col=0)
    pval_df = pd.read_csv(pval_df_path, index_col=0)
    weight = list(((pval_df.values < 0.05) * interaction_strength.values).sum(axis=1).flatten())
    from_list, to_list = [], []
    for item in interaction_strength.index:
        item = item.split('|')
        from_list.append(item[0])
        to_list.append(item[1])
    result = np.c_[from_list, to_list, weight]
    edge_df = pd.DataFrame(result, columns=['from', 'to', 'weight'])
    edge_df['weight'] = np.float64(edge_df['weight'])
    edge_df = edge_df[edge_df.weight > 0]

    node_size_dict = {}
    for i in range(edge_df.shape[0]):
        c1 = edge_df.iloc[i,0]
        c2 = edge_df.iloc[i,1]
        score = edge_df.iloc[i,2]
        if c1 not in node_size_dict:
            node_size_dict[c1] = 0
        if c2 not in node_size_dict:
            node_size_dict[c2] = 0
        node_size_dict[c1] += score
        if c1 != c2:
            node_size_dict[c2] += score
    result = []
    for key in sorted(node_size_dict):
        result.append((key, node_size_dict[key]))
    node_df = pd.DataFrame(result, columns=['id', 'all.size'])
    return node_df, edge_df
