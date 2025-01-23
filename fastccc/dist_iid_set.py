from .distrib import Distribution, get_distribution_from_samples
from .distrib import get_quantile_pmf_for_n_iid_distribution
import pandas as pd
from collections import Counter


def calculate_cluster_mean_distribution(counts_df, labels_df):
    meta_dict = Counter(labels_df.cell_type)
    
    ####### 
    gene_sum_pmf_dict = {}
    for gene in counts_df.columns:
        samples = counts_df[gene].values
        gene_sum_pmf_dict[gene] = {1: get_distribution_from_samples(samples)}
        for item in range(2,30):
            gene_sum_pmf_dict[gene][item] = gene_sum_pmf_dict[gene][item-1] + gene_sum_pmf_dict[gene][1]
            
    for gene in counts_df.columns:
        for item in range(2,30):
            gene_sum_pmf_dict[gene][item] = gene_sum_pmf_dict[gene][item] / item
            
    ####### clusters_mean #######
    clusters_mean_dict = {}
    for celltype in sorted(meta_dict):
        clusters_mean_dict[celltype]  = {}
        n_sum = meta_dict[celltype]
        if n_sum < 30:
            for gene in counts_df.columns:
                clusters_mean_dict[celltype][gene] = gene_sum_pmf_dict[gene][n_sum]
        else:
            for gene in counts_df.columns:
                clusters_mean_dict[celltype][gene] = gene_sum_pmf_dict[gene][1] ** n_sum / n_sum
    mean_pmf = pd.DataFrame(clusters_mean_dict).T
    return mean_pmf


def calculate_cluster_quantile_distribution(counts_df, labels_df, quantile):
    meta_dict = Counter(labels_df.cell_type)
    
    ####### 
    gene_distribution_dict = {}
    for gene in counts_df.columns:
        samples = counts_df[gene].values
        gene_distribution_dict[gene] = get_distribution_from_samples(samples)
        
    ####### clusters_quantile #######
    clusters_quantile_dict = {}
    for celltype in sorted(meta_dict):
        clusters_quantile_dict[celltype]  = {}
        celltype_count = meta_dict[celltype]
        for gene in counts_df.columns:
            distribution = gene_distribution_dict[gene]
            quantile_pmf = \
            get_quantile_pmf_for_n_iid_distribution(distribution, celltype_count, quantile)
            clusters_quantile_dict[celltype][gene] = Distribution(
                dtype = 'other',
                pmf_array = quantile_pmf,
                is_align = True
            )
    quantile_pmf = pd.DataFrame(clusters_quantile_dict).T
    return quantile_pmf