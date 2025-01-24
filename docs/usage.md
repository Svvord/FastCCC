---
title: Basic usage
layout: default
nav_order: 3
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# FastCCC for CCC Analysis

## Introduction
This section provides a tutorial of how to use FastCCC for conducting cell-cell communication (CCC) analyses. We introduce the workflow and core functions of FastCCC to help you get familiar with the [FastCCC] package and start analyzing your own datasets. Please check the [Human CCC reference]({{ site.baseurl }}/cccref.html) sections for examples on how to use human CCC reference panel.

## How FastCCC calculates a $$CS$$ statistic
To help users correctly utilize FastCCC, we first explain how the communication score ($$CS$$) is calculated. The $$CS$$ consists of two components: the ligand score and the receptor score. These scores are derived from the expression levels of the ligand and receptor in the cell types of interest.

$$
\begin{equation}
    {CS}_{c_a,c_b,l,r} = h(s(x_{c_a,l}),s(x_{c_b,r})), \label{eq1}
\end{equation}
$$

where $$s(x_{c_a,l})$$ is a scalar and represents the gene expression summary of ligand $$l$$ across cells of cell type $$c_a$$ ; $$s(x_{c_b,r})$$ is also a scalar and represents the gene expression summary of receptor $$r$$ across cells of cell type $$c_b$$ ; and $$h(\cdot,\cdot)$$ is a function that measures the coordinated expression between $$s(x_{c_a,l})$$ and $$s(x_{c_b,r})$$.

FastCCC supports the use of mean or any $$k^{th}$$ order statistic (e.g., median, 3<sup>rd</sup> quartile, 90<sup>th</sup> percentile) for the $$s(\cdot)$$ function. It also supports arithmetic or geometric mean for the $$h(\cdot)$$ function. For ligand or receptor complexes consisting of multiple genes, the expression summaries for each gene in the complex can be aggregated. Currently, FastCCC supports aggregation methods such as averaging or taking the minimum value.

<p align="center">
  <img src="{{ site.baseurl }}/images/basic_usage_illustration1.svg" width="700">
</p>

Users can choose their preferred method for `single-unit expression summaries`, `multi-unit aggregation`, and `ligand-receptor integration functions`(i.e. $$h(\cdot)$$). Notably, FastCCC can easily perform the same calculation as [CellPhoneDB], making it a fast alternative. A comparison of results with CellPhoneDB is provided [here]({{ site.baseurl }}/sim.html). In the following [section]({{ site.baseurl }}/usage.html#basic-usages-with-single-cs), we will use the same $$CS$$ calculation functions as CellPhoneDB to demonstrate how to use FastCCC for analyzing single-cell transcriptomic data.

Users can also combine multiple approaches. We provide the Cauchy combination test to integrate the results across all branches. For the multi-branch version of FastCCC, please refer [here]({{ site.baseurl }}/usage.html#basic-usages-with-multiple-cs-combination).

Although the calculation rules are straightforward, computing the null distribution of these statistical measures without permutations remains a challenging task. For the specific algorithm, please refer to the Methods section of our paper. The preprint version will be made available shortly.

## Basic usages with single $$CS$$


```python
import fastccc.core as core

cpdb_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../../data/examples/metadata.tsv'
counts_file_path = '../../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

interactions_strength, pvals, percents_analysis = core.statistical_analysis_method(
    cpdb_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type
)
```

## Basic usages with multiple $$CS$$ combination

```python
import fastccc.core as core

cpdb_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../../data/examples/metadata.tsv'
counts_file_path = '../../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

core.Cauchy_combination_of_statistical_analysis_methods(
    cpdb_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type,
    cluster_distrib_method_list = ['Mean', 'Quantile_0.9'],
    complex_distrib_method_list = ['Minimum'],
    LR_distrib_method_list = ['Arithmetic', 'Geometric'],
    save_path = '../../results/temp/',
    use_DEG = True
)
```

## Build your own LRIs database
To be continued



[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[Poetry]: https://python-poetry.org/
[CellPhoneDB]: https://www.nature.com/articles/s41596-020-0292-x