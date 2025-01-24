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

Users can choose their preferred method for `single-unit expression summary`, `multi-unit aggregation`, and `ligand-receptor integration function`(i.e. $$h(\cdot)$$). Notably, FastCCC can easily perform the same calculation as [CellPhoneDB], making it a fast alternative. A comparison of results with CellPhoneDB is provided [here]({{ site.baseurl }}/sim.html). In the following [section]({{ site.baseurl }}/usage.html#basic-usages-with-single-cs), we will use the same $$CS$$ calculation functions as CellPhoneDB to demonstrate how to use FastCCC for analyzing single-cell transcriptomic data.

Users can also combine multiple approaches. We provide the Cauchy combination test to integrate the results across all branches. For the multi-branch version of FastCCC, please refer [here]({{ site.baseurl }}/usage.html#basic-usages-with-multiple-cs-combination).

Although the calculation rules are straightforward, computing the null distribution of these statistical measures without permutations remains a challenging task. For the specific algorithm, please refer to the Methods section of our paper. The preprint version will be made available shortly.

## Basic usages with a single $$CS$$ statistic
In this section, we demonstrate how to leverage the algebraic computation framework provided by FastCCC to easily build a fast alternative to CellPhoneDB. For the single-unit expression summary, we choose the mean function; for multi-unit aggregation, we use the minimum function; and for the $$h(\cdot)$$, we adopt the arithmetic mean of the ligand and receptor summaries.

We use a sample dataset as an example, which can be downloaded from [here](https://github.com/ventolab/CellphoneDB/blob/v5.0.1/notebooks/data_tutorial.zip). Extract the files, and use only `normalised_log_counts.h5ad` and `metadata.tsv`. As the file names suggest, **FastCCC requires the input to be a normalized log1p-transformed count matrix**. The metadata file contains two columns: the first column, `barcode_sample`, represents the cell ID, and the second column, `cell_type`, corresponds to the cell type labels.

```python
import fastccc.core as core

cpdb_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../../data/examples/metadata.tsv'
counts_file_path = '../../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'

interactions_strength, pvals, percents_analysis = core.statistical_analysis_method(
    cpdb_file_path, # format: normalized log1p transformed count matrix
    meta_file_path, 
    counts_file_path,
    convert_type,
    cluster_distrib_method = 'Mean',
    complex_distrib_method = 'Minimum',
    LR_distrib_method = 'Arithmetic'
)
```
Here, we use `core.statistical_analysis_method` for a single $$CS$$ statistic. 

To run the code, you need to specify the following key parameters:
- `cpdb_file_path` refers to the location of the curated LRI database. We have uploaded several commonly used databases from various sources to our [GitHub repository] for convenient use. If you prefer to use your own curated database, detailed [instructions]({{ site.baseurl }}/usage.html#build-your-own-lris-database) are provided on how to easily convert your data into a FastCCC-compatible format using the provided tools.

- `meta_file_path` specifies the location of the metadata file containing cell type annotations. The file should be in TSV format (tab-separated text file), and you can refer to the provided `metadata.tsv` as an example. Since the H5AD format natively supports storing cell type labels, if your cell type information is already stored in `adata.obs[key]`, you can set `meta_file_path = None` and use `meta_key = key` to avoid the need for a separate metadata file.

- `counts_file_path` refers to the normalized log1p transformed count matrix. Please do not pre-filter highly variable genes (HVGs), as this may lead to a decrease in the number of candidate LRIs. Typically, after transformation, the maximum value of any gene in a single cell should not exceed 10. Our code sets a default upper limit of 14, meaning that the raw count before transformation should satisfy $$\text{raw count} < e^{14} \approx 1,202,604$$. This threshold accommodates both CPM normalization and [Seurat] style total count normalization (e.g., total counts = $$1 \times 10^5$$). Therefore, please ensure that your data has been processed correctly.

- `convert_type` specifies the format of `adata.var_names` used for gene identifiers. FastCCC supports both `hgnc_symbol`, which are the standard official gene names, and `ensembl` (Ensembl IDs). Ensure that the gene names in your dataset match one of these formats to guarantee compatibility with FastCCC’s analysis pipeline.

- `cluster_distrib_method` specifies the statistical metric used to construct the null distribution for cell clusters in the single-unit expression summary. FastCCC supports various statistical options, including `Mean`, `Median` or `Q2`, kth order statistics (e.g., `Q3` for 3<sup>rd</sup> quartile, `Quantile_0.9` for 90<sup>th</sup> percentile, `Quantile_x` for $$(x\times100)$$<sup>th</sup> percentile), allowing users to select the most appropriate measure based on their dataset and analysis goals. Due to the sparsity of single-cell data, we recommend selecting a relatively high percentile, such as `Quantile_0.8` or higher. This is because a significant proportion of the values are zeros, and choosing a lower percentile may not provide meaningful insights.

- `complex_distrib_method` specifies how the previously selected single-unit expression summaries are aggregated to construct the null distribution for multi-unit protein complexes. FastCCC currently supports two aggregation methods: `Minimum` and `Average`.

- `LR_distrib_method` sepcifies how to combine the ligand and receptor scores to construct the null distribution for the CS statistic. FastCCC currently supports two aggregation methods: `Arithmetic` and `Geometric` average.

For more details, please refer to [this section]({{ site.baseurl }}/usage.html#statistical_analysis_method). 

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

## Functions

### `statistical_analysis_method`

#### Description
The `statistical_analysis_method` function performs statistical analysis on biological data, supporting various distribution calculation methods. It provides data transformation and filtering options, allowing users to customize analysis parameters and optionally save the results.

#### Function Signature
```python
def statistical_analysis_method(
    database_file_path,
    celltype_file_path,
    counts_file_path,
    convert_type='hgnc_symbol',
    cluster_distrib_method='Mean',
    complex_distrib_method='Minimum',
    LR_distrib_method='Arithmetic',
    quantile=0.9,
    min_percentile=0.1,
    style=None,
    meta_key=None, 
    select_list=[], 
    filter_=False,
    use_DEG=False,
    save_path=None
)
```
#### Parameters

| Parameter               | Type            | Default Value       | Description  |
|-------------------------|-----------------|--------------------|--------------|
| `database_file_path`     | `str`            |                | Path to the database directory containing the candidate LRIs. |
| `celltype_file_path`     | `str`            |                | Path to the cell type information file. |
| `counts_file_path`       | `str`            |                | Path to the expression matrix file. |
| `convert_type`           | `str`            | `'hgnc_symbol'`     | Type of gene identifier conversion, such as `'hgnc_symbol'` or `'ensembl'`. |
| `cluster_distrib_method` | `str`            | `'Mean'`            | Method for calculating cluster distribution, options include `'Mean'`, `'Median'`, etc. |
| `complex_distrib_method` | `str`            | `'Minimum'`         | Method for calculating complex distribution, options include `'Minimum'`, `'Maximum'`, etc. |
| `LR_distrib_method`      | `str`            | `'Arithmetic'`      | Method for ligand-receptor distribution calculation, options include `'Arithmetic'`, `'Geometric'`, etc. |
| `quantile`               | `float`          | `0.9`               | Quantile threshold for data processing. |
| `min_percentile`         | `float`          | `0.1`               | Minimum percentile threshold for filtering. |
| `style`                  | `str` or `None`  | `None`              | Optional style for data analysis output. |
| `meta_key`               | `str` or `None`  | `None`              | Metadata key used for specific filtering purposes. |
| `select_list`            | `list`           | `[]`                | List of specific items to select for analysis. |
| `filter_`                | `bool`           | `False`             | Whether to enable data filtering. |
| `use_DEG`                | `bool`           | `False`             | Whether to use differentially expressed genes (DEGs) in the analysis. |
| `save_path`              | `str` or `None`  | `None`              | Path to save the analysis results; if `None`, results will not be saved. |


#### Returns

| Return Value            | Type                | Description  |
|-------------------------|---------------------|--------------|
| `interactions_strength` | `pandas.DataFrame`   | A dataframe containing the calculated $$CS$$ between different LRIs between sender and receiver cell types. |
| `pvals`                 | `pandas.DataFrame`   | A dataframe with $$p$$-values indicating the statistical significance of the interactions. |
| `percents_analysis`     | `pandas.DataFrame`   | A dataframe summarizing the percentage distribution of the analyzed data. |



Example Usage

```python
result = statistical_analysis_method(
    database_file_path='data/database.csv',
    celltype_file_path='data/celltypes.csv',
    counts_file_path='data/counts.csv',
    convert_type='ensembl',
    cluster_distrib_method='Median',
    complex_distrib_method='Maximum',
    LR_distrib_method='Geometric',
    quantile=0.85,
    min_percentile=0.05,
    style='detailed',
    meta_key='sample_id',
    select_list=['gene1', 'gene2'],
    filter_=True,
    use_DEG=True,
    save_path='output/results.csv'
)
```
Notes
	•	Ensure that input file paths are correct to avoid data loading errors.
	•	If filter_ is set to True, appropriate filtering conditions should be provided.
	•	If save_path is not specified, results will not be saved to a file.

### Version Information
- Author: Siyu Hou
- Version: early access
- Last Updated: 2025-01-22

## Build your own LRIs database
To be continued



[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[Poetry]: https://python-poetry.org/
[CellPhoneDB]: https://www.nature.com/articles/s41596-020-0292-x
[GitHub repository]: https://github.com/Svvord/FastCCC/tree/main/db
[Seurat]: https://satijalab.org/seurat/