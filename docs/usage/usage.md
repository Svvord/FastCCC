---
title: Basic usage
layout: default
nav_order: 3
has_children: true
permalink: usage
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# FastCCC for CCC Analysis

## Introduction
This section provides a tutorial of how to use FastCCC for conducting cell-cell communication (CCC) analyses. We introduce the workflow and core functions of FastCCC to help you get familiar with the [FastCCC] package and start analyzing your own datasets. Please check the [Human CCC reference]({{site.baseurl}}/cccref.html) sections for examples on how to use human CCC reference panel.

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

Users can choose their preferred method for `single-unit expression summary`, `multi-unit aggregation`, and `ligand-receptor integration function`(i.e. $$h(\cdot)$$). Notably, FastCCC can easily perform the same calculation as [CellPhoneDB], making it a fast alternative. A comparison of results with CellPhoneDB is provided [here]({{site.baseurl}}/sim.html). In the following [section]({{ site.baseurl }}/usage.html#basic-usages-with-single-cs), we will use the same $$CS$$ calculation functions as CellPhoneDB to demonstrate how to use FastCCC for analyzing single-cell transcriptomic data.

Users can also combine multiple approaches. We provide the Cauchy combination test to integrate the results across all branches. For the multi-branch version of FastCCC, please refer [here]({{ site.baseurl }}/usage.html#basic-usages-with-multiple-cs-combination).

Although the calculation rules are straightforward, computing the null distribution of these statistical measures without permutations remains a challenging task. For the specific algorithm, please refer to the Methods section of our paper. The preprint version will be made available shortly.

## Basic usages with a single $$CS$$ statistic
In this section, we demonstrate how to leverage the algebraic computation framework provided by FastCCC to easily build a fast alternative to CellPhoneDB. For the single-unit expression summary, we choose the mean function; for multi-unit aggregation, we use the minimum function; and for the $$h(\cdot)$$, we adopt the arithmetic mean of the ligand and receptor summaries.

We use a sample dataset as an example, which can be downloaded from [here](https://github.com/ventolab/CellphoneDB/blob/v5.0.1/notebooks/data_tutorial.zip). Extract the files, and use only `normalised_log_counts.h5ad` and `metadata.tsv`. As the file names suggest, **FastCCC requires the input to be a normalized log1p-transformed count matrix**. The metadata file contains two columns: the first column, `barcode_sample`, represents the cell ID, and the second column, `cell_type`, corresponds to the cell type labels.

```python
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/' # use the LRI database provided on github
meta_file_path = '../../data/example/metadata.tsv'
counts_file_path = '../../data/example/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol' # supports 'hgnc_symbol' or 'ensembl'

interactions_strength, pvals, percents_analysis = fastccc.statistical_analysis_method(
    LRI_db_file_path, 
    meta_file_path,
    counts_file_path, # format: normalized log1p transformed count matrix
    convert_type,
    single_unit_summary = 'Mean',
    complex_aggregation = 'Minimum',
    LR_combination = 'Arithmetic'
)
```
Here, we use `statistical_analysis_method` for a single $$CS$$ statistic. 

To run the code, you need to specify the following key parameters:
- `LRI_db_file_path` refers to the location of the curated LRI database. We have uploaded several commonly used databases from various sources to our [GitHub repository] for convenient use. If you prefer to use your own curated database, detailed [instructions]({{site.baseurl}}/toolkits.html#build-your-own-lris-database) are provided on how to easily convert your data into a FastCCC-compatible format using the provided tools.

- `meta_file_path` specifies the location of the metadata file containing cell type annotations. The file should be in TSV format (***tab***-separated text file), and you can refer to the provided `metadata.tsv` as an example. Since the H5AD format natively supports storing cell type labels, if your cell type information is already stored in `adata.obs['your_key']`, you can set `meta_file_path = None` and use `meta_key = 'your_key'` to avoid the need for a separate metadata file.

- `counts_file_path` refers to the normalized log1p transformed count matrix. Please do not pre-filter highly variable genes (HVGs), as this may lead to a decrease in the number of candidate LRIs. Typically, after transformation, the maximum value of any gene in a single cell should not exceed 10. Our code sets a default upper limit of 14, meaning that the raw count before transformation should satisfy $$\text{raw count} < e^{14} \approx 1,202,604$$. This threshold accommodates both CPM normalization and [Seurat] style total count normalization (e.g., total counts = $$1 \times 10^5$$). Therefore, please ensure that your data has been processed correctly.

- `convert_type` specifies the format of `adata.var_names` used for gene identifiers. FastCCC supports both `hgnc_symbol`, which are the standard official gene names, and `ensembl` (Ensembl IDs). Ensure that the gene names in your dataset match one of these formats to guarantee compatibility with FastCCC’s analysis pipeline.

- `single_unit_summary` specifies the statistical metric used to compute the expression summary of single-unit ligands or receptors within each cell type cluster. FastCCC supports various statistical options, including `Mean`, `Median` or `Q2`, kth order statistics (e.g., `Q3` for 3<sup>rd</sup> quartile, `Quantile_0.9` for 90<sup>th</sup> percentile, `Quantile_x` for $$(x\times100)$$<sup>th</sup> percentile), allowing users to select the most appropriate measure based on their dataset and analysis goals. Due to the sparsity of single-cell data, we recommend selecting a relatively high percentile, such as `Quantile_0.75` or higher. This is because a significant proportion of the values are zeros, and choosing a lower percentile may not provide meaningful insights.

- `complex_aggregation` specifies how the previously selected single-unit expression summaries are aggregated to calculate summaries for multi-unit protein complexes. FastCCC currently supports two aggregation methods: `Minimum` and `Average`.

- `LR_combination` sepcifies how to combine the ligand and receptor scores to calculate the $$CS$$ statistic. FastCCC currently supports two combination methods: `Arithmetic` and `Geometric` average.

For more details, please refer to [this section]({{site.baseurl}}/usage.html#statistical_analysis_method). For more code examples, please refer to [Code snippets]({{site.baseurl}}/usage/snippet.html) section

## Basic usages for combining multiple $$CS$$ statistics.
Users may want to apply different $$CS$$ scoring methods and combine the final $$p$$-value results. Here, we use the Cauchy combination test to integrate the results from all branches.

<p align="center">
  <img src="{{ site.baseurl }}/images/cct.svg" width="400">
</p>

For example, if you choose `Mean` and `Quantile_0.9` as single-unit expression summaries, use `Minimum` aggregation for multiple-units, and apply two different $$h(\cdot)$$, you can also optionally use DEG filtering (for usage demonstration), though this step is not mandatory and can be skipped according to your preference. The final results will be saved in a table format at `../../results/temp/` in this example. Other parameters are the same as those in the [`statistical_analysis_method`]({{site.baseurl}}/usage.html#statistical_analysis_method) function.

```python
import fastccc

LRI_db_file_path = './db/CPDBv4.1.0/'
meta_file_path = '../../data/examples/metadata.tsv'
counts_file_path = '../../data/examples/normalised_log_counts.h5ad'
convert_type = 'hgnc_symbol'
save_path = '../../results/temp/'

fastccc.Cauchy_combination_of_statistical_analysis_methods(
    LRI_db_file_path, 
    meta_file_path,
    counts_file_path,
    convert_type,
    single_unit_summary_list = ['Mean', 'Quantile_0.9'],
    complex_aggregation_list = ['Minimum'],
    LR_combination_list = ['Arithmetic', 'Geometric'],
    save_path = '../../results/temp/',
    use_DEG = True
)
```

## Functions

### `statistical_analysis_method`

#### **Description**
The `fastccc.core.statistical_analysis_method` function performs statistical analysis on cell-cell communication, supporting various distribution calculation methods. It offers options for filtering candidate LRIs and constructing communication scores, allowing users to customize analysis parameters and optionally save the results.

#### **Function Signature**
```python
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
    meta_key = None, 
    select_list = [], 
    filter_ = False,
    use_DEG = False,
    save_path = None
)
```
#### **Parameters**

| Parameter               | Type            | Default Value       | Description  |
|-------------------------|-----------------|--------------------|--------------|
| `database_file_path`     | `str`            |                | Path to the database directory containing the candidate LRIs. |
| `celltype_file_path`     | `str`            |                | Path to the cell type annotation file. If the h5ad count file already contains cell type labels, this can be set to `None`, and the `meta_key` parameter should be specified instead. |
| `counts_file_path`       | `str`            |                | Path to the ***normalized log1p-transformed matrix*** file in h5ad format. |
| `convert_type`           | `str`            | `'hgnc_symbol'`     | Type of gene identifier used in your data, such as `'hgnc_symbol'` or `'ensembl'`. |
| `single_unit_summary` | `str`            | `'Mean'`            | Method for calculating single-unit expression summaries, options include `'Mean'`, `'Median'` or `'Q2'`, `'Q3'`, `'Quantile_x'`, etc. |
| `complex_aggregation` | `str`            | `'Minimum'`         | Method for calculating multi-unit complex summaries, options include `'Minimum'`, `'Average'`. |
| `LR_combination`      | `str`            | `'Arithmetic'`      | Method for combining ligand and receptor score to calculate $$CS$$, options include `'Arithmetic'` and `'Geometric'` average. |
| `min_percentile`         | `float`          | `0.1`               | Minimum non-zero expression percentile threshold for filtering genes within each cell type cluster. If `'Quantile_x'` is selected, percentile will be set as `max(min_percentile, 1-x)`. |
| `meta_key`               | `str` or `None`  | `None`              | Metadata key specifying the column in `adata.obs` that contains the cell type labels. |
| `select_list`            | `list`           | `[]`                | List of specific set of LRIs to select for analysis. |
| `filter_`                | `bool`           | `False`             | Whether to enable data filtering. If set to `True`, empty cells and genes will be removed. |
| `use_DEG`                | `bool`           | `False`             | Whether to use differentially expressed genes (DEGs) for further filtering. |
| `save_path`              | `str` or `None`  | `None`              | Path to save the analysis results; if `None`, results are saved to the default path. |


#### **Returns**

| Return Value            | Type                | Description  |
|-------------------------|---------------------|--------------|
| `interactions_strength` | `pandas.DataFrame`   | A dataframe containing the calculated $$CS$$ between different LRIs between sender and receiver cell types. |
| `pvals`                 | `pandas.DataFrame`   | A dataframe with $$p$$-values indicating the statistical significance of the interactions. |
| `percents_analysis`     | `pandas.DataFrame`   | A dataframe summarizing the percentage anaylsis of the interactions. |

All these results will be also saved to the `save_path` folder.

---

### `Cauchy_combination_of_statistical_analysis_methods`

#### **Function Signature**
```python
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
)
```

#### **Description**

The `fastccc.core.Cauchy_combination_of_statistical_analysis_methods` function performs an advanced statistical analysis by combining multiple single-unit summary, complex aggregation, and ligand-receptor integration methods. It processes biological data to assess CCC using different statistical approaches and combines the results using the Cauchy combination method. The function offers flexibility by allowing users to specify various distribution methods and filtering options, making it suitable for comprehensive CCC analysis in scRNA-seq datasets.

#### **Parameters**

| Parameter                    | Type          | Default Value                                     | Description  |
|------------------------------|---------------|--------------------------------------------------|--------------|
| `database_file_path`          | `str`         |                          | Path to the database directory containing the candidate LRIs. |
| `celltype_file_path`          | `str`         |                          | Path to the cell type annotation file. If the h5ad count file already contains cell type labels, this can be set to `None`, and the `meta_key` parameter should be specified instead. |
| `counts_file_path`            | `str`         |                          | Path to the normalized log1p-transformed matrix file in h5ad format. |
| `convert_type`                | `str`         | `'hgnc_symbol'`                                  | Type of gene identifier used in your data, such as `'hgnc_symbol'` or `'ensembl'`. |
| `single_unit_summary_list` | `list[str]`    | `['Mean', 'Median', 'Q3', 'Quantile_0.9']`       | List of methods for calculating single-unit summaries. |
| `complex_aggregation_list` | `list[str]`    | `['Minimum', 'Average']`                         | List of methods for calculating multi-unit complex summaries. |
| `LR_combination_list`      | `list[str]`    | `['Arithmetic', 'Geometric']`                    | List of methods for ligand-receptor combination functions. |
| `min_percentile`              | `float`       | `0.1`                                            | Minimum non-zero expression percentile threshold for filtering genes within each cell type cluster. For each scoring method, if `'Quantile_x'` is selected, percentile will be set as `max(min_percentile, 1-x)`. |
| `save_path`                   | `str` or `None` | `None`                                           | Path to save analysis results; if `None`, results are saved to the default path. |
| `meta_key`                    | `str` or `None` | `None`                                           | Metadata key specifying the column in `adata.obs` that contains the cell type labels. |
| `select_list`            | `list`           | `[]`                | List of specific set of LRIs to select for analysis. |
| `filter_`                     | `bool`        | `False`                                          | Whether to enable data filtering. If set to `True`, empty cells and genes will be removed. |
| `use_DEG`                     | `bool`        | `False`                                          | Whether to use differentially expressed genes (DEGs) in the analysis. |


#### **Returns**

No variables will be returned directly. Instead, all scoring-specific results, along with the final combined results, will be saved to the user-specified folder.

---

### Version Information
- Author: Siyu Hou
- Version: early access
- Last Updated: 2025-01-25


[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[Poetry]: https://python-poetry.org/
[CellPhoneDB]: https://www.nature.com/articles/s41596-020-0292-x
[GitHub repository]: https://github.com/Svvord/FastCCC/tree/main/db
[Seurat]: https://satijalab.org/seurat/