---
title: Functions
layout: default
parent: Basic usage
nav_order: 1
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>



# Functions

- [`statistical_analysis_method`]({{base.site}}/usage/func.html#statistical_analysis_method)
- [`infer_query_workflow`]({{base.site}}/cccref/func.html#infer_query_workflow)


## `statistical_analysis_method`

### Description
The `fastccc.core.statistical_analysis_method` function performs statistical analysis on cell-cell communication, supporting various distribution calculation methods. It offers options for filtering candidate LRIs and constructing communication scores, allowing users to customize analysis parameters and optionally save the results.

### Function Signature

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
### Parameters

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


### Returns

| Return Value            | Type                | Description  |
|-------------------------|---------------------|--------------|
| `interactions_strength` | `pandas.DataFrame`   | A dataframe containing the calculated $$CS$$ between different LRIs between sender and receiver cell types. |
| `pvals`                 | `pandas.DataFrame`   | A dataframe with $$p$$-values indicating the statistical significance of the interactions. |
| `percents_analysis`     | `pandas.DataFrame`   | A dataframe summarizing the percentage anaylsis of the interactions. |

All these results will be also saved to the `save_path` folder.

---

## `Cauchy_combination_of_statistical_analysis_methods`

### Function Signature

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

### Description

The `fastccc.core.Cauchy_combination_of_statistical_analysis_methods` function performs an advanced statistical analysis by combining multiple single-unit summary, complex aggregation, and ligand-receptor integration methods. It processes biological data to assess CCC using different statistical approaches and combines the results using the Cauchy combination method. The function offers flexibility by allowing users to specify various distribution methods and filtering options, making it suitable for comprehensive CCC analysis in scRNA-seq datasets.

### Parameters

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


### Returns

No variables will be returned directly. Instead, all scoring-specific results, along with the final combined results, will be saved to the user-specified folder.

---

## Version Information
- Author: Siyu Hou
- Version: early access
- Last Updated: 2025-02-25


[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[Poetry]: https://python-poetry.org/
[CellPhoneDB]: https://www.nature.com/articles/s41596-020-0292-x
[GitHub repository]: https://github.com/Svvord/FastCCC/tree/main/db
[Seurat]: https://satijalab.org/seurat/