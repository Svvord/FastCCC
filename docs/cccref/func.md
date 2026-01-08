---
title: Functions
layout: default
parent: Human CCC reference
nav_order: 1
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# Functions

- [`build_reference_workflow`]({{site.baseurl}}/cccref/func.html#build_reference_workflow)
- [`infer_query_workflow`]({{site.baseurl}}/cccref/func.html#infer_query_workflow)

## `build_reference_workflow`

### Description

The `build_reference_workflow` function constructs a reference panel for cell-cell communication analysis.
It processes reference count data, quantifies and ranks it, and prepares the necessary inputs for downstream CCC analysis.
Additionally, the function saves the processed reference configuration and relevant files for future use.

{: .note }
>
- Reads and preprocesses reference **raw count** data from an `.h5ad` file.
- Uses **HGNC symbols** as gene names (i.e., `anndata.var_names` are official gene symbols).
- Count data in **CSR sparse format** (i.e. `type(anndata.X) == scipy.sparse.csr_matrix`).
- Extracts interaction information from a given ligand-receptor interaction (LRI) database.
- Configures and stores reference settings for later analyses.

### Function Signature
```python
def build_reference_workflow(
    database_file_path, 
    reference_counts_file_path, 
    celltype_file_path, 
    reference_name, 
    save_path, 
    meta_key=None, 
    min_percentile = 0.1
)
```

### Parameters

| Parameter                  | Type          | Default Value | Description  |
|----------------------------|--------------|---------------|--------------|
| `database_file_path`        | `str`        | None          | Path to the database directory containing the candidate LRIs. |
| `reference_counts_file_path` | `str`       | None          | Path to the reference ***raw count matrix*** file in h5ad format. |
| `celltype_file_path`        | `str`        | None          | Path to the cell type annotation file for reference count file. If the h5ad count file already contains cell type labels, this can be set to `None`, and the `meta_key` parameter should be specified instead. |
| `reference_name`            | `str`        | None          | Name of the reference dataset. A folder with the same name as `reference_name` will be created under the `save_path` to store the reference panel. Please ensure the name is valid for file naming conventions. |
| `save_path`                 | `str`        | None          | Path where the processed reference panel data will be saved. Used together with reference_name, i.e., `save_path/reference_name/`. |
| `meta_key`                  | `str` or `None` | `None`       | Metadata key specifying the column in `adata.obs` that contains the cell type labels. |
| `min_percentile`            | `float`      | `0.1`         | Minimum percentile threshold for filtering interactions. The same parameter will be used during inference, and it is recommended to keep the default value. |




### Returns
This function does not return values directly but generates and saves multiple output files in the specified `save_path`. Users can ignore these details—once the reference panel is built, it can be easily utilized through `infer_query_workflow`.

---

## `infer_query_workflow`

### Description

The `infer_query_workflow` function performs query inference using a pre-built cell-cell communication reference.  
It processes query count data, applies quality control, aligns metadata, and compares the query dataset with the reference to infer cell interactions.  
This function enables researchers to analyze new datasets in the context of a predefined reference.

{: .note }
>
- Reads and preprocesses query **raw count** data from an `.h5ad` file.
- Uses **HGNC symbols** as gene names (i.e., `anndata.var_names` are official gene symbols).
- Count data in **CSR sparse format** (i.e. `type(anndata.X) == scipy.sparse.csr_matrix`).
- The LRI database must be the same as the one used in the reference panel. (We will introduce a new feature to remove this restriction.)

### Function Signature
```python
def infer_query_workflow(
    database_file_path, 
    reference_path, 
    query_counts_file_path, 
    celltype_file_path, 
    save_path, 
    celltype_mapping_dict=None, 
    meta_key=None
)
```


### Parameters

| Parameter              | Type          | Default Value | Description  |
|------------------------|--------------|---------------|--------------|
| `database_file_path`   | `str`        | None          | Path to the database directory containing the candidate LRIs. It should be the same as the one used in `build_reference_workflow` |
| `reference_path`       | `str`        | None          | Path to the pre-built CCC reference dataset. i.e. `your/save/path/reference_name/` |
| `query_counts_file_path` | `str`      | None          | Path to the user-provided query ***raw count matrix*** file in h5ad format. |
| `celltype_file_path`   | `str`        | None          | Path to the cell type annotation file for query count file. If the h5ad count file already contains cell type labels, this can be set to `None`, and the `meta_key` parameter should be specified instead. |
| `save_path`           | `str`        | None          | Path where the inference results will be saved. |
| `celltype_mapping_dict` | `dict` or `None` | `None`     | Dictionary for mapping reference cell types to query cell types. For example, this can be used to merge more granular cell subtype categories in the reference into broader categories, ensuring consistency with the cell type annotations in the query dataset. If `None`, the cell type annotations stored in the reference will be used directly. See [examples]({{site.baseurl}}/cccref/snippet.html#example2-how-to-adjust-the-granularity-of-cell-type-annotations-in-reference-panel) for details. |
| `meta_key`            | `str` or None | `None`       | Metadata key specifying the column in `adata.obs` that contains the cell type labels. |

### Returns
This function does not return values directly but generates and saves multiple output files in the specified `save_path`. These include:

1. `query_infer_results.txt`, results comparing query interactions against the reference dataset. This is a dataframe file with tab-separated values.

2. `query_percents_analysis.txt`, results of the percent analysis in the user-provided query dataset. This is a dataframe file with tab-separated values.

3. `interactions_strength.txt`, the $$CS$$ of each candidate LRI for each sender-receiver cell type. This is a dataframe file with tab-separated values.

---

##### Version Information
- Author: Siyu Hou
- Version: early access
- Last Updated: 2025-02-05