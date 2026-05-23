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
- [`list_reference_panels`]({{site.baseurl}}/cccref/func.html#list_reference_panels)
- [`generate_infer_report`]({{site.baseurl}}/cccref/func.html#generate_infer_report)
- [`generate_reference_report`]({{site.baseurl}}/cccref/func.html#generate_reference_report)

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

1. `query_infer_results.tsv`, results comparing query interactions against the reference dataset. This is a dataframe file with tab-separated values.

2. `query_percents_analysis.tsv`, results of the percent analysis in the user-provided query dataset. This is a dataframe file with tab-separated values.

3. `query_interactions_strength.tsv`, the $$CS$$ of each candidate LRI for each sender-receiver cell type. This is a dataframe file with tab-separated values.

---

## `list_reference_panels`

### Description

`list_reference_panels` lists healthy tissue reference panels available under a
downloaded FastCCC reference root. Use it before selecting `reference_tissue` in
`generate_reference_report`.

### Function Signature
```python
from fastccc.report import list_reference_panels

list_reference_panels(reference_root=None)
```

### Parameters

| Parameter | Type | Default Value | Description |
|-----------|------|---------------|-------------|
| `reference_root` | `str` or `None` | `None` | Folder containing FastCCC reference panel subfolders. If omitted, FastCCC checks `./reference` and the source checkout reference folder. |

### Returns

A list of tissue panel names whose folders contain `config.toml`.

---

## `generate_infer_report`

### Description

`generate_infer_report` renders an HTML report from existing reference-based
inference outputs. Use this when `query_infer_results.tsv` already exists and
you do not want to rerun inference.

### Function Signature
```python
from fastccc.report import generate_infer_report

generate_infer_report(
    infer_result_dir,
    database_path,
    output_dir=None,
    query_name="Query",
    reference_name=None,
    reference_path=None,
    reference_tissue=None,
    reference_root=None,
    dpi=150,
    save_individual_figures=True,
    top_n_lr=25,
    top_n_celltypes=20,
)
```

### Parameters

| Parameter | Type | Default Value | Description |
|-----------|------|---------------|-------------|
| `infer_result_dir` | `str` | required | Directory containing `query_infer_results.tsv` and optionally `query_interactions_strength.tsv`. |
| `database_path` | `str` | required | LRI database directory used by the inference results. |
| `output_dir` | `str` or `None` | `None` | Where to write `infer_report.html`; defaults to `<infer_result_dir>/infer_report/`. |
| `query_name` | `str` | `"Query"` | Query label shown in the report. |
| `reference_name` | `str` or `None` | `None` | Reference label shown in the report. If omitted with `reference_tissue`, the report derives a healthy tissue label. |
| `reference_path` | `str` or `None` | `None` | Optional custom reference panel path for report metadata. Mutually exclusive with `reference_tissue`. |
| `reference_tissue` | `str` or `None` | `None` | Optional healthy tissue panel name under `reference_root`. |
| `reference_root` | `str` or `None` | `None` | Folder containing healthy tissue reference panels. |
| `dpi` | `int` | `150` | Resolution for saved figure PNGs. |
| `save_individual_figures` | `bool` | `True` | Whether to save each figure as a PNG in addition to embedding it in HTML. |
| `top_n_lr` | `int` | `25` | Number of top L-R pairs shown in dotplots. |
| `top_n_celltypes` | `int` | `20` | Number of top cell types shown in heatmaps and summaries. |

---

## `generate_reference_report`

### Description

`generate_reference_report` runs reference inference and renders the HTML
reference report in one call. Users can choose a FastCCC healthy tissue panel by
name, for example `reference_tissue='liver'`, or pass `reference_path` for a
custom reference panel.

### Function Signature
```python
from fastccc.report import generate_reference_report

generate_reference_report(
    database_path,
    query_counts_file_path,
    infer_result_dir,
    celltype_file_path=None,
    output_dir=None,
    reference_path=None,
    reference_tissue=None,
    reference_root=None,
    celltype_mapping_dict=None,
    meta_key=None,
    query_name="Query",
    reference_name=None,
    dpi=150,
    save_individual_figures=True,
    top_n_lr=25,
    top_n_celltypes=20,
    debug_mode=False,
)
```

### Parameters

| Parameter | Type | Default Value | Description |
|-----------|------|---------------|-------------|
| `database_path` | `str` | required | LRI database directory. It must match the database used to build the selected reference panel. |
| `query_counts_file_path` | `str` or `AnnData` | required | Query raw-count dataset. Reference workflows perform rank preprocessing internally. |
| `infer_result_dir` | `str` | required | Directory where reference inference TSV outputs are saved. |
| `celltype_file_path` | `str` or `None` | `None` | Optional query cell-type metadata file. Use `meta_key` when labels are already in `.obs`. |
| `output_dir` | `str` or `None` | `None` | Where to write `infer_report.html`. |
| `reference_path` | `str` or `None` | `None` | Custom FastCCC reference panel directory. Mutually exclusive with `reference_tissue`. |
| `reference_tissue` | `str` or `None` | `None` | Healthy tissue panel name under `reference_root`, for example `liver`. |
| `reference_root` | `str` or `None` | `None` | Folder containing healthy tissue panel folders. |
| `celltype_mapping_dict` | `dict` or `None` | `None` | Mapping from reference cell type name to query cell type name when granularities differ. |
| `meta_key` | `str` or `None` | `None` | Query `.obs` column containing cell type labels. |
| `query_name` | `str` | `"Query"` | Query label shown in the report. |
| `reference_name` | `str` or `None` | `None` | Reference label shown in the report. |
| `dpi` | `int` | `150` | Resolution for saved figure PNGs. |
| `save_individual_figures` | `bool` | `True` | Whether to save each figure as a PNG in addition to embedding it in HTML. |
| `top_n_lr` | `int` | `25` | Number of top L-R pairs shown in dotplots. |
| `top_n_celltypes` | `int` | `20` | Number of top cell types shown in heatmaps and summaries. |
| `debug_mode` | `bool` | `False` | Passed to reference inference for debugging output. |

### Healthy tissue panel example
```python
from fastccc.report import generate_reference_report

report_path = generate_reference_report(
    database_path          = './db/CPDBv5.0.0',
    query_counts_file_path = '../data/clean/liver_query_disease_exp1.h5ad',
    infer_result_dir       = './results/PBC_vs_healthy_liver',
    output_dir             = './report/PBC_vs_healthy_liver',
    reference_tissue       = 'liver',
    reference_root         = './reference',
    meta_key               = 'cell_type',
    query_name             = 'PBC',
)
```

Use `fastccc.report.list_reference_panels('./reference')` to list tissue panel
names available under a downloaded reference root. The selected panel must use
the same LRI database as `database_path`.

### Custom reference panel example

Use `reference_path` when the baseline is a user-built control panel.

```python
from fastccc.report import generate_reference_report

report_path = generate_reference_report(
    database_path          = './db/CPDBv5.0.0',
    query_counts_file_path = './data/disease_raw_counts.h5ad',
    infer_result_dir       = './results/disease_vs_control_reference',
    output_dir             = './report/disease_vs_control_reference',
    reference_path         = './reference/my_control_panel',
    meta_key               = 'cell_type',
    query_name             = 'Disease',
    reference_name         = 'Control',
)
```

### Report contents

The reference report includes global trend summaries, sender-receiver heatmaps,
top L-R dotplots, communication-score distributions, pathway breakdowns, a
cell-type reference explorer, and a figure-generation audit. The cell-type
explorer lets users switch among cell types, sender/receiver roles, L-R pairs,
and annotated pathway profiles.

---

##### Version Information
- Author: Siyu Hou
- Version: early access
- Last Updated: 2026-05-22
