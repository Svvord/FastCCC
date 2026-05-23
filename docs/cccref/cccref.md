---
title: Human CCC reference
layout: default
nav_order: 4
has_children: true
permalink: cccref
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# FastCCC for Reference-based CCC Analysis

## Download our constructed human CCC reference panel

In the latest version, we leverage FastCCC’s speed to shift part of the reference computation to the first-time usage process. This adds only 1–2 minutes during the initial use of the panel but significantly reduces the storage size of the uploaded reference data. The panel size for each tissue has been reduced from 3GB to just 5MB, retaining only essential information such as gene expression distributions. The human CCC reference panel is now available on [GitHub] and can be downloaded [here](https://github.com/Svvord/FastCCC/tree/main/reference).

You can also download our reference panel using the following commands:
```bash
git clone https://github.com/Svvord/FastCCC.git
cp -r ./FastCCC/reference your/save/path/
```

Currently, we have utilized 16 million cells from 19 different general tissues or tissue locations available on the [CellxGene] platform. The technical platforms, cell counts, and the number of cell types included are shown in the figure below. We plan to generate and upload more reference datasets in the future, allowing users to directly access them without the need to download and process large-scale data, thereby reducing resource consumption.

<p align="center">
  <img src="{{ site.baseurl }}/images/reference.svg" width="700">
</p>



{: .important-title }
> Coming soon
>
> We have uploaded the first version of our constructed human reference panel, based on the CPDBv5.0 LRI database. We are exploring a better approach that will allow users to customize and use any LRI database without needing to create a corresponding reference panel. This feature is coming soon.

{: .highlight-title }
> Tips
>
> If you would like to reproduce our database results or construct your own reference panels, you can follow our tutorials on [downloading reference datasets]({{ site.baseurl }}/cccref.html#download-large-scale-normal-reference-dataset-from-cellxgene) from CellxGene and [constructing a reference panel from scratch]({{ site.baseurl }}/cccref.html#how-to-build-a-human-ccc-reference-panel-from-scratch-using-fastccc). By following these steps, you can recreate the same reference panel and achieve identical results.




## How to perform reference-based CCC analysis on a user-collected query dataset

The process of performing comparative analysis using FastCCC is straightforward. When comparing a query dataset with the reference, users must use the same LRI database that was used to construct the reference (we are currently developing a new approach that will allow users to specify their own LRI database in the next version). Additionally, users should select the appropriate reference tissue based on the tissue origin of the query dataset, specify the path to the query file, and define the location for saving the results.

```python
import fastccc.infer_query
## Modify the file path according to the location where you run the code.
database_file_path = 'FastCCC/db/CPDBv5.0.0/' 
reference_path = 'FastCCC/reference/lung/' # Take "lung" reference panel as an example.
tissue_query_file = 'your/save/path/user_collected_query.h5ad'
save_path = 'your/save/path/user_collected_query/'

fastccc.infer_query.infer_query_workflow(
    database_file_path = database_file_path,
    reference_path = reference_path,
    query_counts_file_path = tissue_query_file,
    celltype_file_path = None,
    save_path = save_path,
    meta_key = 'cell_type' # Use meta_key or celltype_file_path based on your query data.
)
```

To run the same healthy tissue reference workflow and generate the HTML
reference report in one call:

```python
from fastccc.report import generate_reference_report, list_reference_panels

print(list_reference_panels('FastCCC/reference'))

report_path = generate_reference_report(
    database_path          = database_file_path,
    query_counts_file_path = tissue_query_file,
    infer_result_dir       = save_path,
    output_dir             = 'your/save/path/user_collected_query_report/',
    reference_tissue       = 'lung',
    reference_root         = 'FastCCC/reference',
    meta_key               = 'cell_type',
    query_name             = 'Query',
)
```

Pass `reference_path` instead of `reference_tissue` when the baseline is a
custom panel built from a user-curated control dataset.

The reference report is self-contained and includes global trend summaries,
sender-receiver heatmaps, top L-R dotplots, pathway breakdowns, a cell-type
reference explorer, and a figure-generation audit. The cell-type explorer lets
users switch among cell types, sender/receiver roles, L-R pairs, and annotated
pathway profiles.
<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-26 20:54:32</span> | INFO     | Start inferring by using CCC reference: lung
<span class="sr">2025-01-26 20:54:32</span> | INFO     | Reference min_percentile = 0.1
<span class="sr">2025-01-26 20:54:32</span> | INFO     | Reference LRI DB = CPDBv5.0.0
<span class="sr">2025-01-26 20:54:35</span> | INFO     | Reading query adata, (your data) cells x (your data) genes.                               
<span class="sr">2025-01-26 20:54:51</span> | <span class="sb">SUCCESS </span> | <span class="sb">Rank preprocess done.</span>
<span class="sr">2025-01-26 20:54:53</span> | INFO     | Loading LRIs database. hgnc_symbol as gene name is requested.
<span class="sr">2025-01-26 20:54:55</span> | <span class="sb">SUCCESS </span> | <span class="sb">Requested data for fastccc is prepared.</span>
<span class="sr">2025-01-26 20:54:55</span> | INFO     | Loading reference data.
<span class="sr">2025-01-26 20:54:57</span> | INFO     | Reference cell types label will be used directly.
<span class="sr">2025-01-26 20:54:58</span> | <span class="sb">SUCCESS </span> | <span class="sb">Reference data is loaded.</span>
<span class="sr">2025-01-26 20:54:58</span> | INFO     | Calculating CS score for query data.
<span class="sr">2025-01-26 20:55:00</span> | INFO     | Filtering reference data.
<span class="sr">2025-01-26 20:55:01</span> | INFO     | Filtering by using reference.
<span class="sr">2025-01-26 20:55:01</span> | INFO     | Inferring sig. boundaries.
<span class="sr">2025-01-26 20:55:03</span> | INFO     | Saving inference results.
<span class="sr">2025-01-26 20:55:04</span> | <span class="sb">SUCCESS </span> | <span class="sb">Inference workflow done.</span>
</code></pre></div>
</blockquote>

## Download large-scale normal reference dataset from CellxGene

If you want to reproduce our results, you can use the following code as an example to download the raw reference count matrix, which we use to construct the human CCC reference panel. For detailed instructions on using the CellxGene-Census API, see [here](https://chanzuckerberg.github.io/cellxgene-census/).

```python
import cellxgene_census
census = cellxgene_census.open_soma()

filter_condition = "tissue_general == 'lung' "
filter_condition += "and disease == 'normal' "
filter_condition += "and cell_type!='unknown' "

adata = cellxgene_census.get_anndata(
    census=census,
    organism="Homo sapiens",
    obs_value_filter=filter_condition
)

adata.write_h5ad('your/save/path/census_download_lung.h5ad')
```

## How to build a human CCC reference panel from scratch using FastCCC
Here, we use a single-cell dataset of normal lung tissue obtained from [CellxGene] as an example to construct a human reference panel for lung tissue. Users can, of course, use their own curated data to build a personalized CCC reference panel and later compare it with their own query dataset.

First, ensure that the reference dataset is in raw count matrix format. Unlike the conventional FastCCC workflow, to minimize batch effects across different datasets, we use raw count data. The query data should also maintain the raw count matrix format. FastCCC will perform basic filtering, removing cells with insufficient gene counts (`min_genes=50` for `scanpy.pp.filter_cells`). If additional QC steps are required, users should delete the corresponding cells or genes from the raw data after performing their own QC.

Additionally, we do not use highly variable genes (HVGs), so it is recommended to retain all genes or transcripts. A reduced gene count may undermine the objectivity and effectiveness of the quantification. Lastly, ensure that the var_names in your AnnData object correspond to standard HGNC gene symbols.

We processed the lung dataset downloaded from CellxGene using the following code and saved it as the lung reference dataset.

```python
## Reference Dataset Preparation.
import scanpy as sc
tissue_reference_file = 'your/save/path/census_download_lung.h5ad'
reference_adata = sc.read_h5ad(tissue_reference_file)
## Notice:
## Ensure that the var_names are using HGNC gene symbols.
## We need to manually convert the IDs, as CellxGene uses Ensembl ID.
## If the format is already correct, you can ignore the following line of code.
reference_adata.var_names = reference_adata.var.feature_name 
## The following line of code is used to ensure no data leakage in our validation. 
## Users can ignore this.
observation_joinid = sorted(set(reference_adata.obs.observation_joinid))
import pickle
with open(f'your/save/path/lung_reference_joinid.pkl', 'wb') as f:
    pickle.dump(observation_joinid, f)
## Notice ends.
reference_adata.write_h5ad('your/save/path/lung_reference.h5ad')
```
Then, one can build a reference panel like the following code
```python
## Build human lung reference panel.
import fastccc.build_reference
## Modify the file path according to the location where you run the code.
database_file_path = 'FastCCC/db/CPDBv5.0.0/' 
tissue_reference_file = 'your/save/path/lung_reference.h5ad'
save_path = 'your/save/path/reference/'
reference_name = 'lung',

fastccc.build_reference.build_reference_workflow(
    database_file_path = database_file_path,
    reference_counts_file_path = tissue_reference_file,
    celltype_file_path = None,
    reference_name = reference_name,
    save_path = save_path,
    meta_key = 'cell_type'
)
```
<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-01-27 10:33:14</span> | INFO     | Start building CCC reference.
<span class="sr">2025-01-27 10:33:14</span> | INFO     | Reference_name = lung
<span class="sr">2025-01-27 10:33:14</span> | INFO     | min_percentile = 0.1
<span class="sr">2025-01-27 10:33:14</span> | INFO     | LRI database = CPDBv5.0.0
<span class="sr">2025-01-27 10:33:14</span> | <span class="sb">SUCCESS </span> | <span class="sb">Reference save dir your/save/path/reference/lung is created.</span>
<span class="sr">2025-01-27 10:34:03</span> | INFO     | Reading reference adata, 1673947 cells x 60530 genes.                  
<span class="sr">2025-01-27 10:36:55</span> | <span class="sb">SUCCESS </span> | <span class="sb">Rank preprocess done.</span>
<span class="sr">2025-01-27 10:37:23</span> | INFO     | Loading LRIs database. hgnc_symbol as gene name is requested.
<span class="sr">2025-01-27 10:38:08</span> | <span class="sb">SUCCESS </span> | <span class="sb">Requested data for fastccc is prepared.</span>
<span class="sr">2025-01-27 10:38:08</span> | INFO     | Running FastCCC.
<span class="sr">2025-01-27 10:38:36</span> | INFO     | Calculating null distributions.
<span class="sr">2025-01-27 10:39:25</span> | INFO     | Saving reference.
<span class="sr">2025-01-27 10:39:30</span> | INFO     | Saving reference config.
<span class="sr">2025-01-27 10:39:30</span> | <span class="sb">SUCCESS </span> | <span class="sb">Reference 'lung' is built.</span>
</code></pre></div>
</blockquote>

## How to use one dataset with control and case groups as a custom reference workflow

If a user-collected dataset contains both a reference/control group and a
query/case group, the control group can be used to build a custom FastCCC
reference panel. The case group can then be compared against that panel using
reference-based inference. Both subsets should use raw counts.

```python
import scanpy as sc
import fastccc.build_reference
from fastccc.report import generate_reference_report

adata = sc.read_h5ad('your/save/path/cohort_raw_counts.h5ad')

control = adata[adata.obs['condition'] == 'control'].copy()
case = adata[adata.obs['condition'] == 'disease'].copy()

database_file_path = 'FastCCC/db/CPDBv5.0.0/'
reference_root = 'your/save/path/reference/'

fastccc.build_reference.build_reference_workflow(
    database_file_path         = database_file_path,
    reference_counts_file_path = control,
    celltype_file_path         = None,
    reference_name             = 'cohort_control',
    save_path                  = reference_root,
    meta_key                   = 'cell_type'
)

report_path = generate_reference_report(
    database_path          = database_file_path,
    query_counts_file_path = case,
    infer_result_dir       = 'your/save/path/results/disease_vs_cohort_control/',
    output_dir             = 'your/save/path/report/disease_vs_cohort_control/',
    reference_path         = 'your/save/path/reference/cohort_control',
    meta_key               = 'cell_type',
    query_name             = 'Disease',
    reference_name         = 'Cohort Control',
)
```

If a local FastCCC version does not accept an in-memory `AnnData` object in one
of these entrypoints, write `control` and `case` to temporary `.h5ad` files and
pass those file paths instead.


[CellxGene]: https://cellxgene.cziscience.com/
[GitHub]: https://github.com/Svvord/FastCCC
