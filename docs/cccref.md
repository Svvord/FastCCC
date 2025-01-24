---
title: Human CCC reference
layout: default
nav_order: 4
---

# FastCCC for Reference-based CCC Analysis

## Download our constructed human CCC reference panel
We are currently preparing to upload the fully constructed CCC reference panel, which will be made publicly available soon.

In the meantime, if you would like to experiment with our reference panel or you want to construct your own reference panels, you can follow our tutorials on [downloading reference datasets]({{ site.baseurl }}/cccref.html#download-large-scale-normal-reference-dataset-from-cellxgene) from CellxGene and [constructing a reference panel from scratch]({{ site.baseurl }}/cccref.html#how-to-build-a-human-ccc-reference-panel-from-scratch-using-fastccc). By following these steps, you can recreate the same reference panel and achieve identical results.

## How to perform reference-based CCC analysis on a user-collected query dataset

```python
import fastcci.infer_query
import scanpy as sc
## Modify the file path according to the location where you run the code.
database_file_path = 'FastCCC/db/v5.0.0/' 
reference_path = 'your/save/path/reference/lung/'
tissue_query_file = 'your/save/path/user_collected_query.h5ad'
save_path = 'your/save/path/user_collected_query/'

fastcci.infer_query.infer_query_workflow(
    database_file_path = database_file_path,
    reference_path = reference_path,
    query_counts_file_path = tissue_query_file,
    celltype_file_path = None,
    save_path = save_path,
    meta_key = 'cell_type'
)
```

## Download large-scale normal reference dataset from CellxGene
```python
import cellxgene_census
census = cellxgene_census.open_soma()

filter_condition = "tissue_general == 'lung' "
filter_condition += "and disease == 'normal' "
filter_condition += "and is_primary_data == False "
filter_condition += "and tissue=='lung' "
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
import fastcci.build_reference
import scanpy as sc
## Modify the file path according to the location where you run the code.
database_file_path = 'FastCCC/db/v5.0.0/' 
tissue_reference_file = 'your/save/path/lung_reference.h5ad'
save_path = 'your/save/path/reference/'
reference_name = 'lung',

fastcci.build_reference.build_reference_workflow(
    database_file_path = database_file_path,
    reference_counts_file_path = tissue_reference_file,
    celltype_file_path = None,
    reference_name = reference_name,
    save_path = save_path,
    meta_key = 'cell_type'
)
```

## Functions
To be continued.


[CellxGene]: https://cellxgene.cziscience.com/