---
title: Code snippets
layout: default
parent: Human CCC reference
nav_order: 3
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# Code Snippets

## Introduction
This section provides multiple templates, starting from the folder and data storage structure, followed by data download, processing, and the various ways to use FastCCC for ***comparative analyses***. It aims to guide users in understanding how to apply our methods and help them reproduce our results effectively.

## Code snippets for comparative analysis using FastCCC

### Example1: using human liver CCC reference panel

Here, we downloaded a liver disease scRNA-seq dataset from [CellxGene](https://cellxgene.cziscience.com/e/02792605-4760-4023-82ad-40fc4458a5db.cxg/) as an example.

#### **Dataset**

##### Download
```bash
# https://cellxgene.cziscience.com/e/02792605-4760-4023-82ad-40fc4458a5db.cxg/
# caudate lobe of liver, 10x 5' v2, Homo sapiens, 
# primary sclerosing cholangitis and primary biliary cholangitis
cd ./data/raw
wget https://datasets.cellxgene.cziscience.com/b09dac61-836f-418a-ad0a-064d2bc5343d.h5ad
```
##### Preprocess
```python
import scanpy as sc

tissue_query_file = '../data/raw/b09dac61-836f-418a-ad0a-064d2bc5343d.h5ad'
disease_liver_adata = sc.read_h5ad(tissue_query_file)
disease_liver_adata.X = disease_liver_adata.raw.X
disease_liver_adata.var_names = disease_liver_adata.var.feature_name.tolist() # Use hgnc_symbol.
disease_liver_adata = disease_liver_adata[disease_liver_adata.obs.disease!='normal']
disease_liver_adata.write_h5ad('../data/clean/liver_query_disease_exp1.h5ad')
```



#### **Directory Structure**
```
┌─ ...
├─ (your files)
├─ FastCCC (download from github)
    ├─ db
        ├─ CPDBv5.0.0 (LRI database provided.)
        └─ ...
    ├─ reference
        ├─ breast (Specified tissue origin.)
        └─ ...
    ├─ ...
    ├─ code_snippets.py (run your code here.)
    └─ ...
├─ data
    ├─ raw
        ├─ b09dac61-836f-418a-ad0a-064d2bc5343d.h5ad
        └─ ...
    ├─ clean
        ├─ liver_query_disease_exp1.h5ad
        └─ ...
    └─ ...
├─ (your files)
└─ ...
```

#### **templates**
```python
# copy following templates to code_snippets.py
import fastccc

database_file_path = './db/CPDBv5.0.0/'
reference_path = './reference/liver/'
tissue_query_file = '../data/clean/liver_query_disease_exp1.h5ad'
save_path = './results/liver_query_disease_exp1/'

fastccc.infer_query.infer_query_workflow(
    database_file_path = database_file_path,
    reference_path = reference_path,
    query_counts_file_path = tissue_query_file,
    celltype_file_path = None,
    save_path = save_path,
    meta_key = 'cell_type',
    debug_mode = True
)
```

<blockquote class="new-title"> <p>output</p>
<div class="highlight"><pre class="highlight"><code><span class="sr">2025-02-02 15:11:50</span> | INFO     | Start inferring by using CCC reference: liver
<span class="sr">2025-02-02 15:11:50</span> | INFO     | Reference min_percentile = 0.1
<span class="sr">2025-02-02 15:11:50</span> | INFO     | Reference LRI DB = CPDBv5.0.0
<span class="sr">2025-02-02 15:11:50</span> | <span class="sb">SUCCESS </span> | <span class="sb">Query save dir is created.</span>
<span class="sr">2025-02-02 15:12:02</span> | INFO     | Reading query adata, 65396 cells x 33363 genes.                                           
<span class="sr">2025-02-02 15:12:18</span> | <span class="sb">SUCCESS </span> | <span class="sb">Rank preprocess done.</span>
<span class="sr">2025-02-02 15:12:38</span> | INFO     | Loading LRIs database. hgnc_symbol as gene name is requested.
<span class="sr">2025-02-02 15:12:41</span> | <span class="sb">SUCCESS </span> | <span class="sb">Requested data for fastccc is prepared.</span>
<span class="sr">2025-02-02 15:12:41</span> | INFO     | Loading reference data.
<span class="sr">2025-02-02 15:12:45</span> | INFO     | Reference cell types label will be used directly.
<span class="sr">2025-02-02 15:12:47</span> | <span class="sb">SUCCESS </span> | <span class="sb">Reference data is loaded.</span>
<span class="sr">2025-02-02 15:12:47</span> | INFO     | Calculating CS score for query data.
<span class="sr">2025-02-02 15:12:49</span> | INFO     | Filtering reference data.
<span class="sr">2025-02-02 15:13:01</span> | INFO     | Filtering by using reference.
<span class="sr">2025-02-02 15:13:01</span> | INFO     | Inferring sig. boundaries.
<span class="sr">2025-02-02 15:13:16</span> | INFO     | Saving inference results.
<span class="sr">2025-02-02 15:13:18</span> | <span class="sb">SUCCESS </span> | <span class="sb">Inference workflow done.</span>
</code></pre></div>
</blockquote>

### Example2: How to adjust the granularity of cell type annotations in reference panel

FastCCC allows manual modification and mapping of cell type annotations in the reference data to accommodate differences in annotation granularity between the reference and query datasets. If users want to align the query dataset’s annotations with the reference, they can simply modify their query data directly, which we will not elaborate on here.

We use a TOML file to record the cell types and their counts in each tissue panel. To check this information, please refer to `reference/tissue/config.toml`.

In this example, we use breast tissue data to demonstrate how to adjust cell type labels within the reference dataset.

#### **templates**
```python
import fastccc

database_file_path = './db/NicheNetv1.1.1/'
reference_path = './reference/breast_nichenetv1.1.1/'
tissue_query_file = '../data/clean/case_study/breast_tumor_atlas.h5ad'
save_path = '../results/case_study/breast_tumor_atlas_nichenetv1.1.1'

celltype_mapping_dict = {
    # cell type in reference : cell type in query
    "endothelial cell of lymphatic vessel": "Endothelial Cells",  
    "naive B cell": "B Cells",  
    "class switched memory B cell": "B Cells",  
    "IgG plasma cell": "Plasma Cells",  
    "IgA plasma cell": "Plasma Cells",  
    "unswitched memory B cell": "B Cells",  
    "CD8-positive, alpha-beta memory T cell": "CD8+ T Cells",  
    "naive thymus-derived CD4-positive, alpha-beta T cell": "CD4+ T Cells",  
    "Tc1 cell": "CD8+ T Cells",  
    "macrophage": "Macrophages",  
    "CD4-positive helper T cell": "CD4+ T Cells",  
    "plasma cell": "Plasma Cells",  
    "effector memory CD8-positive, alpha-beta T cell": "CD8+ T Cells",  
    "natural killer cell": "NK Cells",  
    "dendritic cell": "Dendritic Cells",  
    "classical monocyte": "Monocytes",  
    "myeloid dendritic cell": "Dendritic Cells",  
    "non-classical monocyte": "Monocytes",  
    "conventional dendritic cell": "Dendritic Cells",  
    "alternatively activated macrophage": "Macrophages",  
    "inflammatory macrophage": "Macrophages",  
    "plasmacytoid dendritic cell": "Dendritic Cells",  
    "neutrophil": "Neutrophils",  
    "mast cell": "Mast Cells",  
    "pericyte": "Perivascular-like (PVL) Cells",  
    "vascular associated smooth muscle cell": "Perivascular-like (PVL) Cells",  
    "CD4-positive, alpha-beta T cell": "CD4+ T Cells",  
    "regulatory T cell": "Regulatory T Cells",  
    "activated CD4-positive, alpha-beta T cell": "CD4+ T Cells",  
    "effector memory CD4-positive, alpha-beta T cell": "CD4+ T Cells",  
    "activated CD8-positive, alpha-beta T cell": "CD8+ T Cells",  
    "capillary endothelial cell": "Endothelial Cells",  
    "vein endothelial cell": "Endothelial Cells",  
    "endothelial cell of artery": "Endothelial Cells",  
    "fibroblast": "Fibroblasts",  
    "luminal epithelial cell of mammary gland": "Epithelial Cells",  
    "mammary gland epithelial cell": "Epithelial Cells",  
    "basal cell": "Epithelial Cells",  
    "fibroblast of mammary gland": "Fibroblasts",  
    "endothelial tip cell": "Endothelial Cells",  
    "perivascular cell": "Perivascular-like (PVL) Cells",  
    "luminal hormone-sensing cell of mammary gland": "Epithelial Cells",  
    "myoepithelial cell of mammary gland": "Myoepithelial Cells"
}

fastccc.infer_query.infer_query_workflow(
    database_file_path = database_file_path,
    reference_path = reference_path,
    query_counts_file_path = tissue_query_file,
    celltype_file_path = None,
    save_path = save_path,
    meta_key = 'cell_type',

    ## Use celltype mapping dict!
    celltype_mapping_dict = celltype_mapping_dict,
)
```