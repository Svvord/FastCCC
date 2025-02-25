---
title: Inputs
layout: default
parent: Human CCC reference
nav_order: 2
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# FastCCC's Inputs for Reference-based Analyses

## Introduction

The functions `infer_query_workflow` requires three input datasets:

1. **ScRNA-seq data** (raw counts are stored in **compressed sparse row** (`CSR`) format within an `.h5ad` file.).
2. **Cell type annotations** (provided as a tab-separated table or as a column in obsm of the .h5ad file).
3. **Ligand-receptor interaction (LRI) database** ([precompiled datasets](https://github.com/Svvord/FastCCC/tree/main/db) are available for direct use).
4. **Constructed reference panel** ([precompiled reference](https://github.com/Svvord/FastCCC/tree/main/reference) are available for direct use).

For additional parameters, refer to the [Function] section. For specific usage examples, refer to the [Basic Usage] section and [Code Snippets].

## 1. ScRNA-seq data

An `.h5ad` file readable by [Scanpy], containing raw scRNA-seq counts stored in **[compressed sparse row](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html)** (`CSR`) format. To compare with reference panel, FastCCC requires this format and will preprocess the query dataset using the same workflow as for the reference data.

Please avoid gene filtering, as FastCCC quantifies and ranks the expression levels of all genes in each cell. In this context, no additional preprocessing of scRNA-seq data is necessary.

Ensure that **HGNC symbols** (official gene symbols) are used as gene names to maintain consistency with the reference panel. (In future versions, we plan to introduce a gene ID conversion feature for easier usage—stay tuned.)

The internal matrix structure of the scRNA-seq data appears as follows:

<p align="center">
  <img src="{{site.baseurl}}/images/cccref_usage_inputs1.png" width="700">
</p>

[Here]({{site.baseurl}}/cccref/snippet.html#preprocess) is an example of how to construct a suitable input if your data does not meet the required format.

## 2. Cell type annotations

For ease of use, cell type annotations can be provided in two formats.

1. A TSV (tab-separated) file containing two columns—one for cell IDs (matching those in the .h5ad file) and another for cell type labels. (see [example]({{site.baseurl}}/usage/input.html#2-cell-type-annotations))
2. Directly using the cell type annotations stored in the .h5ad file’s obs field (if available). In this case, specify the corresponding column name using the meta_key parameter. (see [example]({{site.baseurl}}/usage/input.html#2-cell-type-annotations))

**FastCCC can analyze cell types that are not present in the reference panel**. However, to better leverage the accumulated information in the reference dataset, if a cell type in the query dataset corresponds to an existing cell type in the reference but has a different string label, you can use `celltype_mapping_dict` to map them accordingly.

```python
celltype_mapping_dict = {
    # cell type in reference : cell type in query
    "endothelial cell of lymphatic vessel": "Endothelial Cells",  
    "naive B cell": "B Cells",  
    "class switched memory B cell": "B Cells",  
    "IgG plasma cell": "Plasma Cells",  
    "IgA plasma cell": "Plasma Cells"
}
```

The cell type information contained in the reference panel can be found in the file:

```bash
your_path/reference/tissue/config.toml
```

This file provides an overview of the cell types included in the reference panel, helping users align their query dataset accordingly.

<p align="center">
  <img src="{{site.baseurl}}/images/cccref_usage_inputs2.png" width="700">
</p>


Refer to the [Functions]({{site.baseurl}}/cccref/#infer_query_workflow) and [Example2]({{site.baseurl}}/cccref/snippet.html#example2-how-to-adjust-the-granularity-of-cell-type-annotations-in-reference-panel) sections for more details.

## 3. Ligand-receptor interaction database

We provide preprocessed database on [GitHub] for users to download directly.

Refer to the [Basic usage - Inputs]({{site.baseurl}}/usage/input.html#3-ligand-receptor-interaction-database) section for more details.

## 4. Constructed reference panel

We provide preprocessed database on [GitHub] for users to download directly. A complete reference database includes the following components:
<p align="left">
  <img src="{{site.baseurl}}/images/cccref_usage_inputs3.png" width="250">
</p>

To build your own reference panel, see [here]({{site.baseurl}}/cccref/#how-to-build-a-human-ccc-reference-panel-from-scratch-using-fastccc).

For additional parameters, refer to the [Function] section. 


[Scanpy]: https://scanpy.readthedocs.io/en/stable/
[GitHub]: https://github.com/Svvord/FastCCC/tree/main/db
[Function]: {{site.baseurl}}/cccref/#infer_query_workflow
[Basic Usage]: {{site.baseurl}}/usage/#basic-usages-with-a-single-cs-statistic
[Code Snippets]: {{site.baseurl}}/usage/snippet.html