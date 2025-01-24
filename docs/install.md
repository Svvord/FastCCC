---
title: Software tutorial
layout: default
nav_order: 2
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# Software tutorial

## Introduction
In this tutorial, we introduce the workflow and core functions of FastCCC to help you get familiar with the [FastCCC] package and start analyzing your own datasets. Please check the [Human CCC reference]({{ site.baseurl }}/cccref.html) sections for examples on how we build the human CCC reference panel and infer CCC on user-collected query dataset.

Please note that we are currently in the submission stage, and the code is still under development and not yet fully packaged. In the future, to enhance user convenience, many functionalities may be restructured with new function names and improvements. If you find that the previous code is no longer applicable, please update it promptly and refer to our latest tutorial. We will finalize a user-friendly version as soon as possible.

## Installation
FastCCC is implemented as a Python (>= 3.11) package. If you wish to use our development version, which allows downloading reference datasets from the [CellxGene] using the FastCCC environment, please note that the current Python version must be lower than 3.13, as `CellxGene-Census` does not support higher versions. FastCCC depends on a few other Python packages that include `numpy`, `pandas`, `scipy`, `scanpy`, `loguru`, `openpyxl`, and `gseapy`. Please refer to the package [pyproject.toml] file for details. For your convenience, you can follow this tutorial to set up and use the environment.

### Method 1: Installing via conda
You can install the environment using Conda by following the steps:
```bash
conda create -n FastCCC python=3.11
conda activate FastCCC
```
Get FastCCC from github:
```bash
git clone https://github.com/Svvord/FastCCC.git
```
Go to the folder `FastCCC` and install:
```bash
cd ./FastCCC
pip install -e .
```

### Method 2: Installing via pip
We are currently organizing the code and packaging functionalities to enhance user convenience. Once the code is finalized, we will upload it to PyPI to support installation via pip install. At this stage, please use the code available on GitHub and install it using Conda or Poetry.
```bash
pip install # coming soon.
```

### Method 3: Installing developing version via Poetry
For developing, we are using the [Poetry] package manager. To install Poetry, follow the instructions [here](https://python-poetry.org/docs/#installing-with-pipx).
```bash
git clone https://github.com/Svvord/FastCCC.git
cd ./FastCCC
poetry install
```

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



[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[Poetry]: https://python-poetry.org/