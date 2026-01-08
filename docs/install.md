---
title: Installation
layout: default
nav_order: 2
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# FastCCC Deployment

## Introduction
In this section, we introduce how to install FastCCC. For details on using FastCCC, please check the [Basic Usage]({{ site.baseurl }}/usage.html) section for an overview of the core functions and how to run CCC analysis on a single-cell dataset. For examples on how to build the human CCC reference panel and perform CCC analysis on user-collected query datasets, refer to the [Human CCC Reference]({{ site.baseurl }}/cccref.html) section.

Additionally, FastCCC serves as a comprehensive CCC analysis toolkit with various tools that are widely applicable in single-cell RNA-seq studies. Usage of these tool functions can be found in the [Toolkits]({{site.baseurl}}/cccref.html) section.

## Installation
FastCCC is implemented as a [Python] (>= 3.11) package. If you wish to use our development version, which allows downloading reference datasets from the [CellxGene] using the FastCCC environment, please note that the current Python version must be lower than 3.13, as [CellxGene-Census] does not support higher versions. FastCCC depends on a few other Python packages that include `numpy`, `pandas`, `scipy`, `scanpy`, `loguru`, `openpyxl`, and `gseapy`. Please refer to the package [pyproject.toml] file for details. For your convenience, you can choose one of the following methods to set up and use the environment.

### Method 1: Installing via conda
{: .d-inline-block }

recommend
{: .label .label-green }

You can install the environment using [`conda`] by following the steps:
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

{: .highlight-title }
> Tips
>
> We recommend this installation method for general users, as it has been tested across various servers to ensure a smooth setup and execution.


### Method 2: Installing via pip
```bash
pip install fastccc
```

{: .important-title }
> Note
>
> Please ensure that Python version >= 3.11.

### Method 3: Installing developing version via Poetry
For developing, we are using the [Poetry] package manager. To install Poetry, follow the instructions [here](https://python-poetry.org/docs/#installing-with-pipx).

```bash
git clone https://github.com/Svvord/FastCCC.git
cd ./FastCCC
poetry install
```

{: .highlight-title }
> Tips
>
> We recommend creating a dedicated `FastCCC` virtual environment with Python version 3.11 and installing dependencies using [Poetry].


[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[Poetry]: https://python-poetry.org/
[CellxGene-Census]: https://chanzuckerberg.github.io/cellxgene-census/
[`conda`]: https://www.anaconda.com/download
[PyPI]: https://pypi.org/
[Python]: https://www.python.org/