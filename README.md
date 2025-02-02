# FastCCC: A permutation-free framework for scalable, robust, and reference-based cell-cell communication analysis in single cell transcriptomics studies.

[![Preprint](https://img.shields.io/badge/preprint-available-brightgreen)](https://www.biorxiv.org/content/10.1101/2025.01.27.635115v1) &nbsp;
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen)](https://svvord.github.io/FastCCC/) &nbsp;
[![License](https://img.shields.io/badge/license-MIT-blue)](https://github.com/Svvord/FastCCC/blob/main/LICENSE)

**[2025.02.01]** Update: To minimize the size of transmitted panel data, we leverage FastCCC’s speed to compute essential reference data during first-time usage. This process incurs only an additional 1–2 minutes during initial activation. Meanwhile, the storage requirement for uploading the panel data has been significantly reduced (from 3GB to 5MB per tissue panel).

**[2025.01.23]** We have provided a comprehensive [tutorial](https://svvord.github.io/FastCCC/) on the usage of FastCCC, which includes detailed instructions on installation, usage, and more. We highly recommend referring to this [tutorial](https://svvord.github.io/FastCCC/) for a step-by-step guide.

## Overview
![scheme](./docs/images/overview.v2.0.jpg)
<p align="justify"> Detecting cell-cell communications (CCCs) in single-cell transcriptomics studies is fundamental for understanding the function of multicellular organisms. Here, we introduce FastCCC, a permutation-free framework that enables scalable, robust, and reference-based analysis for identifying critical CCCs and uncovering biological insights. FastCCC relies on fast Fourier transformation-based convolution to compute $p$-values analytically without permutations, introduces a modular algebraic operation framework to capture a broad spectrum of CCC patterns, and can leverage atlas-scale single cell references to enhance CCC analysis on user-collected datasets. To support routine reference-based CCC analysis, we constructed the first human CCC reference panel, encompassing 19 distinct tissue types, over 450 unique cell types, and approximately 16 million cells. We demonstrate the advantages of FastCCC across multiple datasets, most of which exceed the analytical capabilities of existing CCC methods. In real datasets, FastCCC reliably captures biologically meaningful CCCs, even in highly complex tissue environments, including differential interactions between endothelial and immune cells linked to COVID-19 severity, dynamic communications in thymic tissue during T-cell development, as well as distinct interactions in reference-based CCC analysis.  </p>

## Installation
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

## How to use `FastCCC`
Check our [vignettes](https://svvord.github.io/FastCCC/).

## Citing the work
If you find the `FastCCC` package or any of the source code in this repository useful for your work, please [cite](https://www.biorxiv.org/content/10.1101/2025.01.27.635115v1):

> Siyu Hou, Wenjing Ma, and Xiang Zhou (2025). FastCCC: A permutation-free framework 
> for scalable, robust, and reference-based cell-cell communication analysis 
> in single cell transcriptomics studies.

```
@article {hou2025fastCCC,
	author = {Hou, Siyu and Ma, Wenjing and Zhou, Xiang},
	title = {FastCCC: A permutation-free framework for scalable, robust, and reference-based cell-cell communication analysis in single cell transcriptomics studies},
	year = {2025},
	publisher = {Cold Spring Harbor Laboratory},
	journal = {bioRxiv}
}
```

Visit our [group website](https://xiangzhou.github.io/) for more statistical 
tools on analyzing genetics, genomics and transcriptomics data.