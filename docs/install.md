---
title: Installation
layout: default
nav_order: 2
description: Install FastCCC with pip, conda or uv; what ships in the package, where reference panels are stored, and how to work offline.
---

# FastCCC Deployment

## Introduction
In this section, we introduce how to install FastCCC. For details on using FastCCC, please check the [Quick start]({{ site.baseurl }}/quickstart.html) for the three commands that take an `.h5ad` file to an interactive result, and the [Guides]({{ site.baseurl }}/guides.html) for the standard analysis, conditions and samples, reference comparison, ligand–receptor databases and the result explorer. The command-line and Python interfaces are documented in the [Reference]({{ site.baseurl }}/reference.html) section.

## Installation
FastCCC is implemented as a [Python] (>= 3.11) package and is distributed on [PyPI]. It depends on a few other Python packages that include `numpy`, `pandas`, `scipy`, `scanpy`, `loguru`, `typer` and `gseapy`; please refer to the package [pyproject.toml] file for details. If you wish to use the development version to build reference panels from [CellxGene], please note that the Python version must be lower than 3.13, as [CellxGene-Census] does not support higher versions. For your convenience, you can choose one of the following methods to set up and use the environment.

### Method 1: Installing via pip
{: .d-inline-block }

recommend
{: .label .label-green }

FastCCC is on [PyPI]; one command installs it into your current environment:
```bash
pip install fastccc
```
Check the installation:
```bash
fastccc --version
```

{: .highlight-title }
> Tips
>
> We recommend this installation method for general users. Python 3.11, 3.12 and 3.13 are supported; everything needed for a standard analysis, including the ligand–receptor databases, ships inside the package.


### Method 2: Installing via conda
If you prefer a dedicated environment, create one with [`conda`] and install FastCCC from PyPI inside it:
```bash
conda create -n FastCCC python=3.11
conda activate FastCCC
pip install fastccc
```

{: .important-title }
> Note
>
> Please ensure that Python version >= 3.11. Inside a conda environment, install FastCCC with `pip install fastccc`.


### Method 3: Installing via uv
With the [uv] package manager, install into the current environment or as a standalone command:
```bash
uv pip install fastccc
uv tool install fastccc
```
For development, clone the repository and let uv create the environment. To install uv, follow the
instructions [here](https://docs.astral.sh/uv/getting-started/installation/).
```bash
git clone https://github.com/Svvord/FastCCC.git
cd ./FastCCC
uv sync --group test   # creates .venv with FastCCC and pytest
uv run pytest          # unit tests; add --integration for the full workflows
```

{: .highlight-title }
> Tips
>
> We recommend creating a dedicated `FastCCC` virtual environment with Python version 3.11 and installing dependencies using [uv]. `uv sync --group dev` adds the notebook and atlas-building dependencies.


## What is in the package
{: .d-inline-block }

no extra downloads
{: .label .label-blue }

Everything needed for a standard analysis ships inside the wheel, so nothing has to be cloned:

- **22 ligand–receptor databases**: CellPhoneDB v5 (default), CellPhoneDB v4, CellChatDB, NicheNet v1 and v2, the LIANA+ consensus and its 15 source resources, and a mouse consensus. `fastccc db list` prints them; see the [Databases]({{ site.baseurl }}/guides/databases.html) guide.
- **The result explorer**: D3 is bundled, and the `explorer.html` it writes needs no server and no internet.
- **The agent skill** for Claude Code, Codex and other coding agents (`fastccc skill install`).

{: .note-title }
> Reference panels
>
> The 19 healthy tissue panels are **not** in the package. The first `fastccc ref data.h5ad -r liver` downloads the liver panel (3–5 MB) from GitHub into `~/.fastccc/reference/liver` and activates it (computes its null distributions, under a minute). Set `FASTCCC_HOME` to move that folder, or pass `--reference-root` to use panels stored elsewhere, for example a checkout of the repository's `reference/` directory.


## Working offline
- Download the panels you need once with network access: `fastccc panels -d liver` (or `-d all`), or copy `~/.fastccc/reference` to the offline machine.
- The explorer's pathway-enrichment libraries (Hallmark, KEGG, GO) are fetched from Enrichr the first time `explorer.html` is built and cached in `~/.fastccc/genesets`; without network the panel is simply omitted (`fastccc explore --no-genesets` skips it explicitly).
- Everything else, including the explorer page itself, works without internet.


## Checking the installation
```bash
fastccc --help
fastccc db list            # the 22 ligand–receptor databases
fastccc panels             # the 19 reference panels and which are downloaded
fastccc skill show | head  # the agent skill shipped with the package
```

{: .highlight-title }
> Tips
>
> If `fastccc` is not found after `pip install`, the environment's `bin` directory is not on your `PATH`; `python -m fastccc.cli --help` always works.


[FastCCC]: https://github.com/Svvord/FastCCC
[CellxGene]: https://cellxgene.cziscience.com/
[pyproject.toml]: https://github.com/Svvord/FastCCC/blob/main/pyproject.toml
[uv]: https://docs.astral.sh/uv/
[CellxGene-Census]: https://chanzuckerberg.github.io/cellxgene-census/
[`conda`]: https://www.anaconda.com/download
[PyPI]: https://pypi.org/project/fastccc/
[Python]: https://www.python.org/
