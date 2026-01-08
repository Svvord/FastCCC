---
title: Welcome
layout: home
nav_order: 1
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

<div style="margin: 0 auto; text-align: center;"> 
  <img src="{{ site.baseurl }}/images/logo.png" width="500" />
</div>


## Welcome to the FastCCC website!
This website maintains FastCCC tutorial and code for reproducing simulation and real data application results described in the [FastCCC paper](https://doi.org/10.1038/s41467-025-66272-z). 

[Get started now]({{ site.baseurl }}/install.html){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View it on GitHub][FastCCC]{: .btn .btn-blue .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Full text link][Full text]{: .btn .btn-green .fs-5 .mb-4 .mb-md-0 }

## Introduction to FastCCC
FastCCC is a highly scalable, permutation-free statistical toolkit tailored to identify critical cell-cell communications (CCCs) in the form of ligand-receptor interactions (LRIs) and uncover novel biological insights in single-cell transcriptomics studies. FastCCC presents an analytic solution for computing p-values in CCC analysis, enabling scalable analysis without the need for computationally intensive permutations. It introduces a modular communication score computation framework that calculates various communication scores through a range of algebraic operations between ligand and receptor expression levels, capturing a broad spectrum of CCC patterns and ensuring robust analysis. Additionally, FastCCC not only enables the analysis of large-scale datasets containing millions of cells, but also introduces reference-based CCC analysis, where large-scale datasets are treated as reference panels to substantially improve CCC analysis on user-collected datasets.

## FastCCC overview

<p align="center">
  <img src="{{ site.baseurl }}/images/overview.v3.0.jpg" width="700">
</p>

FastCCC introduces unique features in CCC analyses. 
1. ### Ultra-high scalability:
It offers a novel approach to computing $$p$$-values using convolution techniques through Fast Fourier Transform (FFT), eliminating the need for computationally expensive permutations. This approach makes FastCCC significantly faster and more scalable than existing methods, especially for large datasets. 
2. ### Modular algebraic operation framework:
FastCCC features a modular communication score ($$CS$$) computation framework that captures interaction strength through algebraic operations on ligand-receptor expression levels. This framework can handle multi-subunit protein complexes, increasing the power and robustness of CCC detection.
3. ### Reference-based CCC comparison analysis
FastCCC enables reference-based CCC analysis by allowing large-scale datasets to serve as reference panels for more comprehensive analysis of smaller user-collected datasets.
4. ### Human CCC reference panel:
We constructed, to the best of our knowledge, the first human CCC reference panel, which includes 19 tissue types and approximately 16 million cells, thereby enhancing the interpretability and consistency of CCC analysis across various biological contexts.

## Citing the work
If you find the `FastCCC` package or any of the source code in our [repository] useful for your work, please cite:

> Hou, S., Ma, W. & Zhou, X. FastCCC: a permutation-free framework for scalable, 
> robust, and reference-based cell-cell communication analysis in single cell 
> transcriptomics studies. Nat Commun 16, 11428 (2025). 
> https://doi.org/10.1038/s41467-025-66272-z

Visit our [group website](https://xiangzhou.github.io/) for more statistical 
tools on analyzing genetics, genomics and transcriptomics data.

[FastCCC]: https://github.com/Svvord/FastCCC
[repository]: https://github.com/Svvord/FastCCC
[BioRxiv]: https://www.biorxiv.org/content/10.1101/2025.01.27.635115v1
[Full text]:https://doi.org/10.1038/s41467-025-66272-z