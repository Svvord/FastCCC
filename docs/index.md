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
This website maintains FastCCC tutorial and code for reproducing simulation and real data application results described in the upcoming paper. For the software, please refer to [FastCCC].

Please note that we are currently in the submission stage, and the code is still under development and not yet fully packaged. In the future, to enhance user convenience, many functionalities may be restructured with new function names and improvements. If you find that the previous code is no longer applicable, please update it promptly and refer to our latest tutorial. We will finalize a user-friendly version as soon as possible.

## Introduction to FastCCC
FastCCC is a powerful, permutation-free framework designed to identify critical cell-cell communications (CCCs) in single-cell transcriptomics studies. It provides scalable, robust, and reference-based CCC analysis, enabling researchers to uncover meaningful biological insights from large-scale datasets. FastCCC employs innovative techniques like fast Fourier transformation-based convolution for analytical $$p$$-value computation, modular algebraic operations to capture diverse CCC patterns, and reference panels to enhance analysis on user-collected data.

## FastCCC overview

<p align="center">
  <img src="{{ site.baseurl }}/images/overview.v2.0.jpg" width="700">
</p>


## Citing the work
If you find the `FastCCC` package or any of the source code in our [repository] useful for your work, please cite:

> Hou, S., Ma, W., Zhou, X. FastCCC: A permutation-free framework 
> for scalable, robust, and reference-based cell-cell communication analysis 
> in single cell transcriptomics studies.

Visit our [group website](https://xiangzhou.github.io/) for more statistical 
tools on analyzing genetics, genomics and transcriptomics data.

[FastCCC]: https://github.com/Svvord/FastCCC
[repository]: https://github.com/Svvord/FastCCC