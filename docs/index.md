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

## FastCCC toolkits
FastCCC provides three major advances.
1. FastCCC presents a novel, alternative strategy for computing $$p$$-values in CCC analysis. By leveraging convolution techniques through Fast Fourier Transform (FFT), FastCCC derives the analytic solution for the null distribution of $CS$, enabling scalable analytic computation of $p$-values without requiring computationally intensive permutations. This approach makes FastCCC orders of magnitude faster than existing methods, while enhancing accuracy and robustness. The computational gain of FastCCC further increases with larger datasets, making it particularly suited for large-scale, state-of-the-art single cell studies that existing methods cannot handle. Second, FastCCC introduces a modular $CS$ computation framework that captures interaction strength through algebraic operations between ligand and receptor expression level. This framework enables the development of new $CS$ scores beyond common ones by utilizing a wide range of expression summary statistics to capture a broad spectrum of CCC patterns. It also accommodates multi-subunit protein complexes for ligands and receptors, accounting for distinct interaction strength characterized by different subunits. As such, FastCCC substantially enhances the power and robustness of CCC detection across diverse datasets and biological contexts. Finally, FastCCC not only offers a scalable and effective solution for CCC analysis, enabling the analysis of large-scale datasets containing millions of cells, but also, for the first time, allows these large-scale datasets to serve as reference panels, facilitating more comprehensive and context-aware CCC analysis on user collected datasets, which may be much smaller in size. To support routine reference-based CCC analysis with FastCCC, we constructed the first human CCC reference panel with 19 distinct tissue types and approximately sixteen million cells (Methods, Table~\ref{tabA1}). This reference panel enhances the interpretability and consistency of CCC analysis across diverse biological contexts. 
1. FFT-based probabalistic calculation algorithm.
2. 

## Citing the work
If you find the `FastCCC` package or any of the source code in our [repository] useful for your work, please cite:

> Hou, S., Ma, W., Zhou, X. FastCCC: A permutation-free framework 
> for scalable, robust, and reference-based cell-cell communication analysis 
> in single cell transcriptomics studies.

Visit our [group website](https://xiangzhou.github.io/) for more statistical 
tools on analyzing genetics, genomics and transcriptomics data.

Check out our [Zhou Lab](https://xiangzhou.github.io/) website for more softwares :) 


This is a *bare-minimum* template to create a Jekyll site that uses the [Just the Docs] theme. You can easily set the created site to be published on [GitHub Pages] – the [README] file explains how to do that, along with other details.

If [Jekyll] is installed on your computer, you can also build and preview the created site *locally*. This lets you test changes before committing them, and avoids waiting for GitHub Pages.[^1] And you will be able to deploy your local build to a different platform than GitHub Pages.

More specifically, the created site:

- uses a gem-based approach, i.e. uses a `Gemfile` and loads the `just-the-docs` gem
- uses the [GitHub Pages / Actions workflow] to build and publish the site on GitHub Pages

Other than that, you're free to customize sites that you create with this template, however you like. You can easily change the versions of `just-the-docs` and Jekyll it uses, as well as adding further plugins.

[Browse our documentation][Just the Docs] to learn more about how to use this theme.

To get started with creating a site, simply:

1. click "[use this template]" to create a GitHub repository
2. go to Settings > Pages > Build and deployment > Source, and select GitHub Actions

If you want to maintain your docs in the `docs` directory of an existing project repo, see [Hosting your docs from an existing project repo](https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md#hosting-your-docs-from-an-existing-project-repo) in the template README.

----

[^1]: [It can take up to 10 minutes for changes to your site to publish after you push the changes to GitHub](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll/creating-a-github-pages-site-with-jekyll#creating-your-site).

[Just the Docs]: https://just-the-docs.github.io/just-the-docs/
[GitHub Pages]: https://docs.github.com/en/pages
[README]: https://github.com/just-the-docs/just-the-docs-template/blob/main/README.md
[Jekyll]: https://jekyllrb.com
[GitHub Pages / Actions workflow]: https://github.blog/changelog/2022-07-27-github-pages-custom-github-actions-workflows-beta/
[use this template]: https://github.com/just-the-docs/just-the-docs-template/generate
[FastCCC]: https://github.com/Svvord/FastCCC
[repository]: https://github.com/Svvord/FastCCC