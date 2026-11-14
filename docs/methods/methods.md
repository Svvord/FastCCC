---
title: Methods
layout: default
nav_order: 5
has_children: true
permalink: methods
description: "The statistics behind FastCCC: analytic null distributions instead of permutations, the exact distribution of order statistics of discrete variables, sixteen scoring modules and their Cauchy combination, reference-based calibration; and the distribution toolkit you can use on its own."
---

# Methods

FastCCC replaces permutation testing with analytic null distributions. Two building blocks carry the method: the distribution of the mean of \\(n\\) cells by FFT convolution, and the exact distribution of order statistics of \\(n\\) draws from a *discrete* distribution, ties included, which is what makes median- and quantile-based summaries testable on expression data. On top of them sit the sixteen scoring modules, their Cauchy combination, and the reference-based comparison with its calibration.

- [How FastCCC works]({{ site.baseurl }}/methods/how-it-works.html): communication scores, the null without permutations, the sixteen modules, reference panels and calibration, effect sizes, bootstrap and replication, extreme p-values.
- [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html): the functions behind it, exact for any discrete variable, with a worked example and code you can run for your own problems.

The full derivations are in the paper: Hou, S., Ma, W. & Zhou, X. FastCCC: a permutation-free framework for scalable, robust, and reference-based cell-cell communication analysis in single cell transcriptomics studies. *Nature Communications* 16, 11428 (2025). [doi:10.1038/s41467-025-66272-z](https://doi.org/10.1038/s41467-025-66272-z).
