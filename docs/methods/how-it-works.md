---
title: How FastCCC works
layout: default
parent: Methods
nav_order: 1
description: "The statistics behind FastCCC: how a communication score is built, how its null distribution is derived without permutations (FFT convolution and exact discrete order statistics), the sixteen scoring modules and their Cauchy combination, reference panels and calibration, effect sizes, bootstrap and replication."
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# How FastCCC works

This page explains the method with the figures of the paper; the full derivations and the benchmarks are in Hou, Ma & Zhou, *Nature Communications* 16, 11428 (2025). The functions that implement the statistics can be used on their own, see the [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html).

<p align="center">
  <img src="{{ site.baseurl }}/images/figure1.png" width="860" alt="FastCCC overview: null distributions by convolution, the construction of a communication score from single-unit summaries, complex aggregation and ligand-receptor combination, the Cauchy combination, and the reference workflow">
</p>
<p align="center" style="font-size:0.85rem;color:#666;margin-top:-0.4rem">
Figure 1 of the paper. <b>a</b> The null distribution of a communication score is derived directly from the expression distributions by convolution, without permutations. <b>b</b> A score is built in three steps, each with two choices; the sixteen resulting p-values are combined with the Cauchy combination test. <b>c–e</b> Reference panels are built from large healthy datasets and a query dataset is compared against them.
</p>

## How a communication score is calculated

The communication score (\\(CS\\)) of a ligand \\(l\\) sent by cell type \\(c_a\\) to a receptor \\(r\\) on cell type \\(c_b\\) has two components, a ligand score and a receptor score, each a summary of the gene's expression across the cells of that type:

$$
\begin{equation}
    {CS}_{c_a,c_b,l,r} = h\bigl(s(x_{c_a,l}),\, s(x_{c_b,r})\bigr),
\end{equation}
$$

where \\(s(x_{c_a,l})\\) is the expression summary of ligand \\(l\\) across the cells of \\(c_a\\), \\(s(x_{c_b,r})\\) that of receptor \\(r\\) across the cells of \\(c_b\\), and \\(h(\cdot,\cdot)\\) measures their coordinated expression. Each step offers two kinds of choice:

<p align="center">
  <img src="{{ site.baseurl }}/images/basic_usage_illustration1.svg" width="760" alt="Single-unit summary (mean or k-th order statistic), multi-unit aggregation (average or minimum), communication score (arithmetic or geometric mean)">
</p>

| Step | Choices | Option names |
|---|---|---|
| single-unit expression summary \\(s(\cdot)\\) | the mean, or a \\(k^{th}\\) order statistic: median, 3rd quartile, 90th percentile, any quantile | `Mean`, `Median` (`Q2`), `Q3`, `Quantile_0.9`, `Quantile_<q>` |
| multi-unit aggregation for complexes | the average of the subunit summaries, or their minimum | `Average`, `Minimum` |
| ligand–receptor combination \\(h(\cdot,\cdot)\\) | arithmetic mean \\(\tfrac12(S_l + S_r)\\) (zero if either side is zero), or geometric mean \\(\sqrt{S_l S_r}\\) | `Arithmetic`, `Geometric` |

With `Mean`, `Minimum` and `Arithmetic` the score is exactly CellPhoneDB's, so FastCCC is a drop-in, permutation-free replacement for it. An interaction is tested only when the ligand is expressed in at least 10% of the sender cells and the receptor in at least 10% of the receiver cells (`--min-percentile`), the same filter CellPhoneDB applies.

## The null without permutations

Permutation methods shuffle the cell-type labels hundreds of times and recompute the score to estimate how often random labelling reaches the observed value. FastCCC derives that distribution directly (Fig. 1a). Under random labelling the \\(n_c\\) cells of a cell type are \\(n_c\\) draws from the gene's expression distribution over all cells, so the null of each summary follows from that per-cell distribution:

- **Mean.** The mean of \\(n_c\\) draws has the distribution of their sum rescaled, and the sum is the \\(n_c\\)-fold convolution of the per-cell distribution, \\(f_g^{n_c} = f_g * f_g * \cdots * f_g\\) (eq. 11 of the paper). FastCCC discretises the log-normalised expression on a grid of width 0.01 and evaluates the convolution with the fast Fourier transform; for cell types of more than 30 cells the central limit theorem gives the same result in closed form, \\(\mathcal N(\mu_g, \sigma_g^2/n_c)\\).
- **Order statistics** (median, quartile, percentile). The textbook density of the \\(k^{th}\\) order statistic assumes a continuous variable in which no two draws coincide. Expression values are discrete, with many cells at zero, so ties are the rule and that formula does not apply. FastCCC maps the discrete distribution onto the order statistic of uniform variables through its cumulative distribution function and reads off the exact distribution, for any \\(n_c\\) and any discrete distribution (eqs. 12–18); the construction and a worked example are on the [Distribution toolkit]({{ site.baseurl }}/methods/distributions.html) page.
- **Complexes.** The average of \\(n_g\\) subunit summaries is again a convolution (eq. 19); the minimum has cdf \\(1 - \prod_g (1 - F_{c,g}(x))\\) (eqs. 20–21).
- **Ligand–receptor combination.** The arithmetic score is the convolution of the two sides, the geometric score the distribution of \\(\sqrt{S_l S_r}\\), both computed on the same grid.

The observed score is compared with the resulting null and the p-value is its upper tail, computed from the tail itself so that it stays exact far below \\(10^{-16}\\). Every (sender, receiver, ligand–receptor) triple gets an exact null at a cost that does not depend on a number of permutations or on the number of cells, which is what makes million-cell datasets feasible. The null assumes independence between the ligand and receptor summaries under random labelling; the paper shows that most ligand–receptor pairs in real data are weakly or not correlated.

## Sixteen modules and their combination

Two summaries families (mean; median, Q3, 90th percentile) × two complex aggregations × two combinations give the sixteen default scoring modules, each with its own p-value for every interaction. They are combined with the Cauchy combination test, which is valid under arbitrary dependence between the modules:

<p align="center">
  <img src="{{ site.baseurl }}/images/cct.svg" width="420" alt="Cauchy combination test statistic">
</p>

The pattern of the sixteen p-values is informative on its own. When the median-based modules are significant a typical cell expresses the ligand (a *broad* signal); when only the upper-quantile modules are, a fraction of the cells carries it (*sub-population*); when the average-aggregated modules pass but the minimum-aggregated ones fail, one subunit limits a complex; when arithmetic passes but geometric fails, the two sides are unbalanced. The explorer reports these as the *signal type*, and `--summary / --aggregation / --combination` select other module sets at run time.

## Reference-based analysis

The permutation-free null makes a second kind of question tractable: not "which interactions are significant in my data" but "which interactions differ from healthy tissue". A reference panel stores, for every gene and reference cell type, the distribution of within-cell expression ranks (expression is rank-transformed within each cell, so library size and normalisation do not matter), built once from a large healthy dataset (Fig. 1c–d). For a query (Fig. 1e):

1. the query is rank-transformed the same way and its cell types are matched to the panel's (mapping);
2. a calibration factor \\(k\\) is estimated from housekeeping genes: the spread of the query's gene ranks around the reference's (the ratio of mean rank to root-mean-square deviation over rank bins, floored at 3) widens the significance bands accordingly;
3. for every interaction, the reference-derived null of the query cell-type pair (reference per-gene distributions convolved for the query's cell numbers) gives a significance band; the ligand and receptor summaries of the query are compared with their reference intervals; the interaction is significant in the query if its score exceeds the band or its components exceed their intervals, and significant in the reference by the reference's own p-value;
4. the two calls give `Up`, `Down`, `Both Sig` or `Both NS`.

The built-in panels cover 19 human tissues from healthy data in CELLxGENE Census, about 16 million cells across assays and studies; any dataset or control group can be turned into a panel with `build-ref`.

<p align="center">
  <img src="{{ site.baseurl }}/images/reference.png" width="860" alt="Composition of the 19 reference panels: sequencing assays per tissue, and number of cells and cell types per tissue">
</p>
<p align="center" style="font-size:0.85rem;color:#666;margin-top:-0.4rem">
The 19 built-in panels: assay composition per tissue (left) and number of cells and cell types (right).
</p>

## Effect sizes, bootstrap, replication

- **Effect size.** For the mean/average/arithmetic module the null has closed-form moments: a cell type's mean of a gene is the mean of \\(n\\) independent draws, so its null mean is the gene's overall mean \\(\mu\\) and its variance \\(\sigma^2/n\\); a complex averages its subunits, the score averages the two sides. The standardised effect \\((\mathrm{score}-\mathrm{null\ mean})/\mathrm{null\ SD}\\) ranks interactions whose p-values are at the floor.
- **Bootstrap support.** The analysis is repeated on subsamples keeping 80% of the cells of every cell type; the support of an interaction is the fraction of subsamples in which it stays significant.
- **Replication.** With samples as groups, the fraction of samples in which an interaction found in the pooled data is significant on its own.

## Communication patterns and enrichment in the explorer

Communication patterns are a non-negative matrix factorisation of the cell type × ligand–receptor matrix of significant interactions (multiplicative updates, fixed seed), computed separately for senders and receivers. Pathway enrichment is a hypergeometric test of the ligand and receptor genes of the current selection against all genes of the ligand–receptor database, with Benjamini–Hochberg correction, using Enrichr's Hallmark, KEGG and GO Biological Process libraries pruned to that universe.

## Extreme p-values in the combination

Module p-values are computed as tail probabilities (survival function for the analytic nulls, tail sums for the discretised ones) and clamped to \\([0, 1]\\). In the Cauchy combination a module p-value of exactly 1 enters as \\(1 - 1/16\\) and p-values below \\(10^{-15}\\) enter as \\(1/(\pi p)\\), the tail approximation of the ACAT implementation. `small_p_method='legacy'` (`--small-p legacy`) applies the rule of the published method instead: p-values at or below \\(10^{-16}\\) are treated as 0 and \\(p = 1\\) enters as \\(\tan(-\pi/2)\\).
