---
title: Distribution toolkit
layout: default
parent: Methods
nav_order: 2
description: "Exact distributions of sums, means and order statistics of discrete random variables: the engine behind FastCCC's permutation-free p-values, usable on its own for simulations and other tests."
---
<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js">
</script>

# Distribution toolkit

FastCCC never permutes. For every communication score it derives the null distribution analytically from the per-cell distribution of each gene, and the p-value is read from that distribution. The pieces that make this possible are ordinary functions in `fastccc.distrib` and `fastccc.dist_complex`; they work for any discrete random variable, not only expression, and can be used on their own.

Two building blocks carry the method:

1. **Sums and means of $$n$$ i.i.d. draws**, by FFT convolution of probability mass functions (the mean-based summaries).
2. **Order statistics of $$n$$ i.i.d. draws from a discrete distribution**, exactly, ties included (the median, quartile and percentile summaries). The textbook formula is a continuous-variable result; the construction in the paper (Methods, eqs. 13–18) handles discrete data with ties.

## Order statistics of a discrete variable
{: .d-inline-block }

exact, any n
{: .label .label-green }

FastCCC's median, quartile and percentile summaries are order statistics: with $$n$$ cells in a cell type and a percentile $$q$$, the summary is the $$k^{th}$$ smallest value with $$k = \lfloor q\,(n-1) + 1\rfloor$$ (paper, Methods). Its null distribution is therefore the distribution of the $$k^{th}$$ order statistic of $$n$$ i.i.d. draws from the gene's expression distribution. The textbook result (eq. 12 of the paper) is for **continuous** variables:

$$f_{X_{(k)}}(x') = \frac{n!}{(k-1)!\,(n-k)!}\Bigl(\int_0^{x'} f(x)\,dx\Bigr)^{k-1} f(x') \Bigl(1-\int_0^{x'} f(x)\,dx\Bigr)^{n-k}.$$

It assumes that all sampled values are distinct. Single-cell expression is discrete, with a high proportion of zeros, and repeated values across cells are the norm, so this formula cannot be applied to expression data.

**Example.** A gene takes the values 0, 1, 2 with probabilities 0.5, 0.3, 0.2. Draw $$n = 10$$ cells and take the 5th smallest ($$q = 0.5$$, the median of 10 cells).

| value | 0 | 1 | 2 |
|---|---:|---:|---:|
| exact pmf of the 5th of 10 draws | 0.6230 | 0.3706 | 0.0064 |
| continuous formula (raw values 1.2305, 0.0495, 0; sum 1.28), renormalised | 0.9613 | 0.0387 | 0.0000 |
| Monte Carlo, 2 × 10⁶ trials | 0.6234 | 0.3703 | 0.0063 |
| `order_statistic_pmf` / `get_quantile_pmf_for_n_iid_distribution` | 0.6230 | 0.3706 | 0.0064 |

**The construction (paper, eqs. 13–18).** FastCCC extends the classical formulation to discrete distributions with ties through a reference sample. Let $$Y_1, \dots, Y_n$$ be $$n$$ draws from the uniform distribution $$U(0,1)$$, sorted as $$Y_{(1)} < \dots < Y_{(n)}$$. The $$k^{th}$$ order statistic of the uniform sample has the density

$$g_{Y_{(k)}}(y) = \frac{n!}{(k-1)!\,(n-k)!}\,y^{k-1}(1-y)^{n-k} \qquad (13)$$

and no ties, because $$Y$$ is continuous. The expression variable $$X$$ is connected to $$Y$$ through the cumulative distribution functions $$F_X$$ and $$G_Y$$, both with range $$[0,1]$$, by the map $$F_X(x) = G_Y(y)$$:

$$y = G_Y^{-1}\bigl(F_X(x)\bigr), \qquad f(x) = g(y)\,\frac{dy}{dx}, \qquad (14,\,15)$$

so the cumulative distribution of the $$k^{th}$$ order statistic of $$X$$ is obtained by integrating the uniform order statistic up to the mapped point:

$$\int_0^{x'} f_{X_{(k)}}(x)\,dx = \int_0^{\,G_Y^{-1}(F_X(x'))} g_{Y_{(k)}}(y)\,dy. \qquad (16)$$

Only $$F_X$$ enters, which is why ties cause no trouble: all draws equal to $$x$$ are counted by the jump of $$F_X$$ at $$x$$. For computation the interval $$[0,1]$$ is discretised into $$m = 100{,}000$$ segments and the cdf of the order statistic of the discrete uniform $$\tilde Y$$ is the normalised cumulative sum of $$g_{Y_{(k)}}$$ on the grid,

$$G_{\tilde Y_{(k)}}(i) = P(\tilde y \le i) \simeq \frac{1}{Z}\sum_{j=0}^{i} g_{Y_{(k)}}\!\left(\frac{j}{m}\right), \qquad (17)$$

with $$Z$$ the normalising constant. For a discrete expression variable $$\tilde X$$ with support $$x_1 < x_2 < \dots < x_l$$ and cdf $$F_{\tilde X}$$, the probability mass of the $$k^{th}$$ order statistic at $$x_i$$ is the mass of $$\tilde Y_{(k)}$$ between the mapped images of $$x_{i-1}$$ and $$x_i$$:

$$p\bigl(\tilde X_{(k)} = x_i\bigr) = G_{\tilde Y_{(k)}}\bigl(\lfloor m\,F_{\tilde X}(x_i)\rfloor\bigr) - G_{\tilde Y_{(k)}}\bigl(\lfloor m\,F_{\tilde X}(x_{i-1})\rfloor\bigr). \qquad (18)$$

`get_quantile_pmf_for_n_iid_distribution` implements eq. 16 on FastCCC's 0.01 grid: it takes the cumulative mass $$F_{\tilde X}$$ at every grid value, evaluates the cdf of the order statistic there in closed form (the regularised incomplete beta function, see the remark below) and differences along the support. The cost does not depend on $$n$$, so a cell type of 100,000 cells is as cheap as one of 10, and the input can be any pmf on the grid, including the empirical distribution of a gene over all cells.

*Remark.* The integral in eq. 16 is the regularised incomplete beta function, $$\int_0^{u} g_{Y_{(k)}}(y)\,dy = I_u(k, n-k+1)$$, which also equals $$P\bigl(\mathrm{Binomial}(n, u) \ge k\bigr)$$: at least $$k$$ of the $$n$$ draws fall at or below $$x$$. `order_statistic_pmf` uses the same closed form on the exact support (any real values, any $$n$$), which is convenient for a quick calculation or a simulation study.

```python
from fastccc.distrib import order_statistic_pmf

values, pmf = order_statistic_pmf([0, 1, 2], [0.5, 0.3, 0.2], n=10, k=5)
pmf                                   # array([0.623 , 0.3706, 0.0064])
order_statistic_pmf([0, 1, 2], [0.5, 0.3, 0.2], n=10, quantile=0.5)   # the same, by quantile
order_statistic_pmf([0, 1, 2], [0.5, 0.3, 0.2], n=10, k=10)           # the maximum of 10 draws
```

```python
import numpy as np
from fastccc.distrib import Distribution, distribution_from_values, get_quantile_pmf_for_n_iid_distribution, get_pvalue_from_pmf, precision

d = distribution_from_values([0, 1, 2], [0.5, 0.3, 0.2])       # the same pmf on the grid
pmf = get_quantile_pmf_for_n_iid_distribution(d, n=10, quantile=0.5)
pmf[[0, int(1 / precision), int(2 / precision)]]               # array([0.623 , 0.3706, 0.0064])

# a real gene: one cell's distribution from the data, the 90th percentile of a type of 250 cells
x = np.asarray(adata[:, 'CD74'].X.todense()).ravel()            # log1p expression over all cells
gene = Distribution(samples=x)
q90 = Distribution(pmf_array=get_quantile_pmf_for_n_iid_distribution(gene, n=250, quantile=0.9), is_align=True)
q90.get_mean(), get_pvalue_from_pmf(observed, q90)             # null mean; upper-tail p-value of an observed value
```

**Convention.** The quantile $$q$$ of $$n$$ draws is the $$k^{th}$$ order statistic with $$k = \lfloor q\,(n-1) + 1\rfloor$$ (`order_statistic_index(n, q)`; for $$n = 10$$: $$q = 0.5 \to 5$$, $$q = 0.9 \to 9$$). FastCCC's observed summaries use the same index (`numpy.quantile` with `method='lower'`), so observed value and null refer to the same statistic. `fastccc.dist_iid_set.calculate_cluster_quantile_distribution(counts_df, labels_df, quantile)` returns the null for every gene and cell type of a dataset at once.

## Sums and means of i.i.d. draws

`Distribution` holds a random variable as a pmf on a grid of bin width 0.01 over $$[0, 14]$$ (the range of log1p expression), or analytically when it is normal. The operators implement what FastCCC needs:

| expression | distribution of |
|---|---|
| `Distribution(samples=x)` | one draw from the empirical distribution of `x` (one cell's expression) |
| `Distribution(dtype='normal', loc=m, scale=s)` | a normal variable, kept analytic |
| `d ** n` | the sum of $$n$$ i.i.d. draws; FFT convolution up to $$n = 29$$, the central-limit normal from $$n = 30$$ on |
| `d / n` | the variable divided by $$n$$; `(d ** n) / n` is the mean of $$n$$ cells |
| `d1 + d2` | $$X + Y$$ of two independent variables (FFT convolution) |
| `d1 & d2` | $$X + Y$$ when both are positive, 0 when either is 0 (the arithmetic ligand–receptor rule) |
| `d1 * d2` | $$\sqrt{XY}$$ (the geometric ligand–receptor rule) |
| `get_minimum_distribution(d1, d2, ...)` | the minimum of independent variables (complexes aggregated by their weakest subunit) |
| `get_average_distribution(d1, d2, ...)` | their average |
| `get_pvalue_from_pmf(v, d)` | $$P(X \ge v)$$, computed from the tail itself so that values far below $$10^{-16}$$ stay exact |

```python
import numpy as np
from fastccc.distrib import Distribution, get_pvalue_from_pmf
from fastccc.dist_complex import get_minimum_distribution, get_average_distribution

lig = Distribution(samples=np.asarray(adata[:, 'TGFB1'].X.todense()).ravel())
r1  = Distribution(samples=np.asarray(adata[:, 'TGFBR1'].X.todense()).ravel())
r2  = Distribution(samples=np.asarray(adata[:, 'TGFBR2'].X.todense()).ravel())

L = (lig ** 120) / 120                       # mean of TGFB1 over a sender type of 120 cells, under random labels
R = get_minimum_distribution((r1 ** 80) / 80, (r2 ** 80) / 80)   # receptor complex in a receiver type of 80 cells
cs_arith = (L & R) / 2                        # null of the arithmetic communication score
cs_geom = L * R                               # null of the geometric communication score
get_pvalue_from_pmf(observed_score, cs_arith) # the p-value FastCCC reports for the Mean / Minimum / Arithmetic module
```

This is the computation behind every one of FastCCC's sixteen modules: a summary (mean or an order statistic) of the ligand in the sender and of the receptor in the receiver, a complex aggregation (minimum or average), and a combination (arithmetic or geometric), each step applied to distributions instead of numbers. The observed score is then compared with the resulting null. `get_pvalue_from_pmf` returns the upper tail directly rather than $$1 - F$$, which is why FastCCC's p-values do not collapse to 0 or 1 at the floating-point limit.

## Using it for your own problems

- Any discrete variable fits: counts, scores, grades, bounded measurements. Put its support on the grid with `distribution_from_values(values, probs)` (values in $$[0, 14]$$ at 0.01 resolution; `update_precision` changes the resolution) or work on the exact support with `order_statistic_pmf`.
- The grid routines return pmf arrays; `Distribution(pmf_array=..., is_align=True)` wraps one, `.get_mean()`, `.get_std()`, `.pmf`, `.cdf` read it.
- Monte Carlo is the natural check: 2 million draws reproduce the example above to three decimals, and the exact routine is instantaneous for any $$n$$.

The derivation and its use in the scoring are described in the paper (Hou, Ma & Zhou, *Nat Commun* 16, 11428, 2025) and summarised in [How FastCCC works]({{ site.baseurl }}/methods/how-it-works.html).
