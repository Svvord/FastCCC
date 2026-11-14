"""Module p-values must stay inside [0, 1] even when the convolved pmf carries floating-point noise."""
import numpy as np
import pytest

from fastccc.distrib import Distribution, get_pvalue_from_pmf, precision


def test_pvalue_from_noisy_pmf_is_never_negative():
    # a pmf whose entries sum to slightly more than 1, as FFT convolution can produce
    pmf = np.array([0.2, 0.3, 0.5 + 3e-16])
    dist = Distribution(dtype='other', pmf_array=pmf, is_align=True)
    beyond = len(pmf) * precision * 2  # a value past the whole support: p should be exactly 0, not -3e-16
    p = get_pvalue_from_pmf(beyond, dist)
    assert p == 0.0
    assert 0.0 <= get_pvalue_from_pmf(precision, dist) <= 1.0
    assert get_pvalue_from_pmf(0.0, dist) == 1.0


def test_pvalue_from_normal_pmf_uses_the_survival_function():
    dist = Distribution(dtype='normal', loc=1.0, scale=0.1)
    far = get_pvalue_from_pmf(2.0, dist)        # ten standard deviations: 1 - cdf would round to 0
    assert 0 < far < 1e-20
    assert get_pvalue_from_pmf(-10.0, dist) == 1.0
    assert get_pvalue_from_pmf(1.0, dist) == pytest.approx(0.5)


def test_discretised_tail_is_summed_directly():
    pmf = np.zeros(400); pmf[0] = 1 - 1e-12; pmf[-1] = 1e-12   # almost all mass at 0, a sliver at the end
    dist = Distribution(dtype='other', pmf_array=pmf, is_align=True)
    p = get_pvalue_from_pmf(precision * 300, dist)
    assert p == pytest.approx(1e-12, rel=1e-6)      # 1 - sum(head) would return 1e-12 +- 1e-16 noise


def test_power_of_sample_distribution_without_prior_mean_access():
    """``Distribution(samples=...) ** n`` must not depend on the mean having been read before."""
    import numpy as np
    from fastccc.distrib import Distribution
    rng = np.random.default_rng(0)
    x = np.log1p(rng.poisson(2.0, size=2000).astype(float))
    d = Distribution(samples=x)
    s30 = d ** 30                       # normal approximation branch
    assert abs(s30.loc - 30 * Distribution(samples=x).get_mean()) < 1e-9
    assert abs(s30.scale - np.sqrt(30) * Distribution(samples=x).get_std()) < 1e-9
    m30 = s30 / 30
    assert abs(m30.get_mean() - d.get_mean()) < 0.05


def test_order_statistic_pmf_exact_for_discrete_draws():
    """The discrete order statistic must match brute-force enumeration (ties!) and FastCCC's grid routine."""
    import itertools
    import numpy as np
    from fastccc.distrib import order_statistic_pmf, distribution_from_values, get_quantile_pmf_for_n_iid_distribution, precision
    values, probs = [0, 1, 2], [0.5, 0.3, 0.2]
    # brute force over all 3^n outcomes for n = 6, k = 3
    n, k = 6, 3
    brute = np.zeros(3)
    for draw in itertools.product(range(3), repeat=n):
        brute[sorted(draw)[k - 1]] += np.prod([probs[d] for d in draw])
    v, pmf = order_statistic_pmf(values, probs, n, k=k)
    assert np.allclose(pmf, brute, atol=1e-12) and list(v) == values
    # the documented example: 5th of 10 draws
    _, pmf10 = order_statistic_pmf(values, probs, 10, quantile=0.5)
    assert np.allclose(pmf10, [0.6230, 0.3706, 0.0064], atol=5e-4)
    # FastCCC's grid-based routine gives the same numbers
    d = distribution_from_values(values, probs)
    grid = get_quantile_pmf_for_n_iid_distribution(d, 10, quantile=0.5)
    idx = [int(round(x / precision)) for x in values]
    assert np.allclose(grid[idx], pmf10, atol=1e-6) and abs(grid.sum() - 1) < 1e-9
    # the continuous density formula is wrong here: it does not even sum to one
    F = np.cumsum(probs); cont = 1260 * F**4 * np.array(probs) * (1 - F)**5
    assert cont.sum() > 1.2


def test_quantile_null_closed_form_matches_binomial_reference_at_large_n():
    """The production routine evaluates eq. 16 exactly; the former grid tabulation drifts for large n."""
    import numpy as np
    from scipy.stats import binom
    from fastccc.distrib import (distribution_from_values, get_quantile_pmf_for_n_iid_distribution, _quantile_pmf_grid,
                                 order_statistic_index, precision)
    values, probs = [0, 1, 2], [0.5, 0.3, 0.2]
    d = distribution_from_values(values, probs); idx = [int(round(v / precision)) for v in values]
    for n, q in ((10, 0.5), (1000, 0.9), (50000, 0.5)):
        k = order_statistic_index(n, q)
        ref = np.diff([binom.sf(k - 1, n, f) for f in np.cumsum(probs)], prepend=0)
        new = get_quantile_pmf_for_n_iid_distribution(d, n, q)
        assert np.allclose(new[idx], ref, atol=1e-9) and abs(new.sum() - 1) < 1e-9
        assert np.allclose(_quantile_pmf_grid(d, n, q)[idx], ref, atol=2e-3)       # the table is close for a coarse pmf
    # a dense pmf (many support values with small masses): the table degrades with n, the closed form does not
    rng = np.random.default_rng(0)
    x = np.round(np.log1p(rng.negative_binomial(1, 0.2, size=20000).astype(float)), 2)
    from fastccc.distrib import Distribution
    dense = Distribution(samples=x)
    for n, tol in ((100, 3e-3), (50000, 6e-2)):
        a, b = get_quantile_pmf_for_n_iid_distribution(dense, n, 0.75), _quantile_pmf_grid(dense, n, 0.75)
        L = max(len(a), len(b)); tv = 0.5 * np.abs(np.pad(a, (0, L - len(a))) - np.pad(b, (0, L - len(b)))).sum()
        assert tv < tol
