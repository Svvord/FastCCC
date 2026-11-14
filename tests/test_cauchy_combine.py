"""Cauchy combination of module p-values: clamping, the extreme-p rules and the legacy option."""
import numpy as np
import pytest

from fastccc.cauchy_combine import combine_pvalues


def mat(*rows):
    return np.array(rows, dtype=float)[:, :, None]  # (pairs, modules, 1) as in the result tables


def test_uniform_half_gives_half():
    _, p = combine_pvalues(mat([0.5] * 16))
    assert abs(float(p[0, 0]) - 0.5) < 1e-12


def test_matches_formula_for_ordinary_pvalues():
    ps = np.linspace(0.05, 0.95, 16)
    T, p = combine_pvalues(mat(ps))
    expect = 0.5 - np.arctan(np.mean(np.tan(np.pi * (0.5 - ps)))) / np.pi
    assert abs(float(p[0, 0]) - expect) < 1e-12


def test_p_equal_one_is_weak_evidence_not_a_veto():
    rows = mat([1e-14] * 4 + [1.0] * 12)   # four strongly significant summaries, twelve blind ones
    _, p_acat = combine_pvalues(rows)
    _, p_legacy = combine_pvalues(rows, small_p_method='legacy')
    assert float(p_acat[0, 0]) < 1e-10
    assert float(p_legacy[0, 0]) > 0.99


def test_single_weak_module_cannot_carry_the_interaction():
    _, p = combine_pvalues(mat([0.03] + [1.0] * 15))
    assert float(p[0, 0]) > 0.5


def test_exact_zero_gives_zero_and_tiny_p_is_resolved():
    _, p0 = combine_pvalues(mat([0.0] + [0.5] * 15))
    assert float(p0[0, 0]) == 0.0
    _, pa = combine_pvalues(mat([1e-20] + [0.5] * 15))
    _, pl = combine_pvalues(mat([1e-20] + [0.5] * 15), small_p_method='legacy')
    assert float(pa[0, 0]) < 1e-17 < float(pl[0, 0]) < 1e-14


def test_negative_input_is_clamped_and_nan_propagates():
    _, a = combine_pvalues(mat([-1e-16] + [0.5] * 15))
    _, b = combine_pvalues(mat([0.0] + [0.5] * 15))
    assert float(a[0, 0]) == float(b[0, 0])
    _, n = combine_pvalues(mat([np.nan] + [0.5] * 15))
    assert np.isnan(n[0, 0])


def test_unknown_method_rejected():
    with pytest.raises(ValueError):
        combine_pvalues(mat([0.5] * 16), small_p_method='other')



def test_single_module_returns_its_own_pvalues():
    """With one module there is nothing to combine. The 1 - 1/d truncation of p = 1 gave 0 for d = 1,
    which declared every unexpressed interaction significant."""
    import numpy as np
    from fastccc.cauchy_combine import combine_pvalues
    p = np.array([1.0, 0.9, 0.5, 0.05, 1e-3, 1e-20, 0.0])
    for method in ('acat', 'legacy'):
        _, combined = combine_pvalues(p[:, None], small_p_method=method)
        assert combined.shape == p.shape
        assert combined[0] == 1.0, (method, combined)
        # an exact 0 stays 0 under 'acat'; the published rule caps it at the floating-point floor of tan(pi/2)
        assert combined[-1] == 0.0 if method == 'acat' else combined[-1] < 1e-15, (method, combined)
        np.testing.assert_allclose(combined[1:5], p[1:5], rtol=1e-9)
    _, acat = combine_pvalues(p[:, None], small_p_method='acat')
    np.testing.assert_allclose(acat[5], 1e-20, rtol=1e-9)
    # two modules keep the documented rule: p = 1 counts as 1 - 1/2
    _, two = combine_pvalues(np.array([[1.0, 1.0], [1.0, 0.01]]), small_p_method='acat')
    assert two[0] == 0.5 and 0.01 < two[1] < 0.05
