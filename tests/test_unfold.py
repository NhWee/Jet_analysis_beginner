"""Hand-checkable tests of jeteec.unfold (D'Agostini for weighted histograms)."""

import numpy as np
import pytest

from jeteec import unfold


def test_identity_response_returns_data_whatever_prior():
    A = np.eye(3)
    d = np.array([5.0, 3.0, 2.0])
    for prior in (np.ones(3), np.array([10.0, 1.0, 0.1])):
        t = unfold.dagostini(d, A, prior, n_iter=1)
        assert np.allclose(t, d)


def test_diagonal_efficiency_is_divided_out():
    # 50% / 80% / 100% "efficiency" (or mean reco/gen weight ratio) per bin
    A = np.diag([0.5, 0.8, 1.0])
    d = np.array([1.0, 4.0, 2.0])
    t = unfold.dagostini(d, A, np.ones(3), n_iter=3)
    assert np.allclose(t, [2.0, 5.0, 2.0])


def test_prior_equal_to_truth_is_a_fixed_point():
    # hand-made 2x2 response with migration; data = A t + fake exactly
    t_true = np.array([4.0, 6.0])
    A = np.array([[0.7, 0.2], [0.1, 0.6]])   # columns sum to 0.8: 20% of gen weight lost
    fake = np.array([0.5, 0.25])
    d = A @ t_true + fake                     # [2.8+1.2+0.5, 0.4+3.6+0.25] = [4.5, 4.25]
    assert np.allclose(d, [4.5, 4.25])
    for n in (1, 5, 20):
        assert np.allclose(unfold.dagostini(d, A, t_true, n, fake=fake), t_true)


def test_converges_to_truth_from_flat_prior():
    t_true = np.array([4.0, 6.0])
    A = np.array([[0.7, 0.2], [0.1, 0.6]])
    d = A @ t_true
    _, its = unfold.dagostini(d, A, np.array([5.0, 5.0]), 200, return_all=True)
    # first iteration, by hand: nu = A @ [5,5] = [4.5, 3.5]; d = [4.0, 4.0]
    # t1_0 = 5 * (0.7*4/4.5 + 0.1*4/3.5) / 0.8 = 5 * (0.62222 + 0.11429) / 0.8 = 4.6032
    assert its[0][0] == pytest.approx(5 * (0.7 * 4 / 4.5 + 0.1 * 4 / 3.5) / 0.8)
    assert np.allclose(its[-1], t_true, rtol=1e-6)


def test_fake_fraction_equivalent_to_absolute_fakes():
    A = np.array([[0.9, 0.1], [0.0, 0.8]])
    t = np.array([2.0, 3.0])
    fake = np.array([0.6, 0.6])
    d = A @ t + fake
    f = fake / d
    a = unfold.dagostini(d, A, t, 2, fake=fake)
    b = unfold.dagostini(d, A, t, 2, fake_fraction=f)
    assert np.allclose(a, b) and np.allclose(a, t)


def test_zero_efficiency_bin_stays_zero():
    A = np.array([[1.0, 0.0], [0.0, 0.0]])
    t = unfold.dagostini([3.0, 1.0], A, [1.0, 1.0], 4)
    assert np.allclose(t, [3.0, 0.0])


def test_response_from_counts_and_fold():
    M = np.array([[2.0, 1.0], [0.0, 3.0]])
    T = np.array([4.0, 0.0])
    A = unfold.response_from_counts(M, T)
    assert np.allclose(A, [[0.5, 0.0], [0.0, 0.0]])
    assert np.allclose(unfold.fold(A, [4.0, 1.0], fake=[1.0, 1.0]), [3.0, 1.0])


def test_bin_by_bin():
    corr, c = unfold.bin_by_bin([2.0, 3.0], [1.0, 2.0], [2.0, 1.0])
    assert np.allclose(c, [2.0, 0.5]) and np.allclose(corr, [4.0, 1.5])


def test_iteration_pvalue():
    p, chi2 = unfold.iteration_pvalue([1.0, 2.0], [1.0, 1.0], [1.0, 0.0])
    assert chi2 == pytest.approx(0.0) and p == pytest.approx(1.0)
    p, chi2 = unfold.iteration_pvalue([3.0, 1.0], [1.0, 1.0], [1.0, 1.0])
    assert chi2 == pytest.approx(4.0)
    assert p == pytest.approx(np.exp(-2.0))   # chi2 sf with ndf 2 is exp(-chi2/2)


def test_smooth_ratio_exact_without_noise_and_unbiased_with_noise():
    import unfold_mod as U

    x = np.log(U.CENTRES)
    true = np.exp(0.1 * x + 0.01 * x**2)              # a smooth "model change"
    den, err = np.ones_like(x), np.full(len(x), 0.05)
    assert np.allclose(U.smooth_ratio(true, den, err), true, rtol=1e-12)
    rng = np.random.default_rng(0)
    fits = np.array([U.smooth_ratio(true * (1 + 0.05 * rng.standard_normal(len(x))), den, err)
                     for _ in range(200)])
    rel = fits / true - 1
    assert np.max(np.abs(rel.mean(axis=0))) < 0.01      # no bias (noise 5% per bin)
    # scatter of the fit is smaller than the 5% per-bin noise, except near the two ends
    assert np.max(rel.std(axis=0)[5:-5]) < 0.025
