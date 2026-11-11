"""The estimators used for closure tests must be unbiased themselves."""

import numpy as np
import pytest

from jeteec import stats


def sample_power_law(rng, n, size, lo, hi):
    u = rng.uniform(0, 1, size)
    a = 1.0 - n
    return (lo**a + u * (hi**a - lo**a)) ** (1.0 / a)


def test_power_law_index_recovers_the_input():
    rng = np.random.default_rng(3)
    x = sample_power_law(rng, 5.0, 20000, 100.0, 1000.0)
    n, sigma = stats.power_law_index(x, 100.0, 1000.0)
    assert abs(n - 5.0) < 4 * sigma
    assert sigma == pytest.approx(4.0 / np.sqrt(len(x)), rel=0.1)


def test_power_law_index_is_unbiased_with_calibrated_errors():
    # the unweighted least-squares fit this replaces returned 4.85 for an input of 5.0
    rng = np.random.default_rng(4)
    pulls = []
    for _ in range(200):
        x = sample_power_law(rng, 5.0, 3000, 100.0, 1000.0)
        n, sigma = stats.power_law_index(x, 100.0, 1000.0)
        pulls.append((n - 5.0) / sigma)
    pulls = np.array(pulls)
    assert abs(pulls.mean()) < 0.25
    assert 0.8 < pulls.std() < 1.2


def test_untruncated_limit_is_the_hill_estimator():
    rng = np.random.default_rng(5)
    x = sample_power_law(rng, 4.0, 5000, 50.0, 1e12)
    n, _ = stats.power_law_index(x, 50.0)
    assert n == pytest.approx(1.0 + len(x) / np.sum(np.log(x / 50.0)))


def test_power_law_integral():
    assert stats.power_law_integral(5.0, 100.0, 1000.0) == pytest.approx(
        (100.0**-4 - 1000.0**-4) / 4)
    assert stats.power_law_integral(1.0, 1.0, np.e) == pytest.approx(1.0)


def test_binomial_ratio():
    r, err = stats.binomial_ratio(287, 1902)
    assert r == pytest.approx(287 / 1902)
    assert err == pytest.approx(np.sqrt(r * (1 - r) / 1902))
    assert np.isnan(stats.binomial_ratio(0, 0)[0])
