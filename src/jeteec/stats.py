"""Statistics helpers for closure tests: estimators whose own bias has been checked.

The point of a toy sample is that the answer is known. That only helps if the estimator
used to extract the answer is unbiased itself -- the original spectrum fit here was not
(an unweighted least-squares fit to log(counts) returned n = 4.85 for an input of 5.0).
"""

from __future__ import annotations

import numpy as np


def _mean_var_u(t: float, span: float) -> tuple[float, float]:
    """Mean and variance of u in [0, span] with density ~ exp(-t u)."""
    if np.isinf(span):
        return 1.0 / t, 1.0 / t**2
    x = t * span
    if abs(x) < 1e-6:  # series expansion around the flat limit
        return span / 2 - t * span**2 / 12, span**2 / 12
    mean = 1.0 / t - span / np.expm1(x)
    var = 1.0 / t**2 - (span / (2 * np.sinh(x / 2))) ** 2
    return mean, var


def power_law_index(x, lo: float, hi: float = np.inf) -> tuple[float, float]:
    """Maximum-likelihood index ``n`` of ``dN/dx ~ x**-n`` for samples truncated to [lo, hi).

    Unbinned, so there is no histogram, no "enough entries" bin cut and no fit to
    log(counts) -- the ingredients that biased the least-squares estimate. With
    ``u = ln(x / lo)`` the model is an exponential in ``u`` with slope ``t = n - 1``, and
    the likelihood equation reduces to matching the mean of ``u``. For ``hi = inf`` this is
    the Hill estimator, ``n = 1 + N / sum(ln(x / lo))``.

    Returns
    -------
    n, sigma : float
        The estimate and its standard error from the Fisher information.
    """
    x = np.asarray(x, dtype=float)
    x = x[(x >= lo) & (x < hi)]
    n_obs = len(x)
    if n_obs < 2:
        raise ValueError("need at least two values inside [lo, hi)")
    u_bar = float(np.mean(np.log(x / lo)))
    span = float(np.log(hi / lo))

    if np.isinf(span):
        t = 1.0 / u_bar
    else:
        if not 0.0 < u_bar < span:
            raise ValueError("sample mean of ln(x/lo) outside (0, ln(hi/lo))")
        # E_t[u] decreases monotonically from span (t -> -inf) to 0 (t -> +inf): bisect
        a, b = -700.0 / span, 700.0 / span
        for _ in range(200):
            t = 0.5 * (a + b)
            if _mean_var_u(t, span)[0] > u_bar:
                a = t
            else:
                b = t
        t = 0.5 * (a + b)

    var_u = _mean_var_u(t, span)[1]
    return 1.0 + t, 1.0 / np.sqrt(n_obs * var_u)


def power_law_integral(n: float, lo: float, hi: float) -> float:
    """``integral_lo^hi x**-n dx``, e.g. to normalise a fitted spectrum for plotting."""
    if np.isclose(n, 1.0):
        return float(np.log(hi / lo))
    return float((lo ** (1 - n) - hi ** (1 - n)) / (n - 1))


def binomial_ratio(k: int, n: int) -> tuple[float, float]:
    """``k / n`` and its binomial standard error, for a subset ratio such as R32."""
    if n <= 0:
        return float("nan"), float("nan")
    r = k / n
    return r, float(np.sqrt(r * (1 - r) / n))
