"""Generalised Hurst exponent H(q) by the structure-function method.

Reference: T. Di Matteo, "Multi-scaling in finance", Quantitative Finance 7(1), 2007.

For a series X(t) (log prices), the q-th order structure function is

    S_q(tau) = < |X(t + tau) - X(t)|^q >_t

and is assumed to scale as S_q(tau) ~ tau^{q H(q)}.  H(q) is read off the slope
of log S_q(tau) against log tau, divided by q.

    H(q) constant in q  -> uniscaling (monofractal); H = 0.5 is a random walk
    H(q) decreasing in q -> multiscaling
"""
from __future__ import annotations

import numpy as np


def structure_functions(X: np.ndarray, qs: np.ndarray, taus: np.ndarray) -> np.ndarray:
    """S_q(tau) for every q at once: one pass over taus.  Shape (len(qs), len(taus))."""
    X = np.asarray(X, dtype=float)
    qs = np.asarray(qs, dtype=float)
    out = np.empty((len(qs), len(taus)), dtype=float)
    for i, tau in enumerate(taus):
        d = np.abs(X[tau:] - X[:-tau])
        logd = np.log(d, out=np.full_like(d, -np.inf), where=d > 0)
        out[:, i] = np.mean(np.exp(qs[:, None] * logd[None, :]), axis=1)
    return out


def _ols_many(x: np.ndarray, Y: np.ndarray):
    """Least squares of each row of Y on x.  Returns (slopes, se, r2)."""
    n = len(x)
    xm = x.mean()
    sxx = np.sum((x - xm) ** 2)
    ym = Y.mean(axis=1, keepdims=True)
    b = ((Y - ym) @ (x - xm)) / sxx
    a = ym[:, 0] - b * xm
    resid = Y - (a[:, None] + b[:, None] * x[None, :])
    rss = np.sum(resid ** 2, axis=1)
    se_b = np.sqrt((rss / (n - 2)) / sxx)
    ss_tot = np.sum((Y - ym) ** 2, axis=1)
    r2 = np.where(ss_tot > 0, 1.0 - rss / ss_tot, np.nan)
    return b, se_b, r2


def _ols(x: np.ndarray, y: np.ndarray):
    """Least squares y = a + b x. Returns (b, se_b, r2)."""
    n = len(x)
    xm, ym = x.mean(), y.mean()
    sxx = np.sum((x - xm) ** 2)
    b = np.sum((x - xm) * (y - ym)) / sxx
    a = ym - b * xm
    resid = y - (a + b * x)
    dof = n - 2
    s2 = np.sum(resid ** 2) / dof
    se_b = np.sqrt(s2 / sxx)
    ss_tot = np.sum((y - ym) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan
    return b, se_b, r2


def generalised_hurst(X: np.ndarray, qs, tau_max: int = 19, tau_min: int = 1):
    """H(q) for each q in qs.

    Returns dict with arrays: q, H, se (from the OLS slope only), r2.
    The standard error is the regression standard error and is *not* a valid
    confidence interval for H(q) -- the points are strongly dependent.  Use the
    block bootstrap in `block_bootstrap` for that.
    """
    qs = np.atleast_1d(np.asarray(qs, dtype=float))
    taus = np.arange(tau_min, tau_max + 1)
    S = structure_functions(X, qs, taus)
    b, se_b, r2 = _ols_many(np.log(taus), np.log(S))
    return {"q": qs, "H": b / qs, "se": se_b / qs, "r2": r2}


def multiscaling_strength(qs: np.ndarray, H: np.ndarray) -> float:
    """B in the linear fit H(q) = A - B q.  Larger B = stronger multiscaling."""
    b, _, _ = _ols(np.asarray(qs, float), np.asarray(H, float))
    return -b


def h1_minus_h4(qs: np.ndarray, H: np.ndarray) -> float:
    """H(1) - H(4), the simple multiscaling proxy used in Di Matteo et al.

    Raises if qs does not span [1, 4]: np.interp would otherwise clip
    silently to the nearest boundary value instead of extrapolating,
    returning a mislabelled H(1)-H(2)-or-similar with no warning.
    """
    qs = np.asarray(qs, float)
    if qs.min() > 1.0 or qs.max() < 4.0:
        raise ValueError(
            f"h1_minus_h4 needs qs to span [1, 4]; got [{qs.min()}, {qs.max()}]"
        )
    return float(np.interp(1.0, qs, H) - np.interp(4.0, qs, H))


def singularity_spectrum(qs: np.ndarray, H: np.ndarray):
    """f(alpha) via the Legendre transform of the mass exponent tau(q).

        tau(q) = q H(q) - 1
        alpha(q) = d(tau)/dq          (the Holder/singularity exponent)
        f(alpha) = q * alpha(q) - tau(q)

    alpha is a numerical derivative over the q-grid, so it is noisier than H
    itself and that noise is worst exactly where H(q) is already least
    reliable -- the high-q end (THEORY.md SS6-7).  Treat `width` the same way
    `B` is treated: bootstrap/shuffle it before reading anything into it, not
    as a drop-in replacement for B.

    Returns dict with q, alpha, f, and width = alpha.max() - alpha.min().
    """
    qs = np.asarray(qs, float)
    H = np.asarray(H, float)
    tau = qs * H - 1.0
    alpha = np.gradient(tau, qs)
    f = qs * alpha - tau
    return {"q": qs, "alpha": alpha, "f": f, "width": float(alpha.max() - alpha.min())}


# --------------------------------------------------------------------------
# uncertainty
# --------------------------------------------------------------------------

def circular_block_bootstrap(increments: np.ndarray, block: int, rng) -> np.ndarray:
    """One circular-block-bootstrap resample of an increment series."""
    n = len(increments)
    nblocks = int(np.ceil(n / block))
    starts = rng.integers(0, n, size=nblocks)
    idx = (starts[:, None] + np.arange(block)[None, :]).ravel()[:n] % n
    return increments[idx]


def block_bootstrap(X: np.ndarray, qs, tau_max: int = 19, n_boot: int = 500,
                    block: int | None = None, seed: int = 0):
    """Block-bootstrap CIs for H(q) and for the multiscaling summaries.

    Resamples blocks of *increments* (returns), rebuilds a price path by
    cumulative sum, and re-estimates.  Block length defaults to
    max(20, round(n^{1/3})) -- the floor dominates at the sample sizes used
    here (n~7000 gives n^{1/3}~19.1, so every market in this project is
    bootstrapped with block=20) -- long enough to preserve the volatility
    clustering that generates the multiscaling.
    """
    X = np.asarray(X, float)
    d = np.diff(X)
    n = len(d)
    if block is None:
        block = max(20, int(round(n ** (1 / 3))))
    rng = np.random.default_rng(seed)

    qs = np.atleast_1d(np.asarray(qs, float))
    Hs = np.empty((n_boot, len(qs)))
    B = np.empty(n_boot)
    D = np.empty(n_boot)
    for b in range(n_boot):
        dd = circular_block_bootstrap(d, block, rng)
        Xb = np.concatenate([[0.0], np.cumsum(dd)])
        res = generalised_hurst(Xb, qs, tau_max=tau_max)
        Hs[b] = res["H"]
        B[b] = multiscaling_strength(qs, res["H"])
        D[b] = h1_minus_h4(qs, res["H"])
    return {"H": Hs, "B": B, "D": D, "block": block}


def shuffle_surrogates(X: np.ndarray, qs, tau_max: int = 19, n_surr: int = 500,
                       seed: int = 0):
    """Shuffle the increments (destroys all temporal structure, keeps the
    marginal distribution exactly).  Any multiscaling left is an artefact of
    fat tails and finite samples, not of the dynamics."""
    X = np.asarray(X, float)
    d = np.diff(X)
    rng = np.random.default_rng(seed)
    qs = np.atleast_1d(np.asarray(qs, float))
    B = np.empty(n_surr)
    D = np.empty(n_surr)
    Hs = np.empty((n_surr, len(qs)))
    for s in range(n_surr):
        dd = rng.permutation(d)
        Xs = np.concatenate([[0.0], np.cumsum(dd)])
        res = generalised_hurst(Xs, qs, tau_max=tau_max)
        Hs[s] = res["H"]
        B[s] = multiscaling_strength(qs, res["H"])
        D[s] = h1_minus_h4(qs, res["H"])
    return {"H": Hs, "B": B, "D": D}
