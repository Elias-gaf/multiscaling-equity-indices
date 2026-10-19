"""Multifractal Detrended Fluctuation Analysis (MFDFA).

Reference: Kantelhardt et al., "Multifractal detrended fluctuation analysis
of nonstationary time series", Physica A 316, 2002.

A structurally different estimator of the same h(q) that `ghe.py` gets from
structure functions.  Where the structure function averages |X(t+tau)-X(t)|^q
directly, MFDFA builds the cumulative profile of the *increments*, splits it
into non-overlapping segments at each scale s, removes a local polynomial
trend of order `order` within every segment, and averages the q-th moment of
the detrended local variance.  The two methods should agree closely here:
our increments are already stationary, so there's no slow drift for the
detrending step to usefully remove -- which is exactly what makes agreement
between them a check on the estimator rather than a foregone conclusion.

    Y(i) = sum_{k<=i} (d_k - mean(d))      profile of the increments d
    F^2(v, s) = local variance of Y in segment v of length s, after removing
                a degree-`order` polynomial fit to that segment
    F_q(s) = { mean_v [F^2(v,s)]^(q/2) }^(1/q)
    F_q(s) ~ s^h(q)                        h(q) read off by log-log OLS

For q=2 and order>=1 this is standard DFA; h(2) is the Hurst exponent.  As in
`ghe.py`, tau(q) = q h(q) - 1, so `multiscaling_strength`, `h1_minus_h4` and
`singularity_spectrum` from `ghe.py` apply unchanged to MFDFA's h(q).
"""
from __future__ import annotations

import numpy as np

from ghe import _ols_many, circular_block_bootstrap, multiscaling_strength, h1_minus_h4


def profile(X: np.ndarray) -> np.ndarray:
    """Cumulative sum of the mean-centred increments of X (a log-price series)."""
    d = np.diff(np.asarray(X, float))
    return np.cumsum(d - d.mean())


def _segment_variances(Y: np.ndarray, s: int, order: int) -> np.ndarray:
    """F^2(v, s) for every non-overlapping segment of length s, counted from
    both ends of Y (standard MFDFA practice: uses all of a series that isn't
    an exact multiple of s, instead of discarding the remainder)."""
    n = len(Y)
    n_seg = n // s
    t = np.arange(s, dtype=float)
    A = np.vstack([t ** k for k in range(order + 1)]).T
    A_pinv = np.linalg.pinv(A)

    F2 = np.empty(2 * n_seg)
    for v in range(n_seg):
        seg = Y[v * s:(v + 1) * s]
        trend = A @ (A_pinv @ seg)
        F2[v] = np.mean((seg - trend) ** 2)
    for v in range(n_seg):
        seg = Y[n - (v + 1) * s: n - v * s]
        trend = A @ (A_pinv @ seg)
        F2[n_seg + v] = np.mean((seg - trend) ** 2)
    return F2


def fluctuation_function(X: np.ndarray, qs, scales, order: int = 1) -> np.ndarray:
    """F_q(s) for every q and scale at once.  Shape (len(qs), len(scales)).

    q=0 is not handled (it needs a log-average rather than a power mean) --
    unneeded here since the project's q-grid never includes it.
    """
    qs = np.atleast_1d(np.asarray(qs, dtype=float))
    if np.any(qs == 0):
        raise ValueError("fluctuation_function does not support q=0")
    scales = np.atleast_1d(np.asarray(scales, dtype=int))
    Y = profile(X)

    out = np.empty((len(qs), len(scales)))
    for j, s in enumerate(scales):
        F2 = np.clip(_segment_variances(Y, s, order), 1e-300, None)
        out[:, j] = np.mean(F2[None, :] ** (qs[:, None] / 2.0), axis=1) ** (1.0 / qs)
    return out


def mfdfa(X: np.ndarray, qs, scales, order: int = 1):
    """h(q) for each q in qs, by MFDFA.

    Returns dict with arrays: q, H (= h(q)), se, r2, and the scales used --
    same shape as `ghe.generalised_hurst`'s return value, so every downstream
    function written against that (multiscaling_strength, h1_minus_h4,
    singularity_spectrum, plotting code) works unchanged on this one.
    """
    qs = np.atleast_1d(np.asarray(qs, dtype=float))
    scales = np.atleast_1d(np.asarray(scales, dtype=int))
    Fq = fluctuation_function(X, qs, scales, order=order)
    b, se_b, r2 = _ols_many(np.log(scales), np.log(Fq))
    return {"q": qs, "H": b, "se": se_b, "r2": r2, "scales": scales}


# --------------------------------------------------------------------------
# uncertainty -- mirrors ghe.block_bootstrap / ghe.shuffle_surrogates exactly,
# just with mfdfa() in place of generalised_hurst()
# --------------------------------------------------------------------------

def block_bootstrap(X: np.ndarray, qs, scales, order: int = 1, n_boot: int = 500,
                    block: int | None = None, seed: int = 0):
    """Block-bootstrap CIs for MFDFA's H(q) and the multiscaling summaries."""
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
        res = mfdfa(Xb, qs, scales, order=order)
        Hs[b] = res["H"]
        B[b] = multiscaling_strength(qs, res["H"])
        D[b] = h1_minus_h4(qs, res["H"])
    return {"H": Hs, "B": B, "D": D, "block": block}


def shuffle_surrogates(X: np.ndarray, qs, scales, order: int = 1, n_surr: int = 500,
                       seed: int = 0):
    """Shuffle the increments (destroys temporal structure, keeps the marginal
    exactly) and re-estimate with MFDFA -- the same null-distribution check
    as `ghe.shuffle_surrogates`, for the second estimator."""
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
        res = mfdfa(Xs, qs, scales, order=order)
        Hs[s] = res["H"]
        B[s] = multiscaling_strength(qs, res["H"])
        D[s] = h1_minus_h4(qs, res["H"])
    return {"H": Hs, "B": B, "D": D}
