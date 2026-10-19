"""Synthetic processes with known scaling, used to validate the estimator.

    brownian_motion   H(q) = 0.5 for all q          (negative control)
    fractional_bm     H(q) = H  for all q           (negative control, H != 0.5)
    mrw               H(q) = 0.5 + lam2 - lam2*q/2  (positive control)

The MRW is the Bacry-Delour-Muzy multifractal random walk; its structure
function exponent is zeta(q) = q H(q) = q(1/2 + lam2) - lam2 q^2 / 2.
"""
from __future__ import annotations

import numpy as np


def brownian_motion(n: int, rng) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(rng.standard_normal(n))])


def _fgn_davies_harte(n: int, H: float, rng) -> np.ndarray:
    """Fractional Gaussian noise by circulant embedding (exact for fGn)."""
    k = np.arange(0, n + 1)
    g = 0.5 * (np.abs(k - 1) ** (2 * H) - 2 * np.abs(k) ** (2 * H)
               + np.abs(k + 1) ** (2 * H))
    c = np.concatenate([g, g[-2:0:-1]])          # length 2n
    lam = np.fft.fft(c).real
    if lam.min() < -1e-8:
        raise ValueError("circulant embedding not positive definite")
    lam = np.clip(lam, 0.0, None)
    m = len(c)
    V = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    Y = np.fft.fft(np.sqrt(lam / m) * V)
    return Y[:n].real


def fractional_bm(n: int, H: float, rng) -> np.ndarray:
    fgn = _fgn_davies_harte(n, H, rng)
    fgn = fgn / fgn.std()
    return np.concatenate([[0.0], np.cumsum(fgn)])


def _log_correlated_gaussian(n: int, lam2: float, L: int, rng) -> np.ndarray:
    """Gaussian field with Cov(w_i, w_j) = lam2 * ln(L / (|i-j| + 1)) for
    |i-j| < L, by circulant embedding.  Negative circulant eigenvalues are
    clipped -- standard practice, and the reason the recovered exponents carry
    a small bias.  Mean is set to -lam2*ln(L) so that E[exp(2w)] = 1."""
    m = 1
    while m < 4 * n:
        m *= 2
    k = np.arange(m)
    dist = np.minimum(k, m - k)
    cov = np.where(dist < L, lam2 * np.log(L / (dist + 1.0)), 0.0)
    eig = np.fft.fft(cov).real
    eig = np.clip(eig, 0.0, None)
    V = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    w = np.fft.fft(np.sqrt(eig / m) * V).real[:n]
    return w - lam2 * np.log(L)


def mrw(n: int, lam2: float = 0.02, L: int = 2048, sigma: float = 1.0,
        rng=None) -> np.ndarray:
    """Bacry-Delour-Muzy multifractal random walk."""
    rng = rng if rng is not None else np.random.default_rng(0)
    w = _log_correlated_gaussian(n, lam2, L, rng)
    eps = rng.standard_normal(n) * sigma
    return np.concatenate([[0.0], np.cumsum(eps * np.exp(w))])


def mrw_theory_H(q, lam2: float = 0.02) -> np.ndarray:
    q = np.asarray(q, float)
    return 0.5 + lam2 - lam2 * q / 2.0
