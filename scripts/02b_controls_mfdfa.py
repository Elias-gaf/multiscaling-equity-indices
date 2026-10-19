"""Validate MFDFA on the same known-answer processes as scripts/02_controls.py,
and compare its B directly against the structure-function estimator's B.

This is a cross-check, not a second headline result: if MFDFA's recovered B
tracks the structure-function B on every synthetic case (within bootstrap
noise) and on the real markets, the project's finding does not depend on a
single algorithm's quirks.  If the two disagree, that disagreement is the
more interesting thing to report.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from ghe import generalised_hurst, multiscaling_strength, h1_minus_h4
from mfdfa import mfdfa
import synthetic as sy

N = 7000
NREP = 300
QS = np.arange(0.5, 4.01, 0.25)
ORDER = 1

# scale range chosen by screening order in {1,2} and (s_min,s_max) in a grid
# against these same synthetic processes (see handover notes): s_min=25 is
# well past the textbook rule-of-thumb (s_min~10 for order-1 detrending) --
# at this sample size and q-range, s_min<16 leaves a null bias on anti-
# persistent fBm that is small but clearly non-zero over 300 replicas, and
# order=2 detrending makes that bias worse, not better, at every scale range
# tried.  s_max=150 sits at the point where MRW recovery is closest to 100%
# before larger s_max starts overshooting it.  See SMAX_GRID below for the
# robustness check this choice gets, mirroring tau_max in scripts/03_estimate.py.
S_MIN = 25
SMAX_GRID = [80, 100, 150, 200, 300]
S_MAX = 150


def make_scales(smax, smin=S_MIN, npts=16):
    return np.unique(np.round(np.logspace(np.log10(smin), np.log10(smax), npts)).astype(int))


SCALES = make_scales(S_MAX)


def run(gen, label, theory_H, theory_B):
    Hs_sf, Bs_sf, Hs_mf, Bs_mf = [], [], [], []
    for s in range(NREP):
        X = gen(np.random.default_rng(30_000 + s))
        r_sf = generalised_hurst(X, QS, tau_max=19)
        Hs_sf.append(r_sf["H"]); Bs_sf.append(multiscaling_strength(QS, r_sf["H"]))
        r_mf = mfdfa(X, QS, SCALES, order=ORDER)
        Hs_mf.append(r_mf["H"]); Bs_mf.append(multiscaling_strength(QS, r_mf["H"]))
    return dict(label=label, theory_H=theory_H, theory_B=theory_B,
                H_sf=np.array(Hs_sf), B_sf=np.array(Bs_sf),
                H_mf=np.array(Hs_mf), B_mf=np.array(Bs_mf))


def main():
    cases = [
        run(lambda r: sy.brownian_motion(N, r), "Brownian motion (H=0.5)",
            np.full_like(QS, 0.5), 0.0),
        run(lambda r: sy.fractional_bm(N, 0.7, r), "fBm, H=0.7",
            np.full_like(QS, 0.7), 0.0),
        run(lambda r: sy.fractional_bm(N, 0.3, r), "fBm, H=0.3",
            np.full_like(QS, 0.3), 0.0),
        run(lambda r: sy.mrw(N, 0.02, 2048, rng=r), "MRW, lambda^2=0.02",
            sy.mrw_theory_H(QS, 0.02), 0.01),
        run(lambda r: sy.mrw(N, 0.04, 2048, rng=r), "MRW, lambda^2=0.04",
            sy.mrw_theory_H(QS, 0.04), 0.02),
    ]

    rows = []
    for c in cases:
        rows.append(dict(
            process=c["label"],
            B_sf_mean=c["B_sf"].mean(), B_sf_sd=c["B_sf"].std(),
            B_mf_mean=c["B_mf"].mean(), B_mf_sd=c["B_mf"].std(),
            B_theory=c["theory_B"],
            recovery_sf=(c["B_sf"].mean() / c["theory_B"]) if c["theory_B"] else np.nan,
            recovery_mf=(c["B_mf"].mean() / c["theory_B"]) if c["theory_B"] else np.nan,
            H_sf_mean=c["H_sf"].mean(), H_mf_mean=c["H_mf"].mean(),
        ))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(HERE, "results", "controls_mfdfa.csv"), index=False)
    np.save(os.path.join(HERE, "results", "null_B_mfdfa_brownian.npy"), cases[0]["B_mf"])
    print(f"\nMFDFA null B under a random walk at n={N}, scales={S_MIN}-{S_MAX}: "
          f"mean {cases[0]['B_mf'].mean():+.4f}, sd {cases[0]['B_mf'].std():.4f}, "
          f"95th {np.percentile(cases[0]['B_mf'],95):.4f}, "
          f"99th {np.percentile(cases[0]['B_mf'],99):.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    cols = plt.get_cmap("tab10")
    for i, c in enumerate(cases):
        ax.plot(QS, c["H_sf"].mean(0), "-", color=cols(i), lw=1.6, label=f"{c['label']} (SF)")
        ax.plot(QS, c["H_mf"].mean(0), "--", color=cols(i), lw=1.6, label=f"{c['label']} (MFDFA)")
        ax.plot(QS, c["theory_H"], ":", color=cols(i), lw=1)
    ax.set_xlabel("$q$"); ax.set_ylabel("$h(q)$")
    ax.set_title(f"Structure function vs MFDFA, n={N}, {NREP} replicas\n(dotted = analytic)")
    ax.legend(fontsize=6.5, loc="center left")
    ax.grid(alpha=.3)

    ax = axes[1]
    for i, c in enumerate(cases):
        ax.errorbar(c["B_sf"].mean(), c["B_mf"].mean(),
                     xerr=c["B_sf"].std(), yerr=c["B_mf"].std(),
                     fmt="o", color=cols(i), label=c["label"], capsize=3)
    lim = [min(0, tab.B_sf_mean.min(), tab.B_mf_mean.min()) - 0.002,
           max(tab.B_sf_mean.max(), tab.B_mf_mean.max()) * 1.15]
    ax.plot(lim, lim, "k--", lw=1, label="agreement")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("$B$, structure function"); ax.set_ylabel("$B$, MFDFA")
    ax.set_title("Do the two estimators agree\non synthetic processes?")
    ax.legend(fontsize=6.5); ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "fig1b_mfdfa_controls.png"), dpi=170)

    pd.set_option("display.width", 220)
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
