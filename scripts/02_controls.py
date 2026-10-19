"""Validate the estimator on processes whose H(q) is known analytically.

Sample length is matched to the real data (n = 7000, about 28 years of daily
closes) so the biases measured here are the biases that apply to the market
estimates.

  Brownian motion      H(q) = 0.5            (flat)
  fBm, H = 0.3 / 0.7   H(q) = H              (flat, and not 0.5)
  MRW, lam2            H(q) = 0.5 + l - l q/2 (falling, known slope l/2)

The Brownian replicas also give the NULL DISTRIBUTION of the multiscaling
statistic B at this sample size: a market's B only means something if it lies
outside that band.
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
import synthetic as sy

N = 7000
NREP = 300
TAU_MAX = 19
QS = np.arange(0.5, 4.01, 0.25)
QREP = np.array([1., 2., 3., 4.])


def run(gen, label, theory_H, theory_B):
    Hs, Bs, Ds = [], [], []
    for s in range(NREP):
        X = gen(np.random.default_rng(10_000 + s))
        r = generalised_hurst(X, QS, tau_max=TAU_MAX)
        Hs.append(r["H"]); Bs.append(multiscaling_strength(QS, r["H"]))
        Ds.append(h1_minus_h4(QS, r["H"]))
    Hs = np.array(Hs)
    return dict(label=label, H_mean=Hs.mean(0), H_sd=Hs.std(0), H_all=Hs,
                B=np.array(Bs), D=np.array(Ds),
                theory_H=theory_H, theory_B=theory_B)


def main():
    cases = [
        run(lambda r: sy.brownian_motion(N, r), "Brownian motion (H=0.5)",
            np.full_like(QS, 0.5), 0.0),
        run(lambda r: sy.fractional_bm(N, 0.7, r), "fBm, H=0.7",
            np.full_like(QS, 0.7), 0.0),
        run(lambda r: sy.fractional_bm(N, 0.3, r), "fBm, H=0.3",
            np.full_like(QS, 0.3), 0.0),
        run(lambda r: sy.mrw(N, 0.02, 2048, rng=r), "MRW, $\\lambda^2$=0.02",
            sy.mrw_theory_H(QS, 0.02), 0.01),
        run(lambda r: sy.mrw(N, 0.04, 2048, rng=r), "MRW, $\\lambda^2$=0.04",
            sy.mrw_theory_H(QS, 0.04), 0.02),
    ]

    rows = []
    for c in cases:
        rows.append(dict(
            process=c["label"].replace("$\\lambda^2$", "lambda^2"),
            H1=np.interp(1, QS, c["H_mean"]), H2=np.interp(2, QS, c["H_mean"]),
            H4=np.interp(4, QS, c["H_mean"]),
            B_mean=c["B"].mean(), B_sd=c["B"].std(),
            B_theory=c["theory_B"],
            recovery=(c["B"].mean() / c["theory_B"]) if c["theory_B"] else np.nan,
            B_null_p95=np.percentile(c["B"], 95)))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(HERE, "results", "controls.csv"), index=False)

    # iid Student-t: fat tails, ZERO temporal structure.  This is the control
    # that matters most -- it sets how much B the marginal alone can produce.
    tt = []
    for df in (5, 4, 3):
        Bs = [multiscaling_strength(QS, generalised_hurst(
                  np.concatenate([[0.0], np.cumsum(
                      np.random.default_rng(20_000 + s_).standard_t(df, N))]),
                  QS, tau_max=TAU_MAX)["H"]) for s_ in range(NREP)]
        tt.append(dict(dof=df, B_mean=float(np.mean(Bs)), B_sd=float(np.std(Bs)),
                       B_p95=float(np.percentile(Bs, 95))))
    pd.DataFrame(tt).to_csv(os.path.join(HERE, "results", "controls_iid_t.csv"), index=False)
    print("\niid Student-t random walks (no dynamics whatsoever):")
    print(pd.DataFrame(tt).to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    null_B = cases[0]["B"]
    np.save(os.path.join(HERE, "results", "null_B_brownian.npy"), null_B)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    cols = plt.get_cmap("tab10")
    for i, c in enumerate(cases):
        ax.plot(QS, c["H_mean"], "-", color=cols(i), label=c["label"])
        ax.fill_between(QS, c["H_mean"] - 2 * c["H_sd"], c["H_mean"] + 2 * c["H_sd"],
                        color=cols(i), alpha=0.15, lw=0)
        ax.plot(QS, c["theory_H"], "--", color=cols(i), lw=1)
    ax.set_xlabel("$q$"); ax.set_ylabel("$H(q)$")
    ax.set_title(f"Estimator vs theory, n={N}, {NREP} replicas\n(dashed = analytic)")
    ax.legend(fontsize=7.5, loc="center left")
    ax.grid(alpha=.3)

    ax = axes[1]
    ax.hist(null_B, bins=45, color="0.35", alpha=.85)
    for pc, st in [(95, "--"), (99, ":")]:
        v = np.percentile(null_B, pc)
        ax.axvline(v, color="crimson", ls=st, lw=1.4, label=f"{pc}th pct = {v:.4f}")
    ax.set_xlabel("$B$  (slope of $H(q)$, sign-flipped)")
    ax.set_ylabel("count")
    ax.set_title("Null distribution of the multiscaling\nstatistic under a random walk")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "figures", "fig1_controls.png"), dpi=170)

    pd.set_option("display.width", 200)
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(f"\nnull B under a random walk at n={N}: "
          f"mean {null_B.mean():+.4f}, sd {null_B.std():.4f}, "
          f"95th {np.percentile(null_B,95):.4f}, 99th {np.percentile(null_B,99):.4f}")


if __name__ == "__main__":
    main()
