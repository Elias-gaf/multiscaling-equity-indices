"""How much of the developed/emerging gap is dynamics, and how much is tails?

The shuffle surrogates in script 03 destroy all temporal structure while
keeping the marginal distribution EXACTLY.  They come back with *more*
multiscaling than the real series -- and so does a shuffled MRW, a process
whose multiscaling is known to be real.  So the plain shuffle is not a valid
one-sided test; what it gives is a marginal-matched benchmark.  Three things
are done here:

  1. excess multiscaling, dB = B_obs - E[B_shuffled], which nets out whatever
     the marginal distribution alone can manufacture;
  2. a regression of B on excess kurtosis and the developed/emerging dummy;
  3. a low-moment re-estimate using only q <= 2, where moment estimation is
     far less tail-sensitive.

The question in each case is the same: does emerging > developed survive?
"""
import os, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from ghe import generalised_hurst, multiscaling_strength

QS_FULL = np.arange(0.5, 4.01, 0.25)
QS_LOW = np.arange(0.5, 2.01, 0.25)


def ols_t(X, y, names):
    X = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    r2 = 1 - (resid @ resid) / ((y - y.mean()) @ (y - y.mean()))
    return pd.DataFrame({"coef": beta, "se": se, "t": t, "p": p},
                        index=["const"] + names), r2


def main():
    s = pd.read_csv(os.path.join(HERE, "results", "summary.csv"))
    lp = pd.read_csv(os.path.join(HERE, "data", "clean", "logprices.csv"))

    kurt, blow = {}, {}
    for sym, g in lp.groupby("symbol"):
        X = g["logprice"].to_numpy()
        kurt[sym] = float(pd.Series(np.diff(X)).kurtosis())
        blow[sym] = multiscaling_strength(QS_LOW, generalised_hurst(X, QS_LOW)["H"])
    s["exkurt"] = s.symbol.map(kurt)
    s["B_lowq"] = s.symbol.map(blow)
    s["dB"] = s.B - s.B_shuffled
    s["emerging"] = (s.tier == "Emerging").astype(float)
    s.to_csv(os.path.join(HERE, "results", "summary.csv"), index=False)

    def gap(col, label, drop_klse=False):
        d = s[s.symbol != "^KLSE"] if drop_klse else s
        e = d.loc[d.tier == "Emerging", col]
        v = d.loc[d.tier == "Developed", col]
        u, p = stats.mannwhitneyu(e, v, alternative="greater")
        tag = " (ex-Malaysia)" if drop_klse else ""
        print(f"{label+tag:38s} emerging {e.mean():+.4f}  developed {v.mean():+.4f}"
              f"   gap {e.mean()-v.mean():+.4f}   MW p={p:.4f}")

    print("=== does emerging > developed survive? ===")
    gap("B", "raw multiscaling B")
    gap("B", "raw multiscaling B", True)
    gap("dB", "excess over marginal-matched null")
    gap("dB", "excess over marginal-matched null", True)
    gap("B_lowq", "low-moment B (q<=2 only)")
    gap("B_lowq", "low-moment B (q<=2 only)", True)

    print("\n=== is B just a kurtosis ranking? ===")
    r = stats.spearmanr(s.B, s["exkurt"])
    print(f"Spearman(B, excess kurtosis) = {r.statistic:+.3f}  p={r.pvalue:.4f}")
    r2_ = stats.spearmanr(s.B[s.symbol != "^KLSE"], s["exkurt"][s.symbol != "^KLSE"])
    print(f"  ex-Malaysia                = {r2_.statistic:+.3f}  p={r2_.pvalue:.4f}")

    for col in ["B", "B_lowq"]:
        tab, r2 = ols_t(np.column_stack([s["exkurt"], s["emerging"]]), s[col].to_numpy(),
                        ["kurtosis", "emerging"])
        print(f"\nOLS {col} ~ kurtosis + emerging   (R2={r2:.3f}, n={len(s)})")
        print(tab.to_string(float_format=lambda v: f"{v:.5f}"))

    print("\n=== tau_max robustness (raw B) ===")
    tcols = [c for c in s.columns if c.startswith("B_taumax")]
    print(s.groupby("tier")[tcols].mean().to_string(float_format=lambda v: f"{v:.4f}"))
    for c in tcols:
        e = s.loc[s.tier == "Emerging", c]; v = s.loc[s.tier == "Developed", c]
        _, p = stats.mannwhitneyu(e, v, alternative="greater")
        print(f"  {c}: gap {e.mean()-v.mean():+.4f}, p={p:.4f}")


if __name__ == "__main__":
    main()
