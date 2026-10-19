"""Estimate h(q) for every index by MFDFA -- the structurally independent
estimator cross-checked against scripts/03_estimate.py's structure-function B.

Mirrors 03_estimate.py exactly (bootstrap CIs, shuffle surrogates, a scale-
range robustness grid in place of the tau_max grid, Mann-Whitney on the
developed/emerging split) so the two estimators are held to the same bar
before anything is concluded from comparing them. Scale range and order
chosen in scripts/02b_controls_mfdfa.py by screening against known-answer
synthetic processes -- see that file's comments for why s_min=25, s_max=150,
order=1.
"""
import os, sys, warnings
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from ghe import multiscaling_strength, h1_minus_h4
from mfdfa import mfdfa, block_bootstrap, shuffle_surrogates

QS = np.arange(0.5, 4.01, 0.25)
ORDER = 1
S_MIN = 25
S_MAX = 150
SMAX_GRID = [80, 100, 150, 200, 300]
NBOOT = 500
NSURR = 500


def make_scales(smax, smin=S_MIN, npts=16):
    return np.unique(np.round(np.logspace(np.log10(smin), np.log10(smax), npts)).astype(int))


SCALES = make_scales(S_MAX)

META = {r[0]: r[1:] for r in [
    ("^GSPC","S&P 500","United States","Developed"),
    ("^FTSE","FTSE 100","United Kingdom","Developed"),
    ("^GDAXI","DAX","Germany","Developed"),
    ("^FCHI","CAC 40","France","Developed"),
    ("^N225","Nikkei 225","Japan","Developed"),
    ("^AXJO","ASX 200","Australia","Developed"),
    ("^GSPTSE","S&P/TSX","Canada","Developed"),
    ("^SSMI","SMI","Switzerland","Developed"),
    ("^AEX","AEX","Netherlands","Developed"),
    ("^IBEX","IBEX 35","Spain","Developed"),
    ("^BFX","BEL 20","Belgium","Developed"),
    ("^ATX","ATX","Austria","Developed"),
    ("^HSI","Hang Seng","Hong Kong","Developed"),
    ("^BVSP","Bovespa","Brazil","Emerging"),
    ("XU100.IS","BIST 100","Turkey","Emerging"),
    ("^BSESN","BSE Sensex","India","Emerging"),
    ("^JKSE","Jakarta Composite","Indonesia","Emerging"),
    ("^KS11","KOSPI","South Korea","Emerging"),
    ("^TWII","TAIEX","Taiwan","Emerging"),
    ("^MXX","IPC","Mexico","Emerging"),
    ("^KLSE","FTSE Bursa KLCI","Malaysia","Emerging"),
    ("^SET.BK","SET","Thailand","Emerging"),
]}


def analyse(args):
    sym, X = args
    X = np.asarray(X, float)
    point = mfdfa(X, QS, SCALES, order=ORDER)
    B = multiscaling_strength(QS, point["H"])
    D = h1_minus_h4(QS, point["H"])

    boot = block_bootstrap(X, QS, SCALES, order=ORDER, n_boot=NBOOT, seed=abs(hash(sym)) % 2**31)
    surr = shuffle_surrogates(X, QS, SCALES, order=ORDER, n_surr=NSURR, seed=7)
    p = float((surr["B"] >= B).mean())

    smax_sens = {f"B_smax{s}": multiscaling_strength(
                     QS, mfdfa(X, QS, make_scales(s), order=ORDER)["H"])
                 for s in SMAX_GRID}

    name, country, tier = META[sym]
    row = dict(symbol=sym, index=name, country=country, tier=tier, n=len(X),
               H1=float(np.interp(1, QS, point["H"])),
               H2=float(np.interp(2, QS, point["H"])),
               H4=float(np.interp(4, QS, point["H"])),
               B=B, B_lo=float(np.percentile(boot["B"], 2.5)),
               B_hi=float(np.percentile(boot["B"], 97.5)),
               D=D, D_lo=float(np.percentile(boot["D"], 2.5)),
               D_hi=float(np.percentile(boot["D"], 97.5)),
               B_shuffled=float(surr["B"].mean()),
               B_shuffled_p95=float(np.percentile(surr["B"], 95)),
               p_vs_shuffle=p, block=boot["block"],
               r2_min=float(np.min(point["r2"])), **smax_sens)
    return row


def main():
    lp = pd.read_csv(os.path.join(HERE, "data", "clean", "logprices.csv"))
    series = [(s, g["logprice"].to_numpy()) for s, g in lp.groupby("symbol", sort=False)]
    print(f"estimating {len(series)} series by MFDFA ...", flush=True)

    with ProcessPoolExecutor() as ex:
        rows = pd.DataFrame(list(ex.map(analyse, series)))

    rows = rows.sort_values("B", ascending=False).reset_index(drop=True)
    rows.to_csv(os.path.join(HERE, "results", "summary_mfdfa.csv"), index=False)

    sf = pd.read_csv(os.path.join(HERE, "results", "summary.csv"))[
        ["symbol", "B", "B_lo", "B_hi", "p_vs_shuffle"]
    ].rename(columns={"B": "B_sf", "B_lo": "B_sf_lo", "B_hi": "B_sf_hi",
                       "p_vs_shuffle": "p_vs_shuffle_sf"})
    cmp = rows.merge(sf, on="symbol").rename(columns={"B": "B_mf", "B_lo": "B_mf_lo",
                                                       "B_hi": "B_mf_hi",
                                                       "p_vs_shuffle": "p_vs_shuffle_mf"})
    cmp.to_csv(os.path.join(HERE, "results", "sf_vs_mfdfa.csv"), index=False)

    null_mf = np.load(os.path.join(HERE, "results", "null_B_mfdfa_brownian.npy"))
    null95 = float(np.percentile(null_mf, 95))

    cols = ["index", "country", "tier", "H1", "H2", "H4", "B", "B_lo", "B_hi",
            "p_vs_shuffle", "r2_min"]
    pd.set_option("display.width", 220)
    print(rows[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    print(f"\nMFDFA null 95th pct (random walk, n=7000): {null95:.4f}")
    n_beyond = int((rows.B > null95).sum())
    print(f"{n_beyond}/{len(rows)} markets exceed the MFDFA null at B > {null95:.4f}")

    dev = rows.loc[rows.tier == "Developed", "B"]
    eme = rows.loc[rows.tier == "Emerging", "B"]
    try:
        from scipy.stats import mannwhitneyu
        u, pmw = mannwhitneyu(eme, dev, alternative="greater")
    except Exception:
        u, pmw = np.nan, np.nan
    print(f"\nMFDFA B: developed mean {dev.mean():.4f} (n={len(dev)}), "
          f"emerging mean {eme.mean():.4f} (n={len(eme)})")
    print(f"Mann-Whitney U (emerging > developed), MFDFA: U={u}, p={pmw:.4f}")

    print("\nscale-range robustness (B at each s_max, s_min=25 fixed), "
          "developed vs emerging means:")
    grid_cols = [f"B_smax{s}" for s in SMAX_GRID]
    print(rows.groupby("tier")[grid_cols].mean().to_string(float_format=lambda v: f"{v:.4f}"))

    print("\nSF vs MFDFA agreement:")
    corr = cmp[["B_sf", "B_mf"]].corr().iloc[0, 1]
    print(f"Pearson correlation of B across 22 markets: {corr:.3f}")
    disagree = cmp.loc[(cmp.B_sf - cmp.B_sf.median()) * (cmp.B_mf - cmp.B_mf.median()) < 0]
    print(f"markets on opposite sides of their method's own median: "
          f"{len(disagree)}/{len(cmp)} -> {list(disagree.symbol)}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
