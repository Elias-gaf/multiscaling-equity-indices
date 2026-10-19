"""Estimate H(q) for every index, with block-bootstrap CIs, shuffle surrogates
and a tau_max robustness check.  Writes results/ tables and figures/."""
import os, sys, warnings
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from ghe import (generalised_hurst, multiscaling_strength, h1_minus_h4,
                 block_bootstrap, shuffle_surrogates)

QS = np.arange(0.5, 4.01, 0.25)
TAU_MAX = 19
NBOOT = 500
NSURR = 500
TAU_GRID = [10, 15, 19, 25, 30]

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
    point = generalised_hurst(X, QS, tau_max=TAU_MAX)
    B = multiscaling_strength(QS, point["H"])
    D = h1_minus_h4(QS, point["H"])

    boot = block_bootstrap(X, QS, tau_max=TAU_MAX, n_boot=NBOOT, seed=abs(hash(sym)) % 2**31)
    surr = shuffle_surrogates(X, QS, tau_max=TAU_MAX, n_surr=NSURR, seed=7)

    # p-value: how often does destroying the time ordering reproduce this much
    # multiscaling?  Tests whether the effect is dynamics or just fat tails.
    p = float((surr["B"] >= B).mean())

    tau_sens = {f"B_taumax{t}": multiscaling_strength(QS, generalised_hurst(X, QS, tau_max=t)["H"])
                for t in TAU_GRID}

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
               r2_min=float(np.min(point["r2"])), **tau_sens)
    curves = pd.DataFrame(dict(symbol=sym, tier=tier, q=QS, H=point["H"],
                               H_lo=np.percentile(boot["H"], 2.5, axis=0),
                               H_hi=np.percentile(boot["H"], 97.5, axis=0)))
    return row, curves, surr["B"]


def main():
    lp = pd.read_csv(os.path.join(HERE, "data", "clean", "logprices.csv"))
    series = [(s, g["logprice"].to_numpy()) for s, g in lp.groupby("symbol", sort=False)]
    print(f"estimating {len(series)} series ...", flush=True)

    with ProcessPoolExecutor() as ex:
        out = list(ex.map(analyse, series))

    rows = pd.DataFrame([o[0] for o in out])
    curves = pd.concat([o[1] for o in out], ignore_index=True)
    surrB = {o[0]["symbol"]: o[2] for o in out}

    rows = rows.sort_values("B", ascending=False).reset_index(drop=True)
    rows.to_csv(os.path.join(HERE, "results", "summary.csv"), index=False)
    curves.to_csv(os.path.join(HERE, "results", "hq_curves.csv"), index=False)
    np.savez(os.path.join(HERE, "results", "shuffle_B.npz"), **surrB)

    cols = ["index", "country", "tier", "H1", "H2", "H4", "B", "B_lo", "B_hi",
            "D", "p_vs_shuffle", "r2_min"]
    pd.set_option("display.width", 220)
    print(rows[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    dev = rows.loc[rows.tier == "Developed", "B"]
    eme = rows.loc[rows.tier == "Emerging", "B"]
    try:
        from scipy.stats import mannwhitneyu
        u, pmw = mannwhitneyu(eme, dev, alternative="greater")
    except Exception:
        u, pmw = np.nan, np.nan
    print(f"\nB: developed mean {dev.mean():.4f} (n={len(dev)}), "
          f"emerging mean {eme.mean():.4f} (n={len(eme)})")
    print(f"Mann-Whitney U (emerging > developed): U={u}, p={pmw:.4f}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
