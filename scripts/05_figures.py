"""Figures 2-5.  (Figure 1, the controls, is produced by scripts/02_controls.py.)"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "src"))
from ghe import structure_functions, generalised_hurst

FIG = os.path.join(HERE, "figures")
DEV, EME = "#2b6cb0", "#c05621"
plt.rcParams.update({"figure.dpi": 170, "font.size": 9})

s = pd.read_csv(os.path.join(HERE, "results", "summary.csv"))
curves = pd.read_csv(os.path.join(HERE, "results", "hq_curves.csv"))
lp = pd.read_csv(os.path.join(HERE, "data", "clean", "logprices.csv"))
X = {k: g["logprice"].to_numpy() for k, g in lp.groupby("symbol")}

# ---- fig 2: the raw scaling, so the reader can see the power laws ----------
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
taus = np.arange(1, 20)
qs = np.array([1., 2., 3., 4.])
for ax, sym in zip(axes, ["^GSPC", "XU100.IS"]):
    S = structure_functions(X[sym], qs, taus)
    r = generalised_hurst(X[sym], qs)
    for i, q in enumerate(qs):
        y = np.log(S[i]) - np.log(S[i][0])
        ax.plot(np.log(taus), y, "o", ms=3.5, color=plt.get_cmap("viridis")(i / 3.4))
        b = np.polyfit(np.log(taus), np.log(S[i]), 1)
        ax.plot(np.log(taus), np.polyval(b, np.log(taus)) - np.log(S[i][0]), "-",
                lw=1.1, color=plt.get_cmap("viridis")(i / 3.4),
                label=f"$q$={q:.0f}: $qH(q)$={q*r['H'][i]:.3f}, $H$={r['H'][i]:.3f}")
    name = s.loc[s.symbol == sym, "index"].iloc[0]
    ax.set_title(f"{name} ({s.loc[s.symbol==sym,'tier'].iloc[0]})")
    ax.set_xlabel(r"$\log \tau$"); ax.set_ylabel(r"$\log S_q(\tau)$  (offset)")
    ax.legend(fontsize=7); ax.grid(alpha=.3)
fig.suptitle("Structure functions are straight in log-log: the scaling assumption holds", y=1.0)
fig.tight_layout(); fig.savefig(f"{FIG}/fig2_scaling_fits.png"); plt.close(fig)

# ---- fig 3: H(q) curves ----------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
for ax, tier, col in [(axes[0], "Developed", DEV), (axes[1], "Emerging", EME)]:
    sub = curves[curves.tier == tier]
    for sym, g in sub.groupby("symbol"):
        ax.plot(g.q, g.H, color=col, alpha=.55, lw=1.1)
    m = sub.groupby("q")["H"].mean()
    ax.plot(m.index, m.values, color="k", lw=2.2, label=f"{tier} mean")
    ax.axhline(0.5, color="0.4", ls="--", lw=1, label="random walk, $H$=0.5")
    ax.set_xlabel("$q$"); ax.set_title(f"{tier}  (n={sub.symbol.nunique()})")
    ax.legend(fontsize=8); ax.grid(alpha=.3)
axes[0].set_ylabel("$H(q)$")
fig.suptitle("Generalised Hurst exponent, 22 indices, 1998-2026 daily closes", y=1.0)
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_hq_curves.png"); plt.close(fig)

# ---- fig 4: ranking with bootstrap CIs -------------------------------------
d = s.sort_values("B")
fig, ax = plt.subplots(figsize=(7.2, 6.4))
y = np.arange(len(d))
cols = [EME if t == "Emerging" else DEV for t in d.tier]
ax.errorbar(d.B, y, xerr=[d.B - d.B_lo, d.B_hi - d.B], fmt="none",
            ecolor="0.65", lw=1.1, capsize=2)
ax.scatter(d.B, y, c=cols, s=34, zorder=3)
null95 = float(np.percentile(np.load(os.path.join(HERE, "results", "null_B_brownian.npy")), 95))
ax.axvline(0, color="0.3", lw=1)
ax.axvline(null95, color="crimson", ls="--", lw=1.2,
           label=f"95th pct under a random walk ({null95:.4f})")
ax.set_yticks(y); ax.set_yticklabels([f"{r['index']} ({r.country})" for _, r in d.iterrows()],
                                     fontsize=7.5)
ax.set_xlabel("multiscaling strength $B$  (slope of $H(q)$), 95% block-bootstrap CI")
ax.set_title("Emerging markets multiscale more strongly\n"
             f"gap {d[d.tier=='Emerging'].B.mean()-d[d.tier=='Developed'].B.mean():+.4f}, "
             "Mann-Whitney $p$=0.0013")
h = [plt.Line2D([], [], marker="o", ls="", color=DEV, label="Developed"),
     plt.Line2D([], [], marker="o", ls="", color=EME, label="Emerging")]
ax.legend(handles=h + ax.get_legend_handles_labels()[0], fontsize=8, loc="lower right")
ax.grid(alpha=.3, axis="x")
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_ranking.png"); plt.close(fig)

# ---- fig 5: how much of this is tails? -------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(13, 4))
ax = axes[0]
tdf = pd.read_csv(os.path.join(HERE, "results", "controls_iid_t.csv"))
ctl = pd.read_csv(os.path.join(HERE, "results", "controls.csv"))
bm_B = float(ctl.loc[ctl.process.str.startswith("Brownian"), "B_mean"].iloc[0])
ax.bar(["Gaussian\nrandom walk"] + [f"iid $t_{{{d}}}$" for d in tdf.dof],
       [bm_B, *tdf.B_mean], yerr=[0, *tdf.B_sd], capsize=3,
       color=["0.6"] + ["0.35"] * len(tdf))
ax.axhspan(s[s.tier == "Developed"].B.min(), s[s.tier == "Developed"].B.max(),
           color=DEV, alpha=.2, label="range, developed markets")
ax.axhspan(s[s.tier == "Emerging"].B.min(), s[s.tier == "Emerging"].B.max(),
           color=EME, alpha=.2, label="range, emerging markets")
ax.set_ylabel("$B$"); ax.legend(fontsize=7.5)
ax.set_title("iid fat tails alone manufacture $B$\nof the size seen in real markets")
ax.grid(alpha=.3, axis="y")

ax = axes[1]
ax.scatter(s.B_shuffled, s["B"], c=[EME if t == "Emerging" else DEV for t in s.tier], s=32)
lim = [0, max(s.B_shuffled.max(), s["B"].max()) * 1.08]
ax.plot(lim, lim, "k--", lw=1, label="observed = shuffled")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("$B$ of shuffled returns (same marginal, no dynamics)")
ax.set_ylabel("$B$ observed")
ax.set_title("Every market falls BELOW the line:\nthe shuffle is not a one-sided null")
ax.legend(fontsize=8); ax.grid(alpha=.3)

ax = axes[2]
w = 0.35
for i, (tier, col) in enumerate([("Developed", DEV), ("Emerging", EME)]):
    v = s[s.tier == tier]
    ax.bar([0 + (i - .5) * w, 1 + (i - .5) * w], [v["B"].mean(), v.B_lowq.mean()],
           width=w, color=col, label=tier)
ax.set_xticks([0, 1]); ax.set_xticklabels(["all $q \\leq 4$\n$p$=0.0013", "$q \\leq 2$ only\n$p$=0.055"])
ax.set_ylabel("mean $B$")
ax.set_title("The gap lives in the high moments,\nwhich are the least reliably estimated")
ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")
fig.tight_layout(); fig.savefig(f"{FIG}/fig5_tails.png"); plt.close(fig)
print("wrote fig2-fig5")
