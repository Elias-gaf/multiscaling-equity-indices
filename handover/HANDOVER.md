# Handover — Multiscaling in International Equity Indices

**Status: analysis complete. Not yet published. ~20–30 minutes of work left.**

## What this project is

Estimates the generalised Hurst exponent H(q) for 22 national equity indices
(1998–2026, ~7,000 daily closes each), using T. Di Matteo's structure-function
method, to test whether emerging markets multiscale more strongly than
developed ones — and how much of that difference is real dynamics versus an
artefact of fat tails. A second, structurally independent estimator (MFDFA)
is then run as a cross-check on that exact question.

**Headline result:** the textbook finding replicates (multiscaling strength
more than 3x larger in emerging markets, p = 0.0013) and survives controlling
for excess kurtosis. **But** the separation is carried entirely by moments
q > 2 — the moments a 7,000-point sample estimates least reliably. Restricted
to q ≤ 2, the gap shrinks by two-thirds and stops being significant. An iid
Student-t random walk with zero temporal structure produces multiscaling of
the same magnitude as the real markets.

**Updated — a second estimator (MFDFA) was cross-checked against the
headline and does not replicate it: it finds the gap running the other way**
(developed B=0.0213 vs emerging B=0.0104, p=0.99), robustly across a scale-
range grid, and the reason is traceable — MFDFA is demonstrably less
susceptible to the same fat-tail artefact §7 already flagged as the main
threat. This is not a new concern: Barunik, Aste, Di Matteo and Liu (*Physica
A*, 2012) reached the same "mostly fat tails, not dynamics" conclusion on
different data, by the same shuffle test, years before this project
independently reran it. See `docs/THEORY.md` §14 for the full chain of
reasoning. The headline is reported as-is, under its own method, with this
caveat now stated as sharply as the evidence supports — not softened, not
hidden.

Full write-up with method, validation, robustness checks and limitations is
in `../README.md`. A five-level walkthrough, simplest to most technical, is
in `../docs/EXPLAINER.md`. A critical evaluation (precedent comparison,
economic reasoning for/against, data sourcing, concept explainers with
corrected figures) is in `../docs/EVALUATION.md`. Ready-to-paste lines for
the CV / supporting statement / an email to Prof. Di Matteo are in
`../note/for_application.md`.

## Why this exists

Written as evidence of mathematical/statistical depth for a 2027 reapplication
to KCL's MSc Complex Systems and Data Modelling, after a 2026 rejection
attributed to a weak undergraduate pure-mathematics transcript. It follows
Prof. Tiziana Di Matteo's own published method (she is faculty on the target
programme), and is deliberately built like a mini-paper: a validated
estimator, a replication, and one honest negative result — not a from-scratch
novel claim.

## What's actually left before this counts as "done"

1. **Push to GitHub.** No git repo exists yet in this folder (`git init` was
   never run). `.gitignore` already excludes `__pycache__/`, `*.pyc` and
   `.DS_Store`, so it's clean to commit as-is.
2. **Fix the CV.** `Elias_Watts_CV_KCL_ComplexSystems.html` (parent of this
   project's parent folder) still has the placeholder: *"[ONE-LINE FINDING —
   write this once the analysis is done. Add repository URL.]"* — the finding
   exists; swap in the CV line from `note/for_application.md` and the repo
   URL once it's live.
3. **Optional — email Di Matteo.** Her 11 Sep reply implicitly invited this:
   come back once there's a result. The question drafted in
   `note/for_application.md` is now built around the MFDFA cross-check and
   the 2012 Barunik/Aste/Di Matteo/Liu paper it lands on, rather than the
   original q≤2 fragility question alone.

## File map

```
multifractal/
├── README.md              full write-up — read this first
├── requirements.txt       numpy, pandas, scipy, matplotlib, yfinance
├── note/
│   └── for_application.md ready-to-paste CV line, statement lines, Di Matteo question
├── docs/
│   ├── THEORY.md           full derivation, including §14 (the MFDFA cross-check)
│   ├── EXPLAINER.md        the same project, five levels of depth, ELI5 to full formalism
│   └── EVALUATION.md       precedent comparison, economic reasoning, data source, concept explainers
├── src/
│   ├── ghe.py              the structure-function estimator (~150 lines) — core of the project
│   ├── mfdfa.py             the second, independent estimator (§14)
│   └── synthetic.py        validation processes (Brownian motion, fBm, multifractal random walk)
├── scripts/
│   ├── 01_download.py      fetch + clean 22 indices → data/clean/
│   ├── 02_controls.py      validate the estimator against known processes → figures/fig1
│   ├── 02b_controls_mfdfa.py validate MFDFA against the same processes, choose its scale range
│   ├── 03_estimate.py      H(q) for all 22 markets, bootstrap + shuffle surrogates
│   ├── 03b_estimate_mfdfa.py same, by MFDFA, plus the s_max robustness grid
│   ├── 04_tails.py         tests whether fat tails explain the multiscaling gap
│   └── 05_figures.py       figures 2–5
├── data/
│   ├── raw/                22 index CSVs as downloaded
│   └── clean/              cleaned log-prices + a data-quality report
├── results/                controls.csv, hq_curves.csv, summary.csv, bootstrap outputs,
│                           controls_mfdfa.csv, controls_iid_t_mfdfa.csv, summary_mfdfa.csv,
│                           sf_vs_mfdfa.csv, null_B_mfdfa_brownian.npy
└── figures/                fig1_controls … fig5_tails, fig1b_mfdfa_controls, fig6_mfdfa_crosscheck,
                           fig7_Sq_vs_lntau, fig8_legendre_transform
```

## Reproducing from scratch

```bash
pip install -r requirements.txt
python scripts/01_download.py        # fetch + clean, writes data/clean/ and a quality report
python scripts/02_controls.py        # validate the estimator, figure 1
python scripts/02b_controls_mfdfa.py # validate MFDFA, choose its scale range, figure 1b
python scripts/03_estimate.py        # H(q) for 22 markets, bootstrap + surrogates (~1 min)
python scripts/03b_estimate_mfdfa.py # same, by MFDFA, plus the s_max robustness grid
python scripts/04_tails.py           # tails vs dynamics
python scripts/05_figures.py         # figures 2-5
```

## Known limitations (stated honestly in README, repeated here for convenience)

- Bootstrap intervals (circular block, block ≈ 19 days) are honest about
  width but shifted upward, not centred — same mechanism that inflates the
  shuffle-test comparison.
- "Developed vs emerging" is MSCI's classification — a coarse, partly
  political proxy. Hong Kong (classified developed) behaves like the
  emerging group; Malaysia is an extreme outlier (1998 capital controls) and
  every result is reported with and without it.
- The synthetic multifractal random walk's covariance is biased 1–3% low
  (negative eigenvalue clipping in `_log_correlated_gaussian`), which is part
  of why the estimator recovers 66% rather than 100% of the true signal —
  reported B values are lower bounds, not point estimates.
- MFDFA's own null distribution is ~65% wider than the structure function's
  at this sample size (B's 95th percentile under a random walk: 0.0086 vs
  0.0046), and its real-market bootstrap CIs are correspondingly ~46% wider.
  It is not underpowered overall (it distinguishes slightly more individual
  markets from its null than the structure function does), but any single
  market's MFDFA B should be read with that wider uncertainty in mind.
- §14's MFDFA cross-check does not, on its own, prove the structure-function
  headline wrong — the two estimators probe different scale regimes and are
  not required to agree for both to be real. What it adds is independent
  evidence, via a mechanism this project can name (differential sensitivity
  to fat tails), for the same conclusion Barunik, Aste, Di Matteo and Liu
  (2012) already reached on different data.
