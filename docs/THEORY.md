# Multiscaling, Explained

The theory, the code, and the reasoning behind every step of
[`../README.md`](../README.md) — for anyone who wants the derivation, not just
the headline. A companion interactive version (with charts rendered from the
real `results/*.csv` files) exists alongside this document.

## 0. What this is for

The project asks one question: **do emerging stock markets behave differently
from developed ones at a structural level** — not just "more volatile", but
differently shaped across time scales — and if so, is that difference real,
or a statistical illusion produced by a handful of extreme days?

It follows a method published by **Prof. Tiziana Di Matteo** (who teaches on
the KCL MSc Complex Systems and Data Modelling programme this project targets)
and then does something most replications don't: it tries hard to break its
own headline result, and reports honestly where it breaks.

**Why it's useful, beyond the specific finding:**

- It's a worked example of a general, transferable skill — not trusting your
  own measurement until you've tried to fool it. See [§13](#13-where-this-builds-on-earlier-work).
- The scaffolding (synthetic processes with known answers, a null
  distribution, a deliberately adversarial check) is reusable far outside
  finance.
- The specific finding is genuinely interesting: markets do look more
  "multifractal" in emerging economies, but the entire effect turns out to
  live in the four or five largest price moves in each series — a fragile
  place for a finding to live. [§14](#14-a-second-estimator-and-the-di-matteo-precedent-it-lands-on)
  pushes on exactly that fragility with a structurally different estimator,
  and lands on a 2012 paper co-authored by Di Matteo herself that already
  made the same point on different data.

## 1. The raw material: prices, logs, and increments

Everything starts from a daily closing price for each index (~7,000 trading
days, back to 1998 where available). The project works with the **log
price**, `X(t) = ln P(t)`, not the raw price.

**Why logs, not prices.** A move from 100 to 110 and a move from 1000 to 1100
are "the same" 10% change but are 10 and 100 as raw numbers — not comparable.
Logs turn multiplicative change into additive change, so increments from very
different price levels and eras mean the same thing statistically.

The object everything else is a function of is the **increment over a gap of
τ days**:

$$\Delta X_\tau(t) = X(t+\tau) - X(t)$$

τ is a *time scale* — one day, five days, nineteen days. The whole project
reduces to: does the statistical spread of $\Delta X_\tau$ grow with τ the
same way, regardless of how "spread" is measured?

## 2. The structure function, and what `q` does

The central object, for every market and scale τ, is the **q-th order
structure function**:

$$S_q(\tau) = \langle |\Delta X_\tau(t)|^q \rangle_t$$

Take every τ-day move in the whole history, take its absolute size, raise it
to the power `q`, and average. `q` is an exponent — but changing it changes
*whose* behaviour the average reports on.

**Low q (e.g. 0.5) asks about the typical day. High q (e.g. 4) asks almost
exclusively about the single worst day** — raising a small number to the 4th
power shrinks it further, while raising a large number to the 4th power
inflates it, so a handful of extreme moves can dominate the whole average at
high q. This is why the project scans fifteen values of q from 0.5 to 4
(`QS = np.arange(0.5, 4.01, 0.25)`) rather than picking one — a single q only
tells you about one slice of "how extreme a move you've decided to care
about."

## 3. The Hurst exponent, and the generalised H(q)

The method assumes $S_q(\tau)$ scales as a power law in τ:

$$S_q(\tau) \sim \tau^{\,q \cdot H(q)}$$

Taking logs linearises this into exactly what the code fits:

$$\log S_q(\tau) = \text{const} + q \cdot H(q) \cdot \log \tau$$

So `H(q)` is the slope of a straight line on a log–log plot, divided by `q`.
Two reference points:

- **H = 0.5** — a plain random walk (Brownian motion). Moves at different
  scales are statistically independent — the "nothing interesting" baseline.
- **H > 0.5** — persistent (trending); **H < 0.5** — anti-persistent
  (mean-reverting).

**The shape that matters: monofractal vs. multifractal.** A *flat* `H(q)`
curve means one exponent describes price moves at every scale of extremity —
**monofractal**. A *falling* `H(q)` curve means ordinary days and extreme days
are governed by different exponents — **multiscaling**. See
`figures/fig1_controls.png` for the exact theoretical shapes (flat for
Brownian motion and fractional Brownian motion, gently falling for the
multifractal random walk).

## 4. Multiscaling strength, B

A whole curve is unwieldy to rank 22 countries by, so the project fits a
straight line through `H(q)` itself:

$$H(q) = A - B \cdot q$$

`B` is the sign-flipped slope. `B = 0` is flat — monofractal. **Larger B means
a more steeply falling H(q) — stronger multiscaling.** This single number is
what gets ranked, bootstrapped, and compared between developed and emerging
markets.

A second, simpler proxy — `D = H(1) − H(4)` — follows Di Matteo's own paper
directly and agrees with `B` on every market in this dataset; `B` is reported
because it uses all fifteen q-values rather than two.

## 5. Why the estimator is validated before it's trusted

Before `H(q)` is estimated on a single real index, it's run on **three
synthetic processes whose true H(q) is known exactly** (`scripts/02_controls.py`),
because a believable-looking estimate at n≈7,000 is not the same thing as a
correct one.

| Process | True H(q) | Tests | Recovered B | True B |
|---|---|---|---:|---:|
| Brownian motion | 0.5, flat | hallucinated multiscaling? | −0.0002 ± 0.0029 | 0 |
| fBm, H=0.7 | 0.7, flat | same, away from 0.5 | +0.0004 ± 0.0035 | 0 |
| fBm, H=0.3 | 0.3, flat | same, other side | +0.0002 ± 0.0020 | 0 |
| MRW, λ²=0.02 | falling, slope 0.01 | finds real multiscaling? | 0.0066 | 0.01 (66%) |
| MRW, λ²=0.04 | falling, slope 0.02 | same, double strength | 0.0132 | 0.02 (66%) |

This buys two things: **a null band** (on a genuine random walk at this
sample size, `B` exceeds `0.0046` only 5% of the time by chance), and **a
known, constant bias** — when multiscaling is really there, the estimator
only recovers ~66% of the true strength, at both tested levels, so every `B`
on real markets is a lower bound, not a point estimate.

*Where the 66% comes from:* the multifractal random walk is built from a
log-correlated Gaussian field via circulant embedding (an FFT trick), which
requires clipping a handful of negative eigenvalues that appear only from
floating-point rounding — standard practice, but it slightly throttles the
simulated volatility correlations, the traceable cause of the recovery gap.

## 6. The fitting range: τ_max and the q-range

**τ (1–19 trading days).** The log–log line is only straight over a finite
window. Too narrow and the slope is noisy; too wide and the fit may run
through a regime where the power law has broken down. Rather than defend
`τ_max=19` by assertion, the project re-runs the whole estimation at
`τ_max ∈ {10,15,19,25,30}` and checks the developed/emerging gap survives at
every one (it does, p between 0.0007 and 0.026 across the grid).

**q (0.5–4).** The sharper reason: **sample moments of a heavy-tailed
distribution become unreliable as q grows**, because a moment of order q is
dominated by the single largest observation once q exceeds roughly the number
of finite moments the underlying distribution actually has. This is exactly
why §7 exists, and why the project separately re-estimates everything
restricted to q ≤ 2.

## 7. Kurtosis, and "temporal structure" — dynamics, or just a fat distribution?

This is the conceptual fork the whole project turns on.

**Kurtosis** is the standardised fourth moment — how much of the variance
comes from rare, extreme values rather than everyday ones:

$$\text{Kurt}[X] = \frac{E[(X-\mu)^4]}{\sigma^4} - 3$$

A Normal distribution has excess kurtosis 0. The 22 markets here range from
**4.3** (Taiwan) to **85.9** (Malaysia — one historical event, the 1998
capital controls, dominates the entire 28-year sample).

**Why this threatens the headline finding.** The structure function at high q
is, by construction, dominated by the single most extreme move in the sample.
A market needs *no* real time-dependent structure — no volatility clustering,
no memory — to produce a falling H(q) curve, if its one-day return
*distribution* is fat-tailed enough on its own.

**The control that settles it.** The project builds a random walk from an
**iid Student-t distribution** — independent, identically distributed, zero
memory by construction — with a "degrees of freedom" parameter controlling
tail fatness. At `df=5` it produces `B ≈ 0.026`; at `df=3`, `B ≈ 0.069` —
squarely inside the range measured on the real markets. **A process with no
temporal structure at all reproduces the headline effect just from having a
fat-tailed distribution.** This is why "dynamics or tails?" can't be answered
from B alone.

## 8. Two different resampling tricks

**Block bootstrap (confidence interval).** Cuts returns into contiguous
blocks (length ≈ n^(1/3) ≈ 19 days), resamples whole blocks with
replacement, rebuilds a price path, re-estimates B — 500 times. Resampling
*blocks* rather than single days preserves short-range volatility clustering
within a block, which is what makes the spread an honest confidence interval.

**Shuffle surrogates (null distribution).** Randomly permutes individual
days, destroying all order while keeping the marginal distribution of
returns **exactly**. Whatever multiscaling survives this is definitionally
not about dynamics.

**The twist.** The standard sanity check in this literature is: shuffle the
data, and if multiscaling survives, it's a tail artefact. Here, **every one
of the 22 markets comes back with *more* multiscaling after shuffling.** That
looks damning — until the same test is run on the multifractal random walk,
whose multiscaling is real by construction. Shuffling *that* also increases
its measured B (0.0119 → 0.0366 at λ²=0.04). Destroying volatility clustering
isolates large moves from their neighbours while leaving the heavy marginal
distribution intact, inflating this statistic regardless of whether the
original structure was real. **The plain shuffle is a marginal-matched
benchmark, not a one-sided null** — which is why the project reports the
*excess* over the shuffled value, `dB = B − B_shuffled`, as a separate,
honest statistic.

## 9. The Mann–Whitney U test, and what p = 0.0013 actually says

With 13 developed and 9 emerging markets, "is emerging B bigger?" is a
two-group comparison. **Mann–Whitney** is used instead of a t-test because a
t-test assumes rough normality within each group — a bad assumption here, with
Malaysia as an extreme outlier. Mann–Whitney asks a purely rank-based
question: *pick one emerging B and one developed B at random — how often is
the emerging one larger?* Under the null (no real difference), that
probability is 50%.

**Reading p = 0.0013 correctly:** it means that *if* emerging and developed
markets had identically distributed multiscaling, a split this lopsided would
appear about 13 times in 10,000 random re-splits. It does **not** mean "0.13%
chance the null is true" — a common misreading — and it says nothing about
*why* the groups differ, which §7–8 exist to chase separately.

| Comparison | Emerging mean | Developed mean | p-value |
|---|---:|---:|---:|
| Raw B | 0.0281 | 0.0085 | 0.0013 |
| Raw B, excluding Malaysia | 0.0218 | 0.0085 | 0.0017 |
| Excess over shuffled null (dB) | −0.0220 | −0.0343 | 0.0097 |
| Low moments only, q ≤ 2 | 0.0210 | 0.0169 | 0.0545 |
| Low moments, excluding Malaysia | 0.0191 | 0.0169 | 0.1051 |

Read down the p-value column: decisive at the full q-range, still holds
(more weakly) net of the shuffle benchmark, and **stops being significant the
moment the analysis is restricted to reliably-estimable moments.** That
pattern across five deliberately adversarial versions — not any single
p-value — is the actual finding.

## 10. Controlling for kurtosis: a two-variable regression

$$B = \beta_0 + \beta_1 \cdot (\text{excess kurtosis}) + \beta_2 \cdot (\text{1 if emerging else 0}) + \varepsilon$$

If β₂ is still clearly non-zero after kurtosis is accounted for, the
developed/emerging split carries information kurtosis alone doesn't. It
survives at **t = 3.6, p = 0.002, R² = 0.80** — the t-statistic is the
coefficient divided by its own standard error (|t| above ~2 is the
conventional "probably not noise" bar here); R² = 0.80 means these two
predictors jointly explain 80% of the variation in B across the 22 markets.
Kurtosis does real work, but the emerging/developed distinction adds
something on top rather than being swallowed by it.

## 11. The results, assembled

See `figures/fig4_ranking.png` for all 22 markets ranked by B with
bootstrap intervals and the random-walk null line; `figures/fig3_hq_curves.png`
for measured H(q) curves per market; `figures/fig5_tails.png` for the
tails-vs-dynamics analysis; `figures/fig1_controls.png` for the estimator
validation. All numbers above are read directly from `results/summary.csv`,
`results/hq_curves.csv` and `results/controls.csv`.

## 12. Every function, mapped

Estimator in `src/ghe.py` (~150 lines); validation processes in
`src/synthetic.py`.

- **`structure_functions`** — computes $S_q(\tau)$ for every q at once (§2).
  Everything downstream calls this repeatedly on resampled or permuted data.
- **`_ols` / `_ols_many`** — plain straight-line least squares, by closed
  form rather than a library call, because the slope is needed for many rows
  at once. This is the literal mechanism behind both nested fits in §3–4.
- **`generalised_hurst`** — runs the structure functions, takes logs, calls
  `_ols_many`, divides by q to get H(q) (§3). The estimator validated in §5.
- **`multiscaling_strength` / `h1_minus_h4`** — fits `H(q) = A − Bq` and
  returns B; reads off H(1)−H(4) directly (§4).
- **`circular_block_bootstrap` / `block_bootstrap`** — resamples whole
  blocks of increments with replacement, rebuilds a price path, re-estimates
  B 500 times — the confidence intervals in §8 and §11.
- **`shuffle_surrogates`** — permutes individual days, rebuilds a price
  path, re-estimates B 500 times — the null distribution behind §7–9, and
  the function behind the "shuffling increases B" finding in §8.
- **`brownian_motion`, `fractional_bm` (via `_fgn_davies_harte`), `mrw`
  (via `_log_correlated_gaussian`)** — three generators of fake price series
  with mathematically known H(q): plain cumulative Gaussian noise;
  FFT-generated correlated Gaussian noise (Davies–Harte circulant embedding)
  for long-range dependence; Gaussian noise multiplied by a log-correlated
  random volatility field for genuine multifractality. These are the
  known-answer test cases in §5.

Second estimator in `src/mfdfa.py` (§14):

- **`profile` / `_segment_variances`** — builds the cumulative profile of the
  increments and, for a given scale s, the locally-detrended variance of
  every non-overlapping segment (both directions).
- **`fluctuation_function` / `mfdfa`** — F_q(s) for every q and scale, then
  the log-log OLS slope that gives h(q). Same return shape as
  `generalised_hurst`, so `multiscaling_strength`, `h1_minus_h4` and
  `singularity_spectrum` all work on it unchanged.
- **`block_bootstrap` / `shuffle_surrogates`** — the same two resampling
  checks as `ghe.py`, with `mfdfa()` in place of `generalised_hurst()`.

## 13. Where this builds on earlier work

The method is new; the *moves* are not.

- **Exeter group project (CMIP6 sea-ice modelling).** 35 climate models
  filtered to 20 by a confidence-interval test, combined by model selection
  and an emergent-constraint regression, rather than trusting one model's
  projection. Same move as §6: don't trust one estimate — build a band, filter
  by a statistical criterion, check the answer survives across choices.
- **Actuarial work (stochastic pension modelling).** The object of interest
  in multi-decade liability projection is the *tail* of the propagated
  distribution, not the central forecast. Same move as §7: the headline
  result lives in the tail moments (q>2) of the return distribution, which is
  exactly why it needs the same suspicion an actuary gives a tail-driven
  number.
- **Lincoln City FC (tracking & recruitment data).** With 76% of opposition
  ratings missing, the recurring problem was distinguishing real structure in
  high-frequency, many-agent data from measurement artefact. Same move as
  §7–8: the iid Student-t control and the shuffle surrogate exist to answer
  exactly that question for markets instead of match data.

Three different domains, one discipline: build something with a known
answer, check the method recovers it, then ask what would make the pattern
fake before asking what it means if it's real.

## 14. A second estimator, and the Di Matteo precedent it lands on

Every check in §5–§10 stays inside one method: the same structure-function
estimator, re-run on resampled data, restricted q-ranges, different τ_max.
That tests the estimator's *stability*, not whether its answer is an
artefact of the estimator itself. A sharper test is a structurally different
estimator, with different failure modes by construction, applied to the same
22 series. If it agrees, that is real corroboration. If it doesn't, the
disagreement needs a reason, not a shrug — so this section chases the reason.

### 14.1 MFDFA, and choosing its scale range the way τ_max was chosen

[Multifractal Detrended Fluctuation Analysis](../src/mfdfa.py) (Kantelhardt
et al., 2002) estimates the same h(q) a structurally different way: build the
cumulative profile of the increments, cut it into non-overlapping segments at
scale s, remove a local polynomial trend from each segment, and take the
q-th moment of the detrended local variance; h(q) is the log-log slope of
that against s. Where the structure function asks "how does a lag-τ
difference scale?", MFDFA asks "how does a detrended τ-day *window's*
variance scale?" — related, but not the same question.

Polynomial order and scale range (`order`, `s_min`, `s_max`) were screened
against the same synthetic processes as §5, not guessed:

- **Order.** Order-2 (quadratic) detrending made the null bias on
  anti-persistent fBm *worse* at every scale range tried, not better —
  short segments with an extra free parameter overfit. Order 1 (linear
  detrending) was kept.
- **s_min.** The textbook rule of thumb is s_min ≈ 10 for linear detrending.
  Empirically, at this sample size the null bias on fBm(H=0.3) keeps
  shrinking well past that, down to s_min = 25.
- **s_max.** s_max = 150 is where the MRW positive controls' recovered B is
  closest to the analytic value (≈108%, vs the structure function's
  documented 66–72%); larger s_max keeps shrinking the null band but starts
  overshooting the true signal.

Final choice: **order=1, s_min=25, s_max=150** (16 log-spaced scales),
re-validated at 300 replicas each:

| process | B, structure fn | B, MFDFA | true B |
|---|---:|---:|---:|
| Brownian motion | 0.0004 ± 0.0029 | 0.0012 ± 0.0048 | 0 |
| fBm, H=0.7 | 0.0002 ± 0.0035 | 0.0009 ± 0.0059 | 0 |
| fBm, H=0.3 | 0.0000 ± 0.0021 | 0.0020 ± 0.0027 | 0 |
| MRW, λ²=0.02 | 0.0072 (72%) | 0.0113 (113%) | 0.0100 |
| MRW, λ²=0.04 | 0.0137 (69%) | 0.0204 (102%) | 0.0200 |

MFDFA's own null band is wider than the structure function's (95th
percentile of B under a random walk: **0.0086** vs 0.0046) — a real,
quantified property of the estimator at this sample size, not a guess, and
the reason a market-level comparison needs MFDFA's own null, not a reused
one. `scripts/02b_controls_mfdfa.py` reproduces this table and saves the
null band to `results/null_B_mfdfa_brownian.npy`.

### 14.2 It does not replicate the headline, anywhere on the grid

Run on all 22 markets with the same rigor as §8–§11 — block bootstrap,
shuffle surrogates, and an `s_max` grid in place of the τ_max grid
(`scripts/03b_estimate_mfdfa.py`, `results/summary_mfdfa.csv`):

| | developed mean B | emerging mean B | Mann–Whitney p (emerging > developed) |
|---|---:|---:|---:|
| structure function | 0.0085 | 0.0281 | 0.0013 |
| MFDFA | **0.0213** | **0.0104** | **0.9903** |

The gap is not just smaller under MFDFA — it points the other way, and does
so at **every** point on the robustness grid (`s_max ∈ {80,100,150,200,300}`,
developed mean always above emerging mean). Across the 22 markets, `B`
(structure function) and `B` (MFDFA) correlate at **r = −0.13**: 12 of 22
markets land on opposite sides of their own method's median. Malaysia, the
single most extreme structure-function outlier (B=0.0786, 1998 capital
controls, rank 1 of 22), is unremarkable under MFDFA (B=0.0161, rank 12 of
22). See `results/sf_vs_mfdfa.csv` and `figures/fig6_mfdfa_crosscheck.png`.

### 14.3 Two tempting explanations, both ruled out

**Not simple outlier sensitivity.** Replacing Malaysia's single worst day
(a −24% move) with its series' median |return| drops the structure-function
B to 65% of its original value — consistent with that estimator being
outlier-driven, as expected. The *same* edit moves the MFDFA B to **306%**
of its original value. Whatever is driving the disagreement, it is not "MFDFA
just dilutes one bad day" — that story predicts the wrong sign.

**Not a power shortfall.** MFDFA's bootstrap CIs are ~46% wider than the
structure function's (median width 0.0499 vs 0.0341), but it distinguishes
*slightly more* individual markets from its own zero-multiscaling null (16
of 22, vs 15 of 22 for the structure function), with nearly identical median
|B| (0.0172 vs 0.0174). MFDFA is not simply a noisier version of the same
measurement that washes the signal out — it is picking out a different
pattern with comparable confidence.

### 14.4 What actually explains it — and where this was already found

§7 already showed the structure-function method manufactures large B from an
iid Student-t random walk with **zero temporal structure** — B≈0.0264 (df=5)
to B≈0.0689 (df=3), purely from fat tails. That control had never been run
through MFDFA. Doing so now: at the same degrees of freedom, **MFDFA
recovers only 30–50% of that spurious signal** (df=5: 0.0078 vs the
structure function's 0.0264 [30%]; df=3: 0.0346 vs 0.0689 [50%] —
`results/controls_iid_t_mfdfa.csv`). MFDFA is demonstrably, measurably less
fooled by fat tails alone than the structure function is.

That one property closes the loop. §7–9 already established that much of
the structure-function B is a tail artefact rather than real dynamics. An
estimator that is specifically less susceptible to that artefact, applied to
markets whose tail-fatness varies enormously (excess kurtosis 4.3 to 85.9)
and does not track the developed/emerging split cleanly, should be expected
to find a different ranking — and does. This is not two estimators randomly
disagreeing; it is the same mechanism §7 already named, now demonstrated by
its absence under a method built differently.

**This also isn't a new worry — Di Matteo found it first.** Barunik, Aste,
Di Matteo and Liu, "Understanding the source of multifractality in financial
markets" (*Physica A* 391(17), 2012), run exactly the shuffle test of §8 on
stock indices, FX rates and bond yields, and reach the same two conclusions
this project reached independently on a different, more recent dataset:
shuffling (which destroys dynamics, keeps the marginal) *increases* the
measured multiscaling almost everywhere — their words, "empirical data show
even larger degree of multifractality than original ones" — and their
closing section, titled "*Does all multifractality come from the
distribution?*", concludes "most of the multifractality is caused solely by
the broad, fat-tailed distribution of the returns." §8 of this document
rediscovered their headline finding on 22 indices spanning 1998–2026, using
Di Matteo's own 2007 structure-function method, without having read that
2012 paper first. What this section adds on top of it — an independent
estimator (MFDFA, which that 2012 paper only cites in passing as a
generalisation of DFA, never runs) and a direct test showing *why* it
disagrees — is the part that goes beyond what was already published.

### 14.5 What this changes, and what it doesn't

It does not overturn the headline as stated: under T. Di Matteo's own
structure-function method, extensively validated on its own terms, emerging
markets multiscale more strongly than developed ones in this dataset,
B=0.0281 vs 0.0085, p=0.0013. That remains true and is reported as such. Two
methods probing the data differently — different scale regimes (1–19 days vs
25–150-day detrended windows), different functionals (pointwise |Δ|^q vs
segment-averaged local variance) — are not required to agree for both
measurements to be real.

What changes is how the fragility already flagged in §9 (the gap stops being
significant once high moments are excluded) should be read: it is not an
isolated quirk of one estimator's treatment of q>2 moments. A second,
independent estimator — chosen specifically because it is less fooled by the
exact artefact §7 already identified — does not find this effect at all, in
either direction a real but weaker dynamics-driven gap would predict. Read
together with §14.4's precedent, the most defensible position is the one
Barunik, Aste, Di Matteo and Liu already reached in 2012: the multiscaling
measured here is predominantly a signature of the return distribution, not
of time-dependent structure, and the developed/emerging separation in B is
better understood as a difference in tail behaviour than in market dynamics.

## References

- T. Di Matteo, "Multi-scaling in finance", *Quantitative Finance* 7(1), 2007.
- T. Di Matteo, T. Aste, M. M. Dacorogna, "Long-term memories of developed and
  emerging markets", *Journal of Banking & Finance* 29, 2005.
- E. Bacry, J. Delour, J.-F. Muzy, "Multifractal random walk",
  *Physical Review E* 64, 026103, 2001.
- J. Baruník, L. Kristoufek, "On Hurst exponent estimation under heavy-tailed
  distributions", *Physica A* 389, 2010.
- J. W. Kantelhardt, S. A. Zschiegner, E. Koscielny-Bunde, S. Havlin,
  A. Bunde, H. E. Stanley, "Multifractal detrended fluctuation analysis of
  nonstationary time series", *Physica A* 316(1), 87–114, 2002.
- J. Baruník, T. Aste, T. Di Matteo, R. Liu, "Understanding the source of
  multifractality in financial markets", *Physica A* 391(17), 4234–4251,
  2012. — the direct precedent for §14: the same shuffle-increases-B
  paradox and the same "mostly fat tails" conclusion, on different data.
