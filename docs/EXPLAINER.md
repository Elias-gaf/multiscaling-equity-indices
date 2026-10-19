# This project, five times

The same project — what it asks, what it does, what it found — explained five
times at increasing depth. Level 1 assumes nothing. Level 5 assumes you could
reproduce the code from the description alone. Skip to whichever level
matches what you already know; each one is self-contained, so you don't have
to read the ones before it, though each does assume slightly more than the
last.

## Does this align with Di Matteo's own findings, and does it build on them?

Short answer: **yes, to both of her relevant papers, and the "building on"
part turned out to be the most interesting thing in the project.**

- **T. Di Matteo, "Multi-scaling in finance"** (*Quantitative Finance*, 2007)
  and **Di Matteo, Aste, Dacorogna, "Long-term memories of developed and
  emerging markets"** (*Journal of Banking & Finance*, 2005) established the
  method used here (the generalised Hurst exponent `H(q)` from structure
  functions) and the headline pattern it finds: markets at different stages
  of development scale differently, with less-developed markets showing
  stronger multiscaling. This project replicates that pattern on a newer,
  independent dataset (22 indices, 1998–2026) — emerging markets multiscale
  more than three times as strongly as developed ones here too (`B`=0.0281
  vs 0.0085, p=0.0013).
- **Barunik, Aste, Di Matteo, Liu, "Understanding the source of
  multifractality in financial markets"** (*Physica A*, 2012) asked whether
  that multiscaling is really about market *dynamics* (volatility clustering,
  memory) or just about markets having fat-tailed return distributions, with
  no dynamics required. Their answer, from shuffling the data and watching
  multiscaling go *up* rather than down: mostly fat tails. **This project
  reran that exact test on the new dataset without having read their paper
  first, and got the same answer** — every one of the 22 markets' shuffled
  surrogates comes back more "multifractal," not less.
- What goes beyond both papers: this project then built a second,
  structurally different estimator (MFDFA) specifically to pressure-test
  that conclusion, validated it independently, and showed it is measurably
  less fooled by the fat-tail artefact than the original method — and that
  under it, the developed/emerging gap disappears (and weakly reverses)
  entirely. As far as this project's own literature search found, that
  specific cross-estimator test hasn't been run for this specific claim
  before. See [§14 of THEORY.md](THEORY.md#14-a-second-estimator-and-the-di-matteo-precedent-it-lands-on)
  for the full chain.

So: the project doesn't just replicate Di Matteo's headline — it
independently rediscovers her own later, more skeptical finding, and then
pushes one step past it.

---

## Level 1 — No background assumed

Imagine a chocolate-chip cookie. If the chips are spread out perfectly
evenly, a tiny crumb and a huge slice both have "the same kind" of
chocolate-to-dough ratio — one simple rule describes the whole cookie. Now
imagine a messier cookie: dense chocolate clumps in some spots, bare dough in
others. A tiny crumb from the clump and a tiny crumb from the bare patch look
totally different. One rule can't describe this cookie — you need to know
*how the pattern changes depending on where, and how closely, you look.*

Stock markets have prices that move up and down every day. This project asks
a cookie-style question about those price moves: **does a stock market's
"roughness" look the same whether you're looking at an ordinary day or a
historic crash — or does it change?**

It studied 22 countries' stock markets (the US, UK, Germany, Brazil, India,
Turkey and 16 others) over about 28 years, and asked whether markets in
*developing* economies behave more like the messy cookie than markets in
*developed* economies.

**What it found:** yes — at first. Emerging markets did look much "messier"
than developed ones, by a clear and statistically convincing margin. But
digging in, that messiness turned out to come almost entirely from a
handful of the single biggest crashes in each market's history — not from
some everyday, persistent difference in how these markets behave. And when
the question was asked a second way, using a completely different
measuring method built specifically to not be fooled by a few huge days,
**the difference between developed and emerging markets vanished.**

That's the whole project: a result that looks solid, a project that tried
hard to break its own result, and an honest report of exactly where it
broke.

---

## Level 2 — Casual technical

Every stock index here is turned into a **log price** series — day-by-day
levels, not raw dollar values, so a 10% move means the same thing whether
the index is at 100 or 10,000.

The key idea is **scale**: you can measure how much a price typically moves
over 1 day, 5 days, or 19 days. In a "boring," efficient market, the pattern
of those moves is self-similar — zoom in or out and the statistical shape
looks the same, just stretched. That single stretch factor is called the
**Hurst exponent, H**. H=0.5 is a plain random walk (tomorrow is unrelated
to today); H>0.5 means moves tend to keep going the same direction;
H<0.5 means they tend to reverse.

Now add a dial called **q**. Instead of asking "how much do prices move on
average," q lets you ask "how much do prices move if I only care about the
*k*-th most extreme kind of day." Low q (like 0.5) describes ordinary days.
High q (like 4) is dominated almost entirely by the single worst days in the
whole sample, because raising a small number to a high power shrinks it
further while raising a large number to a high power blows it up.

If a market's H stays the same no matter what q you pick, it's
**monofractal** — one rule, like the even cookie. If H *falls* as q rises —
meaning ordinary days and extreme days behave differently — that's
**multiscaling**, like the messy cookie. The project boils the whole falling
curve down to one number, `B` (how steeply it falls), and ranks all 22
markets by it.

**Headline:** emerging markets' `B` averaged 0.0281; developed markets'
averaged 0.0085 — more than three times larger, a gap unlikely to be chance
(about a 1-in-770 shot, if there were really no difference).

**The catch:** that gap disappears almost entirely if you only use the
"ordinary day" end of the q dial (q≤2) and ignore the crash-dominated end.
And a totally artificial process with *zero* memory — just random numbers
drawn from a distribution with fat tails, no connection from one day to the
next — produces a `B` just as large as the real markets. So a lot of what
looked like "emerging markets behave differently over time" might really be
"emerging markets have had a few more enormous one-off days."

To check that properly, the project built a *second* way of measuring the
same thing (called MFDFA) that is specifically less swayed by a handful of
huge days, and reran the whole comparison. Under that method, the
developed/emerging gap isn't just smaller — it's gone, and if anything
points the other way. That's the honest ending: a textbook finding that
replicates cleanly, and a second look that takes much of the air out of it.

---

## Level 3 — Formal (undergraduate level)

**The object of study.** For log price `X(t)`, the project works with the
increment over a gap of τ days, `ΔX_τ(t) = X(t+τ) − X(t)`, and the q-th
order **structure function**:

```
S_q(τ) = ⟨ |ΔX_τ(t)|^q ⟩_t
```

Assuming `S_q(τ) ~ τ^(qH(q))`, the **generalised Hurst exponent** `H(q)` is
the slope of `log S_q(τ)` against `log τ`, divided by `q`, fit by OLS over
τ ∈ [1, 19] days across q ∈ [0.5, 4]. `H(q)` constant in q is monofractal;
`H(q)` falling in q is multiscaling. The project fits a straight line
`H(q) = A − Bq` and reports the slope `B` (sign-flipped) as the headline
statistic — larger `B` means a more steeply falling curve, i.e. stronger
multiscaling. This is T. Di Matteo's own structure-function method
(*Quantitative Finance*, 2007), applied here to 22 national indices,
1998–2026, ~7,000 daily closes each.

**Why validate before trusting.** An estimator can look reasonable on real
data and still be biased or noisy in ways you can't see from the data alone.
So before touching a single real index, the estimator is run 300 times on
three processes with *analytically known* `H(q)`: Brownian motion (`H`=0.5,
flat — the "nothing interesting" case), fractional Brownian motion at two
other H values (also flat, just shifted), and the Bacry–Delour–Muzy
multifractal random walk (falling `H(q)`, by construction, at two known
strengths). This buys a **null band** (how large can `B` get on a pure
random walk, just from noise, at this sample size?) and a **known bias**
(the estimator recovers about 66% of the true multiscaling strength on the
one process built to have real multiscaling — so every reported `B` is a
lower bound, not a point estimate).

**The pipeline**, matching the actual scripts:

1. `01_download.py` — fetch and clean 22 indices into log prices.
2. `02_controls.py` — the validation above; also gives the null band.
3. `03_estimate.py` — `H(q)`, `B`, bootstrap confidence intervals (resample
   blocks of returns, not individual days, to preserve volatility
   clustering), and shuffle surrogates (permute individual days, destroying
   all time-order while keeping the return distribution exactly) for all 22
   markets.
4. `04_tails.py` — tests whether fat tails alone (an iid Student-t random
   walk — literally no time structure) can reproduce the observed `B`.
5. `05_figures.py` — the write-up's figures.

**The replication.** The developed/emerging gap in `B` (0.0085 vs 0.0281,
p=0.0013 by Mann–Whitney) matches the qualitative pattern in Di Matteo's own
2005/2007 papers, which linked scaling behaviour to a market's stage of
development.

**The complication.** Shuffling destroys time-order but keeps the
distribution of returns exactly. If multiscaling survived shuffling, that
would normally be read as "it's just the distribution, not real dynamics."
Here, *every* market's shuffled version has *higher* `B` than the original —
which sounds like the opposite conclusion, except the same thing happens
when you shuffle the synthetic MRW process, whose multiscaling is real by
construction. So the shuffle isn't a clean one-sided null here; it's a
marginal-matched benchmark. This exact pattern — and the same interpretation
— was already published by Barunik, Aste, Di Matteo and Liu in 2012.

---

## Level 4 — Technical depth

**Bootstrap mechanics.** `block_bootstrap` resamples whole blocks of
increments (circular block bootstrap, block length ≈ max(20, n^(1/3)) ≈ 20
days for n≈7,000), rebuilds a price path by cumulative sum, and
re-estimates `B` 500 times. Block resampling (not single-day resampling)
preserves short-range volatility clustering *within* a block, which is what
makes the resulting interval meaningful — but it also breaks correlation
*across* block boundaries, which is documented as inflating `B` somewhat
(the project reports the intervals as honest about width, not centring).

**The singularity spectrum.** `H(q)` has a geometric dual via the Legendre
transform: with `τ(q) = qH(q) − 1` (the mass exponent), the Hölder exponent
`α(q) = dτ/dq` and the singularity spectrum `f(α) = qα(q) − τ(q)` give the
multifractal spectrum directly (`ghe.singularity_spectrum`). Its width
`Δα = α_max − α_min` is the standard alternative multiscaling-strength
statistic in the physics literature. It is *not* used as the project's
headline statistic, because `α` is a numerical derivative over only 15
q-points, concentrating estimation noise exactly where `H(q)` is already
least reliable (the high-q end) — it is offered as a secondary diagnostic,
not a replacement for `B`.

**Why the tails-vs-dynamics question needed a second estimator.** The
iid-Student-t control (zero memory, by construction, with a tail-fatness
knob via degrees of freedom) shows that fat tails *alone* can produce `B` in
the real markets' range (B≈0.026 at df=5, B≈0.069 at df=3). That is
evidence the structure-function method is *susceptible* to this artefact,
but not proof that the real headline gap *is* this artefact. The cleanest
further test is an estimator that is less susceptible to it by construction.

**MFDFA.** Build the cumulative profile of the mean-centred increments,
`Y(i) = Σ_{k≤i}(d_k − d̄)`; cut it into non-overlapping segments of length
`s` (from both ends of the series, to use all the data); remove a
degree-`order` polynomial trend from each segment; the q-th order
fluctuation function is `F_q(s) = {mean_v[F²(v,s)]^(q/2)}^(1/q)`, and
`F_q(s) ~ s^h(q)` gives `h(q)` by log-log OLS — same downstream algebra
(`τ(q) = qh(q) − 1`, same `multiscaling_strength`) as the structure
function.

Its hyperparameters — polynomial order and the scale range `(s_min,
s_max)` — were chosen the same way the structure function's `τ_max` was:
screened against the five synthetic processes, not guessed. Order 2
(quadratic detrending) made the null bias on anti-persistent fBm *worse* at
every scale range tried (short segments with an extra free parameter
overfit their local trend). `s_min=25` is past the textbook rule of thumb
of `s_min≈10` — empirically, at n=7,000, the null bias keeps shrinking well
past that point. `s_max=150` is where the known-answer MRW positive
controls' recovered `B` is closest to its analytic value (≈108% vs the
structure function's documented 66–72%). The resulting null band (95th
percentile of `B` under a random walk) is 0.0086 — about 1.9× the structure
function's 0.0046, a real and quantified difference in the two estimators'
noise floor at this sample size.

**The cross-check result.** Run on all 22 markets with the same rigor
(bootstrap, shuffle, and an `s_max ∈ {80,100,150,200,300}` grid in place of
`τ_max`): developed mean `B`=0.0213, emerging mean `B`=0.0104, Mann–Whitney
p=0.99 — the opposite direction from the structure function, and stable
across the whole grid. Across the 22 markets, the two methods' `B` values
correlate at `r`=−0.13.

Two tempting shortcuts were checked and ruled out, not assumed:

- *Outlier sensitivity?* Replacing Malaysia's single worst day (a −24% move)
  with its own median |return| drops the structure-function `B` to 65% of
  its original value (as expected for an outlier-driven statistic) but
  *raises* the MFDFA `B` to 306% of its value — the wrong sign for a simple
  "MFDFA just dilutes bad days" story.
- *Power shortfall?* MFDFA's bootstrap intervals are ~46% wider, but it
  distinguishes slightly *more* individual markets from its own null (16/22
  vs 15/22) at nearly identical median |B| (0.0172 vs 0.0174) — it isn't
  simply a noisier version of the same signal.

What does explain it: rerunning the iid-Student-t control through MFDFA
shows it recovers only 30–50% of the spurious `B` the structure function
gets from fat tails alone (df=5: 0.008 vs 0.026; df=3: 0.035 vs 0.069).
MFDFA is measurably, specifically less fooled by exactly the mechanism
already implicated in §14.4.

---

## Level 5 — Full picture

**What this project actually establishes, stated precisely.** Under T. Di
Matteo's own structure-function method, with the estimator pre-validated
against three processes of known `H(q)`, emerging-market equity indices
(MSCI classification, 9 markets) show significantly stronger multiscaling
than developed ones (13 markets) over 1998–2026: `B`=0.0281 vs 0.0085,
Mann–Whitney p=0.0013, surviving a joint regression against excess kurtosis
(t=3.6, p=0.002, R²=0.80) and stable across `τ_max ∈ {10,15,19,25,30}`
(p between 0.0007 and 0.026). This is a genuine replication of the
qualitative pattern reported in Di Matteo (2007) and Di Matteo, Aste &
Dacorogna (2005), on an independent, more recent dataset.

**Where the replication gets complicated, precisely.** Restricting to q≤2 —
the moments a ~7,000-point sample estimates reliably — shrinks the gap by
two thirds and drops it below significance (p=0.055, or p=0.105 excluding
Malaysia). An iid Student-t random walk with literally no temporal
structure reproduces `B` of the same magnitude as the real data. Shuffling
the real returns (destroying time order, preserving the marginal exactly)
*increases* measured `B` on every one of the 22 markets — which initially
reads as evidence against a tail explanation, until the same shuffle is
applied to the synthetic MRW process (whose multiscaling is real by
construction) and produces the identical qualitative effect (0.0119 →
0.0366 at λ²=0.04), meaning the shuffle is a marginal-matched benchmark, not
a one-sided null. All of this — the shuffle-increases-multiscaling paradox,
on real indices, FX rates and bond yields, concluding that observed
multifractality is "mainly a consequence of the characteristic fat-tailed
distribution of the returns" — is the 2012 finding of Barunik, Aste, Di
Matteo and Liu (*Physica A* 391(17), 4234–4251), reached independently here
on different data before that paper was found in the literature search.

**What the MFDFA cross-check adds beyond both prior results.** Having
established (following 2012 Di Matteo et al.) that the structure-function
method is susceptible to a fat-tail artefact, the natural next question —
not addressed in either the 2007/2005 or 2012 papers as far as this
project's search found — is whether the specific developed/emerging claim
survives under an estimator chosen for being less susceptible to that exact
artefact. It does not: MFDFA, independently validated (order=1, s_min=25,
s_max=150, chosen by the same synthetic-process screening methodology as
the structure function, not reused parameters), finds the gap running in
the opposite direction (developed 0.0213 > emerging 0.0104, p=0.99),
robustly across a scale-range grid, with a mechanism that is directly
verifiable rather than asserted (its iid-Student-t recovery is 30–50% of
the structure function's, matching the direction and rough magnitude needed
to explain the reversal) and with two simpler alternative explanations
(outlier sensitivity, insufficient power) explicitly tested and rejected.

**What is, and isn't, concluded from this.** The structure-function headline
is not retracted — it is true, as measured, under its own method, which
remains separately and extensively validated. MFDFA disagreeing does not
make it false: the two estimators are not measuring identical objects (one
samples overlapping lag-τ differences directly; the other averages detrended
variance within non-overlapping windows at a different, wider range of
scales), so disagreement is informative rather than contradictory. What
changes is the epistemic status of the finding: a result that replicates a
published pattern, survives several within-method robustness checks, and
then fails to survive a structurally different, independently validated
estimator chosen specifically to probe its weakest documented point, is
better described as "this dataset's manifestation of the tail-driven effect
Di Matteo's own 2012 paper already characterised" than as "emerging markets
have different dynamics." That is a stronger, more specific, and more
falsifiable claim than the one the project started with — and arriving at
it is the actual content of the work.

**For verification.** Every number above is read from
`results/summary.csv`, `results/controls.csv`, `results/controls_iid_t.csv`
(structure function) and `results/summary_mfdfa.csv`,
`results/controls_mfdfa.csv`, `results/controls_iid_t_mfdfa.csv`,
`results/sf_vs_mfdfa.csv` (MFDFA), generated respectively by
`scripts/0{2,3,4}_*.py` and `scripts/0{2b,3b}_*_mfdfa.py`. The estimators
themselves are `src/ghe.py` and `src/mfdfa.py` (~150 lines each); the three
validation processes are shared, in `src/synthetic.py`. `docs/THEORY.md`
§§1–13 derive the structure-function side in full; §14 derives the MFDFA
cross-check and the connection to Barunik, Aste, Di Matteo & Liu (2012).
