# Lines fit for use elsewhere

Every number below is checked against `results/summary.csv` and, for the
MFDFA cross-check, `results/summary_mfdfa.csv` / `results/sf_vs_mfdfa.csv`.
Do not round them differently in different documents.

## CV — one line

> **Multiscaling in international equity indices** (2026) — estimated the
> generalised Hurst exponent H(q) for 22 national equity indices over
> 1998–2026 by two structurally independent methods, validated both against
> known synthetic processes, and traced a sign-reversing disagreement between
> them to the same fat-tail mechanism Barunik, Aste, Di Matteo and Liu (2012)
> identified on different data. Python, ~1,600 lines. [github.com/…]

## Supporting statement — two sentences

> Working through Di Matteo's structure-function method on 22 equity indices, I
> replicated the finding that emerging markets multiscale more strongly than
> developed ones (Mann–Whitney p = 0.0013), then found the separation
> disappears when restricted to moments q ≤ 2, and that an iid Student-t
> random walk with no temporal structure produces multiscaling of the same
> magnitude as the real markets. Cross-checking with a second estimator
> (MFDFA) chosen for being less susceptible to exactly that artefact, the gap
> reverses entirely (p = 0.99) — which, I later found, is the same conclusion
> Di Matteo herself reached by a different route in a 2012 paper I hadn't read
> when I first saw it in my own numbers.

## The question for Prof. Di Matteo

The brief says: write to her only once there is a result, and ask something
specific. The result changed after finding her own 2012 paper with Barunik,
Aste and Liu, so the question is now built around that, not the original q≤2
fragility alone.

> Estimating H(q) by the structure-function method on 22 indices over
> 1998–2026, I found the expected developed/emerging separation in the
> multiscaling slope (0.0085 vs 0.0281, p = 0.0013), surviving a kurtosis
> control. I also found every market's shuffled surrogate returns a *larger*
> slope than the market itself — and only afterwards found your 2012 paper
> with Barunik, Aste and Liu reporting the identical pattern on indices, FX
> rates and bond yields, and concluding it's mostly the fat-tailed marginal,
> not dynamics. I pushed on that by running MFDFA as a second estimator,
> validated the same way as the structure function against known synthetic
> processes, and confirmed it is measurably less susceptible to the same
> fat-tail artefact (it recovers only 30–50% of the spurious multiscaling an
> iid Student-t walk gives the structure function). Run on the real 22
> markets, MFDFA doesn't just weaken the developed/emerging gap — it reverses
> it (p = 0.99), consistently across a scale-range robustness grid. Does that
> match what you'd expect from combining your 2007 and 2012 results, or is
> there a reason the two estimators should disagree this sharply that I'm
> missing?

## What is not claimed

- Not "I discovered that multifractality is an artefact." The literature
  already knows this — Barunik and Kristoufek (2010) is one reference, and
  Barunik, Aste, Di Matteo and Liu (2012) reached the same conclusion on
  different data years before this project. The contribution here is
  independently rediscovering it on a new dataset (22 indices, 1998–2026)
  with a calibrated estimator, then extending it with a second, structurally
  different method and naming the mechanism — not asserting it first.
- Not a finished paper. It is two validated estimators, a replication, and
  one honest negative result that happens to match the literature.
