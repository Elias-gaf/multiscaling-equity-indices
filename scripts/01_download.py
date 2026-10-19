"""Download daily index closes, clean them, and write a data-quality report.

Design decisions, all deliberate:

* Each index keeps its OWN trading calendar.  Aligning 22 exchanges onto a
  common calendar and forward-filling would insert runs of exactly-zero
  returns on other markets' holidays, which suppresses small-tau structure
  functions and biases H(q).  Series are trimmed to a common DATE RANGE
  instead, not a common date index.
* Closes are raw index levels (auto_adjust=False).  Price indices carry no
  dividend adjustment, so there is nothing to adjust.
* "Stale" days (close identical to the previous close) are counted and any
  series above 1% is dropped -- those are exchange feeds that repeat the last
  print through holidays.
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import yfinance as yf

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(HERE, "data", "raw")
CLEAN = os.path.join(HERE, "data", "clean")
os.makedirs(RAW, exist_ok=True); os.makedirs(CLEAN, exist_ok=True)

# MSCI market classification as of 2026.
UNIVERSE = [
    ("^GSPC",    "S&P 500",        "United States", "Developed"),
    ("^FTSE",    "FTSE 100",       "United Kingdom","Developed"),
    ("^GDAXI",   "DAX",            "Germany",       "Developed"),
    ("^FCHI",    "CAC 40",         "France",        "Developed"),
    ("^N225",    "Nikkei 225",     "Japan",         "Developed"),
    ("^AXJO",    "ASX 200",        "Australia",     "Developed"),
    ("^GSPTSE",  "S&P/TSX",        "Canada",        "Developed"),
    ("^SSMI",    "SMI",            "Switzerland",   "Developed"),
    ("^AEX",     "AEX",            "Netherlands",   "Developed"),
    ("^IBEX",    "IBEX 35",        "Spain",         "Developed"),
    ("^BFX",     "BEL 20",         "Belgium",       "Developed"),
    ("^ATX",     "ATX",            "Austria",       "Developed"),
    ("^HSI",     "Hang Seng",      "Hong Kong",     "Developed"),
    ("^BVSP",    "Bovespa",        "Brazil",        "Emerging"),
    ("XU100.IS", "BIST 100",       "Turkey",        "Emerging"),
    ("^BSESN",   "BSE Sensex",     "India",         "Emerging"),
    ("^JKSE",    "Jakarta Composite","Indonesia",   "Emerging"),
    ("^KS11",    "KOSPI",          "South Korea",   "Emerging"),
    ("^TWII",    "TAIEX",          "Taiwan",        "Emerging"),
    ("^MXX",     "IPC",            "Mexico",        "Emerging"),
    ("^KLSE",    "FTSE Bursa KLCI","Malaysia",      "Emerging"),
    ("^SET.BK",  "SET",            "Thailand",      "Emerging"),
]

START, END = "1998-01-01", "2026-07-17"   # common range; SET's feed ends 2026-07-17
MAX_STALE = 0.01


def fetch(sym: str) -> pd.Series:
    cache = os.path.join(RAW, sym.replace("^", "").replace(".", "_") + ".csv")
    if os.path.exists(cache):
        s = pd.read_csv(cache, index_col=0, parse_dates=True).iloc[:, 0]
    else:
        df = yf.download(sym, start="1990-01-01", end="2026-09-09",
                         progress=False, auto_adjust=False, threads=False)
        s = df["Close"].squeeze()
        s.to_csv(cache)
    s.index = pd.to_datetime(s.index)
    return s.sort_index()


def main():
    report, kept = [], {}
    for sym, name, country, tier in UNIVERSE:
        s = fetch(sym)
        n_raw = len(s)
        s = s.dropna()
        s = s[s > 0]
        s = s[~s.index.duplicated(keep="first")]
        s = s.loc[START:END]
        stale = float((s.diff() == 0).mean())
        gap = s.index.to_series().diff().dt.days
        maxgap = int(gap.max()) if len(gap.dropna()) else 0
        drop = stale > MAX_STALE
        report.append(dict(symbol=sym, index=name, country=country, tier=tier,
                           n_raw=n_raw, n_used=len(s),
                           start=str(s.index[0].date()), end=str(s.index[-1].date()),
                           stale_frac=stale, max_gap_days=maxgap,
                           status="dropped: stale>1%" if drop else "kept"))
        if not drop:
            kept[sym] = np.log(s.astype(float))

    rep = pd.DataFrame(report)
    rep.to_csv(os.path.join(CLEAN, "data_quality.csv"), index=False)

    # one long-format file: each series on its own calendar, no fills
    frames = [pd.DataFrame({"symbol": k, "date": v.index, "logprice": v.values})
              for k, v in kept.items()]
    pd.concat(frames).to_csv(os.path.join(CLEAN, "logprices.csv"), index=False)

    pd.set_option("display.width", 200)
    print(rep[["symbol", "index", "tier", "n_used", "start", "end",
               "stale_frac", "max_gap_days", "status"]].to_string(index=False))
    print(f"\nkept {len(kept)} of {len(UNIVERSE)} series -> data/clean/logprices.csv")


if __name__ == "__main__":
    main()
