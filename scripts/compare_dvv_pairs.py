"""dv/v'li model (turkey_v5_dvvpairs) ile dv/v'siz taban (turkey_v5_base) arasında olay başına fark.

Aynı katalog örnekleri (r00 en iyi tahmin + r01..r10) eşleştirilir; fark = (zaman+konum log-olabilirlik) dv/v − taban.
Güven aralığı: ay blokları üzerinden bootstrap, katalog örnekleri birlikte.

    python scripts/compare_dvv_pairs.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE, NEW = ROOT / "runs" / "turkey_v5_base", ROOT / "runs" / "turkey_v5_dvvpairs"
# pilot hücreleri: istasyon kutusu (eval_dvv_pairs ile aynı)
pairs = pd.read_csv(ROOT / "data" / "obs" / "dvv_pairs_series.csv.gz", usecols=["lat", "lon"])
BOX = (pairs["lat"].min() - 0.3, pairs["lat"].max() + 0.3, pairs["lon"].min() - 0.4, pairs["lon"].max() + 0.4)
START = 2012  # dv/v özelliği ancak ~2011 sonundan itibaren dolu (2 yıllık referans + 6 ay)


def load(run, r):
    f = run / "longwalk" / f"m4_r{r:02d}_events.csv.gz"
    d = pd.read_csv(f)
    d["ll"] = d["time_ll_model"] + d["space_ll_model"]
    return d


def boot(diff, month, rng, n=2000):
    months = np.unique(month)
    idx = {m: np.flatnonzero(month == m) for m in months}
    out = np.empty(n)
    for i in range(n):
        pick = rng.choice(months, len(months))
        sel = np.concatenate([idx[m] for m in pick])
        out[i] = diff[sel].mean()
    return np.percentile(out, [2.5, 97.5])


def main():
    rng = np.random.default_rng(0)
    rows = []
    for r in range(11):
        b, n = load(BASE, r), load(NEW, r)
        m = b.merge(n[["time", "ll", "time_ll_model", "space_ll_model"]], on="time", suffixes=("_b", "_n"))
        m["d"] = m["ll_n"] - m["ll_b"]
        m["d_time"] = m["time_ll_model_n"] - m["time_ll_model_b"]
        m["d_space"] = m["space_ll_model_n"] - m["space_ll_model_b"]
        m["r"] = r
        rows.append(m)
    a = pd.concat(rows, ignore_index=True)
    a["month"] = a["time"].str[:7]
    a["pilot"] = a["lat"].between(BOX[0], BOX[1]) & a["lon"].between(BOX[2], BOX[3])
    res = []
    for name, sel in [("Türkiye geneli 1995-2025", a["year"] >= 0),
                      (f"Türkiye geneli {START}-2025", a["year"] >= START),
                      (f"Doğu Marmara pilot kutusu {START}-2025", (a["year"] >= START) & a["pilot"])]:
        s = a[sel]
        lo, hi = boot(s["d"].to_numpy(), s["month"].to_numpy(), rng)
        res.append({"kapsam": name, "olay (örnek başına)": int(len(s) / 11),
                    "olay başına fark (nat)": round(s["d"].mean(), 4), "%95 GA": f"[{lo:.4f}, {hi:.4f}]",
                    "zaman payı": round(s["d_time"].mean(), 4), "konum payı": round(s["d_space"].mean(), 4)})
    res = pd.DataFrame(res)
    c = pd.read_csv(NEW / "longwalk" / "m4_r00_coefs.csv")  # katsayılar yalnız en iyi tahmin kataloğu için yazılıyor
    cy = c[c["özellik"] == "dvv_pair_change_pct"][["yıl", "katsayı", "sd"]]
    (ROOT / "data" / "obs" / "dvv_pairs_compare.json").write_text(json.dumps(
        {"fark": res.to_dict("records"), "katsayı": cy.to_dict("records")}, indent=2, ensure_ascii=False))
    print(res.to_markdown(index=False))
    print()
    print(cy[cy["yıl"] >= START].round(3).to_markdown(index=False))


if __name__ == "__main__":
    main()
