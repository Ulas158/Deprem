"""ETAS varyantlarını ileriye kayan şekilde eski ETAS'a karşı karşılaştır.

Soru: ETAS büyük depremli yıllardan sonra (ve sakin yıllarda) deprem sayısını fazla
tahmin ediyor; bunu düzelten bir ayar, tahmin gücünü bozmadan mümkün mü?

Her yıl Y (2015-2025) için ETAS yalnız Y öncesiyle uydurulur ve Y tahmin edilir:
* yıllık sayı: 200 benzetimden negatif binom öngörü dağılımıyla log-puan, %95 aralık kapsaması
* deprem deprem: koşullu yoğunluğun zaman + konum log-puanı (eski ETAS'a göre fark,
  ay blokları bootstrap'ıyla %95 GA)

Varyantlar:
  A  eski ETAS (benzetimde büyüklükler GR, b ile M8'e kadar; α serbest)
  B  aynı parametreler, benzetimde büyüklükler eğitim kataloğundan yeniden örneklenir
  C  B + α = b·ln10 sabit (küçük depremler daha az, büyükler daha çok üretken)
  D  B + α = 2,0 sabit

Kullanım:
  python scripts/compare_etas.py            # hesapla (~30 dk, 4 çekirdek) ve raporla
  python scripts/compare_etas.py --report   # yalnız kayıtlı sonuçlardan rapor
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from deprem.forecast.etas import ETAS  # noqa: E402
from deprem.forecast.eventwalk import block_ci, score_year  # noqa: E402
from deprem.pipeline import Run, ts  # noqa: E402

CONFIG = ROOT / "configs/train_obs.toml"
OUT = ROOT / "runs/etas_variants/results.pkl"
REPORT = ROOT / "reports/etas_karsilastirma.md"
YEARS = range(2015, 2026)
MAGS = (3.5, 4.0, 5.0)
VARIANTS = {
    "A": dict(label="eski ETAS", samp="gr", fixed=None),
    "B": dict(label="gerçek büyüklük dağılımıyla benzetim", samp="empirical", fixed=None),
    "C": dict(label="B + α = b·ln10", samp="empirical", fixed={"alpha": "b"}),
    "D": dict(label="B + α = 2,0", samp="empirical", fixed={"alpha": 2.0}),
}
_R = {}


def _init():
    r = Run(CONFIG, ROOT)
    _R.update(run=r, cat=r.catalog, c=r.cfg["etas"],
              pj=json.loads((r.dir / "walk_etas_params.json").read_text()))


def job(args):
    v, Y = args
    r, cat, c, pj = _R["run"], _R["cat"], _R["c"], _R["pj"]
    spec = VARIANTS[v]
    t0, t1 = r.window(Y)
    e = ETAS(m0=c["m0"], bandwidth_km=c["bandwidth_km"], r_max_km=c["r_max_km"],
             dt_max_days=c["dt_max_days"], mag_sampling=spec["samp"])
    if spec["fixed"] is None:  # ileriye kayan eğitimdeki parametreler
        tr = cat[(cat.mag >= e.m0) & (cat.time >= ts(c["fit_start"])) & (cat.time < t0)]
        e.b = pj[str(Y)]["b"]
        e._background(tr)
        e.set_mag_pool(tr.mag)
        e.params = pj[str(Y)]["params"]
        ll = None
    else:
        e.fit(cat, ts(c["fit_start"]), t0, init=pj[str(Y)]["params"], fixed=spec["fixed"])
        ll = e.fit_info["loglik"]
    sims = e.simulate(t0, t1, cat, n_sim=200, seed=Y)
    n = np.array([s[0].size for s in sims])
    ev = {mg: score_year(e, cat, t0, t1, r.grid, None, mg) for mg in MAGS}
    return v, Y, n, ev, ll, dict(e.params), e.b


def nb_logp(sims, k):
    m, s2 = sims.mean(), sims.var()
    if s2 <= m:
        return stats.poisson.logpmf(k, m)
    p = m / s2
    return stats.nbinom.logpmf(k, m * p / (1 - p), p)


def report(res) -> str:
    R = {(v, Y): (n, ev) for v, Y, n, ev, *_ in res}
    V = list(VARIANTS)
    Ys = sorted({Y for _, Y in R})
    L = ["# ETAS'ın büyük deprem sonrası fazla tahmini: düzeltme denemeleri", "",
         "> Deneysel araştırma çıktısıdır, deprem tahmini değildir.", "",
         "Her yıl ETAS yalnız o yıldan önceki veriyle uydurulup o yıl tahmin edildi (2015-2025, M≥3,5). "
         "Varyantlar eski ETAS'a (A) karşı karşılaştırıldı.", ""]
    L += ["| varyant | açıklama |", "|---|---|"] + [f"| {k} | {s['label']} |" for k, s in VARIANTS.items()] + [""]
    rows = []
    for Y in Ys:
        obs = len(R[("A", Y)][1][3.5])
        row = {"yıl": Y, "gözlenen": obs}
        for v in V:
            n = R[(v, Y)][0]
            row[f"{v} medyan (%95)"] = f"{int(np.median(n))} ({int(np.quantile(n, .025))}-{int(np.quantile(n, .975))})"
            row[f"{v} log-puan"] = round(float(nb_logp(n, obs)), 2)
        rows.append(row)
    df = pd.DataFrame(rows)
    L += ["## Yıllık deprem sayısı (M≥3,5)", "",
          "Medyan ve %95 aralık 200 benzetimden. Log-puan: gözlenen sayının öngörü dağılımındaki "
          "olasılığının logaritması (büyük = iyi).", "", df.to_markdown(index=False), ""]
    summ = []
    rng = np.random.default_rng(0)
    for v in V:
        sc = df[f"{v} log-puan"].to_numpy()
        d = sc - df["A log-puan"].to_numpy()
        bs = d[rng.integers(0, len(d), (4000, len(d)))].sum(1)
        med = np.array([int(s.split()[0]) for s in df[f"{v} medyan (%95)"]])
        cov = sum(int(s.split("(")[1].split("-")[0]) <= o <= int(s.split("-")[1].rstrip(")"))
                  for s, o in zip(df[f"{v} medyan (%95)"], df["gözlenen"]))
        summ.append({"varyant": v, "toplam log-puan": round(float(sc.sum()), 1),
                     "A'ya göre fark (%95 GA, yıl bootstrap)": "—" if v == "A" else
                     f"{d.sum():+.1f} [{np.quantile(bs, .025):+.1f}, {np.quantile(bs, .975):+.1f}]",
                     "A'dan iyi yıl": "—" if v == "A" else f"{int((d > 0).sum())}/{len(d)}",
                     "medyan / gözlenen (geometrik ort.)": round(float(np.exp(np.mean(np.log(med / df["gözlenen"])))), 2),
                     "%95 aralık gözleneni kapsadı": f"{cov}/{len(df)}"})
    L += ["### Özet", "", pd.DataFrame(summ).to_markdown(index=False), ""]
    ev_rows = []
    for mg in MAGS:
        for v in V[1:]:
            d = []
            for Y in Ys:
                a, b = R[("A", Y)][1][mg], R[(v, Y)][1][mg]
                d.append(pd.DataFrame({"time": a["time"],
                                       "t": b.time_ll_etas - a.time_ll_etas,
                                       "s": b.space_ll_etas - a.space_ll_etas}))
            d = pd.concat(d)
            tot = (d.t + d.s).to_numpy()
            lo, hi = block_ci(tot, d.time.dt.strftime("%Y-%m").to_numpy())
            ev_rows.append({"hedef": f"M≥{mg:g}", "varyant": v, "olay": len(d),
                            "zaman": round(float(d.t.mean()), 3), "konum": round(float(d.s.mean()), 3),
                            "toplam (olay başına)": round(float(tot.mean()), 3), "%95 GA": f"[{lo:.3f}, {hi:.3f}]"})
    L += ["## Deprem deprem (eski ETAS'a göre olay başına log-puan farkı)", "",
          "B, A ile aynı parametreleri kullanır; yalnız benzetimi değiştirir, bu yüzden deprem deprem "
          "puanı aynıdır.", "", pd.DataFrame(ev_rows).to_markdown(index=False), ""]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="yalnız kayıtlı sonuçlardan rapor")
    args = ap.parse_args()
    if not args.report:
        jobs = [(v, Y) for v in VARIANTS for Y in YEARS]
        with Pool(4, initializer=_init) as pool:
            res = pool.map(job, jobs)
        OUT.parent.mkdir(parents=True, exist_ok=True)
        pickle.dump(res, open(OUT, "wb"))
    res = pickle.load(open(OUT, "rb"))
    REPORT.write_text(report(res))
    print(REPORT)


if __name__ == "__main__":
    main()
