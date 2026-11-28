"""Fayın tamamını boşluksuz, veri yoğunluğuna göre değişen boyda pencerelerle döşe.

Ulaş: "boşluk olmadan yap", "daha çok veri olanları daha kısa aralıklarla böl",
"az veri olanları daha uzun aralıklarla böl".

Kural: her kol batıdan doğuya yürünür. Pencere en az L_MIN km'dir; eğitim dönemindeki (1998-2017)
deprem sayısı N_HEDEF'e ulaşana dek 1 km'lik adımlarla uzatılır, en çok L_MAX km. Sona kalan L_MIN'den
kısa parça son pencereye eklenir. Böylece her pencerede yaklaşık eşit eğitim depremi olur; sessiz
yerler L_MAX'ta durur ve "zayıf" işaretlenir.

Deprem sayımı pencerelerin gerçekte kullandığı kayıtlarla aynıdır: faya 10 km içinde, GFZ M>=2
(2006-2017) + EMSC/USGS M>=2,5 (1998-2005); GFZ kapsamı dışında (30,65°D doğusu) yalnız EMSC/USGS.
Sınav sayısı: GFZ 2018-2023 + EMSC M>=2,5 2024-2026.

Hücre boyu: hücre sayısı (L/h)·(14/h) ~4900'ü aşmasın diye h = max(0,2; sqrt(14·L/4900)) km.
Çıktı: reports/tables/pencere_dose.csv, src/deprem/mekanik/pencereler_doseme.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from deprem.mekanik import ozumseme as oz  # noqa: E402

L_MIN, L_MAX = 7.0, 40.0
N_HEDEF = int(sys.argv[1]) if len(sys.argv) > 1 else 100
N_HUCRE = 4900
KX = 111.195 * np.cos(np.radians(40.7))


def sayimlar():
    ag = oz.Ag(ROOT, "mta_kaba")
    y = ag.cells[ag.cells.z_ust == 0]
    g = oz.gfz_katalog(ROOT)
    g = g[g.mw >= 2.0]
    raw = pd.read_parquet(ROOT / "data" / "mekanik" / "katalog_marmara.parquet")
    raw = raw[(raw.time >= "1998-01-01") & (raw.mw >= 2.5) & raw.source.isin(["emsc", "usgs"])]
    raw = raw.assign(t=raw.time.map(oz._ondalik), lon=raw.longitude.astype(float), lat=raw.latitude.astype(float))
    raw = raw.sort_values("t")
    # iki kaynakta aynı deprem: 20 s ve ~20 km içinde olanlardan biri
    tt, la, lo = raw.t.to_numpy(), raw.lat.to_numpy(), raw.lon.to_numpy()
    tut = np.ones(len(raw), bool)
    for i in range(1, len(raw)):
        j = i - 1
        while j >= 0 and (tt[i] - tt[j]) * 365.25 * 86400 < 20:
            if tut[j] and abs(la[i] - la[j]) < 0.2 and abs(lo[i] - lo[j]) < 0.25:
                tut[i] = False
                break
            j -= 1
    raw = raw[tut]
    out = {}
    for ad, d in [("gfz", g), ("emsc", raw)]:
        x, dd, k = oz.izdusum(ag, d.lon.to_numpy(), d.lat.to_numpy())
        out[ad] = pd.DataFrame({"x": x, "d": dd, "kol": k, "t": d.t.to_numpy(), "lon": d.lon.to_numpy()})
    return y, out


def main():
    y, S = sayimlar()
    g, e = S["gfz"], S["emsc"]
    g, e = g[g.d <= 10], e[e.d <= 10]
    rows = []
    for kol in sorted(y.kol.unique()):
        yy = y[y.kol == kol].sort_values("x_km")
        s, x0, x1 = yy.x_km.to_numpy(), yy.x_km.min(), yy.x_km.max()

        def kaba_mi(a, b):
            return bool(kol == 0 and np.interp(0.5 * (a + b), s, yy.lon.to_numpy()) > 30.65)

        def say(a, b, sinav=False):
            gg, ee = g[(g.kol == kol) & (g.x >= a) & (g.x < b)], e[(e.kol == kol) & (e.x >= a) & (e.x < b)]
            if sinav:
                if kaba_mi(a, b):
                    return int(((ee.t >= 2018) & (ee.t < 2026.9)).sum())
                return int(((gg.t >= 2018) & (gg.t < 2024)).sum() + ((ee.t >= 2024) & (ee.t < 2026.9)).sum())
            if kaba_mi(a, b):
                return int((ee.t < 2018).sum())
            return int(((gg.t < 2018)).sum() + (ee.t < 2006).sum())

        a, parca = x0, []
        while a < x1 - 1e-6:
            b = min(a + L_MIN, x1)
            while b < x1 and say(a, b) < N_HEDEF and b - a < L_MAX:
                b = min(b + 1.0, x1)
            parca.append([a, b])
            a = b
        if len(parca) > 1 and parca[-1][1] - parca[-1][0] < L_MIN:
            son = parca.pop()
            parca[-1][1] = son[1]
            a, b = parca[-1]
            if b - a >= 2 * L_MIN and say(a, b) >= 2 * N_HEDEF:  # yoğun yerde birleşen kuyruk ikiye bölünür
                n2, c = say(a, b) / 2, a + L_MIN
                while c < b - L_MIN and say(a, c) < n2:
                    c += 1.0
                parca[-1] = [a, c]
                parca.append([c, b])
        for a, b in parca:
            L = b - a
            xc = 0.5 * (a + b)
            lo = np.interp([a, xc, b], s, yy.lon.to_numpy())
            la = np.interp([a, xc, b], s, yy.lat.to_numpy())
            st = (np.degrees(np.arctan2((lo[2] - lo[0]) * KX, (la[2] - la[0]) * 111.195)) + 360) % 360
            # izin pencere doğrusundan en büyük sapması (km)
            m = (s >= a) & (s <= b)
            px, py = (yy.lon.to_numpy()[m] - lo[0]) * KX, (yy.lat.to_numpy()[m] - la[0]) * 111.195
            ux, uy = (lo[2] - lo[0]) * KX, (la[2] - la[0]) * 111.195
            sap = float(np.max(np.abs(px * uy - py * ux)) / max(np.hypot(ux, uy), 1e-9)) if m.any() else 0.0
            kes = yy.kesim.to_numpy()[np.argmin(np.abs(s - xc))] or {2: "İmralı", 3: "Armutlu"}.get(kol, "")
            n_e, n_s = say(a, b), say(a, b, True)
            h = round(max(0.2, float(np.sqrt(14.0 * L / N_HUCRE))), 3)
            rows.append(dict(ad=f"d{kol}_{int(round(a)):03d}", kol=int(kol), kesim=kes, x0=round(a, 1), x1=round(b, 1),
                             L=round(L, 1), lon=round(lo[1], 3), lat=round(la[1], 3), lon0=round(lo[0], 3), lon1=round(lo[2], 3),
                             strike=round(st, 1), h=h, iz_sapma_km=round(sap, 2), kaba=kaba_mi(a, b),
                             n_egitim=n_e, n_sinav=n_s, zayif=bool(n_e < N_HEDEF)))
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "reports" / "tables" / "pencere_dose.csv", index=False)
    js = {r.ad: {"v": [r.lon, r.lat, r.strike, r.L, 14.0, r.h, r.kol], "kaba": r.kaba, "kesim": r.kesim,
                 "zayif": r.zayif, "n_egitim": r.n_egitim, "n_sinav": r.n_sinav} for r in d.itertuples()}
    (ROOT / "src" / "deprem" / "mekanik" / "pencereler_doseme.json").write_text(json.dumps(js, indent=1, ensure_ascii=False))
    pd.set_option("display.width", 250)
    print(d.drop(columns=["lon0", "lon1"]).to_string())
    print("pencere", len(d), "zayif", int(d.zayif.sum()), "toplam km", round(d.L.sum(), 1),
          "egitim", int(d.n_egitim.sum()), "sinav", int(d.n_sinav.sum()))


if __name__ == "__main__":
    main()
