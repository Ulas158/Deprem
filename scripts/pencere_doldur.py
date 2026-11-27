"""Fayın tamamını 14 km'lik pencerelerle döşe (Ulaş: "boşlukları pencerelerle doldur").

Her kol (kuzey, orta, İmralı, Armutlu) iz boyunca ardışık 14 km'lik parçalara bölünür. Var olan
pencereyle yarıdan fazla örtüşen parça atlanır. Kalan her parçada 1998-2017 eğitim döneminde faya 5 km
içindeki M>=2 (GFZ) ve M>=2,5 (EMSC/USGS) deprem sayılır; en az N_MIN olan parça pencere olur.
GFZ kapsamı dışında (30,65°D doğusu) pencere "kaba" (yalnız EMSC) işaretlenir.
Çıktı: reports/tables/pencere_doldur.csv ve src/deprem/mekanik/pencereler_ek.json
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

L, N_MIN = 14.0, 30
KX = 111.195 * np.cos(np.radians(40.7))


def main():
    ag = oz.Ag(ROOT, oz.BUYUK_AG)
    y = ag.cells[ag.cells.z_ust == 0]
    g = oz.gfz_katalog(ROOT)
    g = g[(g.mw >= 2.0) & (g.t < 2018.0)]
    raw = pd.read_parquet(ROOT / "data" / "mekanik" / "katalog_marmara.parquet")
    raw = raw[(raw.time >= "1998-01-01") & (raw.time < "2018-01-01") & (raw.mw >= 2.5) & raw.source.isin(["emsc", "usgs"])]
    pr = {}
    for ad, lon, lat, t in [("gfz", g.lon.values, g.lat.values, g.t.values),
                            ("emsc", raw.longitude.values.astype(float), raw.latitude.values.astype(float), None)]:
        x, d, k = oz.izdusum(ag, lon, lat)
        pr[ad] = (x, d, k, lon)
    var = []
    for ad, v in oz.YEREL_AGLAR.items():
        lo, la, st, LL, z, h, *kk = v
        kol = kk[0] if kk else 0
        yy = y[y.kol == kol].sort_values("x_km")
        i = np.argmin((yy.lon - lo) ** 2 * KX ** 2 + (yy.lat - la) ** 2 * 111.195 ** 2)
        var.append((kol, float(yy.x_km.iloc[i])))
    rows = []
    for kol in sorted(y.kol.unique()):
        yy = y[y.kol == kol].sort_values("x_km")
        x0, x1 = yy.x_km.min(), yy.x_km.max()
        for a in np.arange(x0, x1 - L / 2, L):
            b = min(a + L, x1)
            xc = 0.5 * (a + b)
            if any(k == kol and abs(xc - xv) < L / 2 for k, xv in var):
                continue
            n_g = int(((pr["gfz"][2] == kol) & (pr["gfz"][1] <= 5) & (pr["gfz"][0] >= a) & (pr["gfz"][0] < b)).sum())
            n_e = int(((pr["emsc"][2] == kol) & (pr["emsc"][1] <= 5) & (pr["emsc"][0] >= a) & (pr["emsc"][0] < b)).sum())
            s = yy.x_km.values
            lo = np.interp([a, xc, b], s, yy.lon.values)
            la = np.interp([a, xc, b], s, yy.lat.values)
            st = (np.degrees(np.arctan2((lo[2] - lo[0]) * KX, (la[2] - la[0]) * 111.195)) + 360) % 360
            kes = yy.kesim.values[np.argmin(np.abs(s - xc))] or {2: "İmralı", 3: "Armutlu"}.get(kol, "")
            kaba = bool(kol == 0 and lo[1] > 30.65)
            n = n_e if kaba else n_g + n_e
            rows.append(dict(kol=int(kol), kesim=kes, x0=round(a, 1), x1=round(b, 1), lon=round((lo[0] + lo[2]) / 2, 3),
                             lat=round((la[0] + la[2]) / 2, 3), strike=round(st, 1), n_gfz=n_g, n_emsc=n_e,
                             kaba=kaba, n=n, secildi=n >= N_MIN and (b - a) >= 10))
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "reports" / "tables" / "pencere_doldur.csv", index=False)
    ek = {}
    for r in d[d.secildi].itertuples():
        ad = f"p{r.kol}_{int(r.x0):04d}"
        ek[ad] = {"v": [r.lon, r.lat, r.strike, L, 14.0, 0.2, r.kol], "kaba": r.kaba, "kesim": r.kesim}
    (ROOT / "src" / "deprem" / "mekanik" / "pencereler_ek.json").write_text(json.dumps(ek, indent=1, ensure_ascii=False))
    pd.set_option("display.width", 200)
    print(d.to_string())
    print("seçilen", int(d.secildi.sum()), "/", len(d))


if __name__ == "__main__":
    main()
