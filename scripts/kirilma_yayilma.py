"""Zincirleme kırılma sınaması (Ulaş: "küçük fay kırılır, yanındakiler sırayla kırılır, uca doğru gider").

Büyük model geçmiş depremlerle 1999 İzmit'in hemen öncesine kadar koşulur (1999 ve sonrası hiçbir deprem
modele verilmez). Sonra yalnız İzmit odağındaki hücre kırılma eşiğine getirilir ve kırılmanın nereye kadar
kendiliğinden yayıldığına bakılır (alan sınırı yok). Ardından model 1 yıl daha kendi başına ilerletilir:
Düzce kesiminde (31,0-31,4°D) kendiliğinden deprem çıkıyor mu, ne zaman?
Gözlem: İzmit 17.08.1999 Mw 7,6, kırık ~29,3°D (Hersek/Çınarcık) - 30,9°D (Gölyaka), ~145 km;
Düzce 12.11.1999 (87 gün sonra) Mw 7,2, ~31,0-31,4°D.
Kullanım: python scripts/kirilma_yayilma.py [--ag mta_kaba] [--takim literatur]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
import mekanik_ikinci_tur as M  # noqa: E402

T_IZMIT = 1999 + (229 - 1 + 0.04) / 365.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ag", default="mta_kaba")
    ap.add_argument("--takim", default="literatur", help="literatur ya da bir theta json dosyası")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    th = dict(oz.baslangic_theta(), zayif=0.0)
    if a.takim != "literatur":
        th.update(json.loads(Path(a.takim).read_text()))
    ag, ol, _ = M._hazirla("buyuk", a.ag, True)
    ol = [o for o in ol if o.t < T_IZMIT - 1e-4]
    A = M.AYAR["buyuk"]
    m, _ = oz.Kosu(ag, th, ol, a.seed, A["isinma"], A["t_bas"], T_IZMIT - 1e-4, T_IZMIT).calistir()
    m.ilerle(T_IZMIT)
    lon, lat = m.c.lon.to_numpy(), m.c.lat.to_numpy()
    z = m.c.z_km.to_numpy()
    ana = (m.c.kol.to_numpy() == 0) & m.vw
    h0 = oz.HIPOSANTR["1999-08-17"]
    k = np.flatnonzero(ana)[np.argmin(((lon - h0[0]) * 84) ** 2 + ((lat - h0[1]) * 111) ** 2 + (z - h0[2]) ** 2)
                            if False else np.argmin((((lon - h0[0]) * 84) ** 2 + ((lat - h0[1]) * 111) ** 2
                                                     + (z - h0[2]) ** 2)[ana])]
    dol = oz.doluluk(m)
    sec = ana & (m.c.lon.to_numpy() > 29.0) & (m.c.lon.to_numpy() < 31.5)
    prof = {f"{l:.1f}": round(float(np.nanmean(dol[sec & (np.abs(lon - l) < 0.05)])), 3)
            for l in np.arange(29.0, 31.5, 0.1)}
    m.S[k] = max(m.S[k], m.S_fail[k]) + 1e-3
    n0 = len(m.olaylar)
    ev = m.kirilma(int(k))
    out = {"ag": a.ag, "takim": a.takim, "doluluk_oncesi": prof}
    if ev is not None:
        kay = np.flatnonzero(m.slip > 0)
        out["izmit_model"] = {"mw": round(ev["mw"], 2), "lon0": None, "lon1": None}
        e = m.olaylar[-1]
        ix = ana & (m.x >= e["x0"]) & (m.x <= e["x1"])
        out["izmit_model"].update(lon0=round(float(lon[ix].min()), 2), lon1=round(float(lon[ix].max()), 2),
                                  boy_km=round(float(e["x1"] - e["x0"]), 1))
    else:
        out["izmit_model"] = None
    m.ilerle(T_IZMIT + 1.0)
    son = [e for e in m.olaylar[n0 + 1:] if e["mw"] >= 6.0]
    dz = []
    for e in son:
        ix = ana & (m.x >= e["x0"]) & (m.x <= e["x1"])
        if not ix.any():
            continue
        dz.append({"gun": round((e["t"] - T_IZMIT) * 365.25, 1), "mw": round(e["mw"], 2),
                   "lon0": round(float(lon[ix].min()), 2), "lon1": round(float(lon[ix].max()), 2)})
    out["sonraki_1_yil_M6"] = dz
    out["gozlem"] = {"izmit": {"mw": 7.6, "lon0": 29.3, "lon1": 30.9, "boy_km": 145},
                     "duzce": {"gun": 87, "mw": 7.2, "lon0": 31.0, "lon1": 31.4}}
    f = ROOT / "reports" / "tables" / f"kirilma_yayilma_{a.ag}_{Path(a.takim).stem}.json"
    f.write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
