"""İkinci tur: bugünkü (2026 sonu) km km doluluk, 30 yıllık M>=6.8 kırılma payı ve İstanbul özeti.

En iyi 5 ayar × 3 tohum (ana sürüm, düzeltilmiş 30 yıllık bakış: *_ana_v2.pkl) ortalaması.
Çıktı: reports/tables/mekanik2_bugun_km.csv, mekanik2_istanbul.csv, mekanik2_ileri_ozet.json
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
from mekanik_ikinci_tur import ISTANBUL_PARCALARI  # noqa: E402

TAB = ROOT / "reports" / "tables"
N_TEKRAR = 50


def main():
    fs = sorted((oz.KALICI / "sinav").glob("buyuk_kaba_*_ana_v2.pkl"))
    kms, ilks, kosu = [], [], []
    for j, f in enumerate(fs):
        r = pickle.loads(f.read_bytes())
        t = r["km"].copy()
        t["p_30yil"] = t["p_30yil"].fillna(0.0)
        t["kosu"] = j
        kms.append(t)
        i = r["ileri_ilk"].copy()
        if len(i):
            i["kosu"] = j
            ilks.append(i)
        kosu.append({"dosya": f.name, "J_ayar": r["J_ayar_arama"], "m68_30yil": len(i) / N_TEKRAR,
                     "p_en_az_bir": (i.tekrar.nunique() / N_TEKRAR) if len(i) else 0.0})
    K = pd.concat(kms)
    g = K.groupby(["kol", "km", "segment"])
    km = g.agg(boylam=("boylam", "mean"), enlem=("enlem", "mean"), doluluk_ort=("doluluk", "mean"),
               doluluk_min=("doluluk", "min"), doluluk_maks=("doluluk", "max"),
               p30_ort=("p_30yil", "mean"), p30_maks=("p_30yil", "max")).reset_index()
    km.to_csv(TAB / "mekanik2_bugun_km.csv", index=False, float_format="%.4f")
    I = pd.concat(ilks) if ilks else pd.DataFrame(columns=["kosu", "tekrar", "x0", "x1", "mw"])
    ana = km[km.kol == 0]
    rows = []
    for ad, (lo0, lo1) in ISTANBUL_PARCALARI.items():
        s = ana[(ana.boylam >= lo0) & (ana.boylam < lo1)]
        x0, x1 = s.km.min(), s.km.max() + 1
        d = K[(K.kol == 0) & K.km.isin(s.km)].groupby("kosu").doluluk.mean()
        e = I[(I.x1 >= x0) & (I.x0 <= x1)]
        p_k = np.array([e[e.kosu == j].tekrar.nunique() / N_TEKRAR for j in range(len(fs))])
        rows.append({"parca": ad, "boylam": f"{lo0}-{lo1}", "km": f"{x0}-{x1}",
                     "doluluk_ort": d.mean(), "doluluk_min_kosu": d.min(), "doluluk_maks_kosu": d.max(),
                     "p30_ort": p_k.mean(), "p30_min_kosu": p_k.min(), "p30_maks_kosu": p_k.max(),
                     "mw_medyan": float(e.mw.median()) if len(e) else np.nan,
                     "mw_maks": float(e.mw.max()) if len(e) else np.nan})
    ist = pd.DataFrame(rows)
    ist.to_csv(TAB / "mekanik2_istanbul.csv", index=False, float_format="%.3f")
    k = pd.DataFrame(kosu)
    oz_ = {"kosu_sayisi": len(fs), "tekrar": N_TEKRAR, "m68_30yil_ort": k.m68_30yil.mean(),
           "m68_30yil_aralik": [k.m68_30yil.min(), k.m68_30yil.max()],
           "p_en_az_bir_ort": k.p_en_az_bir.mean(), "p_en_az_bir_aralik": [k.p_en_az_bir.min(), k.p_en_az_bir.max()],
           "kosular": kosu}
    (TAB / "mekanik2_ileri_ozet.json").write_text(json.dumps(oz_, ensure_ascii=False, indent=1))
    pd.set_option("display.width", 200)
    print(ist.round(3).to_string())
    print(k.round(3).to_string())


if __name__ == "__main__":
    main()
