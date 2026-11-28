"""İkinci tur 3B görünüm verisi: sınav sonu (bugün) hücre doluluğu ve 30 yıllık M>=6.8 olasılığı.

Birden çok sınav koşusunun (en iyi ayarlar × tohumlar) ortalamasını alır.
Kullanım: python scripts/mekanik2_3b.py [--desen 'buyuk_kaba_*_ana_v2.pkl']
Çıktı: reports/mekanik_3b/veri.json
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
from mekanik_3b import ABARTI, ETIKET, SEHIR, hucre_kose, kiyilar, xy  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--desen", default="buyuk_kaba_*_ana_v2.pkl")
    a = ap.parse_args()
    fs = sorted((oz.KALICI / "sinav").glob(a.desen))
    assert fs, "sınav dosyası yok"
    ag = oz.Ag(ROOT, oz.BUYUK_AG)
    c = ag.cells
    kol = c.kol.to_numpy()
    km = np.floor(c.x_km.to_numpy()).astype(int)
    dols, ps, sinifs = [], [], []
    for f in fs:
        r = pickle.loads(f.read_bytes())
        h = r["hucre"]
        dols.append(np.where(h["vw"], np.clip(h["S"] / np.maximum(h["S_fail"], 1), -0.1, 1.0), np.nan))
        sinifs.append(np.where(h["viskoz"], 2, np.where(h["vs"], 1, 0)))
        p = {(int(x.kol), int(x.km)): x.p_30yil for x in r["ileri_p"].itertuples()}
        ps.append(np.array([p.get((a_, b_), 0.0) for a_, b_ in zip(kol, km)]))
    dol = np.nanmean(np.vstack(dols), 0)
    p30 = np.mean(np.vstack(ps), 0)
    # sınıf: koşuların çoğunluğu
    S = np.vstack(sinifs)
    sinif = np.array([np.bincount(S[:, i], minlength=3).argmax() for i in range(S.shape[1])])
    K = hucre_kose(c)
    kiyi = kiyilar()
    kx = [[round(float(t), 2) for t in xy(k[:, 0], k[:, 1])[0]] for k in kiyi]
    ky = [[round(float(t), 2) for t in xy(k[:, 0], k[:, 1])[1]] for k in kiyi]
    ana = c[(c.kol == 0) & (c.z_ust == 0) & (c.katman == "kabuk")].sort_values("x_km")
    et = []
    for ad, lo in ETIKET:
        r = ana.iloc[np.argmin(abs(ana.lon - lo))]
        ex, ey = xy([r.lon], [r.lat])
        et.append({"ad": ad, "x": round(float(ex[0]), 1), "y": round(float(ey[0]), 1)})
    sh = []
    for ad, lo, la in SEHIR:
        ex, ey = xy([lo], [la])
        sh.append({"ad": ad, "x": round(float(ex[0]), 1), "y": round(float(ey[0]), 1)})
    xk = np.where(kol == 0, c.x_km.to_numpy(), c.x_km.to_numpy() - np.select(
        [kol == 1, kol == 2, kol == 3], [600, 1000, 1150], 0))
    KOL = {0: "kuzey kol", 1: "güney kol", 2: "İmralı kolu", 3: "Armutlu kolu"}
    veri = {"kose": np.round(K, 2).reshape(len(c), -1).tolist(), "sinif": sinif.tolist(),
            "dol": [None if not np.isfinite(d) else round(float(d), 3) for d in dol],
            "p30": [round(float(v), 3) for v in p30],
            "x_km": np.round(xk, 1).tolist(), "z": c.z_ust.tolist(),
            "seg": [f"{s} ({KOL[int(k)]})" for s, k in zip(c.seg_ad, kol)],
            "kiyi_x": kx, "kiyi_y": ky, "etiket": et, "sehir": sh, "abarti": ABARTI,
            "durum": f"2026 sonu, ikinci tur modeli ({len(fs)} koşunun ortalaması)"}
    out = ROOT / "reports" / "mekanik_3b"
    (out / "veri.json").write_text(json.dumps(veri, ensure_ascii=False, separators=(",", ":")))
    print("ok", len(fs), len(c), np.nanmean(dol), p30.max())


if __name__ == "__main__":
    main()
