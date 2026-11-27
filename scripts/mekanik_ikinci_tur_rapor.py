"""İkinci tur sınav sonuçlarını özetle (yalnız 'sinav' olayları; ayar olayları ayrı satırda).

Çıktı: reports/mekanik2_*.csv, reports/mekanik2_ozet.json
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from deprem.mekanik import ozumseme as oz  # noqa: E402

REP = ROOT / "reports"


def sira_testi(u, w):
    """Ağırlıklı ortalama sıra (AUC benzeri) ve şansa (0,5) göre tek yönlü p (U(0,1) altında)."""
    u, w = np.asarray(u, float), np.asarray(w, float)
    ok = np.isfinite(u)
    u, w = u[ok], w[ok]
    if len(u) == 0:
        return np.nan, np.nan, 0
    m = np.sum(w * u) / w.sum()
    var = np.sum(w ** 2) / w.sum() ** 2 / 12.0
    return float(m), float(1 - stats.norm.cdf((m - 0.5) / np.sqrt(var))), int(len(u))


def eslesik_fark(um, us, w, n_boot=2000, seed=0):
    """Model sırası - basit kestirme sırası: ağırlıklı ortalama fark ve önyükleme p (fark <= 0)."""
    um, us, w = map(lambda a: np.asarray(a, float), (um, us, w))
    ok = np.isfinite(um) & np.isfinite(us)
    um, us, w = um[ok], us[ok], w[ok]
    if len(um) < 3:
        return np.nan, np.nan
    d = um - us
    f = np.sum(w * d) / w.sum()
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(d), len(d))
        bs.append(np.sum(w[i] * d[i]) / w[i].sum())
    return float(f), float(np.mean(np.asarray(bs) <= 0))


def kucuk_ozet():
    rows = []
    for f in sorted((oz.KALICI / "sinav").glob("kucuk_*.pkl")):
        r = pickle.loads(f.read_bytes())
        k = pd.DataFrame(r["kayit"])
        # ad: kucuk_{ag}{ek}_{sira}_{tohum}_{surum}; ek boşsa alan düşer, o yüzden sondan oku
        par = f.stem.split("_")
        ag, sira, seed, surum = par[1], par[-3], par[-2], par[-1]
        for bol in ("ayar", "sinav"):
            for kesim in ("hepsi", "<2km"):
                a = k[(k.puan == bol) & (k.tip == "kucuk")]
                if kesim == "<2km":
                    a = a[a.sigma < 2.0]
                w = 1.0 / np.maximum(a.sigma.to_numpy(), 0.1) ** 2
                um, pm, n = sira_testi(a.u_model, w)
                us, ps, _ = sira_testi(a.u_sure, w)
                fark, pf = eslesik_fark(a.u_model, a.u_sure, w)
                okt = np.isfinite(a.p_tuttu.to_numpy())
                rows.append({"ag": ag, "ayar_sirasi": sira, "tohum": seed, "bolum": bol, "kesim": kesim,
                             "n": n, "sira_model": um, "p_sans_model": pm, "sira_sure": us, "p_sans_sure": ps,
                             "model_eksi_sure": fark, "p_fark": pf,
                             "tuttu_beklenen": float(np.sum(a.p_tuttu.to_numpy()[okt])),
                             "tuttu_oran_agirlikli": float(np.sum(w[okt] * a.p_tuttu.to_numpy()[okt]) / max(w[okt].sum(), 1e-9)),
                             "tuttu_90_oran": float(np.nanmean(a.p_tuttu_90)) if len(a) else np.nan,
                             "erken_oran": float(np.nanmean(a.p_erken)) if len(a) else np.nan,
                             "doluluk_ort": float(np.nanmean(a.doluluk)) if len(a) else np.nan,
                             "J_ayar_arama": r["J_ayar_arama"]})
    df = pd.DataFrame(rows)
    df.to_csv(REP / "mekanik2_kucuk_sinav.csv", index=False)
    return df


def buyuk_ozet():
    """Büyük deprem sınavı: her sürüm × ayar × tohum için olay olay sonuçlar ve sahte depremler."""
    rows, sah = [], []
    for f in sorted((oz.KALICI / "sinav").glob("buyuk_kaba_*.pkl")):
        par = f.stem.split("_")
        sira, seed, surum = par[2], par[3], "_".join(par[4:])
        if surum == "ana":  # ana_v2 ile aynı koşu (yalnız 30 yıllık bakışı düzeltilmiş)
            continue
        r = pickle.loads(f.read_bytes())
        k = pd.DataFrame(r["kayit"])
        for _, e in k[k.tip == "buyuk"].iterrows():
            rows.append({"surum": surum.replace("_v2", ""), "ayar_sirasi": sira, "tohum": seed, "olay": e["ad"],
                         "bolum": e["puan"], "doluluk": e["doluluk"], "esik_ok": bool(e["dolu_ok"]),
                         "yer_ok": bool(e["yer_ok"]), "zaman_ok": bool(e["zaman_ok"]), "tuttu": bool(e["tuttu"]),
                         "kendiliginden_ortusme": e["ortusme"], "sira_model": e["u_model"],
                         "sira_sure": e["u_sure"]})
        sah.append({"surum": surum.replace("_v2", ""), "ayar_sirasi": sira, "tohum": seed,
                    "sahte_ayar_1719_1894": len(r["sahte_ayar"]), "sahte_sinav_1895_2026": len(r["sahte_sinav"])})
    d, s = pd.DataFrame(rows), pd.DataFrame(sah)
    d.to_csv(REP / "tables" / "mekanik2_buyuk_sinav.csv", index=False, float_format="%.3f")
    s.to_csv(REP / "tables" / "mekanik2_sahte.csv", index=False)
    return d, s


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    b, s = buyuk_ozet()
    print(b.groupby(["surum", "bolum", "olay"])[["doluluk", "tuttu", "kendiliginden_ortusme", "sira_model",
                                                 "sira_sure"]].mean().round(3).to_string())
    t = b[b.bolum == "sinav"]
    for su, g in t.groupby("surum"):
        m, p, n = sira_testi(g.groupby("olay").sira_model.mean(), np.ones(g.olay.nunique()))
        print(su, "sınav sırası", round(m, 3), "p", round(p, 3), "n", n, "tuttu", g.tuttu.mean())
    print(s.groupby("surum").mean(numeric_only=True))
    d = kucuk_ozet()
    pd.set_option("display.width", 250)
    g = d.groupby(["ag", "bolum", "kesim"])[["n", "sira_model", "sira_sure", "model_eksi_sure", "tuttu_oran_agirlikli",
                                             "erken_oran", "doluluk_ort"]].mean()
    print(g.round(3).to_string())
