"""Model kurulumu, uzun dönem istatistikleri, GPS ayarı, tarihi oynatma ve ileri bakış."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from . import hazirlik, tarih
from .simulator import MU, Model, Parametreler

BUYUK_MW = 6.8
# DEPREM_MEKANIK_SURUM=v2: güney/orta kol ve deniz altı yan fayları içeren ağ
SURUM = __import__("os").environ.get("DEPREM_MEKANIK_SURUM", "")
VERI_DIZINI = "mekanik" + (f"_{SURUM}" if SURUM else "")


class Veri:
    """Bir kez yüklenen ortak veri: hücreler, etkileşim matrisleri, D90 profili."""

    def __init__(self, root: Path):
        d = Path(root) / "data" / VERI_DIZINI
        self.cells = pd.read_parquet(d / "hucreler.parquet")
        self.Kt = np.load(d / "K_tau.npy", mmap_mode="r")
        self.Kn = np.load(d / "K_sn.npy", mmap_mode="r")
        self.d90 = pd.read_csv(d / "d90_profil.csv")
        self._K = {}

    def K(self, mu_eff: float):
        key = round(mu_eff, 3)
        if key not in self._K:
            self._K.clear()
            K = np.asfortranarray(np.asarray(self.Kt) + np.float32(mu_eff) * np.asarray(self.Kn))
            self._K[key] = K.astype(np.float32, copy=False)
        return self._K[key]


def model_kur(v: Veri, p: Parametreler, uye: int | None = None) -> Model:
    """``uye`` verilirse bilinmeyen girdiler (katman kalınlığı, levha hızı, sığ krip) örneklenir."""
    rng = None if uye is None else np.random.default_rng(10_000 + uye)
    q, _ = hazirlik.isi_akisi_profili(v.d90, v.cells, rng=rng)
    vpl = hazirlik.levha_hizi(v.cells, rng=rng)
    ek = hazirlik.sig_krip_maskesi(v.cells, rng=rng)
    p2 = replace(p, isi_akisi=q, v_pl=vpl, seed=p.seed if uye is None else 1000 + uye)
    m = Model(v.cells, None, None, p2, ek_krip=ek, K=v.K(p.mu_eff))
    m.sig_rahatla()
    return m


# ------------------------------------------------------------------ istatistik
def b_degeri(mw, m_min=4.5, m_max=6.5):
    m = np.asarray(mw)
    m = m[(m >= m_min) & (m <= m_max)]
    if len(m) < 20:
        return np.nan
    return np.log10(np.e) / (m.mean() - m_min)


def tekrarlanma(cat: pd.DataFrame, t0: float, t1: float, x_grid, mw_min=BUYUK_MW):
    """Her km için büyük depremlerin ortalama tekrarlanma aralığı ve değişkenlik katsayısı."""
    big = cat[(cat.mw >= mw_min) & (cat.t >= t0) & (cat.t <= t1)]
    rec, cov = [], []
    for x in x_grid:
        ts = np.sort(big[(big.x0 <= x) & (big.x1 >= x)].t.to_numpy())
        if len(ts) >= 3:
            d = np.diff(ts)
            rec.append(d.mean())
            cov.append(d.std() / d.mean())
        else:
            rec.append((t1 - t0) / max(len(ts), 0.5))
            cov.append(np.nan)
    return np.array(rec), np.array(cov)


def ozet(cat: pd.DataFrame, t0: float, t1: float, x_ana=(60.0, 400.0)) -> dict:
    c = cat[(cat.t >= t0) & (cat.t <= t1)]
    big = c[c.mw >= BUYUK_MW]
    xg = np.arange(x_ana[0], x_ana[1], 10.0)
    rec, cov = tekrarlanma(c, t0, t1, xg)
    return {"yil": t1 - t0, "n_olay": len(c), "b": b_degeri(c.mw),
            "n_buyuk": len(big), "mw_maks": float(c.mw.max()) if len(c) else np.nan,
            "uzunluk_ort_buyuk": float((big.x1 - big.x0).mean()) if len(big) else np.nan,
            "uzunluk_maks": float((big.x1 - big.x0).max()) if len(big) else np.nan,
            "tekrar_medyan": float(np.nanmedian(rec)), "cov_medyan": float(np.nanmedian(cov)),
            "m0_orani_buyuk": float(big.m0.sum() / max(c.m0.sum(), 1))}


# ------------------------------------------------------------------ tarihi oynatma
def tarihsel_olaylar(v: Veri, uye: int | None):
    df = tarih.tablo(v.cells)
    if uye is not None:
        rng = np.random.default_rng(20_000 + uye)
        df["x0"] = df["x0"] + rng.uniform(-1, 1, len(df)) * df["uc_km"]
        df["x1"] = df["x1"] + rng.uniform(-1, 1, len(df)) * df["uc_km"]
    return df


def zorla(m: Model, ev, derin_max=None, doldur=False):
    """Tarihsel kırılmayı uygula. ``doldur=True`` (veri özümseme): gerçek deprem, kırılan
    kilitli parçanın o anda kırılma sınırında olduğunu söyler; önce %100 doluya çekilir."""
    mask = (m.x >= ev.x0) & (m.x <= ev.x1)
    if doldur:
        d = mask & m.vw
        m.S[d] = np.maximum(m.S[d], m.S_fail[d])
    return m.kirilma(0, maske=mask, zorla=True, tarih=ev.tarih)


def ileri_olasilik(m: Model, ufuk: float, n_tekrar: int, seed: int, mw_min=BUYUK_MW,
                   x_bins=None, ilk=True):
    """Durumdan ileri koşular: her km için ``ufuk`` yıl içinde büyük kırılmaya girme olasılığı.

    Her tekrar, bilinmeyen küçük ölçekli dayanım farklarını yeniden çeker
    (gerilme düzeyi korunur). ``ilk=True``: yalnız ilk büyük kırılma sayılır.
    Dönen: (olasılık dizisi, ilk büyük kırılmaların tablosu).
    """
    if x_bins is None:
        x_bins = np.arange(0, np.ceil(m.x.max()) + 1, 1.0)
    st = m.durum()
    n_ev0 = len(m.olaylar)
    hits = np.zeros(len(x_bins) - 1)
    ilkler = []
    for r in range(n_tekrar):
        m.yukle(st)
        m.dayanimi_yenile(seed * 1000 + r)
        t_bas = m.t
        k0 = len(m.olaylar)
        bulundu = None
        # büyük deprem çıkana ya da ufuk dolana kadar 1'er yıllık dilimler
        while m.t < t_bas + ufuk - 1e-9:
            m.ilerle(min(m.t + 1.0, t_bas + ufuk))
            yeni = [e for e in m.olaylar[k0:] if e["mw"] >= mw_min]
            if yeni:
                bulundu = yeni[0]
                break
        if bulundu is not None:
            sel = (x_bins[1:] > bulundu["x0"]) & (x_bins[:-1] < bulundu["x1"])
            hits[sel] += 1
            ilkler.append({**bulundu, "dt": bulundu["t"] - t_bas, "tekrar": r})
        del m.olaylar[k0:]
    m.yukle(st)
    del m.olaylar[n_ev0:]
    return hits / n_tekrar, pd.DataFrame(ilkler)


# ------------------------------------------------------------------ 1999 sonrası GPS
GPS_PENCERE = [(1999.9, 2001), (2000, 2002), (2001, 2003), (2002, 2004), (2003, 2006),
               (2006, 2009), (2009, 2012), (2012, 2015), (2015, 2019), (2019, 2023)]
GPS_ISTASYON = {"TUBI": (40.7867, 29.4507), "ISTA": (41.1044, 29.0193)}


def gps_gozlem(root: Path) -> pd.DataFrame:
    """TUBI ve ISTA pencere hızlarının 2015-2023 hızına göre farkı (mm/yıl)."""
    from ..observations import gnss
    raw = Path(root) / "data" / "raw" / "gnss"
    steps = gnss.read_steps(raw / "steps.txt")
    rows = []
    for s in GPS_ISTASYON:
        ts = gnss.read_tenv3(raw / "tenv3" / f"{s}.EU.tenv3")
        vel = {}
        for a, b in GPS_PENCERE:
            d = ts[(ts.t >= a) & (ts.t < b)]
            if len(d) < 100:
                continue
            st = steps[(steps.site == s) & (steps.t > a) & (steps.t < b)].t.values
            for c, sc in (("e", "se"), ("n", "sn")):
                beta, _, _ = gnss._robust_fit(d.t.values, d[c].values, d[sc].values, st, (a + b) / 2)
                vel[(a, b, c)] = beta[1]
        for (a, b, c), val in vel.items():
            ref = 0.5 * (vel[(2015, 2019, c)] + vel[(2019, 2023, c)])
            rows.append({"site": s, "a": a, "b": b, "bil": c, "anomali": val - ref})
    return pd.DataFrame(rows)


def gps_model(m: Model, G: dict, t_olay: float) -> pd.DataFrame:
    """Modelde aynı pencere hız farkları. ``m`` 1999 depremlerinden hemen sonraki durumda olmalı."""
    edges = sorted({e for w in GPS_PENCERE for e in w})
    snap = {}
    for t in edges:
        if t > m.t:
            m.ilerle(t, deprem=True)
        snap[t] = np.where(m.creep, m.slip, 0.0)  # GPS'te deprem adımları ayıklanır: yalnız yavaş kayma
    rows = []
    vel = {}
    for s, g in G.items():
        for a, b in GPS_PENCERE:
            du = g @ (snap[b] - snap[a]) / (b - a) * 1000.0  # mm/yıl (e, n, u)
            vel[(s, a, b, "e")], vel[(s, a, b, "n")] = du[0], du[1]
    for (s, a, b, c), val in vel.items():
        ref = 0.5 * (vel[(s, 2015, 2019, c)] + vel[(s, 2019, 2023, c)])
        rows.append({"site": s, "a": a, "b": b, "bil": c, "anomali": val - ref})
    return pd.DataFrame(rows)
