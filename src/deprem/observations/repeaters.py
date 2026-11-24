"""Tekrarlayan depremler (repeaters): fayın derinlikteki sessiz kayma hızı.

Aynı küçük fay yamasının (asperite) çevresindeki krip tarafından tekrar
tekrar yüklenip kırılmasıyla oluşan, neredeyse aynı dalga biçimli depremler.
Her tekrar, Nadeau & Johnson (1998, BSSA) ölçeğiyle bir kayma miktarına
karşılık gelir: d [cm] = 10^(−2,36) · M0^0,17 (M0 dyn·cm).
Aile içi tekrar sıklığı hızlanırsa çevredeki fay daha hızlı kayıyor demektir
(Uchida vd. 2016, Science: Tohoku öncesi hızlanma). Kuzey Anadolu Fayı
üzerinde Marmara'da gözlenmiştir (Schmittbuhl vd. 2016, GRL; Bohnhoff vd.
2017).

İki aşama:
1. **Katalog adayları** (bu modül, her zaman çalışır): aynı noktada
   (≤ ``radius_km``), benzer büyüklükte (|ΔM| ≤ 0,3), en az 3 olay ve
   tekrar aralıkları düzenli (değişim katsayısı < 1; artçı dizileri > 1).
   Katalog konum hatası km düzeyinde olduğundan bunlar **adaydır**.
2. **Dalga biçimi doğrulaması** (``confirm_with_waveforms``): KOERI/AFAD
   FDSN'den kayıt indirilip çapraz korelasyon ≥ 0,9 aranır. Yavaştır;
   ``fetch_observations.py repeaters --confirm`` ile seçili ailelerde çalışır.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

log = logging.getLogger(__name__)


def _xyz(lat, lon, depth, lat0=39.0):
    return np.column_stack([np.asarray(lon) * 111.195 * np.cos(np.radians(lat0)),
                            np.asarray(lat) * 111.195, np.asarray(depth, dtype=float)])


def nj_slip_cm(mw):
    m0_dyncm = 10 ** (1.5 * np.asarray(mw) + 16.1)
    return 10 ** (-2.36) * m0_dyncm ** 0.17


def find_families(cat: pd.DataFrame, radius_km: float = 1.0, dmag: float = 0.3,
                  min_mag: float = 1.5, max_mag: float = 4.5, min_events: int = 3,
                  min_gap_days: float = 30.0, max_cv: float = 1.0, max_family: int = 30) -> pd.DataFrame:
    ev = cat[(cat["mag"] >= min_mag) & (cat["mag"] <= max_mag) & cat["depth_km"].notna()].copy()
    ev = ev.sort_values("time").reset_index(drop=True)
    if len(ev) < min_events:
        return pd.DataFrame()
    # AFAD derinliklerinin çoğu sabit (ör. 7 km) olduğundan yalnız yatay uzaklık kullanılır;
    # derinlik ayrıca 5 km içinde olmalı.
    xyz = _xyz(ev["latitude"], ev["longitude"], np.zeros(len(ev)))
    pairs = cKDTree(xyz[:, :2]).query_pairs(radius_km, output_type="ndarray")
    dz = np.abs(ev["depth_km"].to_numpy()[pairs[:, 0]] - ev["depth_km"].to_numpy()[pairs[:, 1]])
    pairs = pairs[dz <= 5.0]
    m = ev["mag"].to_numpy()
    t = ev["time"].astype("int64").to_numpy() / 86400e9
    ok = np.abs(m[pairs[:, 0]] - m[pairs[:, 1]]) <= dmag
    pairs = pairs[ok]
    # Zincirleme birleşme yerine tohum tabanlı aileler: her aile bir tohum olayın 1 km
    # içindeki benzer büyüklükteki olaylardır (aile çapı ≤ 2 km). Tohumun komşusu çok
    # fazlaysa orası bir sürü/kümedir; katalogdan tekrarlayan ayırt edilemez, atlanır.
    n = len(ev)
    nb = [[] for _ in range(n)]
    for a_, b_ in pairs:
        nb[a_].append(b_)
        nb[b_].append(a_)
    deg = np.array([len(x) for x in nb])
    used = np.zeros(n, bool)
    fam_id = np.full(n, -1)
    fid = 0
    for i in np.argsort(-deg):
        if used[i] or deg[i] + 1 < min_events or deg[i] + 1 > max_family:
            continue
        mem = [i] + [j for j in nb[i] if not used[j]]
        if len(mem) < min_events:
            continue
        tt = np.sort(t[mem])
        dt = np.diff(tt)
        if dt.min() < min_gap_days:
            continue  # aynı dizi içinde çok yakın olaylar: artçı/sürü
        if dt.std() / dt.mean() > max_cv:
            continue
        used[mem] = True
        fam_id[mem] = fid
        fid += 1
    ev["family"] = fam_id
    keep = list(range(fid))
    ev = ev[ev["family"].isin(keep)].copy()
    mw = ev["mw_best"] if "mw_best" in ev else ev["mag"]
    ev["slip_cm"] = nj_slip_cm(mw)
    log.info("Tekrarlayan deprem adayları: %d aile, %d olay", ev["family"].nunique(), len(ev))
    return ev[["time", "latitude", "longitude", "depth_km", "mag", "family", "slip_cm"]]


def cell_features(rep: pd.DataFrame, cells: pd.DataFrame, t0: pd.Timestamp,
                  radius_km: float = 30.0) -> pd.DataFrame:
    """``t0`` öncesi tekrarlayan deprem kayma hızı ve hızlanması (hücre merkezine ``radius_km``)."""
    idx = cells["cell_id"]
    if rep is None or rep.empty:
        return pd.DataFrame(index=idx)
    r = rep[rep["time"] < t0]
    yrs = (t0 - r["time"]).dt.days.to_numpy() / 365.25
    xy = _xyz(r["latitude"], r["longitude"], np.zeros(len(r)))[:, :2]
    cxy = _xyz(cells["lat_c"], cells["lon_c"], np.zeros(len(cells)))[:, :2]
    tree = cKDTree(xy) if len(xy) else None
    s = r["slip_cm"].to_numpy()
    fam = r["family"].to_numpy()
    rate5 = np.full(len(cells), np.nan)
    accel = np.full(len(cells), np.nan)
    nfam = np.zeros(len(cells))
    for k, c in enumerate(cxy):
        if tree is None:
            break
        ii = np.array(tree.query_ball_point(c, radius_km), dtype=int)
        if ii.size == 0:
            continue
        nfam[k] = len(np.unique(fam[ii]))
        nf = max(nfam[k], 1)
        s5 = s[ii][yrs[ii] <= 5].sum() / 5 / nf
        s1 = s[ii][yrs[ii] <= 1].sum() / nf
        s_prev = s[ii][(yrs[ii] > 1) & (yrs[ii] <= 5)].sum() / 4 / nf
        rate5[k] = s5
        accel[k] = np.log((s1 + 0.05) / (s_prev + 0.05))
    return pd.DataFrame({"rep_slip_rate_cm_yr": rate5, "rep_accel_log": accel, "rep_n_families": nfam},
                        index=idx)
