"""Eleman düzeyinde moment (kayma) bütçesi.

Segment başına tek değer yerine her 1 km × 1 km fay elemanı kendi
bütçesini tutar:

    açık_e(t) = max(0, açık_e(t⁻) + Ṁ0_e · Δt − boşalım_e)

* Birikim: Ṁ0_e = μ · A_e · v_e (eleman alanı ve fayın kayma hızı).
* Boşalım: her depremin momenti M0, Mw'den Wells & Coppersmith (1994)
  ile tahmin edilen kırık uzunluğu L ve genişliği W kadar bir alana
  dağıtılır. Dağılım, episantrdan yatay uzaklıkta Gauss, derinlikte
  hiposantra göre Gauss biçimlidir (kenarlara doğru azalan kayma):

      σ_yatay = sqrt((L/4)² + σ_konum²),  σ_düşey = sqrt((W/4)² + σ_derinlik²)

  σ_konum ve σ_derinlik olayın konum belirsizliğidir (``weights``
  modülündeki dönem varsayılanları); böylece konum hatası da yayılır.
* Bilinen kırık izi (``known_ruptures``) verilen depremlerde boşalım,
  izin 5 km yakınındaki elemanlara uzunluk boyunca düzgün dağıtılır.
* Sonuç "kayma açığı" olarak da verilir: açık / (μ · A) metre. Bu, o
  elemanın biriktirip henüz depremle boşaltmadığı kayma miktarıdır.

Başlangıç durumu bilinmez: açık başlangıç yılında 0 kabul edilir. Bu
yüzden değerler mutlak değil, başlangıç yılına göre görelidir. Tarihsel
katalog eklendikçe başlangıç geriye çekilebilir.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ..catalog.weights import fill_uncertainties
from ..geo import point_to_polyline_km
from .coulomb import rupture_dims_km
from .moment import MU, mw_to_m0

LAT0, LON0 = 39.0, 35.0


def _xy(lat, lon):
    return ((np.asarray(lon, float) - LON0) * 111.195 * np.cos(np.radians(LAT0)),
            (np.asarray(lat, float) - LAT0) * 111.195)


def element_budget(elements: pd.DataFrame, events: pd.DataFrame, start: pd.Timestamp,
                   end: pd.Timestamp, snapshots=(), min_mw: float = 4.0,
                   mw_col: str = "mw_best", known_ruptures: pd.DataFrame | None = None,
                   coupling: float = 1.0) -> dict:
    """Eleman bütçesini zaman içinde ilerlet.

    Dönen sözlük: ``final`` (DataFrame: m0_acc, m0_rel, deficit_Nm, slip_deficit_m,
    last_release_year) ve ``snapshots`` ({zaman: deficit dizisi}).
    """
    n = len(elements)
    ex, ey = _xy(elements["lat"], elements["lon"])
    ez = elements["depth_km"].to_numpy()
    rate = coupling * elements["m0_rate_Nm_yr"].fillna(0).to_numpy() / 365.25  # N·m/gün
    area_m2 = elements["area_km2"].to_numpy() * 1e6
    tree = cKDTree(np.column_stack([ex, ey]))

    ev = events[(events[mw_col] >= min_mw) & (events["time"] >= start) & (events["time"] < end)]
    ev = fill_uncertainties(ev).sort_values("time").reset_index(drop=True)
    t_ev = ((ev["time"] - start).dt.total_seconds() / 86400).to_numpy()
    T_end = (end - start).total_seconds() / 86400
    snaps = sorted(snapshots)
    snap_t = [((s - start).total_seconds() / 86400, s) for s in snaps if start <= s <= end]

    deficit = np.zeros(n)
    t_last = np.zeros(n)
    acc = np.zeros(n)
    rel = np.zeros(n)
    last_rel = np.full(n, np.nan)
    out_snaps = {}
    kr = {} if known_ruptures is None else {
        str(r["event_id"]): r["coords"] for _, r in known_ruptures.iterrows()}

    def advance_all(t):
        dt = t - t_last
        deficit[:] += rate * dt
        acc[:] += rate * dt
        t_last[:] = t

    si = 0
    for k in range(len(ev)):
        while si < len(snap_t) and snap_t[si][0] <= t_ev[k]:
            advance_all(snap_t[si][0])
            out_snaps[snap_t[si][1]] = deficit.copy()
            si += 1
        r = ev.iloc[k]
        mw = float(r[mw_col])
        m0 = mw_to_m0(mw)
        L, W = rupture_dims_km(mw)
        x0, y0 = _xy(r["latitude"], r["longitude"])
        z0 = float(r["depth_km"]) if pd.notna(r["depth_km"]) else 10.0
        sh = float(r["horiz_err_km"]) if pd.notna(r["horiz_err_km"]) else 10.0
        sz = float(r["depth_err_km"]) if pd.notna(r["depth_err_km"]) else 10.0
        eid = str(r.get("event_id", ""))
        if eid in kr:
            coords = kr[eid]
            idx = np.asarray(tree.query_ball_point([float(x0), float(y0)], L + 50))
            if idx.size == 0:
                continue
            d = point_to_polyline_km(elements["lat"].to_numpy()[idx], elements["lon"].to_numpy()[idx],
                                     coords)
            w = (d <= 5.0).astype(float) * area_m2[idx]
        else:
            s_h = np.sqrt((L / 4) ** 2 + sh ** 2)
            s_z = np.sqrt((W / 4) ** 2 + sz ** 2)
            idx = np.asarray(tree.query_ball_point([float(x0), float(y0)], 3 * s_h))
            if idx.size == 0:
                continue  # yakınında fay elemanı yok (faysız bölge / deniz)
            d2 = (ex[idx] - x0) ** 2 + (ey[idx] - y0) ** 2
            w = np.exp(-d2 / (2 * s_h ** 2) - (ez[idx] - z0) ** 2 / (2 * s_z ** 2)) * area_m2[idx]
        if w.sum() <= 0:
            continue
        # etkilenen elemanları bu ana kadar ilerlet, sonra boşalımı uygula
        dt = t_ev[k] - t_last[idx]
        deficit[idx] += rate[idx] * dt
        acc[idx] += rate[idx] * dt
        t_last[idx] = t_ev[k]
        r_e = m0 * w / w.sum()
        deficit[idx] = np.maximum(0.0, deficit[idx] - r_e)
        rel[idx] += r_e
        big = r_e > 0.05 * rate[idx] * 365.25  # en az ~bir aylık birikim kadar boşalım
        last_rel[idx[big]] = r["time"].year + (r["time"].dayofyear - 1) / 365.25
    while si < len(snap_t):
        advance_all(snap_t[si][0])
        out_snaps[snap_t[si][1]] = deficit.copy()
        si += 1
    advance_all(T_end)
    final = pd.DataFrame({
        "m0_acc": acc, "m0_rel": rel, "deficit_Nm": deficit,
        "slip_deficit_m": deficit / (MU * area_m2),
        "last_release_year": last_rel}, index=elements.index)
    return {"final": final, "snapshots": out_snaps, "n_events": len(ev)}


def summarize_along_strike(elements: pd.DataFrame, values: pd.DataFrame, bin_km: float = 3.0):
    """Elemanları doğrultu boyunca ``bin_km`` parçalara, derinlik boyunca ortalayarak özetle."""
    df = pd.concat([elements[["fault_id", "along_km", "lon", "lat", "j_dip", "area_km2"]], values], axis=1)
    df["bin"] = (df["along_km"] // bin_km).astype(int)
    top = df[df["j_dip"] == 0]
    rows = []
    for (fid, b), g in df.groupby(["fault_id", "bin"], sort=False):
        t = top[(top["fault_id"] == fid) & (top["bin"] == b)].sort_values("along_km")
        if len(t) < 1:
            continue
        a = g["area_km2"].to_numpy()
        rec = {"fault_id": fid, "bin": int(b),
               "coords": [[round(x, 4), round(y, 4)] for x, y in zip(t["lon"], t["lat"])]}
        for c in values.columns:
            v = g[c].to_numpy(dtype=float)
            ok = np.isfinite(v)
            rec[c] = float(np.average(v[ok], weights=a[ok])) if ok.any() else None
        rows.append(rec)
    out = pd.DataFrame(rows)
    # komşu parçalar arasında boşluk kalmasın: her parçanın sonuna bir sonrakinin ilk noktasını ekle
    for fid, g in out.groupby("fault_id"):
        idx = g.sort_values("bin").index.tolist()
        for a_, b_ in zip(idx[:-1], idx[1:]):
            out.at[a_, "coords"] = out.at[a_, "coords"] + [out.at[b_, "coords"][0]]
    return out[out["coords"].map(len) >= 2]
