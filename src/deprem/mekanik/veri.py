"""Gözlemden model girdisi: deprem derinlikleri (kırılgan katman), GPS hız profili."""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..geo import project_to_polyline  # noqa: F401  (ileride kullanım)

KM_PER_DEG = 111.195
SABIT_DERINLIK = {0.0, 1.0, 2.0, 5.0, 10.0, 15.0, 33.0}


def iz_uzerine_izdusum(lon, lat, izler, max_uzak_km=12.0):
    """Noktaları birleşik iz üzerine izdüşür: (x_km, uzaklık_km, işaretli dik uzaklık_km)."""
    lat0 = np.radians(40.7)
    pts = []
    x_off = 0.0
    for iz in izler:
        c = iz["coords"]
        X = c[:, 0] * KM_PER_DEG * np.cos(lat0)
        Y = c[:, 1] * KM_PER_DEG
        pts.append((X, Y, x_off))
        x_off += np.sum(np.hypot(np.diff(X), np.diff(Y)))
    px = np.asarray(lon) * KM_PER_DEG * np.cos(lat0)
    py = np.asarray(lat) * KM_PER_DEG
    best = np.full(len(px), np.inf)
    xs = np.full(len(px), np.nan)
    sgn = np.zeros(len(px))
    for X, Y, off in pts:
        s0 = np.concatenate([[0], np.cumsum(np.hypot(np.diff(X), np.diff(Y)))])
        for k in range(len(X) - 1):
            ax, ay = X[k], Y[k]
            dx, dy = X[k + 1] - ax, Y[k + 1] - ay
            L2 = dx * dx + dy * dy
            t = np.clip(((px - ax) * dx + (py - ay) * dy) / L2, 0, 1)
            qx, qy = ax + t * dx, ay + t * dy
            d = np.hypot(px - qx, py - qy)
            m = d < best
            best[m] = d[m]
            xs[m] = off + s0[k] + t[m] * np.sqrt(L2)
            cross = dx * (py - ay) - dy * (px - ax)  # >0: izin solu (kuzey)
            sgn[m] = np.sign(cross[m])
    return xs, best, best * sgn


def kirilgan_derinlik(cat: pd.DataFrame, izler, bin_km=20.0, max_uzak_km=10.0, q=0.9,
                      min_n=25, x_max=455.0):
    """İz boyunca D90 (depremlerin %90'ının üstünde kaldığı derinlik).

    Sabit atanmış derinlikler (0, 2, 5, 10, 15, 33 km) atılır. Kaynak
    önceliği: AFAD > EMSC > USGS (yerel ağ derinlikleri daha güvenilir).
    """
    d = cat.copy()
    d = d[d["depth_km"].notna() & ~d["depth_km"].astype(float).round(2).isin(SABIT_DERINLIK)]
    d = d[d["depth_km"] > 0.5]
    x, dist, _ = iz_uzerine_izdusum(d["longitude"].astype(float), d["latitude"].astype(float), izler)
    d = d.assign(x_km=x, dist=dist)
    d = d[d["dist"] <= max_uzak_km]
    rows = []
    for x0 in np.arange(0, x_max, bin_km):
        sel = d[(d["x_km"] >= x0) & (d["x_km"] < x0 + bin_km)]
        r = {"x_km": x0 + bin_km / 2, "n": len(sel)}
        for src in ("afad", "emsc", "usgs", "hepsi"):
            s = sel if src == "hepsi" else sel[sel["source"] == src]
            r[f"d90_{src}"] = float(np.quantile(s["depth_km"], q)) if len(s) >= min_n else np.nan
            r[f"n_{src}"] = len(s)
        rows.append(r)
    out = pd.DataFrame(rows)
    # tercih: AFAD, yoksa EMSC, yoksa hepsi
    out["d90"] = out["d90_afad"].fillna(out["d90_emsc"]).fillna(out["d90_hepsi"])
    return out, d
