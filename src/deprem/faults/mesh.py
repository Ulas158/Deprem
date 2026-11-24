"""Fayları 2B kabuk elemanlara bölme.

Her fay segmenti, yüzey izi boyunca (doğrultu) ``ds_strike`` km ve fay
düzlemi boyunca aşağı (eğim yönü) ``ds_dip`` km'lik dikdörtgen elemanlara
bölünür; üst ve alt sismojenik derinlik arasında kalır. Varsayılan 1 km ×
1 km'dir. Bölgesel olarak daha sık ağ için ``refine`` listesine
(Region, boyut_km) eklenir.

Neden metre değil kilometre? Deprem konumlarının tipik hatası 1-10 km,
fay izlerinin harita hatası yüzlerce metre ile kilometre arasındadır.
1 km'nin altındaki elemanlar bilgi eklemez, yalnızca aynı veriyi
kopyalar. Eleman boyu yine de parametriktir; daha iyi veri (ör. InSAR,
yerel ağlar) geldiğinde küçültülebilir.

Elemanın konumu, eğim (dip) ve eğim yönü kullanılarak 3B hesaplanır.
Eğim yönü için sağ el kuralı varsayılır (izin yönüne bakıldığında fay
sağa doğru dalar); GEM'de dip_dir verilmeyen segmentlerde bu bir
varsayımdır. Düşey faylarda (dip≈90°) önemsizdir.

Bu ağ, Coulomb gerilme aktarımı (Okada) ve ETAS PR'ında kaynak/alıcı
elemanları olarak kullanılacaktır.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..region import Region
from .moment import MU

KM_PER_DEG = 111.195


def _resample_trace(coords, ds):
    """İzi ``ds`` km aralıklı noktalara böl. Dönen: (lon, lat, s_km) ve toplam uzunluk."""
    c = np.asarray(coords, dtype=float)
    lat0 = np.radians(c[:, 1].mean())
    x = c[:, 0] * KM_PER_DEG * np.cos(lat0)
    y = c[:, 1] * KM_PER_DEG
    seg = np.hypot(np.diff(x), np.diff(y))
    cum = np.concatenate([[0], np.cumsum(seg)])
    L = cum[-1]
    n = max(int(np.ceil(L / ds)), 1)
    s = np.linspace(0, L, n + 1)
    xs, ys = np.interp(s, cum, x), np.interp(s, cum, y)
    return xs, ys, s, L, lat0


def element_size_at(lon, lat, default: float, refine) -> float:
    size = default
    for reg, sz in refine or []:
        if reg.min_lat <= lat <= reg.max_lat and reg.min_lon <= lon <= reg.max_lon:
            size = min(size, sz)
    return size


def mesh_fault(row: pd.Series, ds_strike: float = 1.0, ds_dip: float = 1.0) -> pd.DataFrame:
    xs, ys, s, L, lat0 = _resample_trace(row["coords"], ds_strike)
    if L <= 0:
        return pd.DataFrame()
    dip = float(np.clip(row["dip"], 5, 90))
    top, bot = float(row["upper_depth_km"]), float(row["lower_depth_km"])
    width = (bot - top) / np.sin(np.radians(dip))
    n_dip = max(int(np.ceil(width / ds_dip)), 1)
    w_edges = np.linspace(0, width, n_dip + 1)
    w_c = (w_edges[:-1] + w_edges[1:]) / 2
    dw = width / n_dip

    xm, ym = (xs[:-1] + xs[1:]) / 2, (ys[:-1] + ys[1:]) / 2
    dx, dy = np.diff(xs), np.diff(ys)
    seg_len = np.hypot(dx, dy)
    strike = (np.degrees(np.arctan2(dx, dy)) + 360) % 360  # kuzeyden saat yönü
    s_mid = (s[:-1] + s[1:]) / 2
    dip_az = np.radians(strike + 90)  # sağ el kuralı
    horiz = w_c * np.cos(np.radians(dip))  # eğim yönünde yatay kayma (km)
    depth = top + w_c * np.sin(np.radians(dip))

    X = xm[:, None] + horiz[None, :] * np.sin(dip_az)[:, None]
    Y = ym[:, None] + horiz[None, :] * np.cos(dip_az)[:, None]
    n_s = xm.size
    i_idx, j_idx = np.meshgrid(np.arange(n_s), np.arange(n_dip), indexing="ij")
    area = seg_len[:, None] * dw * np.ones((1, n_dip))
    sr = row.get("slip_rate_mm_yr", np.nan)
    el = pd.DataFrame({
        "fault_id": row["fault_id"],
        "i_strike": i_idx.ravel(), "j_dip": j_idx.ravel(),
        "lon": (X / (KM_PER_DEG * np.cos(lat0))).ravel(), "lat": (Y / KM_PER_DEG).ravel(),
        "depth_km": np.repeat(depth[None, :], n_s, axis=0).ravel(),
        "strike": np.repeat(strike, n_dip), "dip": dip, "rake": row.get("rake", np.nan),
        "slip_type": row.get("slip_type"),
        "elem_len_km": np.repeat(seg_len, n_dip), "elem_wid_km": dw,
        "area_km2": area.ravel(),
        "along_km": np.repeat(s_mid, n_dip), "fault_len_km": L,
        "down_dip_km": np.repeat(w_c[None, :], n_s, axis=0).ravel(), "fault_wid_km": width,
        "upper_depth_km": top, "lower_depth_km": bot,
        "slip_rate_mm_yr": sr,
        "slip_rate_min": row.get("slip_rate_min", np.nan),
        "slip_rate_max": row.get("slip_rate_max", np.nan),
    })
    el["along_frac"] = el["along_km"] / L
    el["end_dist_km"] = np.minimum(el["along_km"], L - el["along_km"])
    el["depth_rel"] = (el["depth_km"] - top) / max(bot - top, 1e-6)
    # doğrultu değişimi (bükülme): komşu elemanlarla doğrultu farkı
    d_strike = np.abs((np.diff(strike, prepend=strike[0]) + 180) % 360 - 180)
    el["strike_bend_deg"] = np.repeat(d_strike, n_dip)
    el["m0_rate_Nm_yr"] = MU * el["area_km2"] * 1e6 * sr * 1e-3
    return el


def mesh_faults(faults: pd.DataFrame, ds: float = 1.0, refine=None) -> pd.DataFrame:
    """Tüm fayları ağla. ``refine``: [(Region, boyut_km), ...] bölgesel sıklaştırma."""
    parts = []
    for _, row in faults.iterrows():
        c = np.asarray(row["coords"])
        size = element_size_at(c[:, 0].mean(), c[:, 1].mean(), ds, refine)
        parts.append(mesh_fault(row, size, size))
    el = pd.concat([p for p in parts if not p.empty], ignore_index=True)
    el.insert(0, "element_id", np.arange(len(el)))
    return el
