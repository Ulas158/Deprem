"""Coulomb gerilme aktarımı (ΔCFF) fay elemanları üzerinde.

Bir deprem çevresindeki faylarda gerilmeyi değiştirir. Alıcı bir fayda,
kendi kayma yönündeki makaslama gerilmesi değişimi Δτ ve normal gerilme
değişimi Δσn (çekme pozitif) ile

    ΔCFF = Δτ + μ' Δσn

pozitifse fay kırılmaya yaklaşmış (depremi "öne çekme"), negatifse
uzaklaşmış (gerilme gölgesi) sayılır (King, Stein & Lin 1994). μ' efektif
sürtünme katsayısıdır ve gizli parametrelerden biridir (``mu_eff``).

Elastik yarı-uzay çözümü ``cutde`` paketinin üçgen dislokasyon
formülasyonundan (Nikkhoo & Walter 2015) gelir. Kaynak depremler,
Mw'den Wells & Coppersmith (1994) ölçeklemeleriyle boyutlandırılan,
episantra merkezlenmiş, en yakın fayın doğrultu/eğim/rake'ini alan düzgün
kaymalı dikdörtgenler olarak temsil edilir. Bu, sonlu-fay çözümleri
yerine geçen kaba bir yaklaşımdır: büyük depremlerin gerçek kırık
geometrisi (ör. 1999 İzmit, 2023 Kahramanmaraş) ileride elle girilmeli.

Kaynak kayma vektörünün cutde'nin yerel üçgen koordinatlarına dönüşümü,
her kaynakta fay düzleminin iki yanındaki yerdeğiştirme farkından sayısal
olarak kalibre edilir; böylece üçgen köşe sırası kurallarına bağımlı
kalınmaz.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

MU_SHEAR = 3.0e10
NU = 0.25
LAT0, LON0 = 39.0, 35.0  # Türkiye için yerel projeksiyon merkezi


def to_local_m(lat, lon, depth_km):
    x = (np.asarray(lon, float) - LON0) * 111195.0 * np.cos(np.radians(LAT0))
    y = (np.asarray(lat, float) - LAT0) * 111195.0
    z = -np.asarray(depth_km, float) * 1000.0
    return np.column_stack([x, y, z])


def fault_vectors(strike, dip, rake):
    """Doğu-kuzey-yukarı eksenlerinde birim vektörler.

    Dönen (n, s, d): n tavan bloğuna doğru normal, s kayma yönü (tavan
    bloğunun taban bloğuna göre hareketi, Aki & Richards rake kuralı),
    d eğim aşağı yönü.
    """
    phi, delta, lam = (np.radians(np.asarray(v, float)) for v in (strike, dip, rake))
    strike_v = np.stack([np.sin(phi), np.cos(phi), np.zeros_like(phi)], -1)
    d = np.stack([np.cos(phi) * np.cos(delta), -np.sin(phi) * np.cos(delta), -np.sin(delta)], -1)
    n = np.stack([np.cos(phi) * np.sin(delta), -np.sin(phi) * np.sin(delta), np.cos(delta)], -1)
    s = np.cos(lam)[..., None] * strike_v - np.sin(lam)[..., None] * d
    return n, s, d


def rupture_dims_km(mw, max_width_km=None, dip=90.0):
    """Wells & Coppersmith (1994) yeraltı kırık uzunluğu ve genişliği (tüm türler)."""
    L = 10 ** (-2.44 + 0.59 * mw)
    W = 10 ** (-1.01 + 0.32 * mw)
    if max_width_km is not None:
        W = min(W, max_width_km / max(np.sin(np.radians(dip)), 0.1))
    return L, W


def source_triangles(center_m, strike, dip, L_m, W_m):
    _, _, d = fault_vectors(strike, dip, 0.0)
    phi = np.radians(strike)
    sv = np.array([np.sin(phi), np.cos(phi), 0.0])
    c = np.asarray(center_m, float).copy()
    top = c - d * W_m / 2
    if top[2] > 0:  # yüzeyin üstüne taşmasın
        c[2] -= top[2]
        top = c - d * W_m / 2
    ts, tn = top - sv * L_m / 2, top + sv * L_m / 2
    bs, bn = ts + d * W_m, tn + d * W_m
    return np.array([[ts, tn, bn], [ts, bn, bs]])


def _calibrate_slip(tris, center_m, strike, dip, rake, eps=1.0):
    """İstenen global kayma vektörünü üreten cutde (strike, dip) bileşenlerini bul."""
    import cutde.halfspace as hs

    n, s, _ = fault_vectors(strike, dip, rake)
    c = np.asarray(center_m, float)
    # merkez yüzeye çok yakınsa bir miktar derine al
    probe = np.array([c + n * eps, c - n * eps])
    probe[:, 2] = np.minimum(probe[:, 2], -10.0)
    g = []
    for basis in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]):
        u = hs.disp_free(probe, tris, np.array([basis, basis]), NU)
        g.append(u[0] - u[1])
    G = np.column_stack(g)  # 3x2
    coef, *_ = np.linalg.lstsq(G, s, rcond=None)
    return coef


def stress_from_sources(obs_m: np.ndarray, sources: pd.DataFrame, chunk: int = 200_000):
    """Kaynak tablosundan (center_x/y/z_m, strike, dip, rake, L_m, W_m, slip_m)
    gözlem noktalarında toplam gerilme tensörü (Pa), sütunlar xx,yy,zz,xy,xz,yz."""
    import cutde.halfspace as hs

    tris_all, slips_all = [], []
    for _, r in sources.iterrows():
        center = np.array([r["x_m"], r["y_m"], r["z_m"]])
        tris = source_triangles(center, r["strike"], r["dip"], r["L_m"], r["W_m"])
        coef = _calibrate_slip(tris, center, r["strike"], r["dip"], r["rake"])
        slip = np.array([coef[0], coef[1], 0.0]) * r["slip_m"]
        tris_all.append(tris)
        slips_all.append(np.array([slip, slip]))
    if not tris_all:
        return np.zeros((len(obs_m), 6))
    tris = np.concatenate(tris_all)
    slips = np.concatenate(slips_all)
    out = np.zeros((len(obs_m), 6))
    for i in range(0, len(obs_m), chunk):
        strain = hs.strain_free(np.ascontiguousarray(obs_m[i:i + chunk]), tris, slips, NU)
        out[i:i + chunk] = hs.strain_to_stress(strain, MU_SHEAR, NU)
    return out


def coulomb_on_receivers(stress: np.ndarray, strike, dip, rake, mu_eff: float = 0.4):
    """ΔCFF (Pa). Δσn çekme pozitif; μ' = mu_eff (gözenek basıncı dahil efektif)."""
    n, s, _ = fault_vectors(strike, dip, rake)
    xx, yy, zz, xy, xz, yz = stress.T
    t = np.stack([xx * n[:, 0] + xy * n[:, 1] + xz * n[:, 2],
                  xy * n[:, 0] + yy * n[:, 1] + yz * n[:, 2],
                  xz * n[:, 0] + yz * n[:, 1] + zz * n[:, 2]], -1)
    tau = (t * s).sum(-1)
    sig_n = (t * n).sum(-1)
    return tau + mu_eff * sig_n, tau, sig_n


def build_sources(events: pd.DataFrame, faults: pd.DataFrame, mw_col: str = "mw_best",
                  min_mw: float = 6.0, max_dist_km: float = 30.0) -> pd.DataFrame:
    """Büyük depremleri en yakın fayın mekanizmasıyla dikdörtgen kaynaklara çevir.

    Yakınında fay olmayan depremler için düşey doğrultu atımlı (doğrultu 90°)
    varsayılır ve ``mech_assumed`` işaretlenir.
    """
    from .moment import assign_events_to_faults, mw_to_m0

    big = events[events[mw_col] >= min_mw].copy()
    if big.empty:
        return pd.DataFrame()
    fid = assign_events_to_faults(big, faults, max_dist_km)
    fmap = faults.set_index("fault_id")
    rows = []
    for (idx, ev), f in zip(big.iterrows(), fid):
        if f is not None and f in fmap.index:
            fr = fmap.loc[f]
            if isinstance(fr, pd.DataFrame):
                fr = fr.iloc[0]
            c = np.asarray(fr["coords"])
            # olayın en yakın olduğu iz parçasının doğrultusu
            seg = np.argmin(np.hypot(c[:-1, 0] - ev["longitude"], c[:-1, 1] - ev["latitude"]))
            dx = (c[seg + 1, 0] - c[seg, 0]) * np.cos(np.radians(ev["latitude"]))
            dy = c[seg + 1, 1] - c[seg, 1]
            strike = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
            dip, rake = float(fr["dip"]), float(fr["rake"]) if pd.notna(fr["rake"]) else 180.0
            wmax = float(fr["lower_depth_km"] - fr["upper_depth_km"])
            assumed = False
        else:
            strike, dip, rake, wmax, assumed = 90.0, 90.0, 180.0, 15.0, True
        mw = float(ev[mw_col])
        L, W = rupture_dims_km(mw, wmax, dip)
        slip = mw_to_m0(mw) / (MU_SHEAR * L * 1e3 * W * 1e3)
        depth = min(max(float(ev["depth_km"]) if pd.notna(ev["depth_km"]) else 10.0, 5.0), 15.0)
        xyz = to_local_m(ev["latitude"], ev["longitude"], depth)[0]
        rows.append({"time": ev["time"], "mw": mw, "latitude": ev["latitude"],
                     "longitude": ev["longitude"], "x_m": xyz[0], "y_m": xyz[1], "z_m": xyz[2],
                     "strike": strike, "dip": dip, "rake": rake, "L_m": L * 1e3, "W_m": W * 1e3,
                     "slip_m": slip, "fault_id": f, "mech_assumed": assumed})
    return pd.DataFrame(rows)


def element_dcff(elements: pd.DataFrame, sources: pd.DataFrame, t_ref: pd.Timestamp,
                 mu_eff: float = 0.4, last_only: bool = False) -> np.ndarray:
    """``t_ref`` öncesindeki kaynaklardan elemanlarda birikimli ΔCFF (bar).

    Kaynağın kendi kırık alanı içindeki elemanlarda tekil değerler oluşur;
    bunlar ±50 bar ile kırpılır.
    """
    src = sources[sources["time"] < t_ref]
    if last_only and len(src):
        src = src.loc[[src["time"].idxmax()]]
    obs = to_local_m(elements["lat"], elements["lon"], elements["depth_km"])
    stress = stress_from_sources(obs, src)
    dcff, _, _ = coulomb_on_receivers(stress, elements["strike"].to_numpy(),
                                      elements["dip"].to_numpy(),
                                      elements["rake"].fillna(180.0).to_numpy(), mu_eff)
    return np.clip(dcff / 1e5, -50, 50)
