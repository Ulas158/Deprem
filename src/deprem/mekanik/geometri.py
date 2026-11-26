"""Fay geometrisi ve elastik etkileşim matrisi.

Fay izi GEM Küresel Aktif Fay veritabanının (Styron & Pagani 2020; EMME,
Danciu vd. 2018) Kuzey Anadolu Fayı kuzey kolu segmentlerinden kurulur:
Saros-Ganos (ME_TRCS996, 26.2°D'den doğusu), Orta Marmara-Çınarcık
(ME_TRCS013), İzmit-Sapanca-Akyazı (ME_TRCS012) ve Düzce (ME_TRCS010).

Hücreler: iz boyunca 1 km, derinlikte 0-20 km arası 1 km. 20-40 km arası
sünek (yavaş akan) kuşak 4 km x 4 km hücrelerle temsil edilir. Faylar
düşey ve sağ yanal doğrultu atımlı kabul edilir; iz kıvrımları normal
gerilme değişimi üretir ve bu Coulomb hesabına girer.

Etkileşim matrisi K[i, j]: j hücresinde 1 m sağ yanal kayma olduğunda i
hücresinde kayma yönündeki makaslama gerilmesi (Pa) ve normal gerilme
(çekme pozitif, Pa). Elastik yarı uzay, üçgen dislokasyon çözümü
(Nikkhoo & Walter 2015, ``cutde``). Her dikdörtgen iki üçgendir; gözlem
noktası iki üçgenin ağırlık merkezlerinin ortalamasıdır (kenar
tekilliğinden kaçınmak için).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..faults import gem
from ..faults.coulomb import MU_SHEAR, NU, _calibrate_slip, fault_vectors

LAT0, LON0 = 40.7, 29.0  # izdüşüm merkezi; hücre boyları da aynı enlemle hesaplanır


def to_local_m(lat, lon, depth_km):
    x = (np.asarray(lon, float) - LON0) * KM_PER_DEG * 1e3 * np.cos(np.radians(LAT0))
    y = (np.asarray(lat, float) - LAT0) * KM_PER_DEG * 1e3
    z = -np.asarray(depth_km, float) * 1000.0
    return np.column_stack([x, y, z])
from ..region import Region

KM_PER_DEG = 111.195
SEGMENTLER = [  # (GEM kimliği, ad, batı sınırı boylam)
    ("ME_TRCS996", "Ganos-Saros", 26.20),
    ("ME_TRCS013", "Orta Marmara-Çınarcık", None),
    ("ME_TRCS012", "İzmit-Sapanca-Akyazı", None),
    ("ME_TRCS010", "Düzce", None),
]


def fay_izleri(gem_path: Path) -> list[dict]:
    """Segment izlerini batıdan doğuya sıralı (lon, lat) dizileri olarak döndür."""
    r = Region("m", 40.0, 41.3, 25.0, 32.0)
    f = gem.load_faults(Path(gem_path), r).set_index("fault_id")
    out = []
    for fid, ad, lon_min in SEGMENTLER:
        c = np.asarray(f.loc[fid, "coords"], float)
        if c[0, 0] > c[-1, 0]:
            c = c[::-1]
        if lon_min is not None:
            keep = c[:, 0] >= lon_min
            i0 = max(np.argmax(keep) - 1, 0)
            # kesme noktasına doğrusal enterpolasyon
            a, b = c[i0], c[i0 + 1]
            t = (lon_min - a[0]) / (b[0] - a[0])
            c = np.vstack([a + t * (b - a), c[i0 + 1:]])
        out.append({"fault_id": fid, "ad": ad, "coords": c,
                    "gem_slip_rate": float(f.loc[fid, "slip_rate_mm_yr"])})
    return out


def _yeniden_ornekle(c, ds_km):
    lat0 = np.radians(40.7)
    x = c[:, 0] * KM_PER_DEG * np.cos(lat0)
    y = c[:, 1] * KM_PER_DEG
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    n = max(int(round(s[-1] / ds_km)), 1)
    edges = np.linspace(0, s[-1], n + 1)
    lon = np.interp(edges, s, c[:, 0])
    lat = np.interp(edges, s, c[:, 1])
    return lon, lat, s[-1]


def hucreler(izler: list[dict], ds=1.0, z_kirilgan=20.0, dz=1.0,
             z_sunek=40.0, ds_derin=4.0, dz_derin=4.0) -> pd.DataFrame:
    """Hücre tablosu. ``x_km`` tüm sistem boyunca (batıdan) iz uzunluğu."""
    rows = []
    x_off = 0.0
    for si, iz in enumerate(izler):
        for (d_s, z0, z1, d_z, kat) in [(ds, 0.0, z_kirilgan, dz, "kabuk"),
                                        (ds_derin, z_kirilgan, z_sunek, dz_derin, "sunek")]:
            lon, lat, L = _yeniden_ornekle(iz["coords"], d_s)
            seg_len = L / (len(lon) - 1)
            for k in range(len(lon) - 1):
                lo0, la0, lo1, la1 = lon[k], lat[k], lon[k + 1], lat[k + 1]
                dx = (lo1 - lo0) * np.cos(np.radians(la0))
                dy = la1 - la0
                strike = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
                for z in np.arange(z0, z1, d_z):
                    rows.append({"seg": si, "seg_ad": iz["ad"], "katman": kat,
                                 "lon": (lo0 + lo1) / 2, "lat": (la0 + la1) / 2,
                                 "x_km": x_off + (k + 0.5) * seg_len,
                                 "z_ust": z, "z_km": z + d_z / 2, "L_km": seg_len, "W_km": d_z,
                                 "strike": strike})
        x_off += _yeniden_ornekle(iz["coords"], ds)[2]
    df = pd.DataFrame(rows)
    df["alan_m2"] = df["L_km"] * df["W_km"] * 1e6
    return df


def _ucgenler(cells: pd.DataFrame):
    """Her hücre için iki üçgen (köşeler yerel metre) ve cutde kayma katsayısı."""
    xyz = to_local_m(cells["lat"], cells["lon"], cells["z_km"])
    tris = np.zeros((len(cells), 2, 3, 3))
    coefs = np.zeros((len(cells), 2))
    for i, r in enumerate(cells.itertuples()):
        phi = np.radians(r.strike)
        sv = np.array([np.sin(phi), np.cos(phi), 0.0])
        top = xyz[i] + np.array([0, 0, r.W_km * 500.0])
        top[2] = min(top[2], 0.0)
        L, W = r.L_km * 1e3, r.W_km * 1e3
        ts, tn = top - sv * L / 2, top + sv * L / 2
        bs, bn = ts - np.array([0, 0, W]), tn - np.array([0, 0, W])
        tris[i] = [[ts, tn, bn], [ts, bn, bs]]
        coefs[i] = _calibrate_slip(tris[i], xyz[i], r.strike, 90.0, 180.0)
    return tris, coefs


def etkilesim_matrisi(cells: pd.DataFrame, blok: int = 96, mu_eff: float = 0.4):
    """(K_tau, K_sn) float32 matrisler, Pa / m. Satır alıcı, sütun kaynak."""
    import cutde.halfspace as hs

    tris, coefs = _ucgenler(cells)
    n = len(cells)
    src = tris.reshape(-1, 3, 3)
    slip_vec = np.zeros((n, 2, 3))
    slip_vec[:, :, 0] = coefs[:, None, 0]
    slip_vec[:, :, 1] = coefs[:, None, 1]
    slip_vec = slip_vec.reshape(-1, 3)
    obs = tris.mean(axis=2)  # (n, 2, 3) üçgen ağırlık merkezleri
    nvec, svec, _ = fault_vectors(cells["strike"].to_numpy(), np.full(n, 90.0), np.full(n, 180.0))
    Kt = np.zeros((n, n), np.float32)
    Kn = np.zeros((n, n), np.float32)
    for i0 in range(0, n, blok):
        i1 = min(i0 + blok, n)
        o = np.ascontiguousarray(obs[i0:i1].reshape(-1, 3))
        M = hs.strain_matrix(o, src, NU)  # (no, 6, ntri, 3)
        E = np.einsum("aitk,tk->ait", M, slip_vec)  # (no, 6, ntri)
        E = E.reshape(len(o), 6, n, 2).sum(-1)  # hücre başına iki üçgen
        E = E.reshape(i1 - i0, 2, 6, n).mean(1)  # iki gözlem noktası ortalaması
        S = hs.strain_to_stress(np.moveaxis(E, 1, -1).reshape(-1, 6), MU_SHEAR, NU)
        S = S.reshape(i1 - i0, n, 6)
        xx, yy, zz, xy, xz, yz = np.moveaxis(S, -1, 0)
        nv, sv = nvec[i0:i1, None, :], svec[i0:i1, None, :]
        t = np.stack([xx * nv[..., 0] + xy * nv[..., 1] + xz * nv[..., 2],
                      xy * nv[..., 0] + yy * nv[..., 1] + yz * nv[..., 2],
                      xz * nv[..., 0] + yz * nv[..., 1] + zz * nv[..., 2]], -1)
        Kt[i0:i1] = (t * sv).sum(-1)
        Kn[i0:i1] = (t * nv).sum(-1)
    return Kt, Kn


def yuzey_yerdegistirme(cells: pd.DataFrame, lat, lon) -> np.ndarray:
    """Yüzey noktalarında birim kaymadan yerdeğiştirme: (n_nokta, 3[e,n,u], n_hucre) m/m."""
    import cutde.halfspace as hs

    tris, coefs = _ucgenler(cells)
    n = len(cells)
    src = tris.reshape(-1, 3, 3)
    slip_vec = np.zeros((n, 2, 3))
    slip_vec[:, :, 0] = coefs[:, None, 0]
    slip_vec[:, :, 1] = coefs[:, None, 1]
    obs = to_local_m(lat, lon, np.zeros(len(lat)))
    M = hs.disp_matrix(np.ascontiguousarray(obs), src, NU)  # (no, 3, ntri, 3)
    G = np.einsum("aitk,tk->ait", M, slip_vec.reshape(-1, 3))
    return G.reshape(len(obs), 3, n, 2).sum(-1)
