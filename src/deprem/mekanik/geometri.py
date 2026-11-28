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

# İkinci sürüm (yan kollar). Her kol kendi x başlangıcından sayılır; kollar
# arasında büyük x boşluğu bırakılır ki kırılma komşuluğu kollar arasında
# karışmasın (kollar arası etkileşim yine elastik matristen gelir).
# (kol adı, x başlangıcı km, [(GEM kimliği, ad)])
YAN_KOLLAR = [
    ("Orta kol", 600.0, [("ME_TRCS024", "Biga-Bandırma"), ("ME_TRCS023", "Bandırma-Gemlik"),
                         ("ME_TRCS022", "İznik-Geyve")]),
    ("İmralı", 1000.0, [("ME_TRCS015", "İmralı")]),
    ("Armutlu-Çınarcık güney", 1150.0, [("ME_TRCS014", "Armutlu-Çınarcık güney")]),
]
KOL_ANA_X_MAKS = 500.0  # x bundan küçükse ana (kuzey) kol


def fay_izleri(gem_path: Path, yan_kollar: bool = False) -> list[dict]:
    """Segment izlerini batıdan doğuya sıralı (lon, lat) dizileri olarak döndür."""
    r = Region("m", 39.9, 41.3, 25.0, 32.0)
    f = gem.load_faults(Path(gem_path), r).set_index("fault_id")
    liste = [(fid, ad, lon_min, 0, None) for fid, ad, lon_min in SEGMENTLER]
    if yan_kollar:
        for k, (_, x_bas, segs) in enumerate(YAN_KOLLAR, start=1):
            for j, (fid, ad) in enumerate(segs):
                liste.append((fid, ad, None, k, x_bas if j == 0 else None))
    out = []
    for fid, ad, lon_min, kol, x_bas in liste:
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
        out.append({"fault_id": fid, "ad": ad, "coords": c, "kol": kol, "x_bas": x_bas,
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
        if iz.get("x_bas") is not None:
            x_off = float(iz["x_bas"])
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
                    rows.append({"seg": si, "seg_ad": iz["ad"], "kol": iz.get("kol", 0), "katman": kat,
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


# ---------------------------------------------------------------- yerel ince ağ
def yerel_ag(lon_c: float, lat_c: float, strike: float, L_km: float, z_max: float, h: float,
             x_merkez: float = 0.0, ad: str = "yerel"):
    """Düz, düşey bir fay parçası üzerinde h km'lik kare hücreler (ör. 100 m).

    Parça düzlemsel olduğundan etkileşim yalnız hücreler arası doğrultu farkına ve iki
    derinliğe bağlıdır (öteleme simetrisi): bir sütun alıcı için hesaplanıp bütün matris
    buradan kurulur. Aynı düzlemde normal gerilme değişimi sıfırdır, K = K_tau.
    Dönen: (hücre tablosu, K_tau float32).
    """
    import cutde.halfspace as hs

    nx = int(round(L_km / h))
    nz = int(round(z_max / h))
    phi = np.radians(strike)
    sv = np.array([np.sin(phi), np.cos(phi), 0.0])
    c0 = to_local_m(np.array([lat_c]), np.array([lon_c]), np.array([0.0]))[0]
    # kaynaklar: -(nx-1)..(nx-1) sütun; alıcı sütunu 0
    cols = np.arange(-(nx - 1), nx)
    H = h * 1e3
    tris = []
    for ci in cols:
        for k in range(nz):
            top = c0 + sv * ci * H + np.array([0, 0, -k * H])
            ts, tn = top - sv * H / 2, top + sv * H / 2
            bs, bn = ts - np.array([0, 0, H]), tn - np.array([0, 0, H])
            tris.append([[ts, tn, bn], [ts, bn, bs]])
    tris = np.asarray(tris)  # (ncol*nz, 2, 3, 3)
    ns = len(tris)
    from ..faults.coulomb import MU_SHEAR, NU, _calibrate_slip, fault_vectors  # noqa: F811
    obs_c = tris[(nx - 1) * nz:(nx - 1) * nz + nz].mean(axis=2)  # alıcı sütunu (nz, 2, 3)
    cen = tris.mean(axis=(1, 2))
    coef = _calibrate_slip(tris[0], cen[0], strike, 90.0, 180.0)
    src = tris.reshape(-1, 3, 3)
    slip_vec = np.zeros((ns * 2, 3))
    slip_vec[:, 0], slip_vec[:, 1] = coef[0], coef[1]
    nvec, svec, _ = fault_vectors(np.array([strike]), np.array([90.0]), np.array([180.0]))
    Kcol = np.zeros((nz, ns), np.float32)
    blok = 4
    for i0 in range(0, nz, blok):
        i1 = min(i0 + blok, nz)
        o = np.ascontiguousarray(obs_c[i0:i1].reshape(-1, 3))
        M = hs.strain_matrix(o, src, NU)
        E = np.einsum("aitk,tk->ait", M, slip_vec).reshape(len(o), 6, ns, 2).sum(-1)
        E = E.reshape(i1 - i0, 2, 6, ns).mean(1)
        S = hs.strain_to_stress(np.moveaxis(E, 1, -1).reshape(-1, 6), MU_SHEAR, NU).reshape(i1 - i0, ns, 6)
        xx, yy, zz, xy, xz, yz = np.moveaxis(S, -1, 0)
        nv, s_ = nvec[0], svec[0]
        t = np.stack([xx * nv[0] + xy * nv[1] + xz * nv[2], xy * nv[0] + yy * nv[1] + yz * nv[2],
                      xz * nv[0] + yz * nv[1] + zz * nv[2]], -1)
        Kcol[i0:i1] = (t * s_).sum(-1)
    Kcol = Kcol.reshape(nz, len(cols), nz)  # (z alıcı, sütun farkı, z kaynak)
    # hücre tablosu
    ix, iz = np.meshgrid(np.arange(nx), np.arange(nz), indexing="ij")
    ix, iz = ix.ravel(), iz.ravel()
    off = (ix - (nx - 1) / 2) * H
    xyz = c0[None, :] + off[:, None] * sv[None, :]
    lon = LON0 + xyz[:, 0] / (KM_PER_DEG * 1e3 * np.cos(np.radians(LAT0)))
    lat = LAT0 + xyz[:, 1] / (KM_PER_DEG * 1e3)
    cells = pd.DataFrame({"seg": 0, "seg_ad": ad, "kol": 0, "katman": "kabuk", "lon": lon, "lat": lat,
                          "x_km": x_merkez + off / 1e3, "z_ust": iz * h, "z_km": iz * h + h / 2,
                          "L_km": h, "W_km": h, "strike": strike})
    cells["alan_m2"] = h * h * 1e6
    n = len(cells)
    K = np.empty((n, n), np.float32)
    for i in range(n):
        d = ix - ix[i] + (nx - 1)
        K[i] = Kcol[iz[i], d, iz]
    return cells, K


def hucreler_kademeli(izler: list[dict], h_fonk, z_ince=16.0, z_kirilgan=20.0, z_sunek=40.0,
                      ds_derin=4.0) -> pd.DataFrame:
    """Tek ağ, kademeli kare hücreler (Ulaş: yoğun yerde 200 m, seyrek yerde 1 km'ye kadar).

    ``h_fonk(kol, x_km) -> h_km``: geçmiş deprem yoğunluğuna göre hücre boyu (bir kez belirlenir, sabit).
    0-z_ince km: h×h kare hücre; z_ince-z_kirilgan: 1 km (küçük depremler bu derinlikte yok denecek kadar
    az); z_kirilgan-z_sunek: sünek kuşak ds_derin. İz boyunca h değişirken hücre sınırları kesintisiz.
    """
    rows = []
    x_off = 0.0
    for si, iz in enumerate(izler):
        if iz.get("x_bas") is not None:
            x_off = float(iz["x_bas"])
        kol = iz.get("kol", 0)
        lon_i, lat_i, L = _yeniden_ornekle(iz["coords"], 0.05)
        s_i = np.linspace(0, L, len(lon_i))
        # iz boyunca kenarlar: her adımda yerel h kadar ilerle
        kenar = [0.0]
        while kenar[-1] < L - 1e-6:
            h = float(h_fonk(kol, x_off + kenar[-1]))
            kenar.append(min(kenar[-1] + h, L))
        if len(kenar) > 2 and kenar[-1] - kenar[-2] < 0.5 * (kenar[-2] - kenar[-3]):
            kenar.pop(-2)
        kenar = np.array(kenar)

        def nokta(s):
            return np.interp(s, s_i, lon_i), np.interp(s, s_i, lat_i)

        def ekle(e, z0, z1, dz, kat, ds_kat=None):
            for k in range(len(e) - 1):
                a, b = e[k], e[k + 1]
                lo0, la0 = nokta(a)
                lo1, la1 = nokta(b)
                strike = (np.degrees(np.arctan2((lo1 - lo0) * np.cos(np.radians(la0)), la1 - la0)) + 360) % 360
                d = dz if dz is not None else (b - a)
                for z in np.arange(z0, z1 - 1e-9, d):
                    rows.append({"seg": si, "seg_ad": iz["ad"], "kol": kol, "katman": kat,
                                 "lon": (lo0 + lo1) / 2, "lat": (la0 + la1) / 2, "x_km": x_off + (a + b) / 2,
                                 "z_ust": z, "z_km": z + d / 2, "L_km": b - a, "W_km": d, "strike": strike})

        # ince kabuk: kare hücre (W = L = h), z_ince'ye tam bölünecek h
        for k in range(len(kenar) - 1):
            a, b = kenar[k], kenar[k + 1]
            h = b - a
            nz = max(int(round(z_ince / h)), 1)
            ekle(np.array([a, b]), 0.0, z_ince, z_ince / nz, "kabuk")
        e1 = np.linspace(0, L, max(int(round(L / 1.0)), 1) + 1)
        ekle(e1, z_ince, z_kirilgan, 1.0, "kabuk")
        e4 = np.linspace(0, L, max(int(round(L / ds_derin)), 1) + 1)
        ekle(e4, z_kirilgan, z_sunek, ds_derin, "sunek")
        x_off += L
    df = pd.DataFrame(rows)
    df["alan_m2"] = df["L_km"] * df["W_km"] * 1e6
    return df


def ucgen_hazirla(cells: pd.DataFrame):
    """Etkileşim blokları için hücre üçgenleri, kayma vektörleri ve fay yönleri (bir kez)."""
    tris, coefs = _ucgenler(cells)
    n = len(cells)
    slip_vec = np.zeros((n, 2, 3))
    slip_vec[:, :, 0] = coefs[:, None, 0]
    slip_vec[:, :, 1] = coefs[:, None, 1]
    nvec, svec, _ = fault_vectors(cells["strike"].to_numpy(), np.full(n, 90.0), np.full(n, 180.0))
    return {"tris": tris, "slip": slip_vec, "obs": tris.mean(axis=2), "nvec": nvec, "svec": svec}


def etkilesim_blok(U: dict, rec: np.ndarray, src: np.ndarray, mu_eff: float, blok: int = 64) -> np.ndarray:
    """K = K_tau + mu_eff K_sn, alıcı ``rec`` × kaynak ``src`` (Pa/m, float32); etkilesim_matrisi ile aynı."""
    import cutde.halfspace as hs
    ns = len(src)
    srct = U["tris"][src].reshape(-1, 3, 3)
    sv_ = U["slip"][src].reshape(-1, 3)
    out = np.zeros((len(rec), ns), np.float32)
    for i0 in range(0, len(rec), blok):
        r = rec[i0:i0 + blok]
        o = np.ascontiguousarray(U["obs"][r].reshape(-1, 3))
        M = hs.strain_matrix(o, srct, NU)
        E = np.einsum("aitk,tk->ait", M, sv_).reshape(len(o), 6, ns, 2).sum(-1)
        E = E.reshape(len(r), 2, 6, ns).mean(1)
        S = hs.strain_to_stress(np.moveaxis(E, 1, -1).reshape(-1, 6), MU_SHEAR, NU).reshape(len(r), ns, 6)
        xx, yy, zz, xy, xz, yz = np.moveaxis(S, -1, 0)
        nv, sv = U["nvec"][r][:, None, :], U["svec"][r][:, None, :]
        t = np.stack([xx * nv[..., 0] + xy * nv[..., 1] + xz * nv[..., 2],
                      xy * nv[..., 0] + yy * nv[..., 1] + yz * nv[..., 2],
                      xz * nv[..., 0] + yz * nv[..., 1] + zz * nv[..., 2]], -1)
        out[i0:i0 + len(r)] = ((t * sv).sum(-1) + mu_eff * (t * nv).sum(-1)).astype(np.float32)
    return out
