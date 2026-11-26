"""İkinci tur: veri özümsemeli, geri saran (sıralı) ayar.

Plan: reports/mekanik_ikinci_tur_plan.md (sonuçlardan önce sabitlendi).

Akış
----
* Bütün depremler (tarihsel büyükler, aletsel M>=5, GFZ hassas katalogları,
  doğu boşluğunda ham AFAD, model dışı büyük depremler) tek bir zaman
  çizelgesinde eskiden yeniye işlenir.
* Her gerçek depremde: (1) depremden hemen önce model durumu puanlanır
  (yer, doluluk, zaman), (2) başlangıç bölgesi konum olasılığıyla ağırlıklı
  olarak sınıra çekilir, S <- S + p (A_bölge/A_hücre) (S_kırılma - S),
  (3) deprem kayması uygulanır ve yükü elastik matrisle bütün hücrelere
  aktarılır.
* Model kendi depremlerini de üretir; gerçek kayıtta olmayan büyükler
  "sahte" sayılır.
* Konum belirsizliği: katalog hatası 1 sigma; adaylar fay boyunca +-2 sigma,
  Gauss öncül; ayar döneminde öncül x doluluk olabilirliği, sınavda yalnız öncül.
"""
from __future__ import annotations

import json
import math
import os
import pickle
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from . import hazirlik, tarih, termal
from .simulator import Model, Parametreler

MU_YENI = 3.3e10           # Pa (6 km altı ölçümleri 32-35 GPa)
MU_MATRIS = 3.0e10         # etkileşim matrisinin hesaplandığı değer
KALICI = Path(os.environ.get("DEPREM_KALICI", "/mnt/project-files/deprem/mekanik_ikinci_tur"))

# ------------------------------------------------------------------ parametreler
# (ad, alt, üst, başlangıç, logaritmik)  -- başlangıçlar yalnız literatürden
PARAMETRELER = [
    ("mu_eff", 0.2, 0.6, 0.4, False),          # King, Stein & Lin 1994
    ("dtau0", 5e6, 50e6, 10e6, True),          # Allmann & Shearer 2009
    ("sigma_h", 0.05, 0.5, 0.3, False),        # Mai & Beroza 2002
    ("eta_dyn", 0.0, 0.6, 0.3, False),         # Di Toro vd. 2011
    ("kalan", -0.2, 0.5, -0.1, False),         # deprem sonu kalan seviye (Madariaga 1976)
    ("kayma_ref", 0.2, 4.0, 1.0, True),        # Tinti vd. 2005
    ("tM_ref", 1.0, 50.0, 10.0, True),         # Hearn vd. 2009; Ergintav vd. 2009
    ("krip_carpan", 0.3, 5.0, 1.0, True),      # Blanpied vd. 1995
    ("sigma_grad", 1e6, 17e6, 5e6, True),      # Rice 1992 (17 = hidrostatik)
    ("kalinlik_u", 0.0, 1.0, 0.5, False),      # bölge aralıkları içinde konum
    ("sig_krip_z", 3.0, 8.0, 5.5, False),      # km, 27,6-28,05°D (Yamamoto vd. 2019; Becker vd. 2023)
    ("v_pl_marmara", 12.5e-3, 17e-3, 15e-3, False),  # m/yıl (Ergintav vd. 2014; Amemoutou vd. 2023)
]
P_AD = [p[0] for p in PARAMETRELER]

# kırılgan katman kalınlığı bölgeleri (ana kol, boylam): (lon0, lon1, D_alt, D_üst) km
KALINLIK_BOLGE = [(27.75, 28.60, 10.0, 12.0),   # Orta Marmara ve Kumburgaz
                  (28.60, 29.10, 10.0, 18.0),   # Adalar
                  (25.0, 27.50, 18.0, 24.0)]    # batı (Ganos-Saros)
SIG_KRIP_BOLGE = (27.60, 28.05)
V_PL_MARMARA_LON = (27.50, 29.50)  # Ganos doğusu - Hersek


def birim_to_theta(u: np.ndarray, adlar=P_AD) -> dict:
    th = {}
    for ad, x in zip(adlar, u):
        _, lo, hi, _, lg = PARAMETRELER[P_AD.index(ad)]
        x = float(np.clip(x, 0, 1))
        th[ad] = math.exp(math.log(lo) + x * (math.log(hi) - math.log(lo))) if lg else lo + x * (hi - lo)
    return th


def theta_to_birim(th: dict, adlar=P_AD) -> np.ndarray:
    out = []
    for ad in adlar:
        _, lo, hi, _, lg = PARAMETRELER[P_AD.index(ad)]
        v = th[ad]
        out.append((math.log(v) - math.log(lo)) / (math.log(hi) - math.log(lo)) if lg else (v - lo) / (hi - lo))
    return np.array(out)


def baslangic_theta() -> dict:
    return {p[0]: p[3] for p in PARAMETRELER}


# ------------------------------------------------------------------ ağlar
class Ag:
    """Hücreler + etkileşim matrisi. ``ad``: 'kaba' (2 km, yan kollar) | 'tam' (1 km, yan kollar)
    | yerel ağ adı ('silivri', 'cinarcik')."""

    def __init__(self, root: Path, ad: str):
        self.root, self.ad = Path(root), ad
        if ad in ("kaba", "tam"):
            d = self.root / "data" / ("mekanik_v2_kaba" if ad == "kaba" else "mekanik_v2")
            self.cells = pd.read_parquet(d / "hucreler.parquet")
            self.Kt = np.load(d / "K_tau.npy", mmap_mode="r")
            self.Kn = np.load(d / "K_sn.npy", mmap_mode="r")
            self.d90 = pd.read_csv(d / "d90_profil.csv")
            self.yerel = None
        else:
            self.cells, Kt, self.yerel = yerel_ag_yukle(self.root, ad)
            self.Kt, self.Kn = Kt, None
            self.d90 = pd.read_csv(self.root / "data" / "mekanik_v2" / "d90_profil.csv")
        self._K = {}
        self.h = float(np.median(self.cells.loc[self.cells.katman == "kabuk", "W_km"]))

    def K(self, mu_eff: float):
        key = 0.0 if self.Kn is None else round(mu_eff, 3)
        if key not in self._K:
            self._K.clear()
            s = np.float32(MU_YENI / MU_MATRIS)
            if self.Kn is None:
                K = self.Kt  # yerel ağ: diskte ölçekli ve sütun öncelikli
            else:
                K = np.asfortranarray((np.asarray(self.Kt) + np.float32(mu_eff) * np.asarray(self.Kn)) * s)
            self._K[key] = K.astype(np.float32, copy=False)
        return self._K[key]


YEREL_AGLAR = {
    # ad: (merkez boylam, enlem, doğrultu, boy km, derinlik km, hücre km). Plan 100 m idi; ayar
    # dönemindeki GFZ olaylarının konum hatası >= 2 km olduğu ve 100 m ağda arama ~48 saat CPU
    # gerektirdiği için 200 m kullanıldı (raporda belirtilir).
    # Silivri önü: 2019 dizisi (Kumburgaz parçası); Çınarcık: Adalar-Çınarcık geçişi
    "silivri": (28.20, 40.865, 82.0, 14.0, 14.0, 0.2),
    "cinarcik": (29.05, 40.78, 110.0, 14.0, 14.0, 0.2),
}


def yerel_ag_yukle(root: Path, ad: str):
    from . import geometri
    d = Path(root) / "data" / f"mekanik_yerel_{ad}"
    lon, lat, strike, L, zmax, h = YEREL_AGLAR[ad]
    if not (d / "K.npy").exists():
        d.mkdir(parents=True, exist_ok=True)
        ana = pd.read_parquet(Path(root) / "data" / "mekanik_v2" / "hucreler.parquet")
        ana = ana[(ana.kol == 0) & (ana.z_ust == 0)]
        xm = float(np.interp(lon, ana.sort_values("lon").lon, ana.sort_values("lon").x_km))
        c, K = geometri.yerel_ag(lon, lat, strike, L, zmax, h, x_merkez=xm, ad=ad)
        c.to_parquet(d / "hucreler.parquet")
        # kayma modülü ölçeği uygulanmış, sütun öncelikli: süreçler arasında paylaşılan bellek eşlemesi
        np.save(d / "K.npy", np.asfortranarray(K * np.float32(MU_YENI / MU_MATRIS)))
        del K
    c = pd.read_parquet(d / "hucreler.parquet")
    K = np.load(d / "K.npy", mmap_mode="r")
    return c, K, {"lon": lon, "lat": lat, "strike": strike, "L": L, "h": h}


# ------------------------------------------------------------------ model kurulumu
def kalinlik(ag: Ag, u: float) -> np.ndarray:
    """Her hücre için kırılgan taban derinliği (D90, km)."""
    c = ag.cells
    d90 = ag.d90
    a, e, hh = d90["d90_afad"].to_numpy(), d90["d90_emsc"].to_numpy(), d90["d90_hepsi"].to_numpy()
    lo = np.where(np.isfinite(a) & np.isfinite(e), np.minimum(a, e), hh)
    hi = np.where(np.isfinite(a) & np.isfinite(e), np.maximum(a, e), hh)
    dd = lo + u * (hi - lo)
    ok = np.isfinite(dd)
    x = hazirlik._ana_x(c) if ag.yerel is None else c["x_km"].to_numpy()
    D = np.interp(x, d90["x_km"].to_numpy()[ok], dd[ok])
    ana = (c["kol"] == 0).to_numpy()
    lon = c["lon"].to_numpy()
    for l0, l1, dlo, dhi in KALINLIK_BOLGE:
        s = ana & (lon >= l0) & (lon < l1)
        D[s] = dlo + u * (dhi - dlo)
    return np.clip(D, 8.0, 24.0)


def levha_hizi(ag: Ag, v_marmara: float) -> np.ndarray:
    c = ag.cells
    if ag.yerel is not None:
        return np.full(len(c), v_marmara)
    v = hazirlik.levha_hizi(c)
    lon = c["lon"].to_numpy()
    s = (c["kol"] == 0).to_numpy() & (lon >= V_PL_MARMARA_LON[0]) & (lon <= V_PL_MARMARA_LON[1])
    v[s] = v_marmara
    return v


def model_kur(ag: Ag, th: dict, seed: int) -> Model:
    c = ag.cells
    D = kalinlik(ag, th["kalinlik_u"])
    q = termal.isi_akisi_derinlikten(D)
    lon = c["lon"].to_numpy()
    ek = (c["kol"] == 0).to_numpy() & (lon >= SIG_KRIP_BOLGE[0]) & (lon <= SIG_KRIP_BOLGE[1]) & \
        (c["z_km"].to_numpy() < th["sig_krip_z"])
    p = Parametreler(mu_eff=th["mu_eff"], dtau0=th["dtau0"], sigma_h=th["sigma_h"], eta_dyn=th["eta_dyn"],
                     asiri=-th["kalan"], kayma_ref=th["kayma_ref"], sigma_grad=th["sigma_grad"],
                     tM_ref=th["tM_ref"], krip_carpan=th["krip_carpan"], isi_akisi=q,
                     v_pl=levha_hizi(ag, th["v_pl_marmara"]), seed=seed, mu_kayma=MU_YENI)
    m = Model(c, None, None, p, ek_krip=ek, K=ag.K(th["mu_eff"]))
    m.tM_adim_alt = 0.05
    m.sig_rahatla()
    return m


# ------------------------------------------------------------------ depremler
def L_kirilma(mw):
    """Kırılma boyu (km), M3 ~0,5 km, M4 ~1,6, M5 ~5, M6 ~16 (Wells & Coppersmith mertebesi)."""
    return 10 ** (0.5 * np.asarray(mw, float) - 1.8)


def gfz_katalog(root: Path) -> pd.DataFrame:
    """GFZ yeniden konumlanmış katalogları: Becker 2006-2023, Bentz Silivri 2019, Martínez-Garzón Armutlu 2019."""
    g = Path(root) / "data" / "raw" / "gfz"
    fr = []
    for f in ("2023-002_Becker-et-al_Marmara_Sea_catalog_2006_2020.csv",
              "2025-003_Becker-et-al_Marmara_Sea_catalog_2021_2023.csv"):
        d = pd.read_csv(g / f, comment="#")
        fr.append(pd.DataFrame({"time": pd.to_datetime(d["origin"], utc=True), "lat": d["latitude"],
                                "lon": d["longitude"], "z": d["depth"], "mw": d["magnitude"],
                                "sigma": d["maj_error"], "kaynak": "becker"}))
    # Bentz: şablon eşleme; hipoDD ile yeniden konumlananlar kendi hatasıyla, diğerleri şablon
    # konumunu taşır (hata 1 km varsayıldı)
    b = pd.read_csv(g / "bentz" / "2020-006_Bentz-et-al_Template_Matching_catalog_CentralMarmaraSea_2019.txt",
                    sep=r"\s+", comment="#", header=None, skiprows=1)
    tb = pd.to_datetime(dict(year=b[1], month=b[2], day=b[3], hour=b[4], minute=b[5]), utc=True) + \
        pd.to_timedelta(b[6], unit="s")
    fr.append(pd.DataFrame({"time": tb, "lat": b[8], "lon": b[9], "z": b[10], "mw": b[11],
                            "sigma": 1.0, "kaynak": "bentz"}))
    ai = pd.read_csv(g / "armutlu" / "2021-004_Martinez-Garcon-et-al_Initial_seismicity_catalog.txt",
                     sep=";", comment="#", header=None)
    ar = pd.read_csv(g / "armutlu" / "2021-004_Martinez-Garcon-et-al_Relocated_seismicity_catalog.txt",
                     sep=";", comment="#", header=None)
    ai = ai.set_index(0)
    ar = ar.set_index(0)
    ta = pd.to_datetime(dict(year=ai[1], month=ai[2], day=ai[3], hour=ai[4], minute=ai[5]), utc=True) + \
        pd.to_timedelta(ai[6], unit="s")
    a = pd.DataFrame({"time": ta, "lat": ai[8], "lon": ai[9], "z": ai[10], "mw": ai[11], "sigma": ai[12],
                      "kaynak": "armutlu"})
    j = ar.index.intersection(a.index)
    a.loc[j, "lat"], a.loc[j, "lon"], a.loc[j, "z"], a.loc[j, "sigma"] = ar.loc[j, 1], ar.loc[j, 2], ar.loc[j, 3], ar.loc[j, 4]
    fr.append(a.reset_index(drop=True))
    df = pd.concat(fr, ignore_index=True)
    df = df[np.isfinite(df.mw) & np.isfinite(df.lat)]
    # çakışanlar: 5 s ve 10 km içinde aynı olay; hatası küçük olan kalır
    df = df.sort_values("time").reset_index(drop=True)
    t = df.time.astype("int64").to_numpy() / 1e9
    keep = np.ones(len(df), bool)
    for i in range(len(df)):
        if not keep[i]:
            continue
        j = i + 1
        while j < len(df) and t[j] - t[i] < 5:
            if keep[j] and abs(df.lat[j] - df.lat[i]) * 111 < 10 and abs(df.lon[j] - df.lon[i]) * 84 < 10:
                if df.sigma[j] < df.sigma[i]:
                    keep[i] = False
                    break
                keep[j] = False
            j += 1
    df = df[keep].reset_index(drop=True)
    df["sigma"] = df["sigma"].clip(lower=0.05)
    df["t"] = df.time.dt.year + (df.time.dt.dayofyear - 0.5 + df.time.dt.hour / 24) / 365.25
    return df


def izdusum(ag: Ag, lon, lat, max_uzak=10.0):
    """Noktaları en yakın fay koluna izdüşür: (x_km, dik uzaklık km, kol)."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    kx = 111.195 * np.cos(np.radians(40.7))
    if ag.yerel is not None:
        y = ag.yerel
        phi = np.radians(y["strike"])
        dx, dy = (lon - y["lon"]) * kx, (lat - y["lat"]) * 111.195
        s = dx * np.sin(phi) + dy * np.cos(phi)
        n = -dx * np.cos(phi) + dy * np.sin(phi)
        x0 = float(ag.cells.x_km.mean())
        return x0 + s, np.abs(n), np.zeros(len(lon), int)
    c = ag.cells[(ag.cells.z_ust == 0) & (ag.cells.katman == "kabuk")]
    cx, cy = c.lon.to_numpy() * kx, c.lat.to_numpy() * 111.195
    phi = np.radians(c.strike.to_numpy())
    px, py = lon * kx, lat * 111.195
    d = np.hypot(px[:, None] - cx[None], py[:, None] - cy[None])
    k = np.argmin(d, axis=1)
    dx, dy = px - cx[k], py - cy[k]
    s = dx * np.sin(phi[k]) + dy * np.cos(phi[k])
    n = -dx * np.cos(phi[k]) + dy * np.sin(phi[k])
    xs = c.x_km.to_numpy()[k] + np.clip(s, -c.L_km.to_numpy()[k], c.L_km.to_numpy()[k])
    return xs, np.abs(n), c.kol.to_numpy()[k]


class _YuzeyAg:
    def __init__(self, cells):
        self.cells, self.yerel = cells, None


@dataclass
class Olay:
    t: float
    ad: str
    tip: str            # 'buyuk' (kırılma boyu bilinen), 'kucuk' (noktasal), 'dis' (model dışı)
    mw: float
    sigma: float        # km, konum hatası (1 sigma)
    kol: int
    x: float            # en olası başlangıç (fay boyu, km)
    z: float
    x0: float = np.nan  # büyük olay kırılma sınırları
    x1: float = np.nan
    uc_km: float = 0.0
    sinif: str = "modern"   # 'tarihsel' | 'erken' | 'modern' : eşik ve zaman toleransı
    puan: str = ""          # '' | 'ayar' | 'sinav'
    dis: dict | None = None  # dış kaynak: {'dtau','dsn'} hücre başına
    kaynak: str = ""


ESIK = {"tarihsel": 0.90, "erken": 0.95, "modern": 0.98}
TOL = {"tarihsel": 20.0, "erken": 10.0, "modern": 5.0}

# aletsel dönemden önce ya da yayınlanmış hiposantrı olan büyük olaylar (lon, lat, derinlik)
HIPOSANTR = {"1999-08-17": (29.97, 40.76, 15.0), "1999-11-12": (31.15, 40.77, 12.0),
             "2022-11-23": (30.98, 40.84, 6.8), "2025-04-23": (28.18, 40.86, 6.9)}

# Model dışı büyük depremler (yaklaşık düzlemler; Ambraseys 2002, Ambraseys & Jackson 1998,
# Barka 1996, Kiratzi 2014 vb. mertebe). (tarih, ad, lon0, lat0, lon1, lat1, Mw, eğim, kayma açısı)
DIS_DEPREMLER = [
    ("1855-02-28", "1855 Bursa", 28.70, 40.17, 29.40, 40.13, 7.0, 90, 180),
    ("1935-01-04", "1935 Erdek-Marmara Adası", 27.30, 40.40, 27.70, 40.42, 6.4, 90, 180),
    ("1944-02-01", "1944 Bolu-Gerede", 31.55, 40.75, 32.80, 40.85, 7.4, 90, 180),
    ("1953-03-18", "1953 Yenice-Gönen", 27.15, 39.92, 27.75, 40.10, 7.2, 90, 180),
    ("1957-05-26", "1957 Abant", 30.85, 40.70, 31.55, 40.65, 7.0, 90, 180),
    ("1964-10-06", "1964 Manyas", 27.95, 40.12, 28.45, 40.15, 6.9, 45, -90),
    ("1967-07-22", "1967 Mudurnu", 30.45, 40.70, 31.05, 40.58, 7.1, 90, 180),
    ("1975-03-27", "1975 Saros", 25.90, 40.42, 26.25, 40.50, 6.6, 90, 180),
    ("2014-05-24", "2014 Kuzey Ege (Gökçeada)", 25.20, 40.25, 26.00, 40.35, 6.9, 90, 180),
]
ETKI_ESIGI = 0.01e6  # Pa (0,1 bar)


def dis_kaynak_etkisi(cells: pd.DataFrame, lon0, lat0, lon1, lat1, mw, dip, rake, mu_eff=0.4):
    """Dikdörtgen kaynağın her hücrede (τ, σn) değişimi, Pa. Kayma = M0 / (μ L W)."""
    import cutde.halfspace as hs
    from . import geometri
    from ..faults.coulomb import NU, fault_vectors
    p0 = geometri.to_local_m(np.array([lat0]), np.array([lon0]), np.array([0.0]))[0]
    p1 = geometri.to_local_m(np.array([lat1]), np.array([lon1]), np.array([0.0]))[0]
    L = np.linalg.norm(p1 - p0)
    W = 15e3 / np.sin(np.radians(dip))
    m0 = 10 ** (1.5 * mw + 9.1)
    slip = m0 / (MU_MATRIS * L * W)
    sv = (p1 - p0) / L
    # eğim yönü: doğrultunun sağı (aşağı)
    dd = np.array([sv[1], -sv[0], 0.0])
    down = dd * np.cos(np.radians(dip)) + np.array([0, 0, -1.0]) * np.sin(np.radians(dip))
    a, b = p0, p1
    c_, d_ = p1 + down * W, p0 + down * W
    tris = np.array([[a, b, c_], [a, c_, d_]])
    # kayma vektörü (doğrultu, eğim, açılma) cutde'de üçgen yerel çerçevesinde: [strike-slip, dip-slip, tensile]
    r = np.radians(rake)
    s = np.array([slip * np.cos(r), slip * np.sin(r), 0.0])
    obs = geometri.to_local_m(cells["lat"], cells["lon"], cells["z_km"])
    M = hs.strain_matrix(np.ascontiguousarray(obs), tris, NU)  # (no, 6, 2, 3)
    E = np.einsum("aitk,k->ai", M, s)
    S = hs.strain_to_stress(E, MU_MATRIS, NU)
    nv, svv, _ = fault_vectors(cells["strike"].to_numpy(), np.full(len(cells), 90.0), np.full(len(cells), 180.0))
    xx, yy, zz, xy, xz, yz = S.T
    t = np.stack([xx * nv[:, 0] + xy * nv[:, 1] + xz * nv[:, 2], xy * nv[:, 0] + yy * nv[:, 1] + yz * nv[:, 2],
                  xz * nv[:, 0] + yz * nv[:, 1] + zz * nv[:, 2]], -1)
    dtau, dsn = (t * svv).sum(-1), (t * nv).sum(-1)
    # yaklaşık düzlemin yakın alanı (kaynak modele birkaç km'den yakınsa) çözülemez: ±1 MPa ile sınırla
    return np.clip(dtau, -1e6, 1e6), np.clip(dsn, -1e6, 1e6), float(slip)


def _ondalik(ts: pd.Timestamp) -> float:
    return ts.year + (ts.dayofyear - 0.5 + ts.hour / 24) / 365.25


def cizelge(ag: Ag, root: Path, t0: float, t1: float, ayar_bitis: float, mw_min_kucuk: float,
            gfz_mw_min=None, puan_tarihsel=("1766-05-22", "1766-08-05", "1894-07-10"),
            sinav_tarihsel=("1912-08-09", "1999-08-17", "1999-11-12", "2022-11-23", "2025-04-23"),
            dis=True) -> tuple[list[Olay], dict]:
    """Zaman çizelgesi. Dönen: (olaylar, rapor sayıları)."""
    ol: list[Olay] = []
    say = {}
    c = ag.cells
    xmin, xmax = c.x_km.min(), c.x_km.max()
    # 1) tarihsel büyükler (yalnız büyük ağda)
    if ag.yerel is None:
        df = tarih.tablo(c[c.kol == 0])
        for r in df.itertuples():
            if not (t0 <= r.t <= t1):
                continue
            sinif = "tarihsel" if r.t < 1900 else ("erken" if r.t < 1990 else "modern")
            pu = "ayar" if r.tarih in puan_tarihsel else ("sinav" if r.tarih in sinav_tarihsel else "")
            if r.tarih in HIPOSANTR:
                lo, la, z = HIPOSANTR[r.tarih]
                x, _, _ = izdusum(ag, [lo], [la])
                x, sig = float(x[0]), 2.0
            else:
                x, sig, z = 0.5 * (r.x0 + r.x1), np.nan, 10.0
            ol.append(Olay(r.t, r.ad, "buyuk", r.mw, sig, 0, x, z, r.x0, r.x1, r.uc_km, sinif, pu,
                           kaynak="tarih"))
        say["tarihsel"] = sum(o.tip == "buyuk" for o in ol)
    # 2) aletsel: USGS M>=5 1905-2005 (GFZ öncesi), ham AFAD/EMSC (GFZ dışı bölge ve 2024 sonrası)
    raw = pd.read_parquet(Path(root) / "data" / "mekanik" / "katalog_marmara.parquet")
    raw = raw.assign(t=raw["time"].map(_ondalik), lon=raw.longitude.astype(float), lat=raw.latitude.astype(float),
                     z=raw.depth_km.astype(float).fillna(10.0), m=raw.mw.astype(float))
    tablo_t = [o.t for o in ol]
    secim = []
    if ag.yerel is None:
        u = raw[(raw.source == "usgs") & (raw.m >= 5.0) & (raw.t < 2006.0)]
        for r in u.itertuples():
            if any(abs(r.t - tt) < 0.1 for tt in tablo_t):
                continue  # tarihsel tablodaki olayın kendisi
            secim.append((r.t, f"USGS M{r.m:.1f} {str(r.time)[:10]}", r.lon, r.lat, r.z, r.m,
                          10.0 if r.t < 1999 else 5.0, "usgs"))
        # 2006 sonrası GFZ kapsamı dışı (30,65°D doğusu) ve 2024 sonrası: ham AFAD, yoksa EMSC
        h = raw[(raw.t >= 2006.0) & (raw.m >= max(mw_min_kucuk, 3.0)) & ((raw.lon > 30.65) | (raw.t >= 2024.0))]
        h = h.sort_values(["t", "source"])  # afad < emsc < usgs (alfabetik); ilk geleni tut
        son = []
        for r in h.itertuples():
            if any(abs(r.t - a) * 365.25 * 86400 < 20 and abs(r.lat - b) < 0.2 and abs(r.lon - cc) < 0.25
                   for a, b, cc in son[-20:]):
                continue
            son.append((r.t, r.lat, r.lon))
            if any(abs(r.t - tt) < 0.01 for tt in tablo_t):
                continue
            secim.append((r.t, f"{r.source.upper()} M{r.m:.1f} {str(r.time)[:10]}", r.lon, r.lat, r.z, r.m,
                          5.0, str(r.source)))
    # 3) GFZ hassas katalogları
    g = gfz_katalog(root)
    if isinstance(gfz_mw_min, dict):  # kaynak başına tamlık büyüklüğü
        g = g[g.mw >= g.kaynak.map(gfz_mw_min).fillna(9.0)]
    else:
        g = g[g.mw >= (mw_min_kucuk if gfz_mw_min is None else gfz_mw_min)]
    g = g[(g.t >= t0) & (g.t <= t1)]
    for r in g.itertuples():
        if any(abs(r.t - tt) < 0.01 for tt in tablo_t):
            continue
        secim.append((r.t, f"{r.kaynak} M{r.mw:.1f} {str(r.time)[:19]}", r.lon, r.lat, r.z, r.mw, r.sigma, r.kaynak))
    if secim:
        s = pd.DataFrame(secim, columns=["t", "ad", "lon", "lat", "z", "mw", "sigma", "kaynak"])
        s = s[(s.t >= t0) & (s.t <= t1)]
        x, dperp, kol = izdusum(ag, s.lon, s.lat)
        s = s.assign(x=x, d=dperp, kol=kol)
        if ag.yerel is not None:
            # yerel ağ yalnız ana kolu temsil eder: büyük ağda en yakın kolu ana kol olmayanlar atılır
            if not hasattr(Ag, "_tam_yuzey"):
                Ag._tam_yuzey = _YuzeyAg(pd.read_parquet(Path(root) / "data" / "mekanik_v2" / "hucreler.parquet"))
            _, d_t, kol_t = izdusum(Ag._tam_yuzey, s.lon, s.lat)
            s = s[(kol_t == 0)]
        uzak = s[s.d > 10.0]
        say["uzak_olay"] = int(len(uzak))
        say["uzak_m0"] = float((10 ** (1.5 * uzak.mw + 9.1)).sum())
        s = s[(s.d <= 10.0) & (s.x >= xmin - 1) & (s.x <= xmax + 1)]
        for r in s.itertuples():
            pu = ("ayar" if r.t < ayar_bitis else "sinav") if r.mw >= mw_min_kucuk else ""
            ol.append(Olay(r.t, r.ad, "kucuk", r.mw, float(r.sigma), int(r.kol), float(r.x),
                           float(np.clip(r.z, 0.5 * ag.h, 19.5)), sinif="modern" if r.t > 1990 else "erken",
                           puan=pu if r.kaynak not in ("usgs",) else "", kaynak=r.kaynak))
        say["kucuk"] = int(len(s))
        say["kucuk_kaynak"] = s.kaynak.value_counts().to_dict()
    # 4) model dışı büyükler (yalnız büyük ağ): 0,1 bar kuralı
    if dis and ag.yerel is None:
        kept, drop = [], []
        for tr, ad, lo0, la0, lo1, la1, mw, dip, rake in DIS_DEPREMLER:
            t = tarih.ondalik_yil(tr)
            if not (t0 <= t <= t1):
                continue
            dt_, dn_, slip = dis_kaynak_etkisi(c, lo0, la0, lo1, la1, mw, dip, rake)
            etki = np.abs(dt_ + 0.4 * dn_).max()
            (kept if etki >= ETKI_ESIGI else drop).append((ad, round(etki / 1e6, 4)))
            if etki >= ETKI_ESIGI:
                ol.append(Olay(t, ad, "dis", mw, np.nan, -1, np.nan, np.nan, sinif="tarihsel",
                               dis={"dtau": dt_.astype(np.float32), "dsn": dn_.astype(np.float32)}, kaynak="dis"))
        say["dis_eklenen"], say["dis_atlanan"] = kept, drop
    ol.sort(key=lambda o: o.t)
    return ol, say


# ------------------------------------------------------------------ adaylar
@dataclass
class Aday:
    p: np.ndarray            # öncül olasılık
    bolge: list              # her aday için (hücre indeksleri, alan oranı) başlangıç bölgesi
    yama: list               # her aday için (hücre indeksleri, alan oranı) kırılma yaması
    merkez: np.ndarray       # her adayın merkez hücresi


def adaylar(ag: Ag, o: Olay, kab_idx_kol: dict) -> Aday:
    """Konum adayları (+-2 sigma, Gauss öncül), başlangıç bölgesi ve kırılma yaması."""
    h = ag.h
    c = ag.cells
    idx = kab_idx_kol[o.kol]
    xs, zs = c.x_km.to_numpy()[idx], c.z_km.to_numpy()[idx]
    A_h = h * h
    if o.tip == "buyuk" and not np.isfinite(o.sigma):
        # tarihsel: başlangıç noktası bilinmiyor -> kırılma boyunca (uçlar dahil) eşit
        pos = np.arange(o.x0 - o.uc_km, o.x1 + o.uc_km + 1e-9, max(h, 1.0))
        p = np.ones(len(pos)) / len(pos)
        z0 = np.full(len(pos), 10.0)
    else:
        sig = max(o.sigma, 1e-6)
        if sig < h / 2:
            pos = np.array([o.x])
            p = np.array([1.0])
        else:
            k = np.arange(-math.floor(2 * sig / h), math.floor(2 * sig / h) + 1)
            pos = o.x + k * h
            p = np.exp(-0.5 * (k * h / sig) ** 2)
            p /= p.sum()
            # %1 kuralı 1 km aralık için yazıldı; ince ağda yoğunluk eşdeğeri. Ağın dışına düşen
            # adaylar atılır ama olasılıkları yeniden dağıtılmaz (olay büyük olasılıkla ağ dışında)
            keep = (p >= 0.01 * min(h, 1.0)) & (pos >= xs.min() - h / 2) & (pos <= xs.max() + h / 2)
            ic_ = p[p >= 0.01 * min(h, 1.0)].sum()
            pos, p = pos[keep], (p[keep] / ic_ if ic_ > 0 else p[keep])
        z0 = np.full(len(pos), o.z)
    s_nuc = float(min(L_kirilma(o.mw), 5.0))
    s_yama = float(L_kirilma(o.mw))
    bolge, yama, merkez = [], [], []
    for xp, zp in zip(pos, z0):
        j0 = int(np.argmin((xs - xp) ** 2 + (zs - zp) ** 2))
        merkez.append(idx[j0])

        def kare(s):
            sel = (np.abs(xs - xp) <= max(s, h) / 2 + 1e-9) & (np.abs(zs - zp) <= max(s, h) / 2 + 1e-9)
            ii = idx[sel] if sel.any() else np.array([idx[j0]])
            f = min(1.0, s * s / (len(ii) * A_h))
            return ii, f
        bolge.append(kare(s_nuc))
        if o.tip != "buyuk":
            yama.append(kare(s_yama))
    return Aday(np.asarray(p), bolge, yama, np.asarray(merkez))


# ------------------------------------------------------------------ özümseme koşusu
def doluluk(m: Model) -> np.ndarray:
    return np.where(m.vw, np.clip(m.S / np.maximum(m.S_fail, 1.0), -1.0, 1.2), np.nan)


def _yama_kaymasi(m: Model, ii: np.ndarray, f: float) -> np.ndarray:
    F = ii[m.vw[ii]]
    target = -m.p.asiri * m.dtau[F]
    F = F[m.S[F] > target]
    if len(F) == 0:
        return F, np.zeros(0)
    target = -m.p.asiri * m.dtau[F]
    Kff = np.asarray(m.K[np.ix_(F, F)], dtype=np.float64)
    x = np.linalg.solve(Kff, -(m.S[F] - target))
    return F, np.maximum(x, 0.0) * f


def uygula(m: Model, o: Olay, ad: Aday, w: np.ndarray, mu_eff: float, tam_yama=False):
    """Gerçek depremi modele uygula: başlangıç bölgesini sınıra çek, kaymayı uygula."""
    if o.tip == "dis":
        m.S += (o.dis["dtau"] + np.float32(mu_eff) * o.dis["dsn"]).astype(float) * (MU_YENI / MU_MATRIS)
        return None
    # (1) başlangıç bölgesi: S <- S + p f (S_kırılma - S)
    for wk, (ii, f) in zip(w, ad.bolge):
        if wk <= 0:
            continue
        v = ii[m.vw[ii]]
        m.S[v] += wk * f * np.maximum(m.S_fail[v] - m.S[v], 0.0)
    if o.tip == "buyuk":
        mask = (m.x >= o.x0) & (m.x <= o.x1) & (m.c["kol"].to_numpy() == 0)
        if tam_yama:  # yan yana koşulan sürüm: bütün parça %100
            d = mask & m.vw
            m.S[d] = np.maximum(m.S[d], m.S_fail[d])
        return m.kirilma(0, maske=mask, zorla=True, tarih=o.ad)
    # (2) küçük olay: adaylar üzerinden beklenen kayma (doğrusal üst üste binme)
    d = np.zeros(m.n)
    for wk, (ii, f) in zip(w, ad.yama):
        if wk <= 0:
            continue
        F, x = _yama_kaymasi(m, ii, f)
        d[F] += wk * x
    idx = np.flatnonzero(d > 0)
    if len(idx) == 0:
        return None
    m.S += m.K[:, idx] @ d[idx].astype(np.float32)
    m.slip[idx] += d[idx]
    m0 = m.p.mu_kayma * float(np.sum(d[idx] * m.alan[idx]))
    ev = {"t": m.t, "mw": (np.log10(m0) - 9.1) / 1.5, "m0": m0, "x0": float(m.x[idx].min()),
          "x1": float(m.x[idx].max()), "zorla": True, "tarih": o.ad}
    m.olaylar.append(ev)
    return ev


def _ortusme(e, o: Olay):
    """Kendiliğinden kırılmanın gerçek kırılma boyuyla ortak oranı (uç belirsizliği tolere)."""
    a0, a1 = o.x0 - o.uc_km, o.x1 + o.uc_km
    L = max(o.x1 - o.x0, 1.0)
    ort = max(0.0, min(e["x1"], a1) - max(e["x0"], a0))
    return min(ort / L, 1.0)


class Kosu:
    """Bir parametre takımı ve bir topluluk üyesiyle baştan sona özümseme koşusu."""

    def __init__(self, ag: Ag, th: dict, olaylar: list[Olay], seed: int, isinma: float,
                 t_bas: float, t_son: float, egitim_bitis: float, sonsal=True, tam_yama=False,
                 rastgele_n=200, karsi_olgu=True, kayit_yolu: Path | None = None):
        self.ag, self.th, self.ol = ag, th, olaylar
        self.seed, self.isinma, self.t_bas, self.t_son = seed, isinma, t_bas, t_son
        self.egitim_bitis, self.sonsal, self.tam_yama = egitim_bitis, sonsal, tam_yama
        self.rastgele_n, self.karsi_olgu = rastgele_n, karsi_olgu
        self.kayit_yolu = kayit_yolu
        c = ag.cells
        kab = (c.katman == "kabuk").to_numpy()
        self.kab_idx_kol = {k: np.flatnonzero(kab & (c.kol == k).to_numpy()) for k in sorted(c.kol.unique())}
        self.ad = [None if o.tip == "dis" else adaylar(ag, o, self.kab_idx_kol) for o in olaylar]
        # aynı yerdeki bir önceki gerçek olay (küçük olaylar için zaman şartı)
        self.onceki = self._onceki_olaylar()

    def _onceki_olaylar(self):
        out = [None] * len(self.ol)
        son = {}
        for i, o in enumerate(self.ol):
            if o.tip == "dis":
                continue
            if o.tip == "buyuk":
                # bir önceki büyük kırılma aynı parçada
                for j in range(i - 1, -1, -1):
                    q = self.ol[j]
                    if q.tip == "buyuk" and min(q.x1, o.x1) - max(q.x0, o.x0) > 0:
                        out[i] = j
                        break
                continue
            r = max(o.sigma, L_kirilma(o.mw) / 2, self.ag.h)
            for j in range(i - 1, max(i - 5000, -1), -1):
                q = self.ol[j]
                if q.tip != "kucuk" or q.kol != o.kol:
                    continue
                rr = max(r, q.sigma, L_kirilma(q.mw) / 2)
                if abs(q.x - o.x) <= rr and abs(q.z - o.z) <= max(2.0, rr):
                    out[i] = j
                    break
        return out

    # ---------------------------------------------------------------- ana döngü
    def calistir(self, durum=None):
        th = self.th
        if durum is None:
            m = model_kur(self.ag, th, self.seed)
            m.t = self.t_bas - self.isinma
            m.ilerle(self.t_bas)
            m.olaylar = []
            i0, kayit = 0, []
        else:
            m, i0, kayit = durum
        self.m = m
        rng = np.random.default_rng(self.seed + 77)
        son_yaz = time.time()
        for i in range(i0, len(self.ol)):
            o = self.ol[i]
            n0 = len(m.olaylar)
            m.ilerle(o.t)
            ad = self.ad[i]
            if o.tip != "dis" and len(ad.p) == 0:
                continue  # bütün konum adayları ağ dışında
            if o.tip == "dis":
                uygula(m, o, ad, None, th["mu_eff"])
                continue
            r = self._puanla(m, i, o, ad, rng) if o.puan else None
            # ağırlık: ayar döneminde öncül x doluluk olabilirliği; sınavda yalnız öncül
            w = ad.p
            if self.sonsal and o.t < self.egitim_bitis and len(w) > 1:
                fz = np.array([np.nanmean(doluluk(m)[ii]) if m.vw[ii].any() else 0.0 for ii, _ in ad.bolge])
                lik = np.exp(-np.clip(1.0 - np.nan_to_num(fz, nan=0.0), 0, None) / 0.1)
                w = ad.p * lik
                w = w / w.sum() if w.sum() > 0 else ad.p
            ev = uygula(m, o, ad, w, th["mu_eff"], self.tam_yama)
            # kıyas modeli için: her hücrede son gerçek depremin zamanı (yalnız katalog)
            if not hasattr(self, "_son_gercek"):
                self._son_gercek = np.full(m.n, -np.inf)
            if o.tip == "buyuk":
                self._son_gercek[(m.x >= o.x0) & (m.x <= o.x1) & (m.c["kol"].to_numpy() == 0)] = o.t
            else:
                self._son_gercek[ad.merkez[w >= 0.05 * w.max()]] = o.t
            if r is not None:
                r["mw_model"] = ev["mw"] if ev else np.nan
                kayit.append(r)
            if self.kayit_yolu is not None and time.time() - son_yaz > 600:
                self.ara_kayit(m, i + 1, kayit)
                son_yaz = time.time()
        m.ilerle(self.t_son)
        self.kayit = kayit
        return m, kayit

    def ara_kayit(self, m, i, kayit):
        tmp = self.kayit_yolu.with_suffix(".tmp")
        tmp.write_bytes(pickle.dumps({"durum": m.durum(), "olaylar": m.olaylar, "i": i, "kayit": kayit,
                                      "th": self.th, "seed": self.seed}))
        tmp.replace(self.kayit_yolu)

    # ---------------------------------------------------------------- puanlama
    def _puanla(self, m: Model, i: int, o: Olay, ad: Aday, rng) -> dict:
        dol = doluluk(m)
        esik, tol = ESIK[o.sinif], TOL[o.sinif]
        j = self.onceki[i]
        t_prev = self.ol[j].t if j is not None else -np.inf
        r = {"i": i, "t": o.t, "ad": o.ad, "tip": o.tip, "mw": o.mw, "sigma": o.sigma, "puan": o.puan,
             "sinif": o.sinif, "kol": o.kol, "x": o.x, "z": o.z, "t_prev": t_prev}
        if o.tip == "buyuk":
            mask = (m.x >= o.x0) & (m.x <= o.x1) & (m.c["kol"].to_numpy() == 0) & m.vw
            r["doluluk"] = float(np.average(np.clip(dol[mask], -1, 1.2), weights=m.alan[mask]))
            # kendiliğinden büyük kırılma: önceki büyük olaydan bu yana, ortak boy >= %50
            erken, zamaninda = None, None
            for e in m.olaylar:
                if e.get("zorla") or e["mw"] < 6.8 or e["t"] <= t_prev:
                    continue
                if _ortusme(e, o) >= 0.5:
                    if e["t"] < o.t - tol:
                        erken = e
                    elif zamaninda is None:
                        zamaninda = e
            dt = None
            if zamaninda is not None:
                dt, ort = zamaninda["t"] - o.t, _ortusme(zamaninda, o)
            elif self.karsi_olgu:
                # karşı olgu: zorlamadan ileri sar (en çok 3 tolerans), ilk uyan kırılma
                st, n0 = m.durum(), len(m.olaylar)
                t_bas = m.t
                ort = 0.0
                while m.t < t_bas + 3 * tol - 1e-9 and dt is None:
                    k0 = len(m.olaylar)
                    m.ilerle(min(m.t + 1.0, t_bas + 3 * tol))
                    for e in m.olaylar[k0:]:
                        if e["mw"] >= 6.8 and _ortusme(e, o) >= 0.5:
                            dt, ort = e["t"] - o.t, _ortusme(e, o)
                            break
                m.yukle(st)
                del m.olaylar[n0:]
            r["erken_kirildi"] = erken is not None
            r["dt"] = np.nan if dt is None else float(dt)
            r["ortusme"] = float(ort) if dt is not None else 0.0
            r["yer_ok"] = dt is not None and r["ortusme"] >= 0.5
            r["zaman_ok"] = (erken is None) and dt is not None and abs(dt) <= tol
            r["dolu_ok"] = r["doluluk"] >= esik
            r["dolu_ok_90"] = r["doluluk"] >= 0.90
            r["tuttu"] = bool(r["yer_ok"] and r["zaman_ok"] and r["dolu_ok"])
            r["tuttu_90"] = bool(r["yer_ok"] and r["zaman_ok"] and r["dolu_ok_90"])
            r["yumusak"] = 0.5 * min(max(r["doluluk"], 0) / esik, 1.0) + \
                0.5 * (0.0 if dt is None or erken is not None else math.exp(-0.5 * (dt / tol) ** 2))
            r["u_model"], r["u_sure"] = self._sira(m, dol, ad, o, rng)
            return r
        # küçük olay
        fz = np.array([dol[k] if m.vw[k] else np.nan for k in ad.merkez])
        p = ad.p
        okf = np.isfinite(fz)
        r["kilitli_oran"] = float(p[okf].sum())
        r["doluluk"] = float(np.nansum(p * np.nan_to_num(fz, nan=0)) / max(p[okf].sum(), 1e-9)) if okf.any() else np.nan
        # zaman: önceki olaydan sonra, toleransın başından önce kendiliğinden kırıldı mı
        if np.isfinite(t_prev):
            dT = o.t - t_prev
            tol = 0.25 * dT
            sk = m.son_kendi[ad.merkez]
            erken = (sk > t_prev) & (sk < o.t - tol)
            zam = (sk >= o.t - tol)
            rate, _ = m._oran()
            gap = np.maximum(m.S_fail[ad.merkez] - m.S[ad.merkez], 0)
            kalan = np.where(rate[ad.merkez] > 0, gap / np.maximum(rate[ad.merkez], 1e-30), np.inf)
            hit = (~erken) & (zam | ((np.nan_to_num(fz, nan=-1) >= esik) & (kalan <= tol)))
            r["p_tuttu"] = float(np.sum(p * hit))
            r["p_tuttu_90"] = float(np.sum(p * ((~erken) & (zam | ((np.nan_to_num(fz, nan=-1) >= 0.9) & (kalan <= tol))))))
            r["p_erken"] = float(np.sum(p * erken))
            r["dT"] = dT
        else:
            r["p_tuttu"] = r["p_tuttu_90"] = r["p_erken"] = np.nan
            r["dT"] = np.nan
        r["u_model"], r["u_sure"] = self._sira(m, dol, ad, o, rng)
        return r

    def _sira(self, m, dol, ad: Aday, o: Olay, rng):
        """Olay yerinin doluluk sırası: aynı kol ve derinlik bandındaki rastgele kilitli hücrelere göre
        (0,5 = şans). Kıyas: son gerçek depremden beri geçen süre (yalnız katalogdan)."""
        idx = self.kab_idx_kol[o.kol if o.kol >= 0 else 0]
        z = m.z[idx]
        zc = float(np.average(m.z[ad.merkez], weights=ad.p))
        band = idx[(np.abs(z - zc) <= 2.0) & m.vw[idx]]
        if len(band) < 10:
            return np.nan, np.nan
        rnd = rng.choice(band, size=min(self.rastgele_n, len(band)), replace=False)
        fv = np.nan_to_num(np.array([dol[k] for k in ad.merkez]), nan=-1.0)
        f_ev = float(np.sum(ad.p * fv))
        u_m = float(np.mean(dol[rnd] < f_ev) + 0.5 * np.mean(dol[rnd] == f_ev))
        # basit kestirme: son gerçek depremden beri geçen süre (her hücre için katalog)
        if not hasattr(self, "_son_gercek"):
            self._son_gercek = np.full(m.n, -np.inf)
        sg = self._son_gercek
        t_ev = float(np.sum(ad.p * np.where(np.isfinite(sg[ad.merkez]), m.t - sg[ad.merkez], 1e4)))
        tr = np.where(np.isfinite(sg[rnd]), m.t - sg[rnd], 1e4)
        u_s = float(np.mean(tr < t_ev) + 0.5 * np.mean(tr == t_ev))
        return u_m, u_s



def ozet_puan(kayit: list[dict], sahte: int, tur: str) -> dict:
    """Ayar puanı. Büyük: önce ardışık tutan sayısı (geri sarma), sonra toplam, sonra yumuşak puan,
    sahte büyük depremler ıska sayılır. Küçük: 1/σ² ağırlıklı tutma olasılığı ve sıra (AUC)."""
    k = pd.DataFrame(kayit)
    if len(k) == 0:
        return {"J": -1e9}
    if tur == "buyuk":
        a = k[(k.puan == "ayar") & (k.tip == "buyuk")].sort_values("t")
        ard = 0
        for t in a.tuttu:
            if not t:
                break
            ard += 1
        J = 10 * ard + 3 * a.tuttu.sum() + 3 * a.yumusak.mean() - 2 * sahte
        return {"J": float(J), "ardisik": int(ard), "tuttu": int(a.tuttu.sum()), "n": int(len(a)),
                "yumusak": float(a.yumusak.mean()), "sahte": int(sahte),
                "doluluk": a.doluluk.round(3).tolist(), "dt": a.dt.round(1).tolist()}
    a = k[(k.puan == "ayar") & (k.tip == "kucuk")]
    w = 1.0 / np.maximum(a.sigma.to_numpy(), 0.1) ** 2
    ok = np.isfinite(a.p_tuttu.to_numpy())
    pt = float(np.sum(w[ok] * a.p_tuttu.to_numpy()[ok]) / max(w[ok].sum(), 1e-9)) if ok.any() else 0.0
    um = a.u_model.to_numpy()
    okm = np.isfinite(um)
    auc = float(np.sum(w[okm] * um[okm]) / w[okm].sum()) if okm.any() else 0.5
    J = 100 * pt + 10 * (auc - 0.5) - 1.0 * sahte
    return {"J": float(J), "p_tuttu": pt, "auc": auc, "n": int(len(a)), "n_zaman": int(ok.sum()),
            "sahte": int(sahte)}


def sahte_say(m: Model, olaylar: list[Olay], t0: float, t1: float, mw=6.8) -> list[dict]:
    """Gerçek kayıtta karşılığı olmayan kendiliğinden büyük depremler."""
    out = []
    for e in m.olaylar:
        if e.get("zorla") or e["mw"] < mw or not (t0 <= e["t"] <= t1):
            continue
        es = False
        for o in olaylar:
            if o.tip == "buyuk" and abs(o.t - e["t"]) <= TOL[o.sinif] and _ortusme(e, o) > 0:
                es = True
                break
        if not es:
            out.append(e)
    return out
