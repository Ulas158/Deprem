"""Fay elemanı parametre kaydı.

Her fay elemanı (bkz. ``mesh.py``) bir parametre satırı taşır. Bu modül
hangi parametrelerin var olduğunu, birimlerini, kaynaklarını ve **durumunu**
tek bir yerde tanımlar:

* ``geometry``   ağdan hesaplanan geometri (her zaman dolu)
* ``fault_db``   fay veritabanından (GEM/EMME) gelen değerler
* ``seismicity`` deprem kataloğundan, bir referans zamandan **önceki** veriyle
                 hesaplanan değerler (``element_seismicity``)
* ``placeholder`` literatürde anlamlı ama henüz veri kaynağı bağlanmamış
                 alanlar; NaN olarak tanımlıdır, veri bağlandıkça dolar

Yeni parametre eklemek için ``register(...)`` çağırmak yeterli; ağ tablosu
``ensure_columns`` ile tüm kayıtlı alanları içerir.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ..stats.gr import LOG10E


@dataclass(frozen=True)
class ParamSpec:
    name: str
    unit: str
    group: str
    status: str  # "geometry" | "fault_db" | "seismicity" | "placeholder"
    description: str
    source: str = ""


REGISTRY: dict[str, ParamSpec] = {}


def register(name, unit, group, status, description, source=""):
    REGISTRY[name] = ParamSpec(name, unit, group, status, description, source)


# ------------------------------------------------------------------ geometri
for n, u, d in [
    ("lon", "°", "Eleman merkezi boylamı"), ("lat", "°", "Eleman merkezi enlemi"),
    ("depth_km", "km", "Eleman merkezi derinliği"),
    ("strike", "°", "Yerel doğrultu (kuzeyden saat yönü)"), ("dip", "°", "Eğim"),
    ("area_km2", "km²", "Eleman alanı"),
    ("along_km", "km", "Segment başından doğrultu boyunca uzaklık"),
    ("along_frac", "-", "Segment boyunca göreli konum (0-1)"),
    ("end_dist_km", "km", "En yakın segment ucuna uzaklık (kırılma durma/başlama bölgeleri)"),
    ("down_dip_km", "km", "Fay düzleminde üst kenardan aşağı uzaklık"),
    ("depth_rel", "-", "Sismojenik zon içindeki göreli derinlik (0 üst, 1 alt)"),
    ("strike_bend_deg", "°", "Komşu elemana göre doğrultu değişimi (bükülme/büklüm)"),
    ("fault_len_km", "km", "Ait olduğu segmentin uzunluğu"),
]:
    register(n, u, "geometri", "geometry", d, "mesh.py")

# ------------------------------------------------------------------ fay veritabanı
for n, u, d in [
    ("rake", "°", "Kayma açısı (mekanizma türü)"),
    ("slip_rate_mm_yr", "mm/yıl", "Uzun dönem kayma hızı (en iyi tahmin)"),
    ("slip_rate_min", "mm/yıl", "Kayma hızı alt sınırı"),
    ("slip_rate_max", "mm/yıl", "Kayma hızı üst sınırı"),
    ("upper_depth_km", "km", "Sismojenik zon üst sınırı"),
    ("lower_depth_km", "km", "Sismojenik zon alt sınırı (kilitli derinlik)"),
    ("m0_rate_Nm_yr", "N·m/yıl", "Moment birikim hızı = μ·A·v"),
]:
    register(n, u, "fay", "fault_db", d, "GEM GAF / EMME")

# ------------------------------------------------------------------ sismisite
SEISMICITY_PARAMS = [
    ("seis_n_all", "-", "Elemana R km içindeki tüm olaylar (referans zamandan önce)"),
    ("seis_n_10y", "-", "Son 10 yıldaki olay sayısı (M≥Mc)"),
    ("seis_n_1y", "-", "Son 1 yıldaki olay sayısı (M≥Mc)"),
    ("seis_rate_z", "-", "Hız değişimi z-değeri: son 1 yıl vs önceki 9 yıl (+ hızlanma, − sessizlik)"),
    ("seis_b_10y", "-", "Son 10 yıl b-değeri (Aki-Utsu; n<30 ise boş)"),
    ("seis_max_mag", "-", "Görülen en büyük büyüklük"),
    ("seis_yrs_since_m4", "yıl", "Son M≥4'ten bu yana geçen süre"),
    ("seis_yrs_since_m5", "yıl", "Son M≥5'ten bu yana geçen süre"),
    ("seis_log_m0_10y", "log10 N·m", "Son 10 yılda boşalan sismik moment"),
    ("seis_mean_depth", "km", "Yakın olayların ortalama derinliği"),
    ("seis_depth_std", "km", "Yakın olayların derinlik yayılımı"),
    ("seis_frac_dependent", "-", "Artçı/öncü (bağımlı) olay oranı: kümelenme/sürü göstergesi"),
    ("seis_m0_release_ratio_10y", "-", "Son 10 yılda boşalan / biriken moment"),
]
for n, u, d in SEISMICITY_PARAMS:
    register(n, u, "sismisite", "seismicity", d, "birleşik katalog")

# ------------------------------------------------------------------ henüz veri yok
for n, u, d, s in [
    ("gnss_strain_rate", "1/yıl", "GNSS'ten ikinci değişmez gerinim hızı", "Kreemer vd. 2014 GSRM; Aktuğ vd."),
    ("gnss_slip_deficit_mm_yr", "mm/yıl", "Blok modelinden kilitlenme × kayma açığı", "Reilinger vd. 2006; Ergintav vd."),
    ("locking_ratio", "-", "Jeodezik kilitlenme oranı (0 krip, 1 tam kilitli)", "GNSS/InSAR ters çözüm"),
    ("creep_rate_mm_yr", "mm/yıl", "Gözlenen aseismik krip hızı", "İsmetpaşa, Destek vb. gözlemleri"),
    ("insar_los_rate_mm_yr", "mm/yıl", "InSAR görüş doğrultusu hız alanı", "Sentinel-1 / COMET-LiCSAR"),
    ("heat_flow_mw_m2", "mW/m²", "Yüzey ısı akısı", "Global Heat Flow Database"),
    ("curie_depth_km", "km", "Curie nokta derinliği (termal yapı)", "Manyetik anomali ters çözümleri"),
    ("crustal_thickness_km", "km", "Kabuk kalınlığı (Moho)", "CRUST1.0 / alıcı fonksiyonları"),
    ("vp_vs_ratio", "-", "Vp/Vs oranı (akışkan/doygunluk göstergesi)", "Yerel tomografi"),
    ("lithology_class", "kategori", "Fay boyunca kayaç türü", "MTA 1:500k jeoloji haritası"),
    ("geothermal_proximity_km", "km", "En yakın jeotermal alana uzaklık", "MTA jeotermal envanteri"),
    ("friction_a_minus_b", "-", "Hız-durum sürtünme (a−b) tahmini (litolojiye dayalı)", "Laboratuvar derlemeleri"),
    ("tidal_stress_sensitivity", "-", "Gelgit gerilmesi fazı ile olay zamanlaması korelasyonu", "Katalogdan türetilecek"),
    ("coulomb_dcff_bar", "bar", "Geçmiş depremlerden birikimli Coulomb gerilme değişimi", "Sonraki PR (Okada)"),
    ("coulomb_dcff_last_bar", "bar", "Son büyük depremden kaynaklı ΔCFF", "Sonraki PR (Okada)"),
    ("last_rupture_year", "yıl", "Elemanın bilinen son yüzey kırığı yılı (tarihsel dahil)", "Tarihsel katalog thread'i"),
    ("paleo_recurrence_yr", "yıl", "Paleosismolojik yinelenme aralığı", "MTA hendek çalışmaları"),
    ("fluid_pressure_index", "-", "Akışkan basıncı göstergesi (b, Vp/Vs, jeotermalden birleşik)", "Türetilecek"),
]:
    register(n, u, "henüz veri yok", "placeholder", d, s)


def ensure_columns(elements: pd.DataFrame) -> pd.DataFrame:
    """Kayıttaki tüm parametrelerin sütun olarak var olmasını sağla (eksik = NaN)."""
    out = elements.copy()
    for name in REGISTRY:
        if name not in out.columns:
            out[name] = np.nan
    return out


def registry_table() -> pd.DataFrame:
    return pd.DataFrame([vars(p) for p in REGISTRY.values()])


# ------------------------------------------------------------------ sismisite hesapları
def _local_xyz(lat, lon, depth, lat0=39.0, lon0=35.0):
    x = (np.asarray(lon) - lon0) * 111.195 * np.cos(np.radians(lat0))
    y = (np.asarray(lat) - lat0) * 111.195
    return np.column_stack([x, y, np.asarray(depth, dtype=float)])


def element_seismicity(elements: pd.DataFrame, events: pd.DataFrame, t_ref: pd.Timestamp,
                       mc: float, radius_km: float = 5.0, dm: float = 0.1,
                       independent: pd.Series | None = None) -> pd.DataFrame:
    """Her eleman için ``t_ref`` öncesi olaylardan sismisite parametreleri.

    Olay-eleman eşleşmesi 3B hiposantral uzaklıkla (``radius_km``) yapılır.
    Konum hataları km düzeyinde olduğundan bir olay birden fazla komşu
    elemana sayılabilir; bu kasıtlıdır (yumuşatma etkisi).
    ``independent``: artçı ayıklama bayrakları (True = bağımsız); verilirse
    ``seis_frac_dependent`` hesaplanır.
    """
    ev = events[events["time"] < t_ref]
    if independent is not None:
        independent = independent.loc[ev.index]
    n_el = len(elements)
    out = pd.DataFrame(index=elements.index)
    if ev.empty:
        for n, _, _ in SEISMICITY_PARAMS:
            out[n] = np.nan
        return out
    el_xyz = _local_xyz(elements["lat"], elements["lon"], elements["depth_km"])
    ev_xyz = _local_xyz(ev["latitude"], ev["longitude"], ev["depth_km"].fillna(10.0))
    pairs = cKDTree(el_xyz).sparse_distance_matrix(cKDTree(ev_xyz), radius_km,
                                                   output_type="coo_matrix")
    r, c = pairs.row, pairs.col
    yrs_before = ((t_ref - ev["time"]).dt.total_seconds() / (365.25 * 86400)).to_numpy()[c]
    mag = ev["mag"].to_numpy()[c]
    mw = (ev["mw_best"] if "mw_best" in ev else ev["mag"]).to_numpy()[c]
    depth = ev["depth_km"].fillna(10.0).to_numpy()[c]
    above = mag >= mc - dm / 2

    def count(mask):
        return np.bincount(r[mask], minlength=n_el).astype(float)

    out["seis_n_all"] = count(np.ones_like(r, dtype=bool))
    n10 = count(above & (yrs_before <= 10))
    n1 = count(above & (yrs_before <= 1))
    n_prev = n10 - n1
    out["seis_n_10y"], out["seis_n_1y"] = n10, n1
    r1, r2 = n1 / 1.0, n_prev / 9.0
    with np.errstate(divide="ignore", invalid="ignore"):
        out["seis_rate_z"] = np.where(n10 > 0, (r1 - r2) / np.sqrt(n1 / 1.0 + n_prev / 81.0 + 1e-9),
                                      np.nan)
    # b-değeri: Aki-Utsu, vektörel
    m10 = above & (yrs_before <= 10)
    mb = np.round(mag / dm) * dm
    s = np.bincount(r[m10], weights=mb[m10], minlength=n_el)
    with np.errstate(divide="ignore", invalid="ignore"):
        b = LOG10E / (s / n10 - (mc - dm / 2))
    out["seis_b_10y"] = np.where(n10 >= 30, b, np.nan)
    mx = np.full(n_el, np.nan)
    np.fmax.at(mx, r, mag)
    out["seis_max_mag"] = mx
    for m in (4, 5):
        since = np.full(n_el, np.inf)
        sel = mag >= m
        np.minimum.at(since, r[sel], yrs_before[sel])
        out[f"seis_yrs_since_m{m}"] = np.where(np.isinf(since), np.nan, since)
    m0 = 10 ** (1.5 * mw + 9.1)
    m0_10 = np.bincount(r[yrs_before <= 10], weights=m0[yrs_before <= 10], minlength=n_el)
    with np.errstate(divide="ignore"):
        out["seis_log_m0_10y"] = np.where(m0_10 > 0, np.log10(m0_10), np.nan)
    n_all = out["seis_n_all"].to_numpy()
    sd = np.bincount(r, weights=depth, minlength=n_el)
    sd2 = np.bincount(r, weights=depth ** 2, minlength=n_el)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_d = sd / n_all
        out["seis_mean_depth"] = np.where(n_all > 0, mean_d, np.nan)
        out["seis_depth_std"] = np.where(n_all > 1, np.sqrt(np.clip(sd2 / n_all - mean_d ** 2, 0, None)),
                                         np.nan)
    if independent is not None:
        dep = ~independent.to_numpy()[c]
        nd = np.bincount(r, weights=dep.astype(float), minlength=n_el)
        with np.errstate(divide="ignore", invalid="ignore"):
            out["seis_frac_dependent"] = np.where(n_all > 0, nd / n_all, np.nan)
    else:
        out["seis_frac_dependent"] = np.nan
    # Not: bir olay R içindeki tüm elemanlara sayıldığından moment oranı
    # eleman düzeyinde abartılı olabilir; göreli karşılaştırma için kullanılmalı.
    acc = elements["m0_rate_Nm_yr"].to_numpy() * 10
    with np.errstate(divide="ignore", invalid="ignore"):
        out["seis_m0_release_ratio_10y"] = np.where(acc > 0, m0_10 / acc, np.nan)
    return out
