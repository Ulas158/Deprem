"""İnsan kaynaklı ve doğal akışkan sahaları (baraj, maden, jeotermal, petrol/gaz).

Veri ``external_data/insan_kaynakli_sahalar.csv`` dosyasındadır. Her satırın
başlangıç/bitiş yılı ve durumu (aktif/pasif/kapalı/bilinmiyor/dogal) vardır,
böylece bir saha yalnızca o tarihte faaliyetteyse özellik üretir (gelecekten
sızıntı olmaz). Konumların çoğu yaklaşık ve **doğrulanmamıştır**
(``verified`` sütunu); ``loc_accuracy_km`` belirsizliği verir.

Mekanizmalar (Foulger vd. 2018 HiQuake derlemesi; Gupta 2002; McGarr 2002;
Ellsworth 2013; Keranen & Weingarten 2018):
* Baraj: su yükü anında elastik gerilme yaratır; gözenek basıncı ise aylar–
  yıllar içinde difüzyonla derine iner. Tetiklemenin çoğu ilk dolumdan
  sonraki ilk yıllarda veya hızlı seviye değişimlerinde görülür (Koyna,
  Hindistan; Kremasta, Yunanistan).
* Madencilik: kaya çıkarımı gerilmeyi yeniden dağıtır; kaya patlamaları ve
  ocak çevresinde küçük depremler (genellikle M<5).
* Akışkan enjeksiyonu/çekimi (jeotermal, petrol, gaz, atık su, gaz depolama):
  gözenek basıncı değişimi; en güçlü kanıtlanmış insan kaynaklı etki
  (Oklahoma, Pohang 2017 M5.5, Groningen).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..geo import haversine_km

ROOT = Path(__file__).resolve().parents[3]
SITES_CSV = ROOT / "external_data" / "insan_kaynakli_sahalar.csv"

CATEGORY_GROUPS = {
    "reservoir": ["baraj"],
    "mining": ["maden_komur", "maden_metal", "maden_endustriyel", "maden_cozelti"],
    "geothermal": ["jeotermal"],
    "hydrocarbon": ["petrol", "dogalgaz", "gaz_depolama"],
    "injection": ["maden_cozelti", "gaz_depolama", "jeotermal", "petrol"],
    "gas_seep": ["dogal_gaz_cikisi"],
}


def load_sites(path: Path | None = None) -> pd.DataFrame:
    df = pd.read_csv(path or SITES_CSV)
    for c in ("start_year", "end_year"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def active_at(sites: pd.DataFrame, year: float) -> pd.Series:
    """O yıl faaliyette mi? Başlangıç yılı bilinmeyen 'aktif' sahalar her zaman aktif sayılır."""
    started = sites["start_year"].isna() | (sites["start_year"] <= year)
    not_ended = sites["end_year"].isna() | (sites["end_year"] >= year)
    known = sites["status"] != "planlanan"
    return started & not_ended & known


def site_features(lat, lon, sites: pd.DataFrame, t0: pd.Timestamp) -> pd.DataFrame:
    """Noktalar (hücre/eleman merkezleri) için en yakın faal saha uzaklıkları.

    Rezervuarlar için ayrıca en yakın barajın (50 km içinde) dolumundan bu
    yana geçen yıl verilir: rezervuar tetiklemesi ilk yıllarda yoğunlaşır.
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    yr = t0.year + (t0.dayofyear - 1) / 365.25
    act = sites[active_at(sites, yr)]
    out = {}
    for key, cats in CATEGORY_GROUPS.items():
        s = act[act["category"].isin(cats)]
        if s.empty:
            out[f"site_{key}_dist_km"] = np.full(len(lat), np.nan)
            continue
        d = haversine_km(lat[:, None], lon[:, None], s["lat"].to_numpy()[None, :],
                         s["lon"].to_numpy()[None, :])
        out[f"site_{key}_dist_km"] = d.min(axis=1)
        if key == "reservoir":
            k = d.argmin(axis=1)
            since = yr - s["start_year"].to_numpy()[k]
            out["reservoir_years_since_fill"] = np.where(d.min(axis=1) <= 50, since, np.nan)
    return pd.DataFrame(out)
