"""Ortak katalog şeması.

Tüm kaynaklardan gelen olaylar bu sütunlara dönüştürülür; analiz ve tahmin
modülleri yalnızca bu şemaya dayanır. Böylece yeni bir kaynak (ör. AFAD,
KOERI arşivi) eklemek analiz kodunu değiştirmez.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

COLUMNS = [
    "time",  # UTC, tz-aware pandas Timestamp
    "latitude",
    "longitude",
    "depth_km",
    "mag",
    "mag_type",  # kaynaktaki büyüklük türü (ml, mw, mb, md ...), küçük harf
    "source",  # verinin çekildiği servis: usgs, emsc, afad
    "agency",  # çözümü üreten kurum (KOERI, AFAD, us ...)
    "event_id",  # kaynağa özgü kimlik
]

# İsteğe bağlı belirsizlik alanları (1σ). Kaynak vermezse NaN kalır ve
# ``deprem.catalog.weights.fill_uncertainties`` kaynak/dönem varsayılanlarıyla
# doldurur. Tarihsel (aletsel öncesi) kataloglar da aynı alanları kullanır.
UNCERTAINTY_COLUMNS = [
    "time_err_s",  # oluş zamanı belirsizliği, saniye (tarihsel: gün/yıl ölçeğinde)
    "horiz_err_km",  # episantr belirsizliği, km
    "depth_err_km",  # derinlik belirsizliği, km
    "mag_err",  # büyüklük belirsizliği, büyüklük birimi
]


def empty_catalog() -> pd.DataFrame:
    df = pd.DataFrame({c: pd.Series(dtype="object") for c in COLUMNS})
    df["time"] = pd.to_datetime(df["time"], utc=True)
    for c in ("latitude", "longitude", "depth_km", "mag"):
        df[c] = df[c].astype(float)
    return df


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Sütun türlerini zorla, eksik sütun ekle, zamana göre sırala."""
    out = df.copy()
    for c in COLUMNS + UNCERTAINTY_COLUMNS:
        if c not in out.columns:
            out[c] = np.nan
    extra = [c for c in out.columns if c not in COLUMNS + UNCERTAINTY_COLUMNS]
    out = out[COLUMNS + UNCERTAINTY_COLUMNS + extra]
    for c in UNCERTAINTY_COLUMNS:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype(float)
    out["time"] = pd.to_datetime(out["time"], utc=True, format="ISO8601")
    for c in ("latitude", "longitude", "depth_km", "mag"):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype(float)
    for c in ("mag_type", "source", "agency", "event_id"):
        out[c] = out[c].astype("string")
    out["mag_type"] = out["mag_type"].str.lower()
    return out.sort_values("time").reset_index(drop=True)


def validate(df: pd.DataFrame) -> None:
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Eksik sütunlar: {missing}")
    if not df["time"].is_monotonic_increasing:
        raise ValueError("Katalog zamana göre sıralı değil")
