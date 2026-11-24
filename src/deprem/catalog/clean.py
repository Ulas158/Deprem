"""Katalog temizleme, kaynaklar arası tekilleştirme ve artçı ayıklama."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..geo import haversine_km
from ..region import Region
from .schema import normalize


def epoch_seconds(times: pd.Series) -> np.ndarray:
    """UTC zamanları Unix saniyesine çevir (pandas zaman çözünürlüğünden bağımsız)."""
    return (times - pd.Timestamp("1970-01-01", tz="UTC")).dt.total_seconds().to_numpy()


def basic_clean(df: pd.DataFrame, region: Region | None = None) -> pd.DataFrame:
    """Eksik değerleri, bölge dışını ve aşırı derinliği at."""
    out = df.dropna(subset=["time", "latitude", "longitude", "mag"]).copy()
    out["depth_km"] = out["depth_km"].fillna(10.0)
    if region is not None:
        out = out[region.contains(out) & (out["depth_km"] <= region.max_depth_km)]
    return normalize(out)


def merge_sources(frames: list[pd.DataFrame], priority=("emsc", "usgs", "afad"),
                  dt_s: float = 16.0, dist_km: float = 50.0) -> pd.DataFrame:
    """Birden fazla katalogdaki aynı olayları eşleyip tek kayıt bırak.

    Eşleşme: zaman farkı ``dt_s`` saniyeden ve episantr farkı ``dist_km``
    km'den küçük. Tutulan kayıt ``priority`` sırasına göre seçilir.
    ``mw`` sütunu, eşleşen kayıtlar arasında moment büyüklüğü (mw*) olan
    varsa onun değerini taşır; moment bütçesi hesaplarında bu tercih edilir.
    ``n_sources`` olayın kaç kaynakta görüldüğünü verir.
    """
    all_ = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    all_ = normalize(all_)
    rank = {s: i for i, s in enumerate(priority)}
    all_["_rank"] = all_["source"].map(rank).fillna(len(rank)).astype(int)
    t = epoch_seconds(all_["time"])
    lat = all_["latitude"].to_numpy()
    lon = all_["longitude"].to_numpy()
    mtype = all_["mag_type"].fillna("").to_numpy()
    mag = all_["mag"].to_numpy()
    src = all_["source"].to_numpy()

    n = len(all_)
    group = -np.ones(n, dtype=int)
    g = 0
    order = np.argsort(t, kind="stable")
    for idx_pos, i in enumerate(order):
        if group[i] >= 0:
            continue
        group[i] = g
        j_pos = idx_pos + 1
        while j_pos < n and t[order[j_pos]] - t[i] <= dt_s:
            j = order[j_pos]
            if group[j] < 0 and src[j] != src[i] and \
                    haversine_km(lat[i], lon[i], lat[j], lon[j]) <= dist_km:
                group[j] = g
            j_pos += 1
        g += 1
    all_["_group"] = group

    is_mw = pd.Series(mtype).str.startswith("mw").to_numpy()
    mw = pd.Series(np.where(is_mw, mag, np.nan)).groupby(group).max()
    n_src = all_.groupby("_group")["source"].nunique()

    best = all_.sort_values(["_group", "_rank"]).drop_duplicates("_group")
    best = best.set_index("_group")
    best["mw"] = mw
    best["n_sources"] = n_src
    best = best.reset_index(drop=True).drop(columns=["_rank"])
    best = best.sort_values("time").reset_index(drop=True)
    return best


def moment_magnitude(df: pd.DataFrame) -> pd.Series:
    """Moment hesapları için en iyi Mw tahmini.

    Kaynakta Mw varsa o kullanılır; yoksa verilen büyüklük Mw kabul edilir.
    Bu kaba bir yaklaşımdır: Türkiye için ML→Mw dönüşüm ilişkileri
    (ör. Kadirioğlu & Kartal 2016) ileride eklenmeli. M≥5 olaylarda
    genellikle USGS/GCMT Mw mevcut olduğundan moment bütçesine etkisi
    küçüktür, çünkü toplam momenti büyük olaylar belirler.
    """
    mw = df["mw"] if "mw" in df.columns else pd.Series(np.nan, index=df.index)
    return mw.fillna(df["mag"])


# ---------------------------------------------------------------- artçı ayıklama
def _gk_window(m):
    """Gardner & Knopoff (1974) uzay (km) ve zaman (gün) pencereleri."""
    d = 10 ** (0.1238 * m + 0.983)
    t = np.where(m >= 6.5, 10 ** (0.032 * m + 2.7389), 10 ** (0.5409 * m - 0.547))
    return d, t


def decluster_gardner_knopoff(df: pd.DataFrame) -> pd.Series:
    """Ana şok olmayan olayları işaretle. Dönen seri True = bağımsız olay.

    Basit pencere yöntemi: büyükten küçüğe her olay, penceresindeki daha
    küçük olayları (öncü ve artçı) bağımlı işaretler. ETAS tabanlı
    stokastik ayıklama ileride bununla karşılaştırılmalı.
    """
    n = len(df)
    t = epoch_seconds(df["time"]) / 86400.0
    lat = df["latitude"].to_numpy()
    lon = df["longitude"].to_numpy()
    m = df["mag"].to_numpy()
    d_win, t_win = _gk_window(m)
    independent = np.ones(n, dtype=bool)
    order = np.argsort(-m, kind="stable")
    sort_t = np.argsort(t)
    ts = t[sort_t]
    for i in order:
        if not independent[i]:
            continue
        lo = np.searchsorted(ts, t[i] - t_win[i])
        hi = np.searchsorted(ts, t[i] + t_win[i], side="right")
        cand = sort_t[lo:hi]
        cand = cand[(cand != i) & (m[cand] <= m[i]) & independent[cand]]
        if cand.size == 0:
            continue
        dist = haversine_km(lat[i], lon[i], lat[cand], lon[cand])
        independent[cand[dist <= d_win[i]]] = False
    return pd.Series(independent, index=df.index, name="independent")
