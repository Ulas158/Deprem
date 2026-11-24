"""Fay başına moment birikimi ve deprem kataloğundaki moment boşalımı.

Temel fikir (Brune 1968; Kostrov 1974): bir fay yılda

    M0_dot = mu * A * v      (mu: makaslama modülü, A: fay alanı, v: kayma hızı)

kadar sismik moment biriktirir. Aynı süre içinde depremlerde boşalan
toplam moment ile karşılaştırıldığında "moment açığı" elde edilir.

Bu hesap bir **tahmin değildir**. Sonuçlar şunlara karşı çok duyarlıdır:
kayma hızı belirsizliği (çoğu segmentte 2-4 kat), sismik bağlaşım oranı
(aseismik krip, ör. İsmetpaşa), sismojenik genişlik, katalog süresi
(1900 sonrası, oysa yinelenme süreleri yüzlerce yıl) ve depremlerin
faylara atanması. Bu yüzden çıktı, test edilecek bir hipotez değişkeni
olarak ele alınmalıdır.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..geo import point_to_polyline_km

MU = 3.0e10  # Pa


def mw_to_m0(mw):
    """Hanks & Kanamori (1979): M0 [N·m]."""
    return 10 ** (1.5 * np.asarray(mw, dtype=float) + 9.1)


def m0_to_mw(m0):
    return (np.log10(np.asarray(m0, dtype=float)) - 9.1) / 1.5


def add_geometry(faults: pd.DataFrame, mu: float = MU) -> pd.DataFrame:
    """Genişlik, alan, moment birikim hızı ve kaba maksimum büyüklük ekle."""
    f = faults.copy()
    sin_dip = np.sin(np.radians(f["dip"].clip(10, 90)))
    f["width_km"] = (f["lower_depth_km"] - f["upper_depth_km"]).clip(lower=1) / sin_dip
    f["area_km2"] = f["length_km"] * f["width_km"]
    v = f["slip_rate_mm_yr"] * 1e-3  # m/yıl
    f["m0_rate_Nm_yr"] = mu * f["area_km2"] * 1e6 * v
    for tag in ("min", "max"):
        f[f"m0_rate_{tag}"] = mu * f["area_km2"] * 1e6 * f[f"slip_rate_{tag}"] * 1e-3
    # Leonard (2010) kabuk içi doğrultu atımlı: Mw = log10(A km2) + 3.99 (alan ölçekleme)
    f["mmax_area"] = np.log10(f["area_km2"].clip(lower=1)) + 3.99
    return f


def assign_events_to_faults(events: pd.DataFrame, faults: pd.DataFrame,
                            max_dist_km: float = 15.0) -> pd.Series:
    """Her olayı en yakın fay segmentine ata (uzaklık sınırı içinde).

    Büyük depremlerin episantrı kırılmanın bir ucunda olabilir; bu yüzden
    atama kabadır ve segment düzeyinden çok fay zonu düzeyinde yorumlanmalı.
    """
    if events.empty or faults.empty:
        return pd.Series(pd.NA, index=events.index, dtype="object")
    lat = events["latitude"].to_numpy()
    lon = events["longitude"].to_numpy()
    dist = np.full((len(events), len(faults)), np.inf)
    for j, coords in enumerate(faults["coords"]):
        dist[:, j] = point_to_polyline_km(lat, lon, coords)
    j_best = dist.argmin(axis=1)
    d_best = dist[np.arange(len(events)), j_best]
    ids = faults["fault_id"].to_numpy()[j_best].astype(object)
    ids[d_best > max_dist_km] = None
    return pd.Series(ids, index=events.index, name="fault_id")


def zone_budget(faults: pd.DataFrame, events: pd.DataFrame, zones: dict,
                start: pd.Timestamp, end: pd.Timestamp, mw_col: str = "mw_best",
                coupling: float = 1.0) -> pd.DataFrame:
    """Alt bölge (zon) başına birikim / boşalım karşılaştırması.

    Bir fay, orta noktası zonun içindeyse o zona sayılır. ``coupling``
    sismik bağlaşım katsayısıdır (1 = tüm kayma depremle boşalır).
    """
    years = (end - start).days / 365.25
    mid_lat = faults["coords"].map(lambda c: np.mean([p[1] for p in c]))
    mid_lon = faults["coords"].map(lambda c: np.mean([p[0] for p in c]))
    mid = pd.DataFrame({"latitude": mid_lat, "longitude": mid_lon})
    ev = events[(events["time"] >= start) & (events["time"] < end)]
    rows = []
    for key, reg in zones.items():
        fz = faults[reg.contains(mid)]
        ez = ev[reg.contains(ev)]
        acc = coupling * fz["m0_rate_Nm_yr"].sum() * years
        acc_lo = coupling * fz["m0_rate_min"].sum() * years
        acc_hi = coupling * fz["m0_rate_max"].sum() * years
        rel = mw_to_m0(ez[mw_col]).sum()
        largest = ez.loc[ez[mw_col].idxmax()] if len(ez) else None
        rows.append({
            "zone": reg.name, "n_faults": len(fz), "years": round(years, 1),
            "m0_accumulated": acc, "m0_acc_min": acc_lo, "m0_acc_max": acc_hi,
            "m0_released": rel,
            "release_ratio": rel / acc if acc > 0 else np.nan,
            "deficit_mw_equiv": float(m0_to_mw(acc - rel)) if acc > rel else np.nan,
            "largest_event_mw": float(largest[mw_col]) if largest is not None else np.nan,
            "largest_event_year": int(largest["time"].year) if largest is not None else None,
        })
    return pd.DataFrame(rows)


def fault_summary(faults: pd.DataFrame, events: pd.DataFrame, start: pd.Timestamp,
                  end: pd.Timestamp, mw_col: str = "mw_best", min_mw: float = 6.0,
                  max_dist_km: float = 15.0) -> pd.DataFrame:
    """Segment başına: birikim hızı, katalogdaki son büyük olay ve geçen süre.

    ``elapsed_x_rate`` = son M>=min_mw olaydan (yoksa katalog başından) bu yana
    biriken moment; ``mw_equiv`` bunun tek depremle boşalması halinde
    karşılık geleceği Mw. Katalog başlangıcından önceki depremler (ör. 1766,
    1894) bilinmediği için bu değerler alt sınır değil, yalnızca
    katalog-içi bir göstergedir.
    """
    big = events[(events[mw_col] >= min_mw) & (events["time"] >= start)
                 & (events["time"] < end)].copy()
    big["fault_id"] = assign_events_to_faults(big, faults, max_dist_km)
    last = big.dropna(subset=["fault_id"]).groupby("fault_id")["time"].max()
    out = faults.drop(columns=["coords"]).copy()
    out["last_big_event"] = out["fault_id"].map(last)
    ref = out["last_big_event"].fillna(start)
    out["years_since_ref"] = (end - ref).dt.days / 365.25
    out["m0_since_ref"] = out["m0_rate_Nm_yr"] * out["years_since_ref"]
    with np.errstate(divide="ignore"):
        out["mw_equiv_since_ref"] = m0_to_mw(out["m0_since_ref"].where(out["m0_since_ref"] > 0))
    return out


def event_fault_geometry(events: pd.DataFrame, faults: pd.DataFrame,
                         max_dist_km: float = 25.0) -> pd.DataFrame:
    """Her olayın en yakın fay segmentine göre konumu ve derinliği.

    Sütunlar:
      fault_id            en yakın segment (max_dist_km dışında boş)
      fault_dist_km       faya dik (yüzey izine) uzaklık
      along_frac          segment boyunca göreli konum (0 = bir uç, 1 = diğer uç)
      end_dist_km         en yakın segment ucuna uzaklık (kırılma başlangıç/durma
                          bölgeleri hipotezi için)
      depth_rel           derinlik / sismojenik alt sınır (1'e yakın = kilitli
                          zonun tabanı; >1 = sismojenik zonun altı)
      fault_slip_rate     segmentin kayma hızı (mm/yıl)
    """
    from ..geo import project_to_polyline

    n = len(events)
    out = pd.DataFrame(index=events.index, data={
        "fault_id": pd.Series([None] * n, index=events.index, dtype="object"),
        "fault_dist_km": np.nan, "along_frac": np.nan, "end_dist_km": np.nan,
        "depth_rel": np.nan, "fault_slip_rate": np.nan})
    if n == 0 or faults.empty:
        return out
    lat = events["latitude"].to_numpy()
    lon = events["longitude"].to_numpy()
    best_d = np.full(n, np.inf)
    best_j = np.full(n, -1)
    best_along = np.full(n, np.nan)
    best_len = np.full(n, np.nan)
    for j, coords in enumerate(faults["coords"]):
        d, along, length = project_to_polyline(lat, lon, coords)
        better = d < best_d
        best_d[better] = d[better]
        best_j[better] = j
        best_along[better] = along[better]
        best_len[better] = length
    ok = best_d <= max_dist_km
    jj = best_j[ok]
    out.loc[ok, "fault_id"] = faults["fault_id"].to_numpy()[jj]
    out.loc[ok, "fault_dist_km"] = best_d[ok]
    L = np.where(best_len[ok] > 0, best_len[ok], np.nan)
    out.loc[ok, "along_frac"] = best_along[ok] / L
    out.loc[ok, "end_dist_km"] = np.minimum(best_along[ok], best_len[ok] - best_along[ok])
    out.loc[ok, "depth_rel"] = (events["depth_km"].to_numpy()[ok]
                                / faults["lower_depth_km"].to_numpy()[jj])
    out.loc[ok, "fault_slip_rate"] = faults["slip_rate_mm_yr"].to_numpy()[jj]
    return out
