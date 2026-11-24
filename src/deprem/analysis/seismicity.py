"""Zaman ve mekânda sismisite betimlemesi."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..geo import haversine_km
from ..region import Region
from ..stats.gr import b_value_mle, mc_maxc


def rate_series(df: pd.DataFrame, freq: str = "MS", min_mag: float | None = None) -> pd.Series:
    d = df if min_mag is None else df[df["mag"] >= min_mag]
    s = d.set_index("time").resample(freq).size()
    return s.rename("count")


def b_value_timeseries(df: pd.DataFrame, mc: float, window: int = 300, step: int = 50,
                       dm: float = 0.1) -> pd.DataFrame:
    """Kayan pencerede b(t): her pencere ``window`` olay (M>=Mc)."""
    d = df[df["mag"] >= mc - dm / 2].sort_values("time")
    rows = []
    for i in range(0, len(d) - window + 1, step):
        w = d.iloc[i:i + window]
        b, s, n = b_value_mle(w["mag"], mc, dm)
        rows.append({"time": w["time"].iloc[-1], "b": b, "sigma_b": s, "n": n})
    return pd.DataFrame(rows)


def grid_edges(region: Region, cell_deg: float):
    lats = np.arange(region.min_lat, region.max_lat + 1e-9, cell_deg)
    lons = np.arange(region.min_lon, region.max_lon + 1e-9, cell_deg)
    return lats, lons


def gridded_counts(df: pd.DataFrame, region: Region, cell_deg: float = 0.25) -> tuple:
    lats, lons = grid_edges(region, cell_deg)
    h, _, _ = np.histogram2d(df["latitude"], df["longitude"], bins=[lats, lons])
    return h, lats, lons


def b_value_map(df: pd.DataFrame, region: Region, cell_deg: float = 0.25,
                n_nearest: int = 150, max_radius_km: float = 75.0,
                min_events: int = 50, dm: float = 0.1) -> tuple:
    """Her ızgara noktasında en yakın N olaydan yerel Mc ve b.

    Wiemer & Wyss (1997) yaklaşımı. Yerel Mc maksimum eğrilikle bulunur;
    yarıçap sınırını aşan veya az olaylı düğümler NaN bırakılır.
    """
    lats, lons = grid_edges(region, cell_deg)
    clat = (lats[:-1] + lats[1:]) / 2
    clon = (lons[:-1] + lons[1:]) / 2
    b = np.full((clat.size, clon.size), np.nan)
    mc_map = np.full_like(b, np.nan)
    la = df["latitude"].to_numpy()
    lo = df["longitude"].to_numpy()
    m = df["mag"].to_numpy()
    for i, y in enumerate(clat):
        for j, x in enumerate(clon):
            d = haversine_km(y, x, la, lo)
            idx = np.argpartition(d, min(n_nearest, d.size - 1))[:n_nearest]
            idx = idx[d[idx] <= max_radius_km]
            if idx.size < min_events:
                continue
            mc = mc_maxc(m[idx], dm)
            bb, _, n = b_value_mle(m[idx], mc, dm)
            if n >= min_events:
                b[i, j] = bb
                mc_map[i, j] = mc
    return b, mc_map, lats, lons
