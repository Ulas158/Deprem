"""Gaz/sıvı ve dış etki parametrelerini her fay elemanına (1 km × 1 km) aktarma.

* Ölçüm istasyonları (radon, CO2, He, yeraltı suyu, kaynak suyu ...):
  ``stations`` tablosu (``lat, lon, param, value[, time]``) verilirse her
  parametre elemana **ters uzaklık ağırlıklı** (IDW, üs 2) aktarılır; yalnızca
  ``max_km`` içindeki istasyonlar kullanılır, yoksa alan boş (NaN) kalır.
  Zaman sütunu varsa yalnızca ``t_ref`` öncesi ölçümler, istasyon başına en
  son değer kullanılır (gelecekten sızıntı yok).
* Katalogdan: gelgit/mevsim Schuster testleri, sürü sayısı ve oranı,
  hidrolik yayınım (``radius_km`` içindeki olaylardan).
* Saha listesinden: baraj/maden/jeotermal/petrol-gaz uzaklıkları.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ..geo import haversine_km
from .params import PARAMS
from .periodic import event_phases
from .sites import load_sites, site_features
from .swarms import cluster_stats, link_clusters

STATION_PARAMS = [n for n, _, _, s, _, _ in PARAMS if s == "placeholder"] + ["heat_flow_mw_m2", "vp_vs_ratio"]


def idw(lat, lon, st_lat, st_lon, values, max_km: float = 30.0, power: float = 2.0):
    """Ters uzaklık ağırlıklı aktarım. Döndürür: (değer, en yakın istasyon uzaklığı)."""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    out = np.full(len(lat), np.nan)
    dmin = np.full(len(lat), np.nan)
    if len(values) == 0:
        return out, dmin
    for i in range(0, len(lat), 5000):
        d = haversine_km(lat[i:i + 5000, None], lon[i:i + 5000, None],
                         np.asarray(st_lat)[None, :], np.asarray(st_lon)[None, :])
        w = np.where(d <= max_km, 1.0 / np.maximum(d, 0.5) ** power, 0.0)
        ws = w.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[i:i + 5000] = np.where(ws > 0, (w * np.asarray(values)[None, :]).sum(axis=1) / ws, np.nan)
        dmin[i:i + 5000] = d.min(axis=1)
    return out, dmin


def stations_to_elements(elements: pd.DataFrame, stations: pd.DataFrame | None,
                         t_ref: pd.Timestamp | None = None, max_km: float = 30.0) -> pd.DataFrame:
    out = pd.DataFrame(index=elements.index)
    if stations is None or stations.empty:
        return out
    st = stations.copy()
    if "time" in st.columns and t_ref is not None:
        st["time"] = pd.to_datetime(st["time"], utc=True)
        st = st[st["time"] < t_ref].sort_values("time").groupby(["param", "lat", "lon"], as_index=False).last()
    dmins = []
    for param, g in st.groupby("param"):
        v, d = idw(elements["lat"], elements["lon"], g["lat"], g["lon"], g["value"], max_km)
        out[param] = v
        dmins.append(d)
    if dmins:
        out["station_dist_km"] = np.nanmin(np.vstack(dmins), axis=0)
    return out


def _local_xyz(lat, lon, depth, lat0=39.0, lon0=35.0):
    x = (np.asarray(lon) - lon0) * 111.195 * np.cos(np.radians(lat0))
    y = (np.asarray(lat) - lat0) * 111.195
    return np.column_stack([x, y, np.asarray(depth, dtype=float)])


def _schuster_by_group(rows, phase_deg, n_groups, min_n=10):
    ph = np.radians(phase_deg)
    ok = np.isfinite(ph)
    c = np.bincount(rows[ok], weights=np.cos(ph[ok]), minlength=n_groups)
    s = np.bincount(rows[ok], weights=np.sin(ph[ok]), minlength=n_groups)
    n = np.bincount(rows[ok], minlength=n_groups).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        lp = -(c * c + s * s) / n / np.log(10)
    return np.where(n >= min_n, lp, np.nan)


def element_external(elements: pd.DataFrame, events: pd.DataFrame, t_ref: pd.Timestamp, mc: float,
                     sites: pd.DataFrame | None = None, stations: pd.DataFrame | None = None,
                     radius_km: float = 20.0, phases: pd.DataFrame | None = None) -> pd.DataFrame:
    """Her eleman için gaz/sıvı ve dış etki parametreleri (``t_ref`` öncesi veriyle)."""
    n_el = len(elements)
    out = pd.DataFrame(index=elements.index)
    ev = events[(events["time"] < t_ref) & (events["time"] >= t_ref - pd.DateOffset(years=10))
                & (events["mag"] >= mc)]
    el_xyz = _local_xyz(elements["lat"], elements["lon"], elements["depth_km"])
    tree = cKDTree(el_xyz[:, :2])
    if len(ev):
        # --- gelgit / mevsim (bağımsız olaylar)
        ind = ev[ev["independent"].astype(bool)] if "independent" in ev.columns else ev
        if phases is None:
            phases = event_phases(ind, tidal=False)
        ph = phases.reindex(ind.index)
        ev_xy = _local_xyz(ind["latitude"], ind["longitude"], np.zeros(len(ind)))[:, :2]
        pairs = tree.sparse_distance_matrix(cKDTree(ev_xy), radius_km, output_type="coo_matrix")
        r, c = pairs.row, pairs.col
        for src, name in (("phase_lunar", "tidal_semidiurnal_log10p"),
                          ("phase_fortnight", "tidal_fortnight_log10p"),
                          ("phase_annual", "seasonal_log10p"),
                          ("phase_hour", "hour_of_day_log10p")):
            out[name] = _schuster_by_group(r, ph[src].to_numpy()[c], n_el)
        out["tidal_stress_sensitivity"] = out["tidal_semidiurnal_log10p"]
        # --- sürüler
        cid = link_clusters(ev)
        cs = cluster_stats(ev, cid)
        if not cs.empty:
            cxy = _local_xyz(cs["lat"], cs["lon"], np.zeros(len(cs)))[:, :2]
            pairs = tree.sparse_distance_matrix(cKDTree(cxy), radius_km, output_type="coo_matrix")
            r, c = pairs.row, pairs.col
            sw = cs["is_swarm"].to_numpy()[c].astype(float)
            n_cl = np.bincount(r, minlength=n_el).astype(float)
            n_sw = np.bincount(r, weights=sw, minlength=n_el)
            out["swarm_n_10y"] = n_sw
            with np.errstate(invalid="ignore", divide="ignore"):
                out["swarm_frac_10y"] = np.where(n_cl > 0, n_sw / n_cl, np.nan)
                D = cs["D_m2_s"].to_numpy()[c]
                okd = (sw > 0) & np.isfinite(D) & (D > 0)
                s = np.bincount(r[okd], weights=np.log10(D[okd]), minlength=n_el)
                k = np.bincount(r[okd], minlength=n_el)
                out["hydraulic_diffusivity_log10"] = np.where(k > 0, s / k, np.nan)
    # --- sahalar
    sites = load_sites() if sites is None else sites
    sf = site_features(elements["lat"].to_numpy(), elements["lon"].to_numpy(), sites, t_ref)
    sf.index = elements.index
    out = pd.concat([out, sf, stations_to_elements(elements, stations, t_ref)], axis=1)
    return out
