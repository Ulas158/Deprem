"""Hücre × pencere özellik tablosu için gaz/sıvı ve dış etki özellikleri.

Hepsi ``deprem.features.table.register_feature`` ile kaydedilir ve yalnızca
pencere başlangıcından (t0) önceki veriyi kullanır. ``build_table`` içine
``ctx_extra`` ile şunlar verilebilir:

* ``remote``     USGS küresel büyük deprem listesi (uzak tetikleme için)
* ``sites``      saha listesi (varsayılan: external_data/insan_kaynakli_sahalar.csv)
* ``heat_flow``  lat, lon, value sütunlu ısı akısı noktaları (IHFC indirilirse)
* ``tidal``      False ise yarı günlük gelgit fazı hesaplanmaz (hız için)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..features.table import register_feature
from ..geo import haversine_km
from .dynamic import remote_triggering
from .periodic import event_phases, schuster
from .sites import load_sites, site_features
from .swarms import cluster_stats, link_clusters


def _phases(past: pd.DataFrame, ctx: dict) -> pd.DataFrame:
    cache = ctx.setdefault("_phase_cache", pd.DataFrame(columns=["phase_annual", "phase_fortnight", "phase_semidiurnal"]))
    missing = past.index.difference(cache.index)
    if len(missing):
        new = event_phases(past.loc[missing], tidal=ctx.get("tidal", True))
        cache = new if cache.empty else pd.concat([cache, new])
        ctx["_phase_cache"] = cache
    return cache.loc[past.index]


@register_feature("periodic_loads")
def f_periodic(past, cells, t0, ctx):
    """Gelgit ve mevsim fazı Schuster testleri (son 10 yıl, M≥Mc, bağımsız olaylar)."""
    w = past[(past["time"] >= t0 - pd.DateOffset(years=10)) & (past["mag"] >= ctx["mc"])]
    if "independent" in w.columns:
        w = w[w["independent"].astype(bool)]
    cols = {"phase_semidiurnal": "tidal_semidiurnal_log10p",
            "phase_fortnight": "tidal_fortnight_log10p", "phase_annual": "seasonal_log10p"}
    if w.empty:
        return pd.DataFrame(index=cells["cell_id"], columns=list(cols.values()), dtype=float)
    ph = _phases(w, ctx)
    out = {}
    for src, name in cols.items():
        out[name] = ph[src].groupby(w["cell_id"]).apply(lambda s: schuster(s.to_numpy())[0])
    return pd.DataFrame(out).reindex(cells["cell_id"])


@register_feature("swarms_fluid")
def f_swarms(past, cells, t0, ctx):
    """Son 10 yıldaki sürüler ve göç hızından hidrolik yayınım."""
    w = past[(past["time"] >= t0 - pd.DateOffset(years=10)) & (past["mag"] >= ctx["mc"])]
    names = ["swarm_n_10y", "swarm_frac_10y", "hydraulic_diffusivity_log10"]
    if len(w) < 5:
        return pd.DataFrame(index=cells["cell_id"], columns=names, dtype=float)
    cid = link_clusters(w)
    cs = cluster_stats(w, cid)
    if cs.empty:
        return pd.DataFrame(index=cells["cell_id"], columns=names, dtype=float).assign(swarm_n_10y=0.0)
    first = w.assign(_c=cid).sort_values("time").groupby("_c")["cell_id"].first()
    cs["cell_id"] = cs["cluster"].map(first)
    g = cs.groupby("cell_id")
    out = pd.DataFrame({
        "swarm_n_10y": g["is_swarm"].sum(),
        "swarm_frac_10y": g["is_swarm"].mean(),
        "hydraulic_diffusivity_log10": np.log10(cs[cs["is_swarm"]].groupby("cell_id")["D_m2_s"].median()),
    }).reindex(cells["cell_id"])
    out["swarm_n_10y"] = out["swarm_n_10y"].fillna(0)
    return out


@register_feature("dynamic_triggering")
def f_dynamic(past, cells, t0, ctx):
    remote = ctx.get("remote")
    names = ["dyn_trigger_beta", "dyn_n_remote"]
    if remote is None or remote.empty:
        return pd.DataFrame(index=cells["cell_id"], columns=names, dtype=float)
    r = remote[(remote["time"] < t0) & (remote["time"] >= t0 - pd.DateOffset(years=10))]
    w = past[(past["time"] >= t0 - pd.DateOffset(years=10) - pd.DateOffset(days=40)) & (past["mag"] >= ctx["mc"])]
    lat_c = cells.set_index("cell_id")["lat_c"]
    lon_c = cells.set_index("cell_id")["lon_c"]
    res = remote_triggering(w, r, lat_c, lon_c, "cell_id")
    return res.rename(columns={"dyn_beta": "dyn_trigger_beta"})[names]


@register_feature("human_sites")
def f_sites(past, cells, t0, ctx):
    if "sites" not in ctx:
        ctx["sites"] = load_sites()
    out = site_features(cells["lat_c"].to_numpy(), cells["lon_c"].to_numpy(), ctx["sites"], t0)
    out.index = cells["cell_id"].to_numpy()
    return out


def nearest_point_value(lat, lon, pts: pd.DataFrame, radius_km: float = 50.0) -> np.ndarray:
    """Noktaların ``radius_km`` içindeki ölçümlerin ortancası (yoksa NaN)."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    out = np.full(len(lat), np.nan)
    if pts is None or pts.empty:
        return out
    for i in range(0, len(lat), 2000):
        d = haversine_km(lat[i:i + 2000, None], lon[i:i + 2000, None],
                         pts["lat"].to_numpy()[None, :], pts["lon"].to_numpy()[None, :])
        v = np.where(d <= radius_km, pts["value"].to_numpy()[None, :], np.nan)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out[i:i + 2000] = np.nanmedian(v, axis=1)
    return out


@register_feature("fluid_static")
def f_fluid_static(past, cells, t0, ctx):
    """Zamandan bağımsız akışkan göstergeleri; veri yoksa NaN (alan tanımlı kalır)."""
    if "_fluid_static" not in ctx:
        hf = ctx.get("heat_flow")
        ctx["_fluid_static"] = pd.DataFrame({
            "heat_flow_mw_m2": nearest_point_value(cells["lat_c"], cells["lon_c"], hf),
            "thermal_spring_dist_km": np.nan,
            "radon_anomaly_z": np.nan,
            "groundwater_anomaly_z": np.nan,
        }, index=cells["cell_id"].to_numpy())
    return ctx["_fluid_static"]
