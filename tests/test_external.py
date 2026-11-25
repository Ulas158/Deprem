import numpy as np
import pandas as pd

from deprem.catalog.schema import normalize
from deprem.external.astro import fortnightly_phase, semidiurnal_phase, tidal_gravity, julian_day
from deprem.external.dynamic import remote_triggering
from deprem.external.elements import element_external, idw, stations_to_elements
from deprem.external.periodic import schuster
from deprem.external.sites import active_at, load_sites, site_features
from deprem.external.swarms import cluster_stats, link_clusters
from deprem.faults.fill import (FILL_ALONG, FILL_DEFAULT, FILL_NEIGHBOR, FILL_NONE, FILL_REAL,
                                fill_missing, gaussian_log_prior, trainable_vector)
from deprem.faults.parameters import REGISTRY, ensure_columns
from deprem.features.table import build_table
from deprem.forecast.latent import LATENT
from deprem.region import Region

T0 = pd.Timestamp("2010-01-01", tz="UTC")


def _ev(times, lat, lon, mag, depth=10.0):
    n = len(times)
    return normalize(pd.DataFrame({
        "time": pd.to_datetime(times, utc=True), "latitude": lat, "longitude": lon,
        "depth_km": depth, "mag": mag, "mag_type": "ml", "source": "emsc", "agency": "X",
        "event_id": [str(i) for i in range(n)]}))


def test_moon_phases_and_tide_oscillates():
    t = pd.to_datetime(["2024-01-11 11:57", "2024-01-25 17:54", "2024-01-18 03:52"], utc=True)
    fp = fortnightly_phase(t)
    assert abs(fp[0]) < 5 and abs(fp[1]) < 5 and abs(abs(fp[2]) - 180) < 5
    g = tidal_gravity(julian_day(pd.date_range("2024-03-01", periods=48, freq="h", tz="UTC")), 39, 35)
    assert g.max() - g.min() > 1.0  # yarı günlük salınım


def test_random_times_have_no_tidal_signal():
    rng = np.random.default_rng(0)
    t = T0 + pd.to_timedelta(rng.uniform(0, 3650, 3000), unit="D")
    ph = semidiurnal_phase(t, 39.0, 35.0)
    assert np.isfinite(ph).mean() > 0.95
    lp, _, _ = schuster(ph)
    assert lp > -3


def test_schuster_detects_concentration():
    rng = np.random.default_rng(1)
    lp, mean, n = schuster(rng.normal(30, 40, 500))
    assert lp < -20 and abs(mean - 30) < 10
    assert np.isnan(schuster([1, 2, 3])[0])


def test_swarm_detection_and_diffusivity():
    rng = np.random.default_rng(2)
    # akışkan sürüsü: D = 1 m²/s, r = sqrt(4πDt) zarfı içinde göç, büyüklükler benzer
    D, n = 1.0, 200
    t_s = np.sort(rng.uniform(3600, 20 * 86400, n))
    r_m = np.sqrt(4 * np.pi * D * t_s) * np.sqrt(rng.uniform(0, 1, n))
    ang = rng.uniform(0, 2 * np.pi, n)
    lat = 39.0 + r_m * np.sin(ang) / 111195
    lon = 35.0 + r_m * np.cos(ang) / (111195 * np.cos(np.radians(39)))
    sw = _ev(np.r_[T0, T0 + pd.to_timedelta(t_s, unit="s")], np.r_[39.0, lat], np.r_[35.0, lon],
             np.r_[2.0, np.round(rng.uniform(2.0, 3.0, n), 1)])
    # ana şok-artçı: büyük olay başta
    t_a = np.sort(rng.exponential(2 * 86400, 100))
    ms = _ev(np.r_[T0 + pd.Timedelta(days=200), T0 + pd.Timedelta(days=200) + pd.to_timedelta(t_a, unit="s")],
             np.r_[40.5, 40.5 + rng.normal(0, 0.02, 100)], np.r_[30.0, 30.0 + rng.normal(0, 0.02, 100)],
             np.r_[5.5, np.round(rng.uniform(2.0, 3.5, 100), 1)])
    ev = pd.concat([sw, ms], ignore_index=True)
    cs = cluster_stats(ev, link_clusters(ev))
    assert len(cs) == 2
    s = cs.sort_values("lat").iloc[0]
    m = cs.sort_values("lat").iloc[1]
    assert s["is_swarm"] and not m["is_swarm"]
    assert 0.3 < s["D_m2_s"] < 3.0


def test_remote_triggering_beta():
    rng = np.random.default_rng(3)
    remote = pd.DataFrame({"time": [T0 + pd.Timedelta(days=100 * k) for k in range(1, 11)],
                           "latitude": 0.0, "longitude": 100.0, "mag": 8.0})
    bg = T0 + pd.to_timedelta(rng.uniform(0, 1100, 1100), unit="D")
    burst = np.concatenate([[r + pd.Timedelta(hours=h) for h in range(1, 30)] for r in remote["time"]])
    local = pd.DataFrame({"time": pd.to_datetime(np.r_[bg, burst], utc=True)}).assign(cell=0)
    res = remote_triggering(local, remote, pd.Series([39.0], index=[0]), pd.Series([35.0], index=[0]), "cell")
    assert res.loc[0, "dyn_n_remote"] == 10 and res.loc[0, "dyn_beta"] > 4
    calm = local.iloc[:1100]
    res0 = remote_triggering(calm, remote, pd.Series([39.0], index=[0]), pd.Series([35.0], index=[0]), "cell")
    assert abs(res0.loc[0, "dyn_beta"]) < 3


def test_sites_active_and_distances():
    s = load_sites()
    assert {"baraj", "jeotermal", "maden_komur", "petrol"} <= set(s["category"])
    assert s[["lat", "lon"]].notna().all().all()
    copler = s[s["id"] == "mine_copler"]
    assert active_at(copler, 2015).iloc[0] and not active_at(copler, 2025).iloc[0]
    f = site_features(np.array([37.95]), np.array([28.84]), s, pd.Timestamp("2020-01-01", tz="UTC"))
    assert f["site_geothermal_dist_km"].iloc[0] < 3
    # Ilısu dolumundan önce rezervuar yok sayılmalı
    f0 = site_features(np.array([37.53]), np.array([41.85]), s, pd.Timestamp("2010-01-01", tz="UTC"))
    f1 = site_features(np.array([37.53]), np.array([41.85]), s, pd.Timestamp("2021-01-01", tz="UTC"))
    assert f0["site_reservoir_dist_km"].iloc[0] > 20 and f1["site_reservoir_dist_km"].iloc[0] < 2


def test_registry_and_latent_extended():
    for n in ("radon_bq_m3", "co2_flux_g_m2_d", "he3_he4_rc_ra", "groundwater_level_m",
              "pore_pressure_lambda", "swarm_frac_10y", "tidal_semidiurnal_log10p",
              "dyn_trigger_beta", "site_reservoir_dist_km", "site_mining_dist_km"):
        assert n in REGISTRY
    assert REGISTRY["tidal_stress_sensitivity"].status == "seismicity"
    for n in ("hydraulic_diffusivity", "tidal_gain", "fluid_gain", "injection_gain"):
        assert n in LATENT
    el = ensure_columns(pd.DataFrame({"lat": [39.0], "lon": [35.0]}))
    assert "radon_bq_m3" in el.columns


def test_idw_and_stations_respect_time():
    v, d = idw([39.0, 41.0], [35.0, 35.0], [39.0, 39.1], [35.0, 35.0], [10.0, 20.0], max_km=30)
    assert 10 < v[0] < 15 and np.isnan(v[1])
    el = pd.DataFrame({"lat": [39.0], "lon": [35.0]})
    st = pd.DataFrame({"lat": [39.0, 39.0], "lon": [35.0, 35.0], "param": "radon_bq_m3",
                       "value": [100.0, 900.0], "time": ["2009-01-01", "2011-01-01"]})
    out = stations_to_elements(el, st, T0)
    assert out["radon_bq_m3"].iloc[0] == 100.0  # gelecekteki ölçüm kullanılmaz


def test_fill_missing_along_fault_then_neighbors():
    el = pd.DataFrame({"fault_id": [1] * 5 + [2] * 3, "along_km": [0, 1, 2, 3, 4, 0, 1, 2],
                       "lat": [39.0] * 5 + [39.1] * 3,
                       "lon": [35, 35.01, 35.02, 35.03, 35.04, 35, 35.01, 35.02],
                       "x": [1, np.nan, np.nan, 4, np.nan, np.nan, np.nan, np.nan],
                       "pore_pressure_lambda": np.nan, "radon_bq_m3": np.nan})
    f = fill_missing(el, ["x", "pore_pressure_lambda", "radon_bq_m3"])
    assert list(f["x__fill"]) == [FILL_REAL, FILL_ALONG, FILL_ALONG, FILL_REAL, FILL_ALONG] + [FILL_NEIGHBOR] * 3
    assert f["x"].iloc[1] == 2.0 and f["x"].iloc[2] == 3.0 and f["x"].iloc[4] == 4.0
    assert (f["x__sd"].iloc[5:] > f["x__sd"].iloc[1]).all()
    assert (f["pore_pressure_lambda__fill"] == FILL_DEFAULT).all()
    assert (f["radon_bq_m3__fill"] == FILL_NONE).all() and f["radon_bq_m3"].isna().all()
    idx, mu, sd = trainable_vector(f, "x")
    assert len(idx) == 6 and (sd > 0).all()
    assert gaussian_log_prior(mu, mu, sd) == 0.0


def test_build_table_with_external_features():
    rng = np.random.default_rng(4)
    n = 1500
    cat = _ev(T0 + pd.to_timedelta(np.sort(rng.uniform(0, 3000, n)), unit="D"),
              rng.uniform(38, 40, n), rng.uniform(30, 32, n), np.round(rng.uniform(2, 4, n), 1))
    cat["mw_best"] = cat["mag"]
    remote = pd.DataFrame({"time": [T0 + pd.Timedelta(days=500)], "latitude": 0.0,
                           "longitude": 100.0, "mag": 8.0})
    starts = [pd.Timestamp("2016-01-01", tz="UTC"), pd.Timestamp("2017-01-01", tz="UTC")]
    tab = build_table(cat, Region("t", 38, 40, 30, 32), 1.0, starts, pd.DateOffset(years=1), [3.0], 2.0,
                      ctx_extra={"remote": remote, "tidal": False})
    for c in ("seasonal_log10p", "swarm_n_10y", "dyn_trigger_beta", "site_reservoir_dist_km",
              "heat_flow_mw_m2"):
        assert c in tab.columns
    assert tab["seasonal_log10p"].notna().any()


def test_element_external_runs():
    rng = np.random.default_rng(5)
    n = 600
    ev = _ev(T0 + pd.to_timedelta(np.sort(rng.uniform(0, 3000, n)), unit="D"),
             39 + rng.normal(0, 0.05, n), 35 + rng.normal(0, 0.05, n), np.round(rng.uniform(2, 4, n), 1))
    el = pd.DataFrame({"lat": [39.0, 39.0], "lon": [35.0, 36.0], "depth_km": [5.0, 5.0]})
    out = element_external(el, ev, pd.Timestamp("2019-01-01", tz="UTC"), 2.0)
    assert out.loc[0, "seasonal_log10p"] <= 0 and np.isnan(out.loc[1, "seasonal_log10p"])
    assert "site_mining_dist_km" in out.columns
