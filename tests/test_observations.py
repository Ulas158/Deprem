import numpy as np
import pandas as pd

from deprem.observations import dvv, fit as pfit, gnss, repeaters, shortterm
from deprem.observations.registry import OBS_FEATURES, apply_transform, prior_for
from deprem.region import TURKEY


def test_registry_priors_and_groups():
    assert {f.evidence for f in OBS_FEATURES.values()} == {"güçlü", "zayıf"}
    # zayıf kanıtlılar yalnız kısa vadeli hipotez testinde, öncülleri temkinli
    for f in OBS_FEATURES.values():
        if f.evidence == "zayıf":
            assert f.horizon == "kısa" and abs(f.prior_mean) <= 0.15
    assert prior_for("gnss_strain_rate_ns")[0] > 0
    assert prior_for("dvv_change_pct")[0] < 0
    assert prior_for("gnss_strain_rate_ns__eksik") is None
    x = apply_transform("gnss_strain_rate_ns", [10.0, 100.0])
    assert np.allclose(x, [1, 2])


def _synthetic_ts(ve=20.0, vn=-5.0, step_t=2014.5, step=8.0, transient_from=None, tr_rate=0.0, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(2008.0, 2016.0, 1 / 365.25)
    e = ve * (t - t[0]) + step * (t >= step_t) + 1.5 * np.sin(2 * np.pi * t) + rng.normal(0, 1.5, len(t))
    n = vn * (t - t[0]) + rng.normal(0, 1.5, len(t))
    if transient_from:
        e += tr_rate * np.clip(t - transient_from, 0, None)
    u = rng.normal(0, 4, len(t))
    return pd.DataFrame({"t": t, "e": e, "n": n, "u": u, "se": 1.5, "sn": 1.5, "su": 4.0})


def test_gnss_velocity_and_step():
    ts = _synthetic_ts()
    f = gnss.station_year_features(ts, np.array([2014.5]), 2016.0)
    assert abs(f["ve"] - 20) < 0.5 and abs(f["vn"] + 5) < 0.5
    assert f["tr_z"] < 3


def test_gnss_transient_detected():
    ts = _synthetic_ts(transient_from=2015.0, tr_rate=15.0, step=0.0)
    f = gnss.station_year_features(ts, np.array([]), 2016.0)
    assert f["tr_e"] > 10 and f["tr_z"] > 3


def test_strain_rate_uniform_shear():
    rng = np.random.default_rng(1)
    lat = 39 + rng.uniform(-0.8, 0.8, 40)
    lon = 33 + rng.uniform(-1, 1, 40)
    x, y = gnss._xy(lat, lon)
    x0, y0 = gnss._xy(39.0, 33.0)
    # ve = 0.1 mm/yıl/km · (y − y0) → exy = 0.05 mm/yıl/km = 50 ns/yıl, 2. değişmez = 70.7
    sta = pd.DataFrame({"lat": lat, "lon": lon, "ve": 0.1 * (y - y0), "vn": np.zeros(40)})
    r = gnss.strain_rate_at([39.0], [33.0], sta)
    assert abs(r["max_shear"][0] - 50) < 2 and abs(r["strain_2inv"][0] - 70.7) < 3


def test_repeaters_found_and_not_aftershocks():
    rng = np.random.default_rng(2)
    t0 = pd.Timestamp("2012-01-01", tz="UTC")
    rep = pd.DataFrame({"time": [t0 + pd.Timedelta(days=200 * k + rng.integers(0, 10)) for k in range(8)],
                        "latitude": 40.7 + rng.normal(0, 0.002, 8), "longitude": 29.0 + rng.normal(0, 0.002, 8),
                        "depth_km": 10 + rng.normal(0, 0.2, 8), "mag": 2.5 + rng.normal(0, 0.05, 8)})
    # artçı dizisi: aynı yerde ama düzensiz, birkaç günde yoğun
    aft = pd.DataFrame({"time": [t0 + pd.Timedelta(days=1000 + float(d)) for d in [0, 0.1, 0.3, 1, 2, 40, 300]],
                        "latitude": 38.0, "longitude": 35.0, "depth_km": 8.0, "mag": 2.4})
    cat = pd.concat([rep, aft], ignore_index=True)
    fam = repeaters.find_families(cat)
    assert fam["family"].nunique() == 1 and len(fam) == 8
    cells = pd.DataFrame({"cell_id": [0, 1], "lat_c": [40.7, 37.0], "lon_c": [29.0, 40.0]})
    cf = repeaters.cell_features(fam, cells, pd.Timestamp("2016-06-01", tz="UTC"))
    assert cf.loc[0, "rep_slip_rate_cm_yr"] > 0 and np.isnan(cf.loc[1, "rep_slip_rate_cm_yr"])


def test_prior_fit_learns_and_falls_back():
    rng = np.random.default_rng(3)
    z = rng.normal(size=(4000, 1))
    y = rng.poisson(np.exp(-1 + 0.5 * z[:, 0]))
    beta, sd, _ = pfit.fit(z, np.zeros(4000), y, np.array([0.1]), np.array([0.15]))
    assert abs(beta[1] - 0.5) < 0.1 and sd[1] < 0.05
    # veri yoksa öncülde kalır
    beta0, sd0, _ = pfit.fit(z[:3], np.zeros(3), np.zeros(3), np.array([0.1]), np.array([0.15]))
    assert abs(beta0[1] - 0.1) < 0.05


def test_shortterm_detects_real_precursor_and_rejects_noise():
    rng = np.random.default_rng(4)
    boxes = shortterm.make_boxes(TURKEY, 2.0)
    wins = pd.date_range("2012-01-01", "2025-12-31", freq="5D")
    panel = pd.DataFrame([(b, w) for w in wins for b in boxes["box"]], columns=["box", "win"])
    sig = rng.normal(size=len(panel))
    noise = rng.normal(size=len(panel))
    panel["n30"], panel["n5"] = rng.poisson(1, len(panel)), rng.poisson(0.2, len(panel))
    panel["y"] = rng.poisson(0.02 * np.exp(0.8 * sig))
    OBS_FEATURES  # tec_z_max ve thermal_dts_z öncülleri kullanılır
    res, coefs = shortterm.evaluate(panel, {"tec_z_max": sig, "thermal_dts_z": noise},
                                    pd.Timestamp("2020-01-01"), pd.Timestamp("2026-01-01"))
    r = res.set_index("hipotez")
    assert r.loc["tec_z_max", "bilgi_kazancı_olay_başına"] > 0.1 and r.loc["tec_z_max", "p_değeri"] < 0.01
    assert abs(r.loc["thermal_dts_z", "bilgi_kazancı_olay_başına"]) < 0.05


def test_dvv_stretching_recovers_velocity_change():
    lag = np.arange(-400, 401) / dvv.FS
    rng = np.random.default_rng(5)
    ref = np.convolve(rng.normal(size=len(lag)), np.hanning(15), "same") * np.exp(-np.abs(lag) / 8)
    true = 0.004  # %0,4 hız artışı → varışlar erken: cur(t) = ref(t(1+dv))
    cur = np.interp(lag * (1 + true), lag, ref)
    est, cc = dvv.stretching(ref, cur, lag)
    assert abs(est - true) < 5e-4
    assert cc > 0.95


def test_eventwalk_model_equals_etas_when_eta_zero():
    from deprem.forecast.base import ForecastGrid
    from deprem.forecast.etas import ETAS
    from deprem.forecast.eventwalk import score_year, summarize
    rng = np.random.default_rng(6)
    n = 400
    t = pd.Timestamp("2010-01-01", tz="UTC") + pd.to_timedelta(np.sort(rng.uniform(0, 6 * 365, n)), unit="D")
    cat = pd.DataFrame({"time": t, "latitude": rng.uniform(37, 41, n), "longitude": rng.uniform(28, 40, n),
                        "mag": 3.5 + rng.exponential(0.43, n), "depth_km": 10.0})
    e = ETAS(m0=3.5)
    train = cat[cat["time"] < pd.Timestamp("2015-01-01", tz="UTC")]
    e._background(train)
    e.params = {"mu": 0.15, "K": 0.02, "alpha": 1.5, "c": 0.01, "p": 1.1, "D": 2.0, "q": 1.7, "gamma": 0.8}
    e.b = 1.0
    grid = ForecastGrid(TURKEY, 0.5, np.array([3.5, 10.0]))
    n_cells = (grid.lat_edges.size - 1) * (grid.lon_edges.size - 1)
    d = score_year(e, cat, pd.Timestamp("2015-01-01", tz="UTC"), pd.Timestamp("2016-01-01", tz="UTC"),
                   grid, np.zeros(n_cells), 3.5)
    assert len(d) > 20
    assert np.allclose(d["time_ll_model"], d["time_ll_etas"]) and np.allclose(d["space_ll_model"], d["space_ll_etas"])
    s = summarize(d)
    assert s["toplam kazanç"] == 0


def test_longwalk_compensator_and_weighted_background():
    from deprem.forecast.base import ForecastGrid
    from deprem.forecast.etas import _omori_integral
    from deprem import longwalk as L
    e = L.WETAS(m0=4.0)
    tr = pd.DataFrame({"time": pd.to_datetime(["2000-01-01", "2001-01-01"], utc=True),
                       "latitude": [39.0, 40.0], "longitude": [35.0, 30.0], "mag": [4.5, 4.5],
                       "depth_km": 10.0, "weight": [1.0, 0.0]})
    e._background(tr)
    from deprem.forecast.etas import _xy_km
    x0, y0 = _xy_km([39.0], [35.0])
    x1, y1 = _xy_km([40.0], [30.0])
    # ağırlığı 0 olan olay arka plana katkı vermez
    assert e._u(x0, y0)[0] > 1e3 * e._u(x1, y1)[0]
    e.params = {"mu": 0.01, "K": 0.05, "alpha": 1.0, "c": 0.01, "p": 1.2, "D": 2.0, "q": 1.7, "gamma": 0.8}
    grid = ForecastGrid(TURKEY, 0.5, np.array([4.0, 10.0]))
    n = (grid.lat_edges.size - 1) * (grid.lon_edges.size - 1)
    U = L.background_share(e, grid)
    assert abs(U.sum() - 1) < 1e-9
    s0, s1 = pd.Timestamp("2000-01-01", tz="UTC"), pd.Timestamp("2001-01-01", tz="UTC")
    par = tr.iloc[[0]].assign(time=pd.Timestamp("2000-06-01", tz="UTC"))
    comp = L.compensator(e, par, s0, s1, U, grid, n)
    T = 366.0
    tp = (par["time"].iloc[0] - s0).days
    expect = 0.05 * np.exp(0.5) * _omori_integral(T - tp, 0.01, 1.2)
    assert abs(comp.sum() - (0.01 * T + expect)) < 1e-6


def test_walkforward_train_offset_default_is_compensator():
    from deprem.pipeline import Run

    class R:
        cfg = {}
    assert Run._train_offset_mode(R()) == "compensator"
    R.cfg = {"walkforward": {"train_offset": "forecast"}}
    assert Run._train_offset_mode(R()) == "forecast"


def test_etas_fixed_alpha_and_empirical_magnitudes():
    from deprem.forecast.etas import ETAS
    rng = np.random.default_rng(7)
    n = 300
    t = pd.Timestamp("2010-01-01", tz="UTC") + pd.to_timedelta(np.sort(rng.uniform(0, 4 * 365, n)), unit="D")
    cat = pd.DataFrame({"time": t, "latitude": rng.uniform(37, 41, n), "longitude": rng.uniform(28, 40, n),
                        "mag": 3.5 + rng.exponential(0.43, n), "depth_km": 10.0})
    e = ETAS(m0=3.5, mag_sampling="empirical").fit(cat, pd.Timestamp("2011-01-01", tz="UTC"),
                                                   pd.Timestamp("2014-01-01", tz="UTC"), fixed={"alpha": 1.3})
    assert abs(e.params["alpha"] - 1.3) < 1e-9
    m = e._sample_mag(1000, rng) + e.m0
    pool = set(np.round(cat.loc[cat["time"] >= pd.Timestamp("2011-01-01", tz="UTC"), "mag"], 9))
    assert set(np.round(m, 9)) <= pool
