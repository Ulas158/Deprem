import numpy as np
import pandas as pd
import pytest

from deprem.catalog.clean import decluster_gardner_knopoff, merge_sources
from deprem.catalog.schema import normalize
from deprem.catalog.sources import parse_afad, parse_emsc, parse_usgs
from deprem.faults.gem import parse_triple
from deprem.faults.moment import add_geometry, m0_to_mw, mw_to_m0, zone_budget
from deprem.features.table import build_table, time_split
from deprem.forecast.base import ForecastGrid, SmoothedSeismicity, UniformPoisson
from deprem.forecast.tests import information_gain, l_test, n_test
from deprem.geo import haversine_km, point_to_polyline_km
from deprem.region import Region
from deprem.stats.gr import b_value_mle, fit_gr, mc_maxc, synthetic_gr

REG = Region("test", 38.0, 40.0, 30.0, 32.0)


def _cat(n=2000, b=1.0, seed=0):
    rng = np.random.default_rng(seed)
    t0 = pd.Timestamp("2000-01-01", tz="UTC")
    return normalize(pd.DataFrame({
        "time": t0 + pd.to_timedelta(np.sort(rng.uniform(0, 20 * 365, n)), unit="D"),
        "latitude": rng.uniform(38, 40, n), "longitude": rng.uniform(30, 32, n),
        "depth_km": 10.0, "mag": synthetic_gr(n, b, 2.0, seed=seed, dm=0.1),
        "mag_type": "ml", "source": "emsc", "agency": "X",
        "event_id": [str(i) for i in range(n)],
    }))


# ------------------------------------------------------------------ GR
@pytest.mark.parametrize("b_true", [0.8, 1.0, 1.2])
def test_b_value_recovers_truth(b_true):
    m = synthetic_gr(20000, b_true, 2.0, seed=1, dm=0.1)
    b, sigma, n = b_value_mle(m, 2.0)
    assert abs(b - b_true) < 3 * sigma + 0.02


def test_mc_maxc_on_incomplete_catalog():
    rng = np.random.default_rng(2)
    m = synthetic_gr(50000, 1.0, 0.5, seed=2)
    # 2.0 altını giderek kaybeden algılama
    keep = rng.uniform(size=m.size) < 1 / (1 + np.exp(-(m - 2.0) / 0.15))
    mc = mc_maxc(m[keep])
    assert 1.8 <= mc <= 2.5


def test_fit_gr_annual_rate():
    m = synthetic_gr(10000, 1.0, 2.0, seed=3, dm=0.1)
    fit = fit_gr(m, years=10, mc=2.0, n_boot=20)
    assert fit.annual_rate(2.0) == pytest.approx(1000, rel=0.05)
    assert fit.annual_rate(4.0) == pytest.approx(10, rel=0.4)


# ------------------------------------------------------------------ parsers
def test_parse_usgs_and_emsc():
    usgs = ("time,latitude,longitude,depth,mag,magType,nst,gap,dmin,rms,net,id,updated,place,"
            "type,horizontalError,depthError,magError,magNst,status,locationSource,magSource\n"
            "2024-01-27T03:17:35.860Z,40.5,28.8,10,3,ml,7,175,1.3,0.67,us,us1,"
            "2024-03-30T22:34:46.040Z,\"x\",earthquake,3,2,0.1,14,reviewed,us,us\n")
    emsc = ("#EventID|Time|Latitude|Longitude|Depth/km|Author|Catalog|Contributor|ContributorID|"
            "MagType|Magnitude|MagAuthor|EventLocationName\n"
            "e1|2024-01-27T03:17:36.0Z|40.52|28.81|9.3|KOERI|EMSC-RTS|KOERI|1|ml|2.9|KOERI|W TURKEY\n")
    a, b = parse_usgs(usgs), parse_emsc(emsc)
    assert len(a) == 1 and len(b) == 1
    assert a["time"].dt.tz is not None and b["mag"].iloc[0] == 2.9
    merged = merge_sources([a, b])
    assert len(merged) == 1 and merged["source"].iloc[0] == "emsc"
    assert merged["n_sources"].iloc[0] == 2


def test_parse_afad():
    js = ('[{"rms":"0.37","eventID":"734705","location":"x","latitude":"38.08","longitude":"37.88",'
          '"depth":"6.93","type":"ML","magnitude":"1.5","date":"2026-11-20T23:40:02"}]')
    df = parse_afad(js)
    assert df["mag"].iloc[0] == 1.5 and df["mag_type"].iloc[0] == "ml"
    assert str(df["time"].dt.tz) == "UTC"


# ------------------------------------------------------------------ declustering
def test_gk_removes_aftershocks():
    t0 = pd.Timestamp("2010-01-01", tz="UTC")
    rows = [{"time": t0, "latitude": 39.0, "longitude": 31.0, "mag": 6.5}]
    rows += [{"time": t0 + pd.Timedelta(days=d), "latitude": 39.05, "longitude": 31.05, "mag": 3.0}
             for d in (1, 5, 20)]
    rows += [{"time": t0 + pd.Timedelta(days=3000), "latitude": 39.0, "longitude": 31.0, "mag": 3.0}]
    df = normalize(pd.DataFrame(rows))
    ind = decluster_gardner_knopoff(df)
    assert ind.tolist() == [True, False, False, False, True]


# ------------------------------------------------------------------ geo / faults
def test_geo_helpers():
    assert haversine_km(0, 0, 0, 1) == pytest.approx(111.19, rel=1e-3)
    d = point_to_polyline_km([40.1], [30.5], [(30.0, 40.0), (31.0, 40.0)])
    assert d[0] == pytest.approx(11.12, rel=1e-2)


def test_parse_triple():
    assert parse_triple("(0.2,0.1,0.4)") == (0.2, 0.1, 0.4)
    s = parse_triple("(15.0,,)")
    assert s[0] == 15.0 and np.isnan(s[1])


def test_moment_budget():
    assert m0_to_mw(mw_to_m0(7.4)) == pytest.approx(7.4)
    faults = pd.DataFrame([{
        "fault_id": "F1", "slip_rate_mm_yr": 20.0, "slip_rate_min": 15.0, "slip_rate_max": 25.0,
        "dip": 90.0, "upper_depth_km": 0.0, "lower_depth_km": 15.0, "length_km": 100.0,
        "coords": [(30.5, 39.0), (31.5, 39.0)]}])
    f = add_geometry(faults)
    # 3e10 * 100e3*15e3 m2 * 0.02 m/yıl = 9e17 Nm/yıl
    assert f["m0_rate_Nm_yr"].iloc[0] == pytest.approx(9e17)
    ev = normalize(pd.DataFrame([{"time": "1950-06-01T00:00:00Z", "latitude": 39.0,
                                   "longitude": 31.0, "mag": 7.0}]))
    ev["mw_best"] = ev["mag"]
    zb = zone_budget(f, ev, {"z": REG}, pd.Timestamp("1900-01-01", tz="UTC"),
                     pd.Timestamp("2000-01-01", tz="UTC"))
    assert zb["release_ratio"].iloc[0] == pytest.approx(mw_to_m0(7.0) / (9e17 * 99.99), rel=1e-3)


# ------------------------------------------------------------------ forecasts
def test_forecast_tests_and_baselines():
    cat = _cat()
    grid = ForecastGrid(REG, 0.5, np.array([2.0, 3.0, 4.0, 9.0]))
    fit_end = pd.Timestamp("2015-01-01", tz="UTC")
    test_end = pd.Timestamp("2020-01-01", tz="UTC")
    u = UniformPoisson(2.0).fit(cat, cat["time"].min(), fit_end)
    s = SmoothedSeismicity(2.0).fit(cat, cat["time"].min(), fit_end)
    obs = grid.count(cat[(cat["time"] >= fit_end) & (cat["time"] < test_end)])
    eu = u.expected(grid, fit_end, test_end)
    assert eu.shape == grid.shape
    assert eu.sum() == pytest.approx(obs.sum(), rel=0.15)
    assert n_test(eu, obs).quantile > 0.01
    assert 0 <= l_test(eu, obs, n_sim=100).quantile <= 1
    ig = information_gain(s.expected(grid, fit_end, test_end), eu, obs)
    assert np.isfinite(ig.statistic)


def test_feature_table_has_no_leakage():
    cat = _cat()
    cat["mw_best"] = cat["mag"]
    starts = [pd.Timestamp(f"{y}-01-01", tz="UTC") for y in range(2005, 2020)]
    tab = build_table(cat, REG, 0.5, starts, pd.DateOffset(years=1), target_mags=[3.0, 3.5], mc=2.0)
    assert {"y_any_m3", "y_count_m3.5"} <= set(tab.columns)
    assert len(tab) == 16 * len(starts)
    # Özellik: geçmiş 1 yıl sayımı; bir sonraki pencerenin etiketiyle aynı veri olmamalı
    row = tab[(tab["window_start"] == starts[3]) & (tab["cell_id"] == 0)].iloc[0]
    past = cat[(cat["time"] < starts[3]) & (cat["time"] >= starts[2])
               & (cat["latitude"] < 38.5) & (cat["longitude"] < 30.5)]
    assert row["n_mc_1y"] == len(past)
    train, test = time_split(tab, pd.Timestamp("2015-01-01", tz="UTC"))
    assert train["window_start"].max() < test["window_start"].min()


def test_reliability_and_transfer():
    from deprem.features.evaluate import reliability_curve, transfer_report
    rng = np.random.default_rng(5)
    lam3 = rng.gamma(0.5, 0.4, size=20000)
    y3 = rng.poisson(lam3)
    y4 = rng.binomial(y3, 0.1)  # b=1: M4 oranı M3'ün %10'u
    tab = pd.DataFrame({"p3": 1 - np.exp(-lam3), "y_any_m3": (y3 > 0).astype(int),
                        "y_count_m3": y3, "y_any_m4": (y4 > 0).astype(int), "y_count_m4": y4})
    rc = reliability_curve(tab["p3"], tab["y_any_m3"])
    assert np.allclose(rc["p_mean"], rc["obs_freq"], atol=0.03)
    rep = transfer_report(tab, "p3", 3.0, [3.0, 4.0], b=1.0)
    assert rep.loc[1, "n_expected"] == pytest.approx(rep.loc[1, "n_events"], rel=0.1)


def test_uncertainty_weights():
    from deprem.catalog.weights import compute_weights
    df = normalize(pd.DataFrame([
        {"time": "1912-08-09T01:29:00Z", "latitude": 40.6, "longitude": 27.2, "mag": 7.3},
        {"time": "2020-01-01T00:00:00Z", "latitude": 39.0, "longitude": 31.0, "mag": 2.1},
        {"time": "2020-01-01T00:00:00Z", "latitude": 39.0, "longitude": 31.0, "mag": 2.1,
         "horiz_err_km": 40.0, "mag_err": 0.5},
        {"time": "1800-01-01T00:00:00Z", "latitude": 39.0, "longitude": 31.0, "mag": 7.0,
         "time_err_s": 86400 * 30.0},
    ]))
    # normalize zamana göre sıralar: [1800, 1912, 2020 kesin, 2020 belirsiz]
    w = compute_weights(df).to_numpy()
    # modern küçük olay > eski büyük olay > tarihsel; belirsiz modern < kesin modern
    assert w[2] > w[1] > w[0] and w[2] > w[3]
    assert w[2] > 0.5 and (w > 0).all()


def test_latent_fit_recovers_mix_and_validates():
    from deprem.forecast.latent import LATENT, fit_latent
    from deprem.forecast.tests import poisson_loglik
    assert len(LATENT) >= 10  # deprem.external yeni gizli parametreler ekler
    rng = np.random.default_rng(0)
    n = 400
    a = rng.dirichlet(np.ones(n))  # "fay moment" uzamsal yoğunluğu
    s = rng.dirichlet(np.ones(n))  # "düzleştirilmiş sismisite" yoğunluğu
    w_true, total = 0.7, 3000

    def lam(p):
        w = p["background_mix"]
        return total * (w * a + (1 - w) * s)

    true = total * (w_true * a + (1 - w_true) * s)
    y_train, y_val = rng.poisson(true), rng.poisson(true)
    res = fit_latent(lambda p: poisson_loglik(lam(p), y_train), ["background_mix"],
                     prior_weight=0.0, val_loglik_fn=lambda p: poisson_loglik(lam(p), y_val))
    assert abs(res.values["background_mix"] - w_true) < 0.1
    assert res.val_loglik > poisson_loglik(total * s, y_val)


def test_fault_mesh_geometry():
    from deprem.faults.mesh import mesh_fault
    row = pd.Series({"fault_id": "F", "coords": [(30.0, 39.0), (30.0, 39.0 + 10 / 111.195)],
                     "dip": 45.0, "upper_depth_km": 0.0, "lower_depth_km": 10.0, "rake": 0.0,
                     "slip_rate_mm_yr": 10.0, "slip_rate_min": 5.0, "slip_rate_max": 15.0})
    el = mesh_fault(row, 1.0, 1.0)
    width = 10 / np.sin(np.radians(45))
    assert len(el) == 10 * int(np.ceil(width))
    assert el["area_km2"].sum() == pytest.approx(10 * width, rel=1e-3)
    assert el["strike"].iloc[0] == pytest.approx(0.0, abs=1e-6)
    # kuzeye doğrultu, sağ el kuralı → doğuya dalar: derin elemanlar daha doğuda
    deep = el[el["j_dip"] == el["j_dip"].max()]["lon"].mean()
    shallow = el[el["j_dip"] == 0]["lon"].mean()
    assert deep > shallow
    assert el["depth_km"].max() < 10.0 and el["depth_rel"].between(0, 1).all()


def test_element_seismicity():
    from deprem.faults.parameters import element_seismicity, ensure_columns, REGISTRY
    el = pd.DataFrame({"lat": [39.0, 40.0], "lon": [30.0, 30.0], "depth_km": [10.0, 10.0],
                       "m0_rate_Nm_yr": [1e15, 1e15]})
    t_ref = pd.Timestamp("2020-01-01", tz="UTC")
    ev = normalize(pd.DataFrame({
        "time": [t_ref - pd.Timedelta(days=d) for d in (100, 200, 3000)],
        "latitude": [39.01, 39.0, 39.0], "longitude": [30.0, 30.01, 30.0],
        "depth_km": [9.0, 11.0, 10.0], "mag": [3.0, 4.2, 5.1]}))
    out = element_seismicity(el, ev, t_ref, mc=2.0, radius_km=5.0)
    assert out["seis_n_all"].tolist() == [3, 0]
    assert out["seis_n_1y"].iloc[0] == 2
    assert out["seis_yrs_since_m5"].iloc[0] == pytest.approx(3000 / 365.25, rel=1e-3)
    assert np.isnan(out["seis_max_mag"].iloc[1])
    full = ensure_columns(el)
    assert set(REGISTRY) <= set(full.columns)


def test_rupture_extent_resets_neighbouring_segments():
    from deprem.faults.moment import fault_summary
    # doğu-batı uzanan üç ardışık 50 km'lik segment; M7.5 batı ucunda
    k = 50 / (111.195 * np.cos(np.radians(40)))
    segs = [[(30 + i * k, 40.0), (30 + (i + 1) * k, 40.0)] for i in range(3)]
    faults = add_geometry(pd.DataFrame({
        "fault_id": ["A", "B", "C"], "coords": segs, "slip_rate_mm_yr": 20.0,
        "slip_rate_min": 10.0, "slip_rate_max": 25.0, "dip": 90.0,
        "upper_depth_km": 0.0, "lower_depth_km": 15.0,
        "length_km": 50.0}))
    ev = normalize(pd.DataFrame([{"time": "1999-08-17T00:01:39Z", "latitude": 40.0,
                                   "longitude": 30.1, "mag": 7.5}]))
    ev["mw_best"] = ev["mag"]
    fs = fault_summary(faults, ev, pd.Timestamp("1900-01-01", tz="UTC"),
                       pd.Timestamp("2026-01-01", tz="UTC"))
    hit = fs.set_index("fault_id")["last_big_event"].notna()
    # L(7.5) ≈ 90 km → yarısı ≈ 45 km; episantr A içinde, B başına ~42 km: A ve B kırılmış, C (≥100 km) değil
    assert hit.to_dict() == {"A": True, "B": True, "C": False}


def test_coulomb_classic_pattern():
    pytest.importorskip("cutde")
    from deprem.faults.coulomb import coulomb_on_receivers, stress_from_sources
    # sağ yanal düşey K-G fay, 30 km × 15 km, 1 m kayma
    src = pd.DataFrame([{"x_m": 0, "y_m": 0, "z_m": -7500, "strike": 0, "dip": 90, "rake": 180,
                         "L_m": 30000, "W_m": 15000, "slip_m": 1.0}])
    obs = np.array([[0, 25000, -7500], [8000, 0, -7500], [0, -25000, -7500]], float)
    st = stress_from_sources(obs, src)
    dcff, _, _ = coulomb_on_receivers(st, np.zeros(3), np.full(3, 90.0), np.full(3, 180.0), 0.4)
    bar = dcff / 1e5
    # uçların ötesi yüklenir (öne çekme), yanlar gölgede kalır
    assert bar[0] > 1 and bar[2] > 1 and bar[1] < -1
    assert bar[0] == pytest.approx(bar[2], rel=1e-6)


def test_element_budget_releases_only_ruptured_part():
    from deprem.faults.budget import element_budget
    from deprem.faults.mesh import mesh_fault
    k = 1 / (111.195 * np.cos(np.radians(39.0)))
    row = pd.Series({"fault_id": "F", "coords": [(35.0, 39.0), (35.0 + 200 * k, 39.0)], "dip": 90.0,
                     "upper_depth_km": 0.0, "lower_depth_km": 15.0, "rake": 180.0,
                     "slip_rate_mm_yr": 20.0, "slip_rate_min": 10.0, "slip_rate_max": 30.0})
    el = mesh_fault(row, 1.0, 1.0)
    ev = normalize(pd.DataFrame([{"time": "1950-01-01T00:00:00Z", "latitude": 39.0,
                                   "longitude": 35.0 + 20 * k, "depth_km": 8.0, "mag": 7.0,
                                   "horiz_err_km": 2.0, "depth_err_km": 2.0}]))
    ev["mw_best"] = ev["mag"]
    res = element_budget(el, ev, pd.Timestamp("1900-01-01", tz="UTC"),
                         pd.Timestamp("2000-01-01", tz="UTC"))
    f = res["final"]
    # toplam boşalım M0'a eşit, birikim 100 yıl × hız
    assert f["m0_rel"].sum() == pytest.approx(10 ** (1.5 * 7 + 9.1), rel=1e-6)
    near = el["along_km"] < 30
    far = el["along_km"] > 120
    assert f.loc[far, "m0_rel"].sum() == 0
    assert f.loc[near, "slip_deficit_m"].mean() < f.loc[far, "slip_deficit_m"].mean()
    assert f.loc[far, "slip_deficit_m"].mean() == pytest.approx(0.020 * 100, rel=0.01)
