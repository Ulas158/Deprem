import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from deprem_tarihsel.build import SCHEMA_COLUMNS, build  # noqa: E402
from deprem_tarihsel.sources import ms_to_mw, parse_intensity  # noqa: E402
from deprem_tarihsel.timeutil import astronomical_year, summarize  # noqa: E402


@pytest.fixture(scope="module")
def catalog(tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    meta = build(out_dir=out)
    allr = pd.read_csv(out / "tarihsel_tum_kayitlar.csv")
    pref = pd.read_csv(out / "tarihsel_katalog.csv")
    return meta, allr, pref


def test_time_summary():
    s = summarize(1766, 5, 22)
    assert s["time"] == "1766-05-22T12:00:00Z"
    assert s["time_precision"] == "day"
    assert s["time_err_s"] == pytest.approx(86400 / 12**0.5, rel=1e-3)
    y = summarize(1509)
    assert y["time_precision"] == "year" and y["time_min"] == "1509-01-01T00:00:00Z"
    assert summarize(1855, 2, 28, 1, 30)["time_precision"] == "minute"


def test_bc_years():
    assert astronomical_year(-282, historical_bc=True) == -281
    s = summarize(-281)
    assert s["time"].startswith("-0281-")
    assert pd.to_datetime(pd.Series([s["time"]]), utc=True).iloc[0].year == -281


def test_conversions():
    assert ms_to_mw(7.0) == pytest.approx(7.01)
    assert ms_to_mw(5.0) == pytest.approx(5.42)
    assert parse_intensity("8-9") == 8.5
    assert parse_intensity("10") == 10.0
    assert parse_intensity("F") is None
    assert parse_intensity(None) is None


def test_schema_and_ranges(catalog):
    _, allr, pref = catalog
    assert list(allr.columns) == SCHEMA_COLUMNS
    assert (allr["source"] == "historical").all()
    assert allr["latitude"].between(34, 44).all() and allr["longitude"].between(24, 46).all()
    assert allr["decimal_year"].max() < 1900
    for c in ["w_mag", "w_loc", "w_time", "w_unc", "training_weight"]:
        assert allr[c].between(0, 1).all(), c
    # Tercih edilmeyen yinelemeler eğitimde sayılmaz.
    assert (allr.loc[~allr["is_preferred"], "training_weight"] == 0).all()
    # Tarihler pandas ile UTC olarak okunabilmeli (aletsel katalogla birleşme için).
    # MÖ yıllar ("-0281-...") nedeniyle format="ISO8601" gerekir.
    t = pd.to_datetime(allr["time"], utc=True, format="ISO8601")
    assert t.notna().all()
    lo = pd.to_datetime(allr["time_min"], utc=True, format="ISO8601")
    hi = pd.to_datetime(allr["time_max"], utc=True, format="ISO8601")
    assert ((lo <= t) & (t < hi)).all()


def test_every_magnitude_is_flagged_estimated(catalog):
    _, allr, _ = catalog
    assert (allr["mag_is_estimated"] == allr["mag"].notna()).all()
    assert allr.loc[allr["mag"].notna(), "mag_err_used"].gt(0).all()


def test_duplicate_groups(catalog):
    _, allr, pref = catalog
    g = allr.groupby("duplicate_group")
    assert (g["agency"].nunique() == g.size()).all(), "bir grupta aynı kataloğun iki kaydı olamaz"
    assert (g["is_preferred"].sum() == 1).all()
    assert len(pref) == allr["duplicate_group"].nunique()
    assert pref["event_id"].is_unique


@pytest.mark.parametrize(
    "date,lat,lon,mmin",
    [
        ("1509-09-10", 41.0, 28.8, 6.8),   # İstanbul "Küçük Kıyamet"
        ("1766-05-22", 40.8, 29.4, 6.5),   # İstanbul / Doğu Marmara
        ("1766-08-05", 40.7, 27.1, 6.8),   # Ganos / Batı Marmara
        ("1894-07-10", 40.9, 29.0, 6.5),   # İstanbul
        ("1855-02-28", 40.2, 29.1, 6.5),   # Bursa
        ("1872-04-03", 36.4, 36.4, 6.8),   # Amik Gölü / Antakya
        ("1822-08-13", 36.7, 36.5, 7.0),   # Antakya-Halep
        ("1840-07-02", 39.6, 44.0, 7.0),   # Ağrı Dağı
        ("1668-08-17", 40.7, 36.0, 6.5),   # Kuzey Anadolu (Amasya)
    ],
)
def test_landmark_events_present(catalog, date, lat, lon, mmin):
    _, _, pref = catalog
    day = pref["time"].str[:10] == date
    near = (pref["latitude"] - lat).abs().lt(1.0) & (pref["longitude"] - lon).abs().lt(1.0)
    hit = pref[day & near]
    assert len(hit) == 1, f"{date} bulunamadı veya birden fazla"
    assert hit["group_mag_max"].iloc[0] >= mmin


def test_ancient_events_kept(catalog):
    _, _, pref = catalog
    assert (pref["decimal_year"] < 1000).sum() > 20
    assert (pref["decimal_year"] < 0).any()


def test_completeness_table(catalog):
    meta, _, _ = catalog
    starts = [meta["completeness"][k]["start_year"] for k in sorted(meta["completeness"], key=float)]
    assert all(np.diff(starts) <= 0), "büyük büyüklükler için tamlık daha erken başlamalı"
