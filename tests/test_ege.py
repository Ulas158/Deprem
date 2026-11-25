import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import ege_canli  # noqa: E402
import ege_onkayit  # noqa: E402


def test_ig_zero_when_equal():
    n = np.array([3, 0, 5.0])
    lam = np.array([2.0, 1.0, 4.0])
    assert abs(ege_onkayit.ig(n, lam, lam)) < 1e-12


def test_ig_positive_when_model_puts_rate_where_events_are():
    n = np.array([10.0, 0.0])
    le = np.array([5.0, 5.0])
    lm = np.array([9.0, 1.0])
    assert ege_onkayit.ig(n, lm, le) > 0
    lo, hi, p = ege_onkayit.cell_bootstrap(np.r_[n, n], np.r_[lm, lm], np.r_[le, le], n_boot=500)
    assert lo <= ege_onkayit.ig(n, lm, le) <= hi


def test_verdict_rules():
    assert ege_onkayit.verdict(-0.01, 0.9, 100) == "kaldı"
    assert ege_onkayit.verdict(0.05, 0.2, 100) == "kararsız"
    assert ege_onkayit.verdict(0.05, 0.01, 100) == "geçti (zayıf)"
    assert ege_onkayit.verdict(0.2, 0.01, 100) == "geçti (güçlü)"
    assert ege_onkayit.verdict(0.2, 0.01, 10).startswith("kararsız")


def test_cell_of_matches_grid_order():
    # Türkiye ızgarası: 14 enlem × 39 boylam satırı, hücre = i * 39 + j
    assert ege_canli.cell_of(35.6, 25.6)[()] == 0
    assert ege_canli.cell_of(36.1, 26.1)[()] == 1 * 39 + 1
    assert ege_canli.cell_of(50.0, 30.0)[()] == -1


def test_daily_etas_background_only():
    U = np.full(14 * 39, 1 / (14 * 39))
    etas = {"params": {"mu": 2.0, "K": 0.0, "alpha": 1.0, "c": 0.01, "p": 1.1, "D": 2.0, "q": 1.7,
                       "gamma": 0.5}, "m0": 3.5, "b": 1.0, "dt_max_days": 3650}
    cat = pd.DataFrame({"time": pd.to_datetime(["2026-11-20"], utc=True), "latitude": [38.0],
                        "longitude": [27.0], "mag": [5.0]})
    d0 = pd.Timestamp("2026-11-26", tz="UTC")
    lam = ege_canli.daily_etas(cat, d0, d0 + pd.Timedelta(days=1), d0, etas, U, seed=0)
    assert np.isclose(lam.sum(), 2.0)
