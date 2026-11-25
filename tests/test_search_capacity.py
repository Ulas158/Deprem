import numpy as np
import pandas as pd

from deprem import search as sr
from deprem.faults import capacity as cp


def test_holm_monotone_and_bounded():
    adj = sr.holm([0.01, 0.04, 0.03, 0.5])
    assert np.all(adj <= 1) and np.all(adj >= [0.01, 0.04, 0.03, 0.5])
    assert adj[0] == 0.04  # 4 × 0.01


def test_column_groups():
    assert sr.column_group("fill_mean") == "doluluk oranı"
    assert sr.column_group("fillc_max") == "doluluk + Coulomb"
    assert sr.column_group("cascade_mw") == "zincirleme kırılma"
    assert sr.column_group("dcff_max_mu0.4") == "Coulomb gerilme aktarımı"
    assert sr.column_group("slip_deficit_mean__eksik") == "kayma açığı (enerji birikimi)"


def test_cascade_stops_at_barrier():
    # 0-1-2-3-4 zincir; 2 boş (bariyer)
    fill = np.array([1.2, 0.9, 0.0, 0.95, 0.95])
    nb = [[1], [0, 2], [1, 3], [2, 4], [3]]
    b = cp.cascade(fill, nb, 0, thr=1.0, alpha=0.5)
    assert b[0] and b[1] and not b[3] and not b[4]


def test_cascade_spreads_when_neighbours_full():
    fill = np.array([1.2, 0.9, 0.9, 0.9])
    nb = [[1], [0, 2], [1, 3], [2]]
    assert cp.cascade(fill, nb, 0, thr=1.0, alpha=0.5).all()


def test_calibration_partial_pooling():
    rng = np.random.default_rng(0)
    rows = []
    for reg, mu in [("a", -2.0), ("b", -3.0)]:
        for k in range(8):
            rows.append({"time": pd.Timestamp(f"{1960 + k}-01-01", tz="UTC"), "mw": 7.0, "bölge": reg,
                         "fay": f"{reg}{k % 2}", "kayma_hızı_ağırlıklı": 10.0,
                         "doluluk_fiziksel": float(np.exp(mu + 0.3 * rng.normal()))})
    c = cp.hierarchical_calibration(pd.DataFrame(rows), 2000)
    assert c["n"] == 16
    assert c["bölge:a"][0] > c["bölge:b"][0]
    assert -3.5 < c["bölge:b"][0] < -2.5


def test_shuffled_table_keeps_labels():
    tab = pd.DataFrame({"year": [1, 1, 1, 2, 2, 2], "cell_id": [0, 1, 2, 0, 1, 2], "f": [1., 2, 3, 4, 5, 6],
                        "y_count_m3": [0, 1, 2, 3, 4, 5]})
    d = sr.SearchData(tab=tab, off_test={}, off_train={}, b={}, mags=[3.0], train_mag=3.0,
                      feature_cols=["f"], dcff_by_mu={0.4: []})
    t = sr.shuffled_table(d, 0)
    assert (t["y_count_m3"] == tab["y_count_m3"]).all()
    assert sorted(t.loc[t.year == 1, "f"]) == [1, 2, 3]
