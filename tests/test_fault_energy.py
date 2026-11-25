import numpy as np
import pandas as pd

from deprem.fault_energy import components, deficit_matrix, length_to_mw, release_matrix, segment_graph


def _faults():
    # üç segment uç uca (0-1-2), dördüncüsü uzakta
    coords = [[(30.0, 40.0), (30.2, 40.0)], [(30.2, 40.0), (30.4, 40.0)], [(30.4, 40.0), (30.6, 40.0)],
              [(35.0, 38.0), (35.2, 38.0)]]
    f = pd.DataFrame({"coords": coords, "length_km": [17.0] * 4, "width_km": [15.0] * 4,
                      "area_km2": [255.0] * 4, "slip_rate_mm_yr": [10.0, 10.0, 10.0, 1.0]})
    return f


def test_graph_and_components():
    f = _faults()
    nb = segment_graph(f)
    assert nb[0] == [1] and sorted(nb[1]) == [0, 2] and nb[3] == []
    L = f["length_km"].to_numpy()
    c = components(np.array([True, True, False, True]), nb, L)
    assert np.allclose(c, [34, 34, 0, 17])
    c = components(np.array([True, True, True, False]), nb, L)
    assert np.allclose(c[:3], 51)
    assert length_to_mw(100) > length_to_mw(10)


def test_release_conserves_moment_and_no_lookahead():
    f = _faults()
    ev = pd.DataFrame({"latitude": [40.0, 38.0], "longitude": [30.3, 35.1], "mw": [6.0, 4.0],
                       "t": [1000.0, 2000.0]})
    I, K, V = release_matrix(ev, f, f["width_km"].to_numpy())
    m0 = 10 ** (1.5 * ev["mw"].to_numpy() + 9.1)
    assert np.isclose(V[I == 0].sum(), m0[0], rtol=1e-6)
    assert np.isclose(V[I == 1].sum(), m0[1], rtol=1e-6)
    tau = np.array([0.0, 999.0, 1000.0, 1500.0, 3000.0])
    D = deficit_matrix(ev, f, tau, I, K, V)
    # adım 1 (999. gün) olaydan önce: boşalım yok
    load = 3e10 * 255e6 * 10e-3
    assert np.isclose(D[1, 1], load * 999 / 365.25)
    # adım 2 olayla başlıyor: boşalım düşülmüş
    assert D[2, 1] < load * 1000 / 365.25
