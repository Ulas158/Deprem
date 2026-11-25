import numpy as np
import pandas as pd

from deprem import regional as rg
from deprem import search as sr
from deprem.features.table import make_cells
from deprem.region import TURKEY


def _toy(n_cells=200, years=range(2000, 2010), seed=0):
    """İki bölge: özellik 'b_10y' bölge 0'da depremi artırır, bölge 1'de azaltır."""
    rng = np.random.default_rng(seed)
    reg = np.r_[np.zeros(n_cells // 2, int), np.ones(n_cells - n_cells // 2, int)]
    rows, off = [], {}
    for Y in years:
        f = rng.normal(size=n_cells)
        eff = np.where(reg == 0, 0.8, -0.8) * f
        lam = np.exp(np.log(2.0) + eff)
        y = rng.poisson(lam)
        rows.append(pd.DataFrame({"year": Y, "cell_id": np.arange(n_cells), "b_10y": f,
                                  "y_count_m3": y, "y_w": y.astype(float)}))
        off[Y] = np.full(n_cells, np.log(2.0))
    tab = pd.concat(rows, ignore_index=True)
    ys = list(years)
    d = sr.SearchData(tab=tab, off_test=off, off_train={(Y, yp): off[yp] for Y in ys for yp in ys if yp < Y},
                      b={Y: 1.0 for Y in ys}, mags=[3.0], train_mag=3.0, feature_cols=["b_10y"],
                      dcff_by_mu={0.4: []}, years=ys, score_mags=(3.0,))
    return d, reg


def test_cell_regions_cover_known_places():
    cells = make_cells(TURKEY, 0.5)
    reg = rg.cell_regions(cells)
    assert len(reg) == len(cells)
    at = lambda lat, lon: reg[((cells.lat_c - lat).abs() + (cells.lon_c - lon).abs()).idxmin()]
    assert rg.REGION_NAMES[at(40.75, 29.25)].startswith("Marmara")      # İzmit
    assert rg.REGION_NAMES[at(37.25, 37.25)].startswith("Doğu Anadolu Fay")  # Kahramanmaraş
    assert rg.REGION_NAMES[at(38.25, 27.25)].startswith("Ege")          # İzmir
    assert rg.REGION_NAMES[at(39.25, 32.75)] == rg.OTHER                # Ankara


def test_partial_pooling_limits():
    d, reg = _toy()
    cfg = dict(sr.default_config(["b-değeri"]), uncertainty_weights=True)
    years = [2007, 2008, 2009]
    glob = sr.evaluate(d, cfg, years)["score"]
    tiny = rg.evaluate(d, cfg, {"tau": 1e-3, "tau0": 1e-3, "tau_r": {}, "vary": ["b-değeri"]}, years, reg)["score"]
    big = rg.evaluate(d, cfg, {"tau": 5.0, "tau0": 1.0, "tau_r": {}, "vary": ["b-değeri"]}, years, reg)["score"]
    assert abs(tiny - glob) < 0.02           # τ → 0: tek ortak ayar
    assert big > glob + 0.1                  # gerçek bölge farkı varsa bölgesel model kazanır
    fake = rg.fake_regions(reg, 0)
    big_fake = rg.evaluate(d, cfg, {"tau": 5.0, "tau0": 1.0, "tau_r": {}, "vary": ["b-değeri"]}, years, fake)["score"]
    assert big_fake < big                    # sahte bölgeler gerçek farkı yakalayamaz
