"""Sıralı (prequential) yürüyüşün geleceği görmediğini sınar.

Sentetik küçük bir veriyle: k. olaydan sonraki olayların yeri/büyüklüğü değiştirildiğinde
ilk k olayın puanları değişmemeli.
"""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("jax")

from deprem.faultwalk import VARIANTS, walk  # noqa: E402
from deprem.faultwalk_data import WalkData  # noqa: E402


def _synthetic(seed=0, perturb_after=None):
    rng = np.random.default_rng(seed)
    N, C = 24, 6
    t = np.sort(rng.uniform(0, 3650, N))
    m = 5.5 + rng.exponential(0.4, N)
    cell = rng.integers(0, C, N)
    xs = np.arange(C) * 50.0
    if perturb_after is not None:
        rng2 = np.random.default_rng(99)
        m[perturb_after:] = 5.5 + rng2.exponential(0.4, N - perturb_after)
        cell[perturb_after:] = rng2.integers(0, C, N - perturb_after)
    ex = xs[cell] + 5.0
    d2 = ((ex[:, None, None] - (xs[None, :, None] + np.array([-15, 0, 15])[None, None, :])) ** 2).repeat(3, -1)
    years = np.arange(10)
    tau = np.unique(np.concatenate([np.arange(11) * 365.0, t]))
    S = tau.size - 1
    step_year = np.clip((tau[:-1] // 365).astype(int), 0, 9)
    ev_bstep = np.searchsorted(tau, t)
    units = pd.DataFrame({"unit": ["a", "b"], "kind": ["fay", "fay dışı"], "system": ["X", "diğer: Y"], "sys": [0, 1]})
    uoc = np.array([0, 0, 0, 1, 1, 1])
    smap = np.full((10, C), 1.0 / C)
    cat = pd.DataFrame({"time": pd.to_datetime("2000-01-01", utc=True) + pd.to_timedelta(t, unit="D"),
                        "latitude": 39.0, "longitude": 35.0 + ex / 85.0, "mw": m})
    return WalkData(cat=cat, lat_c=np.zeros(C), lon_c=np.zeros(C), area=np.full(C, 2500.0), unit_of_cell=uoc,
                    units=units, systems=["X", "diğer: Y"], major=np.array([True, False]), ev_cell=cell,
                    ev_unit=uoc[cell], ev_sys=units["sys"].to_numpy()[uoc[cell]], t_ev=t, m_ev=m,
                    d2=d2.astype(np.float32), cfs=np.zeros((N, C), np.float32), fmap=np.full(C, 1.0 / C),
                    smap=smap, years=years, tau=tau, step_year=step_year, step_era=np.full(S, 3),
                    L=np.zeros((S, 2)), ev_step=ev_bstep - 1, ev_era=np.full(N, 3), ev_bstep=ev_bstep)


def test_no_future_leak():
    k = 12
    a, _ = walk(_synthetic(), VARIANTS["partial"], refit_every=3, maxiter=10, first_maxiter=20, keep_params_every=0)
    b, _ = walk(_synthetic(perturb_after=k), VARIANTS["partial"], refit_every=3, maxiter=10, first_maxiter=20,
                keep_params_every=0)
    np.testing.assert_allclose(a["score"].to_numpy()[:k], b["score"].to_numpy()[:k], rtol=1e-6, atol=1e-8)
    assert not np.allclose(a["score"].to_numpy()[k:], b["score"].to_numpy()[k:])
