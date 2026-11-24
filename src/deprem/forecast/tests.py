"""CSEP tarzı tutarlılık ve karşılaştırma testleri (Zechar vd. 2010).

* N-testi: toplam olay sayısı tahminle tutarlı mı?
* L-testi: gözlemin Poisson birleşik log-olabilirliği, modelin kendi
  benzetimlerine göre makul mü?
* Bilgi kazancı (Rhoades vd. 2011): iki model arasında olay başına
  log-olabilirlik farkı ve t-testi.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats
from scipy.special import gammaln


def poisson_loglik(expected: np.ndarray, observed: np.ndarray) -> float:
    lam = np.clip(expected, 1e-300, None)
    return float((-lam + observed * np.log(lam) - gammaln(observed + 1)).sum())


@dataclass
class TestResult:
    name: str
    statistic: float
    quantile: float  # tek taraflı olasılık(lar)
    detail: dict


def n_test(expected: np.ndarray, observed: np.ndarray) -> TestResult:
    n_fore = float(expected.sum())
    n_obs = int(observed.sum())
    d1 = 1 - stats.poisson.cdf(n_obs - 1, n_fore)  # P(N >= n_obs)
    d2 = stats.poisson.cdf(n_obs, n_fore)  # P(N <= n_obs)
    return TestResult("N", n_obs, min(d1, d2),
                      {"n_forecast": n_fore, "n_observed": n_obs, "delta1": d1, "delta2": d2})


def l_test(expected: np.ndarray, observed: np.ndarray, n_sim: int = 1000,
           seed: int = 0) -> TestResult:
    rng = np.random.default_rng(seed)
    ll_obs = poisson_loglik(expected, observed)
    sims = np.array([poisson_loglik(expected, rng.poisson(expected)) for _ in range(n_sim)])
    gamma = float((sims <= ll_obs).mean())
    return TestResult("L", ll_obs, gamma, {"ll_sim_mean": float(sims.mean())})


def information_gain(expected_a: np.ndarray, expected_b: np.ndarray,
                     observed: np.ndarray) -> TestResult:
    """A'nın B'ye göre olay başına bilgi kazancı (pozitif = A daha iyi)."""
    a = np.clip(expected_a, 1e-300, None).ravel()
    b = np.clip(expected_b, 1e-300, None).ravel()
    obs = observed.ravel()
    n = obs.sum()
    if n == 0:
        return TestResult("IG", np.nan, np.nan, {"n": 0})
    idx = np.repeat(np.arange(obs.size), obs.astype(int))
    per_event = np.log(a[idx]) - np.log(b[idx])
    ig = per_event.mean() - (a.sum() - b.sum()) / n
    if n > 1:
        t, p = stats.ttest_1samp(per_event - (a.sum() - b.sum()) / n, 0.0)
    else:
        t, p = np.nan, np.nan
    return TestResult("IG", float(ig), float(p), {"t": float(t), "n": int(n)})
