"""Olasılık tahminlerinin eşik eşik değerlendirilmesi.

Bir model her (hücre, pencere) için P(en az bir M>=m olay) üretir. Bu modül:

* Güvenilirlik (kalibrasyon) eğrisi: "%10 dediği hücrelerde gerçekten
  ~%10'unda deprem oldu mu?"
* Brier skoru ve referansa göre Brier beceri skoru.
* Olasılığı Poisson beklenen sayıya çevirip (λ = -ln(1-p)) CSEP
  log-olabilirlik testleriyle ETAS / düzleştirilmiş sismisite
  referanslarıyla aynı ölçekte karşılaştırma.

``transfer_report`` küçük eşikte (ör. M3) kalibre edilen olasılıkların GR
ölçeklemesiyle büyük eşiklere (M4, M5) taşındığında ne kadar tutarlı
kaldığını ölçer. Büyük eşiklerde olay sayısı az olduğundan belirsizlik
aralıkları geniştir; bu tabloda her eşiğin olay sayısı da raporlanır.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..forecast.tests import poisson_loglik


def reliability_curve(p: np.ndarray, y: np.ndarray, n_bins: int = 10,
                      strategy: str = "quantile") -> pd.DataFrame:
    p = np.asarray(p, float)
    y = np.asarray(y, float)
    if strategy == "quantile":
        edges = np.unique(np.quantile(p, np.linspace(0, 1, n_bins + 1)))
    else:
        edges = np.linspace(0, 1, n_bins + 1)
    idx = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, len(edges) - 2)
    rows = []
    for k in range(len(edges) - 1):
        sel = idx == k
        if not sel.any():
            continue
        n = sel.sum()
        obs = y[sel].mean()
        # Wilson %95 aralığı
        z = 1.96
        den = 1 + z ** 2 / n
        c = (obs + z ** 2 / (2 * n)) / den
        h = z * np.sqrt(obs * (1 - obs) / n + z ** 2 / (4 * n ** 2)) / den
        rows.append({"p_mean": p[sel].mean(), "obs_freq": obs, "n": int(n),
                     "lo95": c - h, "hi95": c + h})
    return pd.DataFrame(rows)


def brier(p, y) -> float:
    return float(np.mean((np.asarray(p, float) - np.asarray(y, float)) ** 2))


def brier_skill(p, p_ref, y) -> float:
    """1 = mükemmel, 0 = referans kadar, <0 = referanstan kötü."""
    return 1 - brier(p, y) / brier(p_ref, y)


def prob_to_rate(p):
    return -np.log1p(-np.clip(np.asarray(p, float), 0, 1 - 1e-12))


def rate_to_prob(lam):
    return -np.expm1(-np.asarray(lam, float))


def scale_rate_gr(lam, m_from: float, m_to: float, b: float):
    """GR ile λ(M>=m_from) -> λ(M>=m_to)."""
    return np.asarray(lam, float) * 10 ** (-b * (m_to - m_from))


def transfer_report(table: pd.DataFrame, p_col: str, m_cal: float, m_targets,
                    b: float, p_ref_col: str | None = None) -> pd.DataFrame:
    """M=m_cal için verilmiş olasılıkları her hedef eşiğe taşı ve puanla."""
    lam_cal = prob_to_rate(table[p_col])
    rows = []
    for m in m_targets:
        y_any = table[f"y_any_m{m:g}"].to_numpy()
        y_cnt = table[f"y_count_m{m:g}"].to_numpy()
        lam = scale_rate_gr(lam_cal, m_cal, m, b)
        p = rate_to_prob(lam)
        row = {"m": m, "n_cells_with_event": int(y_any.sum()), "n_events": int(y_cnt.sum()),
               "n_expected": float(lam.sum()), "brier": brier(p, y_any),
               "loglik_poisson": poisson_loglik(lam, y_cnt)}
        if p_ref_col is not None:
            p_ref = rate_to_prob(scale_rate_gr(prob_to_rate(table[p_ref_col]), m_cal, m, b))
            row["brier_skill_vs_ref"] = brier_skill(p, p_ref, y_any)
            row["loglik_gain_vs_ref"] = row["loglik_poisson"] - poisson_loglik(
                prob_to_rate(p_ref), y_cnt)
        rows.append(row)
    return pd.DataFrame(rows)
