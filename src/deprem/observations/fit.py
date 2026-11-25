"""Öncüllü log-doğrusal Poisson: literatür başlangıç etkisi + geçmiş depremlerle güncelleme.

λ_i = exp(offset_i + β0 + Σ_k β_k z_ik),   β_k ~ N(m_k, s_k²)

En büyük sonsal (MAP) tahmin ve Laplace yaklaşımıyla sonsal standart sapma
döner. Öncülü olmayan sütunlar için geniş, sıfır merkezli öncül kullanılır.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def fit(Z: np.ndarray, offset: np.ndarray, y: np.ndarray, prior_mean: np.ndarray,
        prior_sd: np.ndarray, weights: np.ndarray | None = None):
    X = np.column_stack([np.ones(len(Z)), Z])
    m = np.r_[0.0, prior_mean]
    prec = np.r_[1e-4, 1 / np.asarray(prior_sd) ** 2]
    w = np.ones(len(y)) if weights is None else weights

    def f(beta):
        eta = np.clip(offset + X @ beta, -30, 30)
        lam = np.exp(eta)
        nll = (w * (lam - y * eta)).sum() + 0.5 * (prec * (beta - m) ** 2).sum()
        grad = X.T @ (w * (lam - y)) + prec * (beta - m)
        return nll, grad

    res = minimize(f, m.copy(), jac=True, method="L-BFGS-B")
    beta = res.x
    lam = np.exp(np.clip(offset + X @ beta, -30, 30))
    H = (X * (w * lam)[:, None]).T @ X + np.diag(prec)
    try:
        cov = np.linalg.inv(H)
        sd = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        sd = np.full(len(beta), np.nan)
    return beta, sd, res


def predict(Z, offset, beta):
    X = np.column_stack([np.ones(len(Z)), Z])
    return np.exp(np.clip(offset + X @ beta, -30, 30))


def standardize(train: np.ndarray, *others):
    mu = np.nanmean(train, axis=0)
    sd = np.nanstd(train, axis=0)
    sd = np.where((sd > 0) & np.isfinite(sd), sd, 1.0)
    mu = np.where(np.isfinite(mu), mu, 0.0)

    def z(a):
        return np.nan_to_num(np.clip((a - mu) / sd, -5, 5), nan=0.0)
    return [z(train)] + [z(o) for o in others]
