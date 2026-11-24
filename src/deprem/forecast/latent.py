"""Veriden öğrenilecek gizli (serbest) parametreler ve en iyileme arayüzü.

Doğrudan ölçülemeyen ama fiziksel anlamı olan 10 parametre tanımlıdır.
Model bunları **zamana göre ayrılmış eğitim verisinde** Poisson
log-olabilirliği en büyütecek şekilde ayarlar; başarı yalnızca eğitimde
görülmemiş doğrulama döneminde ölçülür.

Aşırı uyumu önlemek için:
* her parametrenin fiziksel sınırları ve Gauss öncül dağılımı vardır
  (dönüştürülmüş uzayda; öncül cezası = düzenlileştirme),
* bir model yalnızca kullandığı parametreleri en iyiler (``uses``),
* doğrulama dönemi en iyilemeye hiç girmez.

Bir parametreyi "maksimize etmek" burada log-olabilirliği en büyütecek
değeri bulmak demektir; parametrenin kendisini büyütmek değil.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class LatentParam:
    name: str
    description: str
    init: float
    lower: float
    upper: float
    prior_mean: float
    prior_sd: float  # dönüştürülmüş uzayda (log parametrelerde log10 birimi)
    log: bool = False

    def to_u(self, v):
        return np.log10(v) if self.log else v

    def from_u(self, u):
        return 10 ** u if self.log else u


LATENT: dict[str, LatentParam] = {p.name: p for p in [
    LatentParam("mu_eff", "Efektif sürtünme katsayısı (Coulomb ΔCFF = Δτ + μ'Δσn)",
                0.4, 0.0, 0.9, 0.4, 0.2),
    LatentParam("a_sigma_mpa", "Hız-durum sürtünme Aσ (Dieterich 1994 sismisite yanıtı)",
                0.1, 0.005, 2.0, -1.0, 0.5, log=True),
    LatentParam("pore_pressure_ratio", "Gözenek basıncı / litostatik basınç (λ)",
                0.4, 0.0, 0.95, 0.4, 0.2),
    LatentParam("aftershock_duration_yr", "Artçı süresi t_a = Aσ / gerilme yükleme hızı",
                10.0, 0.1, 300.0, 1.0, 0.6, log=True),
    LatentParam("locking_scale", "Kayma hızına uygulanan kilitlenme/bağlaşım çarpanı",
                0.7, 0.05, 1.0, 0.7, 0.2),
    LatentParam("etas_K", "ETAS üretkenlik katsayısı", 0.02, 1e-4, 1.0, -1.7, 0.7, log=True),
    LatentParam("etas_alpha", "ETAS büyüklük ölçekli üretkenlik üsteli (α)",
                1.5, 0.5, 2.5, 1.5, 0.4),
    LatentParam("coulomb_gain", "ΔCFF'nin oran değişimine etkisi (0 = etkisiz)",
                1.0, 0.0, 10.0, 1.0, 1.0),
    LatentParam("background_mix", "Arka plan oranında fay-moment payı (1-w: düzleştirilmiş sismisite)",
                0.3, 0.0, 1.0, 0.3, 0.3),
    LatentParam("smoothing_km", "Düzleştirme çekirdeği bant genişliği",
                20.0, 2.0, 100.0, 1.3, 0.3, log=True),
]}


@dataclass
class FitResult:
    values: dict
    train_loglik: float
    prior_penalty: float
    val_loglik: float | None = None
    success: bool = True
    message: str = ""
    history: list = field(default_factory=list)


def fit_latent(loglik_fn: Callable[[dict], float], uses: list[str],
               fixed: dict | None = None, prior_weight: float = 1.0,
               val_loglik_fn: Callable[[dict], float] | None = None,
               method: str = "L-BFGS-B") -> FitResult:
    """``uses`` listesindeki gizli parametreleri en iyile.

    ``loglik_fn(params) -> float`` eğitim log-olabilirliğidir (büyük = iyi).
    Amaç: max loglik − prior_weight · Σ ½((u−μ)/σ)².
    """
    fixed = dict(fixed or {})
    ps = [LATENT[n] for n in uses]
    u0 = np.array([p.to_u(p.init) for p in ps])
    bounds = [(p.to_u(p.lower), p.to_u(p.upper)) for p in ps]
    hist = []

    def unpack(u):
        vals = {p.name: float(p.from_u(x)) for p, x in zip(ps, u)}
        vals.update(fixed)
        return vals

    def penalty(u):
        return prior_weight * sum(0.5 * ((x - p.prior_mean) / p.prior_sd) ** 2
                                  for p, x in zip(ps, u))

    def obj(u):
        ll = loglik_fn(unpack(u))
        val = -(ll - penalty(u))
        hist.append(float(-val))
        return val if np.isfinite(val) else 1e30

    res = minimize(obj, u0, method=method, bounds=bounds)
    vals = unpack(res.x)
    out = FitResult(vals, float(loglik_fn(vals)), float(penalty(res.x)),
                    success=bool(res.success), message=str(res.message), history=hist)
    if val_loglik_fn is not None:
        out.val_loglik = float(val_loglik_fn(vals))
    return out
