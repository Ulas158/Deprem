"""Gutenberg-Richter istatistikleri.

log10 N(>=M) = a - b M

* Tamamlılık büyüklüğü Mc: maksimum eğrilik (MAXC, Wiemer & Wyss 2000)
  ve iyi-uyum testi (GFT, %90/%95).
* b-değeri: Aki (1965) / Utsu (1966) en çok olabilirlik kestirimi,
  büyüklük gruplama düzeltmesiyle.
* Belirsizlik: Shi & Bolt (1982) ve önyükleme (bootstrap).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

LOG10E = np.log10(np.e)


def bin_mags(m, dm: float = 0.1):
    return np.round(np.asarray(m, dtype=float) / dm) * dm


def fmd(m, dm: float = 0.1):
    """Frekans-büyüklük dağılımı: (bin merkezleri, artımlı sayı, birikimli sayı)."""
    mb = bin_mags(m, dm)
    if mb.size == 0:
        return np.array([]), np.array([]), np.array([])
    bins = np.round(np.arange(mb.min(), mb.max() + dm / 2, dm), 6)
    inc = np.array([(np.isclose(mb, x)).sum() for x in bins])
    cum = inc[::-1].cumsum()[::-1]
    return bins, inc, cum


def mc_maxc(m, dm: float = 0.1, correction: float = 0.2) -> float:
    """Maksimum eğrilik; Woessner & Wiemer (2005) +0.2 düzeltmesiyle."""
    bins, inc, _ = fmd(m, dm)
    if bins.size == 0:
        return np.nan
    return float(bins[np.argmax(inc)] + correction)


def b_value_mle(m, mc: float, dm: float = 0.1):
    """Aki-Utsu b-değeri ve Shi-Bolt standart hatası. Dönen: (b, sigma_b, n)."""
    mb = bin_mags(m, dm)
    x = mb[mb >= mc - dm / 2]
    n = x.size
    if n < 2:
        return np.nan, np.nan, n
    mean = x.mean()
    b = LOG10E / (mean - (mc - dm / 2))
    sigma = 2.3 * b ** 2 * np.sqrt(((x - mean) ** 2).sum() / (n * (n - 1)))
    return float(b), float(sigma), int(n)


def a_value(n_above_mc: int, b: float, mc: float, years: float = 1.0) -> float:
    """Yıllık a-değeri: log10 N(>=0) / yıl."""
    return float(np.log10(n_above_mc / years) + b * mc)


def mc_gft(m, dm: float = 0.1, level: float = 90.0, min_n: int = 50):
    """İyi-uyum testi. Kalıntı R = 100 - %fark; R >= level olan ilk Mc.

    Hiçbir Mc seviyeyi sağlamazsa en iyi R'li Mc döner. Dönen: (Mc, R).
    """
    bins, _, _ = fmd(m, dm)
    mb = bin_mags(m, dm)
    best = (np.nan, -np.inf)
    for mc in bins:
        x = mb[mb >= mc - dm / 2]
        if x.size < min_n:
            break
        b, _, n = b_value_mle(x, mc, dm)
        a = np.log10(n) + b * mc
        bb, inc, _ = fmd(x, dm)
        cum_obs = inc[::-1].cumsum()[::-1]
        cum_syn = 10 ** (a - b * bb)  # birikimli sentetik GR
        r = 100 - 100 * np.abs(cum_obs - cum_syn).sum() / cum_obs.sum()
        if r > best[1]:
            best = (float(mc), float(r))
        if r >= level:
            return float(mc), float(r)
    return best


@dataclass
class GRFit:
    mc: float
    b: float
    sigma_b: float
    b_boot_std: float
    n: int
    a_annual: float
    years: float

    def annual_rate(self, m: float) -> float:
        """Yılda M>=m olay sayısı (GR ekstrapolasyonu; yalnızca betimsel)."""
        return 10 ** (self.a_annual - self.b * m)


def fit_gr(m, years: float, dm: float = 0.1, mc: float | None = None,
           n_boot: int = 200, seed: int = 0) -> GRFit:
    m = np.asarray(m, dtype=float)
    if mc is None:
        mc = mc_maxc(m, dm)
    b, sigma, n = b_value_mle(m, mc, dm)
    rng = np.random.default_rng(seed)
    boots = []
    if n >= 10:
        for _ in range(n_boot):
            s = rng.choice(m, size=m.size, replace=True)
            boots.append(b_value_mle(s, mc, dm)[0])
    return GRFit(mc=float(mc), b=b, sigma_b=sigma,
                 b_boot_std=float(np.nanstd(boots)) if boots else np.nan,
                 n=n, a_annual=a_value(n, b, mc, years) if n else np.nan, years=years)


def synthetic_gr(n: int, b: float, mmin: float, seed: int = 0, dm: float = 0.0):
    """Test için GR dağılımlı rastgele büyüklükler (sürekli üstel)."""
    rng = np.random.default_rng(seed)
    beta = b * np.log(10)
    # Gruplanacaksa sürekli dağılım kutu alt kenarından (mmin - dm/2) başlar
    m = mmin - dm / 2 + rng.exponential(1 / beta, size=n)
    return bin_mags(m, dm) if dm else m
