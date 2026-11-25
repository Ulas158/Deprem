"""Uzay-zaman ETAS referans modeli (Ogata 1998; Zhuang vd. 2002).

Koşullu yoğunluk (olay/gün/km²):

    λ(t, x, y) = μ u(x, y) + Σ_{t_j < t} κ(m_j) g(t - t_j) f(x - x_j, y - y_j | m_j)

    κ(m) = K exp(α (m - m0))                  üretkenlik
    g(t) = (c + t)^(-p)                       Omori-Utsu (normalize edilmemiş)
    f(r | m) = (q - 1) / (π σ(m)) (1 + r² / σ(m))^(-q),  σ(m) = D² exp(γ (m - m0))

``u(x, y)`` artçıları ayıklanmış eğitim olaylarından Gauss çekirdekle
düzleştirilmiş, bölgede integrali 1 olan arka plan yoğunluğudur; μ
günlük arka plan olay sayısıdır. Parametreler en çok olabilirlikle
(L-BFGS-B, log-dönüşümlü) bulunur. Uzaysal integral bölge sınırında
kaybolan kütle ihmal edilerek 1 alınır (yaygın yaklaşım).

Tahmin (``expected``): pencere başından önceki tüm olaylar bilinir;
pencere içindeki beklenen sayı = arka plan + bilinen olayların doğrudan
artçıları; pencere içinde doğacak olayların kendi artçıları (kademeli
tetiklenme) dallanma oranı n ile 1/(1-n) çarpanı olarak yaklaşıklanır.
Bu, benzetim tabanlı ETAS tahminlerine göre kaba ama hızlı bir yaklaşımdır.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from ..catalog.clean import decluster_gardner_knopoff
from .base import Forecast, ForecastGrid, _gr_mag_fractions
from ..stats.gr import b_value_mle

LAT0, LON0 = 39.0, 35.0
PARAM_NAMES = ["mu", "K", "alpha", "c", "p", "D", "q", "gamma"]
INIT = {"mu": 0.5, "K": 0.02, "alpha": 1.5, "c": 0.01, "p": 1.1, "D": 2.0, "q": 1.7, "gamma": 0.8}
# (alt, üst) sınırlar; mu, K, c, D log uzayında
BOUNDS = {"mu": (1e-4, 1e3), "K": (1e-5, 5.0), "alpha": (0.2, 3.0), "c": (1e-4, 1.0),
          "p": (1.001, 2.0), "D": (0.1, 30.0), "q": (1.05, 3.0), "gamma": (0.0, 2.5)}
LOG_PARAMS = {"mu", "K", "c", "D"}


def _xy_km(lat, lon):
    x = (np.asarray(lon, float) - LON0) * 111.195 * np.cos(np.radians(LAT0))
    y = (np.asarray(lat, float) - LAT0) * 111.195
    return x, y


def _days(t, origin):
    return ((t - origin).dt.total_seconds() / 86400.0).to_numpy()


def _omori_integral(tau, c, p):
    """∫_0^τ (c+s)^-p ds (τ ≥ 0)."""
    tau = np.maximum(tau, 0.0)
    return ((c + tau) ** (1 - p) - c ** (1 - p)) / (1 - p)


@dataclass
class ETAS(Forecast):
    m0: float = 3.0
    bandwidth_km: float = 15.0
    r_max_km: float = 300.0
    dt_max_days: float = 3650.0
    aux_years: float = 5.0  # eğitim başından önceki tetikleyici dönemi
    b: float | None = None
    # benzetimde büyüklük dağılımı: "gr" (b ile M8'e kadar) veya "empirical" (eğitim
    # büyüklükleri; olabilirlikte kullanılan gerçek büyüklüklerle tutarlı)
    mag_sampling: str = "gr"
    name: str = "etas"
    params: dict = field(default_factory=dict)
    fit_info: dict = field(default_factory=dict)

    # ------------------------------------------------------------ arka plan
    def _background(self, train: pd.DataFrame):
        ind = decluster_gardner_knopoff(train)
        bg = train[ind]
        self._bg_xy = np.column_stack(_xy_km(bg["latitude"], bg["longitude"]))
        self._bg_norm = len(bg) * 2 * np.pi * self.bandwidth_km ** 2

    def _u(self, x, y):
        tree = cKDTree(self._bg_xy)
        pts = np.column_stack([x, y])
        out = np.zeros(len(pts))
        h2 = 2 * self.bandwidth_km ** 2
        nb = tree.query_ball_point(pts, 5 * self.bandwidth_km)
        for i, idx in enumerate(nb):
            if idx:
                d2 = ((self._bg_xy[idx] - pts[i]) ** 2).sum(1)
                out[i] = np.exp(-d2 / h2).sum()
        return out / self._bg_norm + 1e-12

    # ------------------------------------------------------------ uyum
    def fit(self, catalog, fit_start, fit_end, init: dict | None = None, fixed: dict | None = None):
        """En çok olabilirlik. ``fixed``: sabit tutulacak parametreler; değer "b" ise
        α = b·ln10 (ör. {"alpha": "b"})."""
        cat = catalog[(catalog["mag"] >= self.m0)].sort_values("time")
        aux_start = fit_start - pd.DateOffset(years=int(self.aux_years))
        hist = cat[(cat["time"] >= aux_start) & (cat["time"] < fit_end)].reset_index(drop=True)
        train = hist[hist["time"] >= fit_start]
        self.origin = fit_start
        self.b = self.b or b_value_mle(train["mag"], self.m0)[0]
        self._background(train)
        self.set_mag_pool(train["mag"])

        t = _days(hist["time"], fit_start)
        x, y = _xy_km(hist["latitude"], hist["longitude"])
        m = hist["mag"].to_numpy() - self.m0
        is_target = t >= 0
        T = (fit_end - fit_start).days
        tgt = np.flatnonzero(is_target)
        u_t = self._u(x[tgt], y[tgt])

        # çiftler: hedef i, tetikleyici j (t_j < t_i, r <= r_max, dt <= dt_max)
        tree = cKDTree(np.column_stack([x, y]))
        pairs = tree.query_pairs(self.r_max_km, output_type="ndarray")
        a, b = pairs[:, 0], pairs[:, 1]
        i = np.where(t[a] > t[b], a, b)
        j = np.where(t[a] > t[b], b, a)
        dt = t[i] - t[j]
        keep = is_target[i] & (dt > 0) & (dt <= self.dt_max_days)
        i, j, dt = i[keep], j[keep], dt[keep]
        r2 = (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2
        pos = np.full(len(t), -1)
        pos[tgt] = np.arange(len(tgt))
        ii = pos[i]
        mj = m[j]
        # integral terimindeki tetikleyiciler: t_j < T
        trig_t, trig_m = t, m

        def unpack(v):
            return {n: (np.exp(val) if n in LOG_PARAMS else val) for n, val in zip(PARAM_NAMES, v)}

        def nll(v):
            P = unpack(v)
            kappa = P["K"] * np.exp(P["alpha"] * mj)
            sig = P["D"] ** 2 * np.exp(P["gamma"] * mj)
            f = (P["q"] - 1) / (np.pi * sig) * (1 + r2 / sig) ** (-P["q"])
            g = (P["c"] + dt) ** (-P["p"])
            trig = np.bincount(ii, weights=kappa * g * f, minlength=len(tgt))
            lam = P["mu"] * u_t + trig
            ll = np.log(lam).sum()
            kap_all = P["K"] * np.exp(P["alpha"] * trig_m)
            lo = np.maximum(0.0, -trig_t)
            hi = np.minimum(T - trig_t, self.dt_max_days)
            integ = P["mu"] * T + (kap_all * (_omori_integral(hi, P["c"], P["p"])
                                              - _omori_integral(lo, P["c"], P["p"]))).sum()
            return -(ll - integ)

        v0 = np.array([np.log(INIT[n]) if n in LOG_PARAMS else INIT[n] for n in PARAM_NAMES])
        v0[0] = np.log(max(len(tgt) / T * 0.4, 1e-3))
        if init:  # sıcak başlangıç (ör. bir önceki yılın parametreleri)
            v0 = np.array([np.log(init[n]) if n in LOG_PARAMS else init[n] for n in PARAM_NAMES])
        bounds = [(np.log(BOUNDS[n][0]), np.log(BOUNDS[n][1])) if n in LOG_PARAMS else BOUNDS[n]
                  for n in PARAM_NAMES]
        for n, val in (fixed or {}).items():
            if val == "b":
                val = self.b * np.log(10)
            k = PARAM_NAMES.index(n)
            v = np.log(val) if n in LOG_PARAMS else float(val)
            v0[k] = v
            bounds[k] = (v, v)
        res = minimize(nll, v0, method="L-BFGS-B", bounds=bounds, options={"maxiter": 500})
        self.params = unpack(res.x)
        self.fit_info = {"loglik": float(-res.fun), "n_targets": int(len(tgt)),
                         "n_pairs": int(len(i)), "success": bool(res.success),
                         "message": str(res.message), "b": float(self.b),
                         "branching_ratio": self.branching_ratio(),
                         "branching_ratio_truncated": self.branching_ratio_truncated()}
        return self

    def branching_ratio(self) -> float:
        P = self.params
        beta = self.b * np.log(10)
        if beta <= P["alpha"]:
            return float("inf")
        return float(P["K"] * P["c"] ** (1 - P["p"]) / (P["p"] - 1) * beta / (beta - P["alpha"]))

    def branching_ratio_truncated(self) -> float:
        """Tetiklenme süresi dt_max ile kesildiğinde dallanma oranı (olabilirlikle tutarlı)."""
        P = self.params
        beta = self.b * np.log(10)
        if beta <= P["alpha"]:
            return float("inf")
        return float(P["K"] * _omori_integral(self.dt_max_days, P["c"], P["p"])
                     * beta / (beta - P["alpha"]))

    # ------------------------------------------------------------ benzetim
    def _sample_bg_xy(self, n, rng):
        k = rng.integers(0, len(self._bg_xy), n)
        return self._bg_xy[k] + rng.normal(0, self.bandwidth_km, (n, 2))

    def _sample_offspring(self, t_par, x_par, y_par, m_par, T, rng):
        """Ebeveynlerin [max(0,t), T] aralığındaki doğrudan artçılarını örnekle."""
        P = self.params
        c, p = P["c"], P["p"]
        lo = np.maximum(0.0, -t_par)  # ebeveyn zamanına göre gecikme alt sınırı
        hi = np.minimum(T - t_par, self.dt_max_days)
        F_lo, F_hi = _omori_integral(lo, c, p), _omori_integral(hi, c, p)
        mean = P["K"] * np.exp(P["alpha"] * m_par) * np.clip(F_hi - F_lo, 0, None)
        n = rng.poisson(mean)
        if n.sum() == 0:
            return np.empty(0), np.empty(0), np.empty(0), np.empty(0)
        par = np.repeat(np.arange(len(n)), n)
        U = rng.uniform(F_lo[par], F_hi[par])
        lag = (U * (1 - p) + c ** (1 - p)) ** (1 / (1 - p)) - c
        t = t_par[par] + lag
        sig = P["D"] ** 2 * np.exp(P["gamma"] * m_par[par])
        r = np.sqrt(sig * ((1 - rng.uniform(size=par.size)) ** (1 / (1 - P["q"])) - 1))
        th = rng.uniform(0, 2 * np.pi, par.size)
        m = self._sample_mag(par.size, rng)
        return t, x_par[par] + r * np.cos(th), y_par[par] + r * np.sin(th), m

    def set_mag_pool(self, mags):
        """Benzetimde büyüklükleri GR yerine eğitim kataloğundan yeniden örnekle."""
        m = np.asarray(mags, float)
        self._mag_pool = m[m >= self.m0] - self.m0

    def _sample_mag(self, n, rng, mmax: float = 8.0):
        pool = getattr(self, "_mag_pool", None)
        if self.mag_sampling == "empirical" and pool is not None and len(pool) >= 50:
            return pool[rng.integers(0, len(pool), n)]
        beta = self.b * np.log(10)
        u = rng.uniform(size=n)
        span = 1 - np.exp(-beta * (mmax - self.m0))
        return -np.log(1 - u * span) / beta  # m - m0

    def simulate(self, start, end, catalog, n_sim: int = 200, seed: int = 0,
                 max_events: int = 200_000):
        """Pencere için ``n_sim`` katalog benzet. Dönen: liste (t_gün, x_km, y_km, m-m0)."""
        rng = np.random.default_rng(seed)
        P = self.params
        T = (end - start).days
        hist = catalog[(catalog["mag"] >= self.m0) & (catalog["time"] < start)
                       & (catalog["time"] >= start - pd.Timedelta(days=self.dt_max_days))]
        th = _days(hist["time"], start)
        xh, yh = _xy_km(hist["latitude"], hist["longitude"])
        mh = hist["mag"].to_numpy() - self.m0
        sims = []
        for _ in range(n_sim):
            nb = rng.poisson(P["mu"] * T)
            xy = self._sample_bg_xy(nb, rng)
            gen = (rng.uniform(0, T, nb), xy[:, 0], xy[:, 1], self._sample_mag(nb, rng))
            hist_off = self._sample_offspring(th, xh, yh, mh, T, rng)
            gen = tuple(np.concatenate([a, b]) for a, b in zip(gen, hist_off))
            allev = [gen]
            total = gen[0].size
            while gen[0].size and total < max_events:
                gen = self._sample_offspring(*gen, T, rng)
                keep = gen[0] < T
                gen = tuple(a[keep] for a in gen)
                allev.append(gen)
                total += gen[0].size
            sims.append(tuple(np.concatenate([g[k] for g in allev]) for k in range(4)))
        return sims

    # ------------------------------------------------------------ tahmin
    def expected(self, grid: ForecastGrid, start, end, catalog: pd.DataFrame | None = None,
                 n_sim: int = 200, seed: int = 0, return_sims: bool = False):
        """Benzetimle beklenen sayılar. ``start`` öncesi olaylar ``catalog``dan alınır.

        Hücre başına M>=m0 ortalaması benzetimden, büyüklük kutularına dağılım
        GR'den gelir (büyük kutulardaki benzetim gürültüsünü önler).
        """
        cat = catalog if catalog is not None else self._catalog
        sims = self.simulate(start, end, cat, n_sim, seed)
        lat_e, lon_e = grid.lat_edges, grid.lon_edges
        counts = np.zeros((n_sim, lat_e.size - 1, lon_e.size - 1))
        for k, (t, x, y, m) in enumerate(sims):
            lat = y / 111.195 + LAT0
            lon = x / (111.195 * np.cos(np.radians(LAT0))) + LON0
            counts[k], _, _ = np.histogram2d(lat, lon, bins=[lat_e, lon_e])
        spatial = counts.mean(0) + 1e-6 / counts[0].size  # boş hücrelere çok küçük taban
        mmin = grid.mag_edges[0]
        scale = 10 ** (-self.b * (mmin - self.m0))
        mag = _gr_mag_fractions(grid.mag_edges, self.b, mmin)
        exp = spatial[:, :, None] * scale * mag[None, None, :]
        if return_sims:
            return exp, counts.sum(axis=(1, 2)) * scale
        return exp
