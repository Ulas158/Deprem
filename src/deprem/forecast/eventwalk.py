"""Deprem deprem ileri değerlendirme: her olaydan önce "bir sonraki deprem" tahmini.

Her hedef deprem i için, i'den hemen önceki ana kadar olan tüm geçmiş bilinir.
Model ve ETAS şu soruya olasılık verir: bir sonraki deprem **ne zaman**, **nerede**
ve **hangi büyüklükte** olacak? Gerçek olay gelince puanlanır, sonra geçmişe eklenir.
Koşullu yoğunluk her olayda güncellenir; parametreler yılda bir (yalnız o yıl
öncesi veriyle) yeniden eğitilir.

Olay başına log-olabilirlik üç parçaya ayrılır (Ogata 1998; Zechar vd. 2013):

* zaman:  log Λ(t_i) − ∫_{t_{i−1}}^{t_i} Λ(t) dt   (bekleme süresi yoğunluğu)
* konum:  log λ(t_i, x_i, y_i) − log Λ(t_i)       (olayın yerindeki uzaysal yoğunluk, km⁻²)
* büyüklük: Gutenberg-Richter log yoğunluğu (b her yıl yeniden)

ETAS: λ_E = μ u(x,y) + Σ κ_j g(t−t_j) f_j(x,y)
Model: λ_M = λ_E · exp(η_k)  (η_k: o yılın ileri eğitilmiş hücre çarpanı; arka plan
ve tetiklenme aynı çarpanla ölçeklenir. Tetikleme çekirdeği birkaç km genişliğinde
olduğundan uzaysal integralde kütle ebeveynin 0,5° hücresine atanır.)

Ayrıca sezgisel ölçüt: olay, o andaki yoğunluğa göre en olası %1 / %5 / %10
hücreden birinde mi gerçekleşti?
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .etas import ETAS, _omori_integral, _xy_km

log = logging.getLogger(__name__)


def _cell_index(lat, lon, grid):
    i = np.searchsorted(grid.lat_edges, lat, side="right") - 1
    j = np.searchsorted(grid.lon_edges, lon, side="right") - 1
    nlon = grid.lon_edges.size - 1
    ok = (i >= 0) & (i < grid.lat_edges.size - 1) & (j >= 0) & (j < nlon)
    return np.where(ok, i * nlon + j, -1)


def score_year(e: ETAS, catalog: pd.DataFrame, t0, t1, grid, eta: np.ndarray | None,
               targets_min_mag: float) -> pd.DataFrame:
    """[t0, t1) içindeki M ≥ targets_min_mag olayları için olay başına puanlar."""
    P = e.params
    m0 = e.m0
    cat = catalog[catalog["mag"] >= m0].sort_values("time")
    hist = cat[(cat["time"] >= t0 - pd.Timedelta(days=e.dt_max_days)) & (cat["time"] < t1)].reset_index(drop=True)
    origin = t0
    t = ((hist["time"] - origin).dt.total_seconds() / 86400.0).to_numpy()
    x, y = _xy_km(hist["latitude"], hist["longitude"])
    m = hist["mag"].to_numpy()
    cell = _cell_index(hist["latitude"].to_numpy(), hist["longitude"].to_numpy(), grid)
    n_cells = (grid.lat_edges.size - 1) * (grid.lon_edges.size - 1)
    kappa = P["K"] * np.exp(P["alpha"] * (m - m0))
    sig = P["D"] ** 2 * np.exp(P["gamma"] * (m - m0))
    # hücre ağırlıkları ve arka plan kütlesi
    w = np.ones(n_cells) if eta is None else np.exp(eta)
    lat_c = 0.5 * (grid.lat_edges[:-1] + grid.lat_edges[1:])
    lon_c = 0.5 * (grid.lon_edges[:-1] + grid.lon_edges[1:])
    LA, LO = np.meshgrid(lat_c, lon_c, indexing="ij")
    cx, cy = _xy_km(LA.ravel(), LO.ravel())
    area = (np.diff(grid.lat_edges)[:, None] * 111.195 * np.diff(grid.lon_edges)[None, :] * 111.195
            * np.cos(np.radians(LA))).ravel()
    U = e._u(cx, cy) * area
    U = U / U.sum()
    bgM = P["mu"] * (w * U).sum()
    wc = np.where(cell >= 0, w[np.clip(cell, 0, None)], 1.0)
    # ebeveyn çiftleri (uzaysal yoğunluk için)
    is_t = (t >= 0) & (m >= targets_min_mag - 1e-9) & (cell >= 0)
    tgt = np.flatnonzero(is_t)
    tree = cKDTree(np.column_stack([x, y]))
    rows = []
    prev_t = None
    prev_idx_all = np.flatnonzero(t < 0)
    # büyüklük-üstü hedefler için önceki hedef zamanı (bekleme süresi)
    tgt_scale = 10 ** (-e.b * (targets_min_mag - m0))
    before = cat[(cat["time"] < t0) & (cat["mag"] >= targets_min_mag)]
    prev_t = ((before["time"].iloc[-1] - origin).total_seconds() / 86400.0) if len(before) else None
    for i in tgt:
        par = np.flatnonzero((t < t[i]) & (t >= t[i] - e.dt_max_days))
        dtp = t[i] - t[par]
        g = (P["c"] + dtp) ** (-P["p"])
        # anlık toplam yoğunluk
        LamE = P["mu"] + (kappa[par] * g).sum()
        LamM = bgM + (wc[par] * kappa[par] * g).sum()
        # uzaysal yoğunluk olay noktasında
        near = np.array(tree.query_ball_point([x[i], y[i]], e.r_max_km), dtype=int)
        near = near[(t[near] < t[i]) & (t[near] >= t[i] - e.dt_max_days)]
        r2 = (x[near] - x[i]) ** 2 + (y[near] - y[i]) ** 2
        f = (P["q"] - 1) / (np.pi * sig[near]) * (1 + r2 / sig[near]) ** (-P["q"])
        lamE = P["mu"] * e._u(np.array([x[i]]), np.array([y[i]]))[0] + (
            kappa[near] * (P["c"] + t[i] - t[near]) ** (-P["p"]) * f).sum()
        lamM = lamE * w[cell[i]]
        # bekleme süresi integrali (önceki hedef büyüklükteki olaydan)
        if prev_t is None:
            prev_t = t[i] - 30.0
        lo_ = np.maximum(prev_t - t[par], 0.0)
        hi_ = t[i] - t[par]
        gi = _omori_integral(hi_, P["c"], P["p"]) - _omori_integral(lo_, P["c"], P["p"])
        # ebeveyn kendisi bekleme aralığında doğduysa lo_=0 doğru
        IE = P["mu"] * (t[i] - prev_t) + (kappa[par] * gi).sum()
        IM = bgM * (t[i] - prev_t) + (wc[par] * kappa[par] * gi).sum()
        # hücre sırası: o anki hücre yoğunlukları
        cellE = P["mu"] * U + np.bincount(np.clip(cell[par], 0, None), weights=kappa[par] * g * (cell[par] >= 0),
                                          minlength=n_cells)
        cellM = cellE * w
        rankE = (cellE > cellE[cell[i]]).mean()
        rankM = (cellM > cellM[cell[i]]).mean()
        s = tgt_scale
        beta = e.b * np.log(10)
        rows.append({
            "time": hist["time"].iloc[i], "lat": hist["latitude"].iloc[i], "lon": hist["longitude"].iloc[i],
            "mag": m[i], "wait_days": t[i] - prev_t,
            "time_ll_etas": np.log(LamE * s) - IE * s, "time_ll_model": np.log(LamM * s) - IM * s,
            "space_ll_etas": np.log(lamE) - np.log(LamE), "space_ll_model": np.log(lamM) - np.log(LamM),
            "mag_ll": np.log(beta) - beta * (m[i] - targets_min_mag),
            "cell_rank_etas": rankE, "cell_rank_model": rankM,
        })
        prev_t = t[i]
    return pd.DataFrame(rows)


def block_ci(values: np.ndarray, blocks: np.ndarray, n_boot: int = 2000, seed: int = 0):
    """Ay blokları üzerinden bootstrap (aynı dizi içindeki olaylar bağımsız değildir)."""
    rng = np.random.default_rng(seed)
    ub = np.unique(blocks)
    if len(ub) < 3:
        return (np.nan, np.nan)
    sums = np.array([values[blocks == b].sum() for b in ub])
    cnts = np.array([(blocks == b).sum() for b in ub])
    pick = rng.integers(0, len(ub), (n_boot, len(ub)))
    bs = sums[pick].sum(1) / cnts[pick].sum(1)
    return float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def summarize(df: pd.DataFrame) -> dict:
    if df.empty:
        return {}
    blocks = df["time"].dt.strftime("%Y-%m").to_numpy()
    out = {"olay": int(len(df))}
    for part, a, b in [("zaman", "time_ll_model", "time_ll_etas"), ("konum", "space_ll_model", "space_ll_etas")]:
        d = (df[a] - df[b]).to_numpy()
        lo, hi = block_ci(d, blocks)
        out[f"{part} kazancı"] = round(float(d.mean()), 4)
        out[f"{part} %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
    tot = (df["time_ll_model"] + df["space_ll_model"] - df["time_ll_etas"] - df["space_ll_etas"]).to_numpy()
    lo, hi = block_ci(tot, blocks)
    out["toplam kazanç"] = round(float(tot.mean()), 4)
    out["toplam %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
    for q in (0.01, 0.05, 0.10):
        out[f"en olası %{int(q * 100)} hücrede (ETAS)"] = f"%{100 * (df['cell_rank_etas'] < q).mean():.0f}"
        out[f"en olası %{int(q * 100)} hücrede (model)"] = f"%{100 * (df['cell_rank_model'] < q).mean():.0f}"
    out["büyüklük log-puanı (ort.)"] = round(float(df["mag_ll"].mean()), 3)
    return out
