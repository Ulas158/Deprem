"""Bölge bazında parametre ayarı: kısmi havuzlama (hiyerarşik öncül).

Soru: tüm Türkiye için tek ortak ayar yerine her tektonik bölgenin (Marmara, Kuzey Anadolu,
Doğu Anadolu, Ege ...) kendi katsayıları ve ayarı olursa, dokunulmamış 2020-2025
yıllarında ETAS'a göre tahmin iyileşir mi?

Model (bölge r'deki hücre i, yıl Y):

    λ_i = ETAS_i · exp(β0 + α_r + Σ_k (β_k + δ_rk) z_ik)
    β_k  ~ N(m_k, s_k²)                   ortak (Türkiye geneli) katsayı, önceki aramadaki gibi
    δ_rk ~ N(0, (τ · τ_r · s_k)²)         bölge sapması, ortak katsayıya doğru çekilir
    α_r  ~ N(0, τ0²)                      bölgenin ETAS'a göre genel düzeltmesi

τ → 0: tek ortak ayar (önceki PR). τ → ∞: her bölge tamamen ayrı (havuzlamasız).
Arada: veri azsa bölge ayarı ortak ayara çekilir (kısmi havuzlama). Bölgenin hangi
parametre gruplarında ortaktan sapabileceği, τ, τ0 ve bölgeye özel τ_r arama
yıllarında aranır; sınama yılları (2020-2025) aramada hiç kullanılmaz.

Şans kontrolü: aynı arama, hücrelerin bölge etiketleri rastgele karıştırılmış
("sahte bölgeler", aynı bölge büyüklükleri) verilerde tekrarlanır. Gerçek bölgeler
sahte bölgelerden daha fazla kazandırmıyorsa bölgesel ayrım anlamsızdır.

Bölge sınırları sonuçlara bakılmadan, ``region.SUBREGIONS`` kutularından (batı ucu
Kuzey Ege dahil 25,5°D'ye uzatılarak) önceden belirlendi.

HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize

from . import search as sr
from .observations import fit as pfit

# (ad, min_lat, max_lat, min_lon, max_lon) — ilk eşleşen kazanır
REGIONS = [
    ("Marmara (KAF batı)", 40.0, 41.3, 25.5, 30.5),
    ("Kuzey Anadolu (KAF orta)", 40.0, 41.2, 30.5, 39.5),
    ("Doğu Anadolu Fay Zonu", 36.0, 39.0, 35.8, 41.0),
    ("Ege / Batı Anadolu", 36.5, 40.0, 25.5, 30.5),
    ("Doğu Anadolu platosu (Van-Kars)", 38.0, 41.5, 41.0, 45.0),
]
OTHER = "diğer (İç Anadolu, Akdeniz, Karadeniz, Kıbrıs yayı)"
REGION_NAMES = [r[0] for r in REGIONS] + [OTHER]


def cell_regions(cells: pd.DataFrame) -> np.ndarray:
    """cell_id sırasında her hücrenin bölge numarası (REGION_NAMES dizini)."""
    out = np.full(len(cells), len(REGIONS), dtype=int)
    lat, lon = cells["lat_c"].to_numpy(), cells["lon_c"].to_numpy()
    for k in range(len(REGIONS) - 1, -1, -1):  # ilk eşleşen kazansın diye tersten yaz
        _, a, b, c, d = REGIONS[k]
        out[(lat >= a) & (lat < b) & (lon >= c) & (lon < d)] = k
    order = np.argsort(cells["cell_id"].to_numpy())
    return out[order]


def default_rcfg(groups) -> dict:
    return {"tau": 0.5, "tau0": 0.5, "tau_r": {}, "vary": sorted(groups)}


# ------------------------------------------------------------------ uydurma
def _region_design(Z: pd.DataFrame, reg: np.ndarray, vary_cols: list, R: int):
    """[1 | Z | Z_vary·1(r) her r için | 1(r) her r için] — seyrek; bölge sapma ve sabit sütunları."""
    Zv = sparse.csr_matrix(Z[vary_cols].to_numpy()) if vary_cols else sparse.csr_matrix((len(Z), 0))
    blocks = [sparse.csr_matrix(np.ones((len(Z), 1))), sparse.csr_matrix(Z.to_numpy())]
    for r in range(R):
        blocks.append(sparse.diags((reg == r).astype(float)) @ Zv)
    blocks.append(sparse.csr_matrix(np.eye(R)[reg]))
    return sparse.hstack(blocks, format="csr")


def _fit_sparse(X, offset, y, m, sd, x0=None):
    """pfit.fit ile aynı MAP amacı (Gauss öncüllü Poisson), seyrek X, Hessian'sız."""
    m = np.r_[0.0, m]
    prec = np.r_[1e-4, 1 / np.asarray(sd) ** 2]
    XT = X.T.tocsr()

    def f(beta):
        eta = np.clip(offset + X @ beta, -30, 30)
        lam = np.exp(eta)
        return (lam - y * eta).sum() + 0.5 * (prec * (beta - m) ** 2).sum(), XT @ (lam - y) + prec * (beta - m)

    b0 = m.copy() if x0 is None else x0
    return minimize(f, b0, jac=True, method="L-BFGS-B", options={"maxiter": 2000}).x


def _predict_sparse(X, offset, beta):
    return np.exp(np.clip(offset + X @ beta, -30, 30))


def fit_predict(data: sr.SearchData, cfg, rcfg, Y: int, region_of_cell: np.ndarray, tab=None,
                fit_year=None, cache=None, return_coef=False):
    """Y yılını bölgesel kısmi havuzlamalı modelle tahmin et (bkz. modül açıklaması)."""
    tab = data.tab if tab is None else tab
    R = len(REGION_NAMES)
    cols = sr.selected_columns(data, cfg)
    FY = Y if fit_year is None else fit_year
    te = tab[tab["year"] == Y]
    reg_te = region_of_cell[te["cell_id"].to_numpy()]
    if cache is not None and FY in cache:
        beta, cols_z, st, vary_cols = cache[FY]
    else:
        tr_mask = (tab["year"] < FY) & tab["year"].isin([yp for (yy, yp) in data.off_train if yy == FY])
        if cfg.get("window"):
            tr_mask &= tab["year"] >= FY - cfg["window"]
        tr = tab[tr_mask]
        off_tr = np.concatenate([data.off_train[(FY, int(yp))] for yp in sorted(tr["year"].unique())])
        ycol = "y_w" if cfg.get("uncertainty_weights", True) else f"y_count_m{data.train_mag:g}"
        y = tr[ycol].to_numpy()
        reg_tr = region_of_cell[tr["cell_id"].to_numpy()]
        if cols:
            Z, st = sr.design(tr, cols)
            pm, ps = sr.priors(Z.columns, cfg)
            if cfg.get("ard"):
                beta0, _, _ = pfit.fit(Z.to_numpy(), off_tr, y, pm, ps)
                mag = np.abs(beta0[1:] - pm)
                rel = np.clip(mag / max(np.median(mag[mag > 0]) if (mag > 0).any() else 1.0, 1e-6), 0.05, 3.0)
                ps = ps * rel
        else:
            Z, st, pm, ps = pd.DataFrame(index=tr.index), None, np.zeros(0), np.zeros(0)
        vary = set(rcfg.get("vary", []))
        vary_cols = [c for c in Z.columns if sr.column_group(c) in vary]
        vidx = [list(Z.columns).index(c) for c in vary_cols]
        tau_r = rcfg.get("tau_r", {})
        m_all, s_all = [pm], [ps]
        for r in range(R):
            mult = rcfg["tau"] * float(tau_r.get(REGION_NAMES[r], 1.0))
            m_all.append(np.zeros(len(vidx)))
            s_all.append(np.maximum(ps[vidx] * mult, 1e-4))
        m_all.append(np.zeros(R))
        s_all.append(np.full(R, max(rcfg["tau0"], 1e-4)))
        X = _region_design(Z, reg_tr, vary_cols, R)
        # başlangıç: ortak model çözümü (bölge sapmaları 0) — hızlı yakınsama
        x0 = None
        if cols:
            b_g, _, _ = pfit.fit(Z.to_numpy(), off_tr, y, pm, ps)
            x0 = np.r_[b_g, np.zeros(X.shape[1] - len(b_g))]
        beta = _fit_sparse(X, off_tr, y, np.concatenate(m_all), np.concatenate(s_all), x0)
        cols_z = list(Z.columns)
        if cache is not None:
            cache[FY] = (beta, cols_z, st, vary_cols)
    if cols_z:
        Zs, _ = sr.design(te, cols, st)
        Zs = Zs.reindex(columns=cols_z, fill_value=0.0)
    else:
        Zs = pd.DataFrame(index=te.index)
    Xs = _region_design(Zs, reg_te, vary_cols, R)
    lam = _predict_sparse(Xs, data.off_test[Y], beta)
    if return_coef:
        names = ["(sabit)"] + cols_z + [f"{c} @ {REGION_NAMES[r]}" for r in range(R) for c in vary_cols] \
            + [f"(bölge sabiti) @ {n}" for n in REGION_NAMES]
        return lam, dict(zip(names, beta.tolist()))
    return lam


def evaluate(data: sr.SearchData, cfg, rcfg, years, region_of_cell, tab=None, score_mags=None) -> dict:
    """sr.evaluate ile aynı puan (olay başına bilgi kazancı, ETAS'a göre); ayrıca bölge bazında."""
    tab = data.tab if tab is None else tab
    score_mags = data.score_mags if score_mags is None else score_mags
    per = {mg: [] for mg in data.mags}
    per_reg = {mg: [] for mg in data.mags}  # her olayın bölgesi (per ile aynı sıra)
    cache = {}
    k = max(int(getattr(data, "refit_every", 1)), 1)
    y0 = min(years)
    for Y in years:
        fy = Y - ((Y - y0) % k)
        if rcfg is None:
            lam = sr.fit_predict(data, cfg, Y, tab=tab, fit_year=fy, cache=cache)
        else:
            lam = fit_predict(data, cfg, rcfg, Y, region_of_cell, tab=tab, fit_year=fy, cache=cache)
        te = tab[tab["year"] == Y]
        reg_te = region_of_cell[te["cell_id"].to_numpy()]
        for mg, v in sr.per_event_gains(data, lam, Y, tab=tab).items():
            per[mg].append(v)
            yv = te[f"y_count_m{mg:g}"].to_numpy().astype(int)
            per_reg[mg].append(np.repeat(reg_te, yv) if yv.sum() else np.array([], int))
    pooled = {mg: (np.concatenate(v).mean() if sum(len(x) for x in v) else np.nan) for mg, v in per.items()}
    score = float(np.nanmean([pooled[m] for m in score_mags if m in pooled]))
    return {"score": score, "pooled": pooled, "per_year": per, "per_year_region": per_reg}


def region_breakdown(res, mags, R=None):
    """Bölge başına olay sayısı ve olay başına kazanç (toplam kazanç her olayda aynı
    beklenen sayı farkı terimini içerdiği için bölgeler arası doğrudan toplanabilir)."""
    R = len(REGION_NAMES) if R is None else R
    rows = []
    for r in range(R):
        row = {"bölge": REGION_NAMES[r]}
        for mg in mags:
            g = np.concatenate(res["per_year"][mg]) if res["per_year"][mg] else np.array([])
            rr = np.concatenate(res["per_year_region"][mg]) if res["per_year_region"][mg] else np.array([])
            sel = g[rr == r] if len(g) else g
            row[f"M≥{mg:g} olay"] = int(len(sel))
            row[f"M≥{mg:g} kazanç"] = float(sel.mean()) if len(sel) else np.nan
        rows.append(row)
    return rows


# ------------------------------------------------------------------ arama
def suggest(trial, groups_on) -> dict:
    vary = [g for g in sorted(groups_on) if trial.suggest_categorical(f"v::{g}", [True, False])]
    return {"tau": trial.suggest_float("tau", 0.01, 5.0, log=True),
            "tau0": trial.suggest_float("tau0", 0.01, 2.0, log=True),
            "tau_r": {n: trial.suggest_float(f"r::{n}", 0.2, 5.0, log=True) for n in REGION_NAMES},
            "vary": vary}


def run_search(data, cfg, years, region_of_cell, n_trials, seed=0, tab=None, extra=()):
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    groups_on = list(cfg["groups"])
    sampler = optuna.samplers.TPESampler(seed=seed, n_startup_trials=min(40, n_trials // 3), multivariate=True)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    rows = []

    def obj(trial):
        rc = suggest(trial, groups_on)
        r = evaluate(data, cfg, rc, years, region_of_cell, tab=tab)
        rows.append({"trial": trial.number, "score": r["score"],
                     **{f"M≥{m:g}": r["pooled"][m] for m in data.mags}, "rcfg": json.dumps(rc, ensure_ascii=False)})
        return r["score"] if np.isfinite(r["score"]) else -9.0

    for rc in extra:
        p = {f"v::{g}": g in rc["vary"] for g in groups_on}
        p.update(tau=rc["tau"], tau0=rc["tau0"], **{f"r::{n}": rc["tau_r"].get(n, 1.0) for n in REGION_NAMES})
        study.enqueue_trial(p)
    study.optimize(obj, n_trials=n_trials)
    return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)


def fake_regions(region_of_cell: np.ndarray, seed: int) -> np.ndarray:
    """Bölge etiketlerini hücreler arasında karıştır (bölge büyüklükleri korunur)."""
    return np.random.default_rng(seed).permutation(region_of_cell)


def paired_ci(per_a, per_b, n_boot=4000, seed=1):
    """İki modelin aynı olaylardaki kazanç farkı (a - b) için yıl blok bootstrap."""
    d = [a - b for a, b in zip(per_a, per_b)]
    return sr.year_block_ci(d, n_boot=n_boot, seed=seed)
