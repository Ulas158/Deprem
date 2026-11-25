"""En eski veriden bugüne deprem deprem ileri eğitim (uzun yürüyüş).

Ulas'ın isteği (24 Kasım 2026): eğitim 2012'den değil, ilk verilerden başlasın;
her adımda bir sonraki depremin zamanı, yeri ve büyüklüğü tahmin edilsin, gerçek
olay gelince puanlansın ve geçmişe eklensin. Eski kayıtların belirsizliği hem
eğitimde hem güven aralıklarında hesaba katılsın.

Veri
----
* Aletsel katalog (1904 → bugün, AFAD/ISC/EMSC/USGS birleşik)
* Tarihsel katalog (``tarihsel/``: EPICA, SHARE, GEM GHEC, NOAA; MÖ 2000 → 1899).
  1700 sonrası tarihli kayıtlar zaman çizelgesine girer; daha eskileri yalnız uzun
  dönem arka plan haritasına (nerede deprem olur) düşük ağırlıkla katkı verir.

Dönemler ve tamlık
------------------
Katalog her dönemde farklı büyüklükten itibaren eksiksizdir. Her koşu kendi alt
büyüklüğünü (m0) kullanır ve yalnız o büyüklükte eksiksiz olan yıllarla eğitilir:

* M≥5,5: 1875 → bugün (tarihsel katalog Stepp analizi: 1875 sonrası eksiksiz)
* M≥4,5: 1970 → bugün
* M≥4,0: 1990 → bugün

Her yıl Y için (yalnız Y öncesi veriyle):
1. ETAS yeniden uydurulur (bir önceki yılın parametrelerinden sıcak başlangıç).
   Arka plan haritası olayların belirsizlik ağırlıklarıyla düzleştirilir.
2. Model λ = λ_ETAS · exp(η_hücre) eğitim yıllarında Poisson olabilirliğiyle eğitilir.
   Ofset, ETAS'ın o hücre-yıldaki koşullu yoğunluk integralidir (kompansatör); bu,
   olay olay puanlamayla aynı olasılık modelidir. Eski yıllar ağırlıkla düşer:
   yılın ağırlığı = o yıldaki olayların ortalama belirsizlik ağırlığı / en iyi yıl.
3. Y yılındaki her M≥m0 deprem sırayla puanlanır (``forecast.eventwalk.score_year``).
4. Y için beklenen deprem sayısı benzetimle %95 tahmin aralığıyla verilir.

Belirsizlik
-----------
Katalog R kez yeniden örneklenir: her olayın zamanı, konumu ve büyüklüğü kendi 1σ
hatasıyla (yoksa dönem varsayılanı: 1900'lerde ±30 km, ±0,4 büyüklük; bugün ±3 km,
±0,15) rastgele kaydırılır. Her örnekte yürüyüş baştan yapılır (ETAS ve model
yeniden eğitilir). Güven aralıkları hem bu katalog örneklerinden hem de ay blokları
bootstrap'ından gelir; eski dönemlerde hatalar büyük olduğundan aralıklar kendiliğinden
genişler. Rapor, aralığın ne kadarının katalog belirsizliğinden geldiğini de gösterir.
Özellik tablosu (fay, Coulomb, geçmiş sismisite) hesap maliyeti nedeniyle en iyi
tahmin katalogundan bir kez hesaplanır.
"""

from __future__ import annotations

import json
import logging
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .catalog.clean import decluster_gardner_knopoff
from .catalog.weights import fill_uncertainties
from .features.table import FEATURE_REGISTRY, assign_cell, build_table, make_cells
from .forecast.base import ForecastGrid
from .forecast.etas import ETAS, LAT0, LON0, _days, _omori_integral, _xy_km
from .forecast.eventwalk import _cell_index, score_year
from .observations import fit as pfit
from .region import TURKEY

log = logging.getLogger("deprem.longwalk")

PERIODS = [(1900, 1959, "1900-1959"), (1960, 1999, "1960-1999"), (2000, 2100, "2000-bugün")]


# ------------------------------------------------------------------ katalog
def load_historical(root: Path, timed_from: float = 1700.0):
    """Tarihsel katalog: (zaman çizelgesine girenler, yalnız arka plana girenler)."""
    h = pd.read_csv(Path(root) / "tarihsel" / "data" / "tarihsel_katalog.csv")
    h = h[h["is_preferred"].astype(bool) & h["mag"].notna() & ~h["suspect"].astype(bool)].copy()
    old = h[h["decimal_year"] < timed_from]
    new = h[(h["decimal_year"] >= timed_from) & (h["decimal_year"] < 1900)].copy()
    new["time"] = pd.to_datetime(new["time"], utc=True, format="ISO8601")
    keep = ["time", "latitude", "longitude", "depth_km", "mag", "mag_type", "time_err_s",
            "horiz_err_km", "depth_err_km", "mag_err"]
    new = new[keep].copy()
    new["source"] = "tarihsel"
    new["agency"] = h.loc[new.index, "agency"]
    new["event_id"] = h.loc[new.index, "event_id"]
    new["mw_best"] = new["mag"]
    new["n_sources"] = h.loc[new.index, "n_sources"]
    bg_old = old[["latitude", "longitude", "mag", "decimal_year", "training_weight"]].rename(
        columns={"training_weight": "weight"})
    return new, bg_old


def build_long_catalog(root: Path, out: Path) -> pd.DataFrame:
    inst = pd.read_csv(Path(root) / "data" / "processed" / "turkey_catalog.csv.gz")
    inst["time"] = pd.to_datetime(inst["time"], utc=True, format="ISO8601")
    hist, bg_old = load_historical(root)
    hist = hist[hist["time"] < inst["time"].min()]
    cat = pd.concat([hist, inst], ignore_index=True).sort_values("time").reset_index(drop=True)
    cat["depth_km"] = cat["depth_km"].fillna(10.0)
    out.parent.mkdir(parents=True, exist_ok=True)
    cat.to_csv(out, index=False)
    bg_old.to_csv(out.with_name("historical_background.csv"), index=False)
    log.info("Uzun katalog: %d olay (%d tarihsel 1700-1899), arka plan için %d eski kayıt",
             len(cat), len(hist), len(bg_old))
    return cat


def perturb(cat: pd.DataFrame, rng) -> pd.DataFrame:
    """Her olayın zamanını, konumunu, büyüklüğünü kendi 1σ hatasıyla kaydır."""
    f = fill_uncertainties(cat)
    out = cat.copy()
    n = len(out)
    dt = rng.normal(0, 1, n) * f["time_err_s"].to_numpy()
    out["time"] = out["time"] + pd.to_timedelta(dt, unit="s")
    h = f["horiz_err_km"].to_numpy()
    out["latitude"] = out["latitude"] + rng.normal(0, 1, n) * h / 111.195
    out["longitude"] = out["longitude"] + rng.normal(0, 1, n) * h / (111.195 * np.cos(np.radians(out["latitude"])))
    out["mag"] = out["mag"] + rng.normal(0, 1, n) * f["mag_err"].to_numpy()
    if "mw_best" in out:
        out["mw_best"] = out["mw_best"] + (out["mag"] - cat["mag"])
    return out.sort_values("time").reset_index(drop=True)


# ------------------------------------------------------------------ ağırlıklı ETAS
class WETAS(ETAS):
    """Arka plan haritası belirsizlik ağırlıklarıyla düzleştirilen ETAS.

    ``bg_extra``: eğitim penceresinden önceki (tarihsel) olaylar; (lat, lon, ağırlık).
    Tetiklenme parametreleri olağan ETAS gibi uydurulur.
    """

    bg_extra: pd.DataFrame | None = None

    def _background(self, train: pd.DataFrame):
        ind = decluster_gardner_knopoff(train)
        bg = train[ind]
        lat, lon = bg["latitude"].to_numpy(), bg["longitude"].to_numpy()
        w = bg["weight"].to_numpy() if "weight" in bg else np.ones(len(bg))
        ex = getattr(self, "bg_extra", None)
        if ex is not None and len(ex):
            ex = ex[ex["mag"] >= self.m0]
            lat, lon = np.r_[lat, ex["latitude"]], np.r_[lon, ex["longitude"]]
            w = np.r_[w, ex["weight"]]
        self._bg_xy = np.column_stack(_xy_km(lat, lon))
        self._bg_w = np.asarray(w, float)
        self._bg_tree = cKDTree(self._bg_xy)
        self._bg_norm = self._bg_w.sum() * 2 * np.pi * self.bandwidth_km ** 2

    def _u(self, x, y):
        pts = np.column_stack([x, y])
        out = np.zeros(len(pts))
        h2 = 2 * self.bandwidth_km ** 2
        for i, idx in enumerate(self._bg_tree.query_ball_point(pts, 5 * self.bandwidth_km)):
            if idx:
                d2 = ((self._bg_xy[idx] - pts[i]) ** 2).sum(1)
                out[i] = (self._bg_w[idx] * np.exp(-d2 / h2)).sum()
        return out / self._bg_norm + 1e-12

    def _sample_bg_xy(self, n, rng):
        k = rng.choice(len(self._bg_xy), n, p=self._bg_w / self._bg_w.sum())
        return self._bg_xy[k] + rng.normal(0, self.bandwidth_km, (n, 2))


# ------------------------------------------------------------------ yardımcılar
def _cell_geometry(grid):
    lat_c = 0.5 * (grid.lat_edges[:-1] + grid.lat_edges[1:])
    lon_c = 0.5 * (grid.lon_edges[:-1] + grid.lon_edges[1:])
    LA, LO = np.meshgrid(lat_c, lon_c, indexing="ij")
    cx, cy = _xy_km(LA.ravel(), LO.ravel())
    area = (np.diff(grid.lat_edges)[:, None] * 111.195 * np.diff(grid.lon_edges)[None, :] * 111.195
            * np.cos(np.radians(LA))).ravel()
    return cx, cy, area


def background_share(e, grid):
    cx, cy, area = _cell_geometry(grid)
    U = e._u(cx, cy) * area
    return U / U.sum()


def compensator(e, cat_m0: pd.DataFrame, s0, s1, U, grid, n_cells) -> np.ndarray:
    """ETAS'ın [s0, s1) içinde hücre başına koşullu yoğunluk integrali (gerçekleşen geçmişle)."""
    P = e.params
    T = (s1 - s0).total_seconds() / 86400.0
    par = cat_m0[(cat_m0["time"] >= s0 - pd.Timedelta(days=e.dt_max_days)) & (cat_m0["time"] < s1)]
    out = P["mu"] * U * T
    if len(par) == 0:
        return out
    tp = _days(par["time"], s0)
    lo = np.maximum(0.0, -tp)
    hi = np.maximum(np.minimum(T - tp, e.dt_max_days), lo)
    gi = P["K"] * np.exp(P["alpha"] * (par["mag"].to_numpy() - e.m0)) * (
        _omori_integral(hi, P["c"], P["p"]) - _omori_integral(lo, P["c"], P["p"]))
    cell = _cell_index(par["latitude"].to_numpy(), par["longitude"].to_numpy(), grid)
    ok = cell >= 0
    return out + np.bincount(cell[ok], weights=gi[ok], minlength=n_cells)


def year_counts(cat_m0, years, grid, n_cells):
    yrs = cat_m0["time"].dt.year.to_numpy()
    cell = _cell_index(cat_m0["latitude"].to_numpy(), cat_m0["longitude"].to_numpy(), grid)
    out = {}
    for y in years:
        s = (yrs == y) & (cell >= 0)
        out[y] = np.bincount(cell[s], minlength=n_cells).astype(float)
    return out


def year_weights(cat_m0, years, w_col="weight"):
    """Yılın eğitim ağırlığı: o yılın (±5 yıl) olaylarının ortalama belirsizlik ağırlığı / en iyisi."""
    yrs = cat_m0["time"].dt.year.to_numpy()
    w = cat_m0[w_col].to_numpy()
    raw = {}
    for y in years:
        s = np.abs(yrs - y) <= 5
        raw[y] = float(w[s].mean()) if s.any() else np.nan
    v = pd.Series(raw).interpolate(limit_direction="both").fillna(1.0)
    return (v / v.max()).clip(0.02, 1.0).to_dict()


def forecast_counts(e, cat, t0, t1, grid, eta, n_sim, seed):
    sims = e.simulate(t0, t1, cat, n_sim=n_sim, seed=seed)
    w = np.exp(eta)
    nE, nM = [], []
    for (t, x, y, m) in sims:
        lat = y / 111.195 + LAT0
        lon = x / (111.195 * np.cos(np.radians(LAT0))) + LON0
        c = _cell_index(lat, lon, grid)
        c = c[c >= 0]
        nE.append(len(c))
        nM.append(float(w[c].sum()))
    return np.array(nE, float), np.array(nM, float)


# ------------------------------------------------------------------ özellikler
def feature_matrix(run, years, m0, root) -> tuple[pd.DataFrame, list[str]]:
    """Hücre × yıl özellik tablosu (yalnız yıl başından önceki veri)."""
    cache = run.dir / f"long_features_m{m0:g}.csv.gz"
    if cache.exists():
        tab = pd.read_csv(cache)
    else:
        starts = [pd.Timestamp(f"{y}-01-01", tz="UTC") for y in years]
        tab = build_table(run.catalog, TURKEY, run.cfg["data"]["cell_deg"], starts,
                          pd.DateOffset(years=1), target_mags=[m0], mc=m0, faults=run.faults,
                          features=list(FEATURE_REGISTRY.values()))
        tab["year"] = pd.to_datetime(tab["window_start"], utc=True).dt.year
        mu = run.cfg["longwalk"].get("coulomb_mu", 0.4)
        cc = pd.read_csv(run.dir / "coulomb_cells.csv")
        cc = cc[np.isclose(cc["mu_eff"], mu)].drop(columns="mu_eff")
        tab = tab.merge(cc, on=["cell_id", "year"], how="left")
        bc = pd.read_csv(run.dir / "budget_cells.csv")
        tab = tab.merge(bc, on=["cell_id", "year"], how="left")
        if (run.dir / "obs_features.csv.gz").exists():
            tab = tab.merge(pd.read_csv(run.dir / "obs_features.csv.gz"), on=["cell_id", "year"], how="left")
        tab.to_csv(cache, index=False)
    drop = {"cell_id", "window_start", "lat_c", "lon_c", "year"}
    cols = [c for c in tab.columns if c not in drop and not c.startswith("y_")]
    return tab, cols


# ------------------------------------------------------------------ tek yürüyüş
_G = {}


def walk_one(spec: dict, rep: int) -> dict:
    """Bir koşu (m0) için bir katalog örneğinde 1. yıldan sona kadar yürü."""
    run, root, cfg = _G["run"], _G["root"], _G["cfg"]
    out_dir = run.dir / "longwalk"
    out_dir.mkdir(exist_ok=True)
    tag = f"m{spec['m0']:g}_r{rep:02d}"
    f_ev, f_fc, f_cf = (out_dir / f"{tag}_events.csv.gz", out_dir / f"{tag}_forecast.csv",
                        out_dir / f"{tag}_coefs.csv")
    if f_ev.exists() and f_fc.exists():
        return {"tag": tag, "cached": True}
    t_start = time.time()
    m0 = float(spec["m0"])
    cat = run.catalog if rep == 0 else perturb(run.catalog, np.random.default_rng(1000 + rep))
    if rep:
        from .catalog.weights import compute_weights
        cat["weight"] = compute_weights(cat)
    cat_m0 = cat[cat["mag"] >= m0].reset_index(drop=True)
    grid = ForecastGrid(TURKEY, cfg["data"]["cell_deg"], np.array([m0, 10.0]))
    n_cells = (grid.lat_edges.size - 1) * (grid.lon_edges.size - 1)
    fit_start = pd.Timestamp(f"{spec['fit_start']}-01-01", tz="UTC")
    y0, y1 = int(spec["score_from"]), int(spec["score_to"])
    all_years = list(range(int(spec["fit_start"]), y1 + 2))
    tab, cols = feature_matrix(run, all_years, m0, root)
    tab = tab.sort_values(["year", "cell_id"]).reset_index(drop=True)
    counts = year_counts(cat_m0, all_years, grid, n_cells)
    yw = year_weights(cat_m0, all_years)
    bg_old = _G["bg_old"].copy()
    if rep:
        rng = np.random.default_rng(5000 + rep)
        bg_old["latitude"] += rng.normal(0, 50 / 111.195, len(bg_old))
        bg_old["longitude"] += rng.normal(0, 50 / 111.195 / np.cos(np.radians(39)), len(bg_old))
        bg_old["mag"] += rng.normal(0, 0.5, len(bg_old))
    lw = cfg["longwalk"]
    prior_sd = lw.get("prior_sd", 0.5)
    ev_rows, fc_rows, cf_rows = [], [], []
    init = None
    for Y in range(y0, y1 + 2):  # son adım (y1+1): gelecek yıl tahmini, puanlama yok
        t0, t1 = pd.Timestamp(f"{Y}-01-01", tz="UTC"), pd.Timestamp(f"{Y + 1}-01-01", tz="UTC")
        e = WETAS(m0=m0, bandwidth_km=cfg["etas"]["bandwidth_km"], r_max_km=cfg["etas"]["r_max_km"],
                  dt_max_days=cfg["etas"]["dt_max_days"])
        e.bg_extra = bg_old[bg_old["decimal_year"] < spec["fit_start"]]
        try:
            e.fit(cat_m0, fit_start, t0, init=init)
        except Exception as exc:  # çok az olay
            log.warning("%s %d ETAS uydurulamadı: %s", tag, Y, exc)
            continue
        init = dict(e.params)
        U = background_share(e, grid)
        # eğitim: fit_start..Y-1, ofset = bu ETAS'ın kompansatörü
        tr_years = list(range(int(spec["fit_start"]), Y))
        Z_rows, off, yv, wv = [], [], [], []
        tr = tab[tab["year"].isin(tr_years)]
        for yp in tr_years:
            s0, s1 = pd.Timestamp(f"{yp}-01-01", tz="UTC"), pd.Timestamp(f"{yp + 1}-01-01", tz="UTC")
            comp = compensator(e, cat_m0, s0, s1, U, grid, n_cells)
            off.append(np.log(comp + 1e-12))
            yv.append(counts[yp])
            wv.append(np.full(n_cells, yw[yp]))
        off, yv, wv = np.concatenate(off), np.concatenate(yv), np.concatenate(wv)
        use = [c for c in cols if tr[c].notna().any() and tr[c].nunique(dropna=True) > 1]
        Zt, st = run._design(tr, use)
        pm, ps = run._priors(Zt.columns, prior_sd)
        beta, sd, _ = pfit.fit(Zt.to_numpy(), off, yv, pm, ps, weights=wv)
        te = tab[tab["year"] == Y]
        if te.empty:
            te = tab[tab["year"] == Y - 1].assign(year=Y)
        Zs, _ = run._design(te, use, st)
        eta = np.zeros(n_cells)
        eta[te["cell_id"].to_numpy()] = beta[0] + Zs.to_numpy() @ beta[1:]
        eta = np.clip(eta, -6, 6)
        if rep == 0:
            for c_, b_, s_ in zip(["sabit"] + list(Zt.columns), beta, sd):
                cf_rows.append({"yıl": Y, "özellik": c_, "katsayı": float(b_), "sd": float(s_)})
        # olay olay puanlama
        if Y <= y1:
            d = score_year(e, cat, t0, t1, grid, eta, m0)
            if len(d):
                d["year"] = Y
                ev_rows.append(d)
            n_obs = len(d)
        else:
            n_obs = np.nan
        nE, nM = forecast_counts(e, cat, t0, t1, grid, eta, lw.get("n_sim", 100), 7 * Y + rep)
        fc_rows.append({"year": Y, "observed": n_obs, "b": float(e.b), "branching": e.branching_ratio_truncated(),
                        **{f"etas_q{q}": float(np.quantile(nE, q / 100)) for q in (2.5, 50, 97.5)},
                        **{f"model_q{q}": float(np.quantile(nM, q / 100)) for q in (2.5, 50, 97.5)},
                        "etas_mean": float(nE.mean()), "model_mean": float(nM.mean()),
                        "etas_sims": " ".join(f"{v:.0f}" for v in nE),
                        "model_sims": " ".join(f"{v:.1f}" for v in nM)})
        log.info("%s %d: %s olay, ETAS %.1f, model %.1f", tag, Y, n_obs, nE.mean(), nM.mean())
    ev = pd.concat(ev_rows, ignore_index=True) if ev_rows else pd.DataFrame()
    ev.to_csv(f_ev, index=False)
    pd.DataFrame(fc_rows).to_csv(f_fc, index=False)
    if cf_rows:
        pd.DataFrame(cf_rows).to_csv(f_cf, index=False)
    return {"tag": tag, "seconds": round(time.time() - t_start, 1)}


def _init_worker(run, root, cfg, bg_old):
    _G.update(run=run, root=root, cfg=cfg, bg_old=bg_old)


def run_all(run, root: Path, workers: int = 3):
    cfg = run.cfg
    lw = cfg["longwalk"]
    bg_path = (Path(root) / cfg["data"]["catalog"]).with_name("historical_background.csv")
    bg_old = pd.read_csv(bg_path)
    _ = run.catalog  # önce yükle (çatallanan süreçler paylaşır)
    # özellik tablolarını önce sırayla üret (süreçler aynı dosyayı yazmasın)
    for spec in lw["runs"]:
        feature_matrix(run, list(range(int(spec["fit_start"]), int(spec["score_to"]) + 2)), float(spec["m0"]), root)
    jobs = [(spec, r) for r in range(lw.get("replicates", 10) + 1) for spec in lw["runs"]]
    _init_worker(run, root, cfg, bg_old)
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        for res in pool.starmap(walk_one, jobs, chunksize=1):
            log.info("bitti: %s", res)


# ------------------------------------------------------------------ özet
def _blocks(df):
    return pd.to_datetime(df["time"], utc=True, format="ISO8601").dt.strftime("%Y-%m").to_numpy()


def pooled_ci(reps: list[tuple[np.ndarray, np.ndarray]], n_boot=2000, seed=0):
    """Katalog örnekleri × ay blokları bootstrap: her çekilişte rastgele bir örnek, sonra bloklar."""
    rng = np.random.default_rng(seed)
    pre = []
    for v, b in reps:
        if len(v) == 0:
            continue
        ub, inv = np.unique(b, return_inverse=True)
        pre.append((np.bincount(inv, weights=v), np.bincount(inv).astype(float)))
    if not pre:
        return (np.nan, np.nan)
    out = np.empty(n_boot)
    for k in range(n_boot):
        s, c = pre[rng.integers(len(pre))]
        pick = rng.integers(0, len(s), len(s))
        out[k] = s[pick].sum() / max(c[pick].sum(), 1)
    return float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975))


def summarize(run_dir: Path, lw: dict) -> dict:
    d = Path(run_dir) / "longwalk"
    res = {"runs": [], "forecast_next": [], "by_year": []}
    for spec in lw["runs"]:
        m0 = float(spec["m0"])
        evs, fcs = [], []
        for f in sorted(d.glob(f"m{m0:g}_r*_events.csv.gz")):
            r = int(f.name.split("_r")[1][:2])
            try:
                ev = pd.read_csv(f)
            except pd.errors.EmptyDataError:
                ev = pd.DataFrame()
            fc = pd.read_csv(d / f"m{m0:g}_r{r:02d}_forecast.csv")
            evs.append((r, ev))
            fcs.append((r, fc))
        if not evs:
            continue
        R = len(evs) - 1
        for p0, p1, pname in [(spec["score_from"], spec["score_to"], "tüm dönem")] + PERIODS:
            p0, p1 = max(p0, spec["score_from"]), min(p1, spec["score_to"])
            if p0 > p1:
                continue
            row = {"eşik": f"M≥{m0:g}", "dönem": f"{p0}-{p1}" + (" (tümü)" if pname == "tüm dönem" else ""),
                   "katalog örneği": R}
            parts = {"toplam": ("time_ll_model", "space_ll_model", "time_ll_etas", "space_ll_etas"),
                     "zaman": ("time_ll_model", None, "time_ll_etas", None),
                     "konum": ("space_ll_model", None, "space_ll_etas", None)}
            for pn, (a, b2, c, d2) in parts.items():
                reps, base, rep_means = [], None, []
                for r, ev in evs:
                    if ev.empty:
                        continue
                    s = ev[(ev["year"] >= p0) & (ev["year"] <= p1)]
                    if s.empty:
                        continue
                    v = (s[a] + (s[b2] if b2 else 0) - s[c] - (s[d2] if d2 else 0)).to_numpy()
                    reps.append((v, _blocks(s)))
                    rep_means.append(float(v.mean()))
                    if r == 0:
                        base = (v, _blocks(s))
                        row["olay"] = int(len(s))
                if not reps or base is None:
                    continue
                lo, hi = pooled_ci(reps)
                lo0, hi0 = pooled_ci([base])
                row[f"{pn} kazanç"] = round(float(base[0].mean()), 3)
                row[f"{pn} %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
                if pn == "toplam":
                    row["katalog örnekleri ortalaması"] = round(float(np.mean(rep_means)), 3)
                    row["yalnız örnekleme GA"] = f"[{lo0:.3f}, {hi0:.3f}]"
                    w_all, w0 = hi - lo, hi0 - lo0
                    row["katalog belirsizliğinin payı"] = (f"%{100 * max(0.0, 1 - (w0 / w_all) ** 2):.0f}"
                                                         if w_all > 0 else "-")
            for q in (0.05,):
                for who in ("etas", "model"):
                    vals = [((ev[(ev["year"] >= p0) & (ev["year"] <= p1)][f"cell_rank_{who}"] < q).mean())
                            for _, ev in evs if not ev.empty]
                    row[f"en olası %5 hücrede ({'ETAS' if who == 'etas' else 'model'})"] = (
                        f"%{100 * vals[0]:.0f} [{100 * np.quantile(vals, 0.025):.0f}-{100 * np.quantile(vals, 0.975):.0f}]")
            # sayı tahmini kapsaması (yıllık %95 aralık gerçekleşeni içeriyor mu)
            for who in ("etas", "model"):
                cov = []
                for _, fc in fcs:
                    s = fc[(fc["year"] >= p0) & (fc["year"] <= p1) & fc["observed"].notna()]
                    if len(s):
                        cov.append(((s["observed"] >= np.floor(s[f"{who}_q2.5"])) &
                                    (s["observed"] <= np.ceil(s[f"{who}_q97.5"]))).mean())
                if cov:
                    row[f"yıllık sayı aralığı tuttu ({'ETAS' if who == 'etas' else 'model'})"] = f"%{100 * np.mean(cov):.0f}"
            res["runs"].append(row)
        # yıl yıl (örnek 0 + katalog aralığı)
        fc0 = fcs[0][1]
        for _, r_ in fc0.iterrows():
            y = int(r_["year"])
            sims_e = np.concatenate([np.array(fc.loc[fc["year"] == y, "etas_sims"].iloc[0].split(), float)
                                     for _, fc in fcs if (fc["year"] == y).any()])
            sims_m = np.concatenate([np.array(fc.loc[fc["year"] == y, "model_sims"].iloc[0].split(), float)
                                     for _, fc in fcs if (fc["year"] == y).any()])
            obs = [fc.loc[fc["year"] == y, "observed"].iloc[0] for _, fc in fcs if (fc["year"] == y).any()]
            item = {"eşik": f"M≥{m0:g}", "yıl": y,
                    "gerçekleşen": None if pd.isna(r_["observed"]) else int(r_["observed"]),
                    "gerçekleşen aralığı (katalog belirsizliği)": (
                        None if pd.isna(r_["observed"]) else
                        f"{np.nanquantile(obs, 0.025):.0f}-{np.nanquantile(obs, 0.975):.0f}"),
                    "ETAS beklenen": round(float(sims_e.mean()), 1),
                    "ETAS %95": f"{np.quantile(sims_e, 0.025):.0f}-{np.quantile(sims_e, 0.975):.0f}",
                    "model beklenen": round(float(sims_m.mean()), 1),
                    "model %95": f"{np.quantile(sims_m, 0.025):.0f}-{np.quantile(sims_m, 0.975):.0f}"}
            if y > spec["score_to"]:
                b = float(r_["b"])
                for mm in (6.0, 7.0):
                    s = 10 ** (-b * (mm - m0))
                    item[f"en az bir M≥{mm:g} olasılığı (model)"] = f"%{100 * np.mean(1 - np.exp(-sims_m * s)):.0f}"
                    item[f"en az bir M≥{mm:g} olasılığı (ETAS)"] = f"%{100 * np.mean(1 - np.exp(-sims_e * s)):.0f}"
                res["forecast_next"].append(item)
            else:
                res["by_year"].append(item)
    (Path(run_dir) / "longwalk.json").write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str))
    return res
