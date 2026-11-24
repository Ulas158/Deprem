"""Tek komutla, yeniden üretilebilir ve kaldığı yerden devam edebilen eğitim hattı.

Aşamalar (her biri ``runs/<run_name>/`` altına çıktısını yazar; çıktı varsa
atlanır):

1. ``etas``       ETAS'ı eğitim döneminde en çok olabilirlikle uydur
2. ``etas_fc``    her yıl penceresi için benzetimle ETAS beklenen sayıları
3. ``coulomb``    büyük depremlerden fay alıcılarında ΔCFF, hücre/yıl özetleri
3b. ``budget``    her 1 km elemanın kendi moment/kayma bütçesi, yıl başı anlık görüntüleri
3c. ``observations`` GPS, tekrarlayan deprem, sismik hız, InSAR, GONAF yıllık özellikleri
                  (``deprem.observations``; yalnız yıl başından önceki veri)
4. ``features``   hücre × yıl özellik tablosu (+ Coulomb + ETAS + gözlem özellikleri)
5. ``fit``        ETAS'ı taban alan log-doğrusal Poisson model, μ' seçimi. Gözlem
                  özelliklerinin öncülü literatür başlangıç etkisidir (registry)
6. ``evaluate``   sınama yıllarında ETAS'a karşı: log-olabilirlik, bilgi kazancı,
                  N-testi, Brier, kalibrasyon; eşik eşik (M3/M4/M5/M6)
6b. ``obs_eval``  gözlemlerin katkısı (gözlemli/gözlemsiz model, grup grup) ve
                  zayıf kanıtlı hipotezlerin kısa vadeli testi
7. ``report``     sade dilli rapor + harita katmanları

Model:   λ(hücre, yıl) = Λ_ETAS(hücre, yıl) · exp(β0 + Σ β_k z_k)
ETAS'tan daha iyi değilse β'lar sıfıra yakın kalır ve model ETAS'a eşitlenir;
yani model ancak özellikler gerçekten bilgi taşıyorsa ETAS'ı geçebilir.
"""

from __future__ import annotations

import json
import logging
import shutil
import time
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy import stats

from .catalog.clean import decluster_gardner_knopoff
from .catalog.weights import compute_weights
from .faults.coulomb import build_sources, coulomb_on_receivers, stress_from_sources, to_local_m
from .faults.gem import load_faults
from .faults.mesh import mesh_faults
from .faults.moment import add_geometry
from .features.evaluate import reliability_curve, brier
from .features.table import FEATURE_REGISTRY, assign_cell, build_table, make_cells
from .forecast.base import ForecastGrid
from .forecast.etas import ETAS
from .forecast.tests import n_test, poisson_loglik
from .observations import fit as pfit
from .observations.registry import OBS_FEATURES, prior_for
from .region import SUBREGIONS, TURKEY

log = logging.getLogger("deprem.pipeline")
STAGES = ["etas", "etas_fc", "coulomb", "budget", "observations", "features", "fit", "evaluate",
          "obs_eval", "report"]


def ts(s):
    return pd.Timestamp(s, tz="UTC")


class Run:
    def __init__(self, config_path: Path, root: Path, force=()):
        self.root = Path(root)
        self.cfg = tomllib.loads(Path(config_path).read_text())
        self.dir = self.root / self.cfg["out_dir"] / self.cfg["run_name"]
        self.dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(config_path, self.dir / "config.toml")
        self.force = set(STAGES) if "all" in force else set(force)
        self.status_path = self.dir / "status.json"
        self.status = json.loads(self.status_path.read_text()) if self.status_path.exists() else {}
        c = self.cfg
        self.train_years = list(range(c["split"]["first_train_year"], c["split"]["first_test_year"]))
        self.test_years = list(range(c["split"]["first_test_year"], c["split"]["last_test_year"] + 1))
        self.years = self.train_years + self.test_years
        self.mags = [float(m) for m in c["model"]["eval_mags"]]
        self.train_mag = float(c["model"]["train_mag"])
        self.grid = ForecastGrid(TURKEY, c["data"]["cell_deg"],
                                 np.array(sorted(set(self.mags + [self.train_mag])) + [10.0]))
        self._cache = {}

    # ------------------------------------------------------------ yardımcılar
    def done(self, stage):
        return self.status.get(stage, {}).get("done") and stage not in self.force

    def mark(self, stage, **info):
        self.status[stage] = {"done": True, "at": pd.Timestamp.now(tz="UTC").isoformat(), **info}
        self.status_path.write_text(json.dumps(self.status, indent=2, default=str))

    @property
    def catalog(self):
        if "cat" not in self._cache:
            cat = pd.read_csv(self.root / self.cfg["data"]["catalog"])
            cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
            cat = cat.sort_values("time").reset_index(drop=True)
            cat["weight"] = compute_weights(cat)
            cat["independent"] = decluster_gardner_knopoff(cat)
            self._cache["cat"] = cat
        return self._cache["cat"]

    @property
    def faults(self):
        if "faults" not in self._cache:
            f = add_geometry(load_faults(self.root / self.cfg["data"]["faults"], TURKEY))
            self._cache["faults"] = f[f["slip_rate_mm_yr"].notna()].reset_index(drop=True)
        return self._cache["faults"]

    def receivers(self):
        if "rec" not in self._cache:
            el = mesh_faults(self.faults, 1.0)
            step = self.cfg["coulomb"]["receiver_strike_step"]
            mid = el.groupby("fault_id")["j_dip"].transform(lambda j: j.max() // 2)
            rec = el[(el["j_dip"] == mid) & (el["i_strike"] % step == 0)].reset_index(drop=True)
            self._cache["rec"] = rec
        return self._cache["rec"]

    def window(self, y):
        return ts(f"{y}-01-01"), ts(f"{y + self.cfg['split']['window_years']}-01-01")

    # ------------------------------------------------------------ aşamalar
    def stage_etas(self):
        c = self.cfg["etas"]
        fit_end = ts(f"{self.cfg['split']['first_test_year']}-01-01")
        e = ETAS(m0=c["m0"], bandwidth_km=c["bandwidth_km"], r_max_km=c["r_max_km"],
                 dt_max_days=c["dt_max_days"]).fit(self.catalog, ts(c["fit_start"]), fit_end)
        out = {"params": {k: float(v) for k, v in e.params.items()}, "fit_info": e.fit_info,
               "fit_start": c["fit_start"], "fit_end": str(fit_end.date())}
        (self.dir / "etas_params.json").write_text(json.dumps(out, indent=2))
        self._cache["etas"] = e
        return out["fit_info"]

    def etas_model(self):
        if "etas" not in self._cache:
            c = self.cfg["etas"]
            fit_end = ts(f"{self.cfg['split']['first_test_year']}-01-01")
            saved = json.loads((self.dir / "etas_params.json").read_text())
            e = ETAS(m0=c["m0"], bandwidth_km=c["bandwidth_km"], r_max_km=c["r_max_km"],
                     dt_max_days=c["dt_max_days"])
            # arka planı ve b'yi aynı eğitim verisinden yeniden kur, parametreleri yükle
            cat = self.catalog
            train = cat[(cat["mag"] >= e.m0) & (cat["time"] >= ts(c["fit_start"]))
                        & (cat["time"] < fit_end)]
            e.b = saved["fit_info"]["b"]
            e._background(train)
            e.params = saved["params"]
            self._cache["etas"] = e
        return self._cache["etas"]

    def stage_etas_fc(self):
        e = self.etas_model()
        out, totals = {}, {}
        for y in self.years:
            s, t = self.window(y)
            exp, tot = e.expected(self.grid, s, t, catalog=self.catalog,
                                  n_sim=self.cfg["etas"]["n_sim"], seed=self.cfg["seed"] + y,
                                  return_sims=True)
            out[str(y)] = exp
            totals[str(y)] = tot
            log.info("ETAS %d: beklenen M≥%.1f = %.1f", y, self.grid.mag_edges[0], exp.sum())
        np.savez_compressed(self.dir / "etas_forecasts.npz", **out)
        np.savez_compressed(self.dir / "etas_sim_totals.npz", **totals)

    def stage_coulomb(self):
        c = self.cfg["coulomb"]
        src = build_sources(self.catalog, self.faults, min_mw=c["min_mw"])
        src = src[src["time"] >= ts(c["source_since"])].reset_index(drop=True)
        rec = self.receivers()
        obs = to_local_m(rec["lat"], rec["lon"], rec["depth_km"])
        per_src = np.zeros((len(src), len(rec), 6), dtype=np.float32)
        for k in range(len(src)):
            per_src[k] = stress_from_sources(obs, src.iloc[[k]])
        np.save(self.dir / "coulomb_stress_per_source.npy", per_src)
        src.to_csv(self.dir / "coulomb_sources.csv", index=False)
        rec.to_csv(self.dir / "coulomb_receivers.csv.gz", index=False)
        # hücre × yıl özetleri her μ' için
        cells = make_cells(TURKEY, self.cfg["data"]["cell_deg"])
        rec_ev = rec.rename(columns={"lat": "latitude", "lon": "longitude"})
        rec_cell = assign_cell(rec_ev, TURKEY, cells)
        rows = []
        times = pd.to_datetime(src["time"], utc=True)
        for mu in c["mu_eff_grid"]:
            for y in self.years + [self.years[-1] + 1]:
                t0, _ = self.window(y)
                sel = (times < t0).to_numpy()
                stress = per_src[sel].sum(0) if sel.any() else np.zeros((len(rec), 6))
                dcff, _, _ = coulomb_on_receivers(stress.astype(float), rec["strike"].to_numpy(),
                                                  rec["dip"].to_numpy(),
                                                  rec["rake"].fillna(180.0).to_numpy(), mu)
                bar = np.clip(dcff / 1e5, -50, 50)
                df = pd.DataFrame({"cell_id": rec_cell, "dcff": bar})
                df = df[df["cell_id"] >= 0]
                g = df.groupby("cell_id")["dcff"]
                agg = pd.DataFrame({"dcff_max": g.max(), "dcff_mean": g.mean(),
                                    "dcff_frac_pos": g.apply(lambda v: (v > 0.1).mean())})
                agg["mu_eff"], agg["year"] = mu, y
                rows.append(agg.reset_index())
        pd.concat(rows).to_csv(self.dir / "coulomb_cells.csv", index=False)
        return {"n_sources": len(src), "n_receivers": len(rec)}

    def elements(self):
        if "el" not in self._cache:
            self._cache["el"] = mesh_faults(self.faults, 1.0)
        return self._cache["el"]

    def stage_budget(self):
        from .faults.budget import element_budget
        c = self.cfg.get("budget", {})
        el = self.elements()
        end = ts(f"{self.cfg['split']['last_test_year'] + 1}-01-01")
        now = self.catalog["time"].max().ceil("D")
        snaps = [self.window(y)[0] for y in self.years]
        res = element_budget(el, self.catalog, ts(c.get("start", "1900-01-01")), max(end, now),
                             snapshots=snaps, min_mw=c.get("min_mw", 4.0))
        fin = pd.concat([el[["element_id", "fault_id", "lon", "lat", "depth_km", "area_km2"]],
                         res["final"]], axis=1)
        fin.to_csv(self.dir / "element_budget_final.csv.gz", index=False)
        cells = make_cells(TURKEY, self.cfg["data"]["cell_deg"])
        cid = assign_cell(el.rename(columns={"lat": "latitude", "lon": "longitude"}), TURKEY, cells)
        area = el["area_km2"].to_numpy() * 1e6
        rows = []
        for t0, d in res["snapshots"].items():
            slip = d / (3.0e10 * area)
            df = pd.DataFrame({"cell_id": cid, "slip": slip, "w": area})
            df = df[df["cell_id"] >= 0]
            g = df.groupby("cell_id")
            agg = pd.DataFrame({
                "slip_deficit_mean": g.apply(lambda x: np.average(x["slip"], weights=x["w"])),
                "slip_deficit_max": g["slip"].max()})
            agg["year"] = t0.year
            rows.append(agg.reset_index())
        pd.concat(rows).to_csv(self.dir / "budget_cells.csv", index=False)
        return {"n_events": res["n_events"], "n_elements": len(el)}

    def obs_enabled(self):
        return self.cfg.get("observations", {}).get("use", False)

    def stage_observations(self):
        if not self.obs_enabled():
            return {"skipped": True}
        from .observations.features import yearly_table
        cells = make_cells(TURKEY, self.cfg["data"]["cell_deg"])
        inc = self.cfg["observations"].get("include", "all")
        tab = yearly_table(self.root, cells, self.years, self.catalog,
                           include=None if inc == "all" else inc)
        tab.to_csv(self.dir / "obs_features.csv.gz", index=False)
        cov = {c: round(float(tab[c].notna().mean()), 3) for c in tab.columns if c in OBS_FEATURES}
        (self.dir / "obs_coverage.json").write_text(json.dumps(
            {"overall": cov, "by_year": {int(y): {c: round(float(g[c].notna().mean()), 3) for c in cov}
                                         for y, g in tab.groupby("year")}}, indent=2))
        return {"columns": list(cov)}

    def stage_features(self):
        cat = self.catalog
        cells = make_cells(TURKEY, self.cfg["data"]["cell_deg"])
        etas_fc = np.load(self.dir / "etas_forecasts.npz")
        starts = [self.window(y)[0] for y in self.years]
        ctx_extra = {"etas": {self.window(int(y))[0]: etas_fc[y].sum(-1).ravel() for y in etas_fc.files}}

        def f_etas(past, cells_, t0, ctx):
            return pd.DataFrame({"etas_log_rate": np.log(ctx_extra["etas"][t0] + 1e-9)},
                                index=cells_["cell_id"])

        groups = self.cfg["model"]["features"]
        fns = list(FEATURE_REGISTRY.values()) if groups == "all" else [FEATURE_REGISTRY[g] for g in groups]
        tab = build_table(cat, TURKEY, self.cfg["data"]["cell_deg"], starts,
                          pd.DateOffset(years=self.cfg["split"]["window_years"]),
                          target_mags=self.mags, mc=self.cfg["data"]["mc"],
                          faults=self.faults, features=fns + [f_etas])
        # ağırlıklı etiket: eğitim büyüklüğünün üstündeki olayların belirsizlik ağırlıkları toplamı
        cat2 = cat.copy()
        cat2["cell_id"] = assign_cell(cat2, TURKEY, cells)
        wl = []
        for t0 in starts:
            t1 = t0 + pd.DateOffset(years=self.cfg["split"]["window_years"])
            sel = cat2[(cat2["time"] >= t0) & (cat2["time"] < t1) & (cat2["mag"] >= self.train_mag)
                       & (cat2["cell_id"] >= 0)]
            w = sel.groupby("cell_id")["weight"].sum().reindex(cells["cell_id"]).fillna(0)
            wl.append(w.to_numpy())
        tab["y_wcount_train"] = np.concatenate(wl)
        # Coulomb özellikleri (her μ' için ayrı sütun grubu)
        cc = pd.read_csv(self.dir / "coulomb_cells.csv")
        tab["year"] = tab["window_start"].dt.year
        for mu, g in cc.groupby("mu_eff"):
            g = g.drop(columns="mu_eff").rename(columns={
                "dcff_max": f"dcff_max_mu{mu:g}", "dcff_mean": f"dcff_mean_mu{mu:g}",
                "dcff_frac_pos": f"dcff_frac_pos_mu{mu:g}"})
            tab = tab.merge(g, on=["cell_id", "year"], how="left")
        bc = pd.read_csv(self.dir / "budget_cells.csv")
        tab = tab.merge(bc, on=["cell_id", "year"], how="left")
        if self.obs_enabled() and (self.dir / "obs_features.csv.gz").exists():
            ob = pd.read_csv(self.dir / "obs_features.csv.gz")
            tab = tab.merge(ob, on=["cell_id", "year"], how="left")
        tab.to_csv(self.dir / "feature_table.csv.gz", index=False)
        return {"rows": len(tab), "cols": tab.shape[1]}

    # ------------------------------------------------------------ model
    @staticmethod
    def _design(tab, cols, stats_=None):
        X = tab[cols].astype(float)
        miss_cols = [c for c in cols if X[c].isna().any()]
        if stats_ is None:
            stats_ = {"mean": X.mean().fillna(0).to_dict(), "std": X.std().replace(0, 1).fillna(1).to_dict(),
                      "miss": miss_cols}
        Z = (X - pd.Series(stats_["mean"])) / pd.Series(stats_["std"])
        Z = Z.clip(-5, 5).fillna(0.0)
        for c in stats_["miss"]:
            Z[f"{c}__eksik"] = X[c].isna().astype(float)
        return Z, stats_

    @staticmethod
    def _priors(columns, prior_sd):
        """Her sütun için öncül (ortalama, sd): gözlem özellikleri literatürden, diğerleri 0 ± prior_sd."""
        m, s = [], []
        for c in columns:
            p = prior_for(c)
            m.append(p[0] if p else 0.0)
            s.append(p[1] if p else prior_sd)
        return np.array(m), np.array(s)

    @classmethod
    def _fit_poisson(cls, Z, offset, y, prior_sd, return_sd=False):
        pm, ps = cls._priors(Z.columns, prior_sd)
        beta, sd, res = pfit.fit(Z.to_numpy(), np.asarray(offset), np.asarray(y), pm, ps)
        return (beta, res, sd) if return_sd else (beta, res)

    def feature_columns(self, tab, mu):
        drop = {"cell_id", "window_start", "lat_c", "lon_c", "year", "y_max_mag", "y_wcount_train",
                "etas_log_rate", "etas_train", "y_fit"}
        excl = set(getattr(self, "_exclude_cols", ()))
        cols = [c for c in tab.columns if c not in drop and not c.startswith("y_")
                and not c.startswith("dcff_") and c not in excl]
        cols += [f"dcff_max_mu{mu:g}", f"dcff_mean_mu{mu:g}", f"dcff_frac_pos_mu{mu:g}"]
        return cols

    def stage_fit(self):
        tab = pd.read_csv(self.dir / "feature_table.csv.gz", parse_dates=["window_start"])
        e = self.etas_model()
        scale = 10 ** (-e.b * (self.train_mag - self.grid.mag_edges[0]))
        tab["etas_train"] = np.exp(tab["etas_log_rate"]) * scale
        train = tab[tab["year"].isin(self.train_years)].copy()
        raw_col = f"y_count_m{self.train_mag:g}"
        if self.cfg["model"]["use_uncertainty_weights"]:
            # Ağırlıklar göreli güvenilirliktir: toplam sayı korunacak şekilde ölçekle,
            # yoksa model oranı sistematik olarak düşük öğrenir.
            k = train[raw_col].sum() / max(train["y_wcount_train"].sum(), 1e-9)
            train["y_fit"] = train["y_wcount_train"] * k
            tab["y_fit"] = tab["y_wcount_train"] * k
        else:
            train["y_fit"] = train[raw_col]
            tab["y_fit"] = tab[raw_col]
        y_col = "y_fit"
        # μ' seçimi: eğitimin son 2 yılı doğrulama
        inner_tr = train[train["year"] < self.train_years[-2]]
        inner_va = train[train["year"] >= self.train_years[-2]]
        sel = []
        for mu in self.cfg["coulomb"]["mu_eff_grid"]:
            cols = self.feature_columns(tab, mu)
            Zt, st = self._design(inner_tr, cols)
            beta, _ = self._fit_poisson(Zt, np.log(inner_tr["etas_train"]), inner_tr[y_col].to_numpy(),
                                        self.cfg["model"]["prior_sd"])
            Zv, _ = self._design(inner_va, cols, st)
            lam = np.exp(np.log(inner_va["etas_train"]) + np.column_stack(
                [np.ones(len(Zv)), Zv.to_numpy()]) @ beta)
            ll = poisson_loglik(lam, inner_va[f"y_count_m{self.train_mag:g}"].to_numpy())
            sel.append({"mu_eff": mu, "val_loglik": ll})
        best_mu = max(sel, key=lambda r: r["val_loglik"])["mu_eff"]
        cols = self.feature_columns(tab, best_mu)
        Z, st = self._design(train, cols)
        beta, res, sd = self._fit_poisson(Z, np.log(train["etas_train"]), train[y_col].to_numpy(),
                                          self.cfg["model"]["prior_sd"], return_sd=True)
        coef = dict(zip(["(sabit)"] + list(Z.columns), beta.tolist()))
        pm, ps = self._priors(Z.columns, self.cfg["model"]["prior_sd"])
        prior_post = [{"özellik": c, "grup": OBS_FEATURES[c].group, "kanıt": OBS_FEATURES[c].evidence,
                       "öncül_ort": float(m_), "öncül_sd": float(s_), "öğrenilen": float(b_),
                       "sonsal_sd": float(d_)}
                      for c, m_, s_, b_, d_ in zip(Z.columns, pm, ps, beta[1:], sd[1:]) if c in OBS_FEATURES]
        out = {"mu_eff_selection": sel, "mu_eff": best_mu, "columns": cols, "design_stats": st,
               "coef": coef, "coef_sd": dict(zip(["(sabit)"] + list(Z.columns), sd.tolist())),
               "obs_prior_posterior": prior_post, "converged": bool(res.success), "etas_b": e.b}
        (self.dir / "model.json").write_text(json.dumps(out, indent=2, default=float))
        return {"mu_eff": best_mu, "n_coef": len(beta)}

    def predict(self, tab):
        m = json.loads((self.dir / "model.json").read_text())
        e = self.etas_model()
        scale = 10 ** (-e.b * (self.train_mag - self.grid.mag_edges[0]))
        base = np.exp(tab["etas_log_rate"]) * scale
        Z, _ = self._design(tab, m["columns"], m["design_stats"])
        beta = np.array([m["coef"]["(sabit)"]] + [m["coef"][c] for c in Z.columns])
        lam = base * np.exp(np.column_stack([np.ones(len(Z)), Z.to_numpy()]) @ beta)
        return base.to_numpy(), lam.to_numpy(), m["etas_b"]

    def stage_evaluate(self):
        tab = pd.read_csv(self.dir / "feature_table.csv.gz", parse_dates=["window_start"])
        base, lam, b = self.predict(tab)
        tab["lam_etas"], tab["lam_model"] = base, lam
        test = tab[tab["year"].isin(self.test_years)]
        rows, rel = [], {}
        for m in self.mags:
            s = 10 ** (-b * (m - self.train_mag))
            le, lm = test["lam_etas"].to_numpy() * s, test["lam_model"].to_numpy() * s
            y = test[f"y_count_m{m:g}"].to_numpy()
            ya = test[f"y_any_m{m:g}"].to_numpy()
            ll_e, ll_m = poisson_loglik(le, y), poisson_loglik(lm, y)
            n = y.sum()
            if n > 0:
                idx = np.repeat(np.arange(len(y)), y.astype(int))
                per = np.log(lm[idx]) - np.log(le[idx]) - (lm.sum() - le.sum()) / n
                ig = per.mean()
                p = stats.ttest_1samp(per, 0).pvalue if n > 1 else np.nan
            else:
                ig, p = np.nan, np.nan
            pe, pm = 1 - np.exp(-le), 1 - np.exp(-lm)
            rows.append({
                "eşik": f"M≥{m:g}", "gözlenen_olay": int(n),
                "ETAS_beklenen": round(le.sum(), 1), "model_beklenen": round(lm.sum(), 1),
                "ETAS_N_testi_q": round(n_test(le, y).quantile, 4),
                "model_N_testi_q": round(n_test(lm, y).quantile, 4),
                "logL_ETAS": round(ll_e, 1), "logL_model": round(ll_m, 1),
                "bilgi_kazancı_olay_başına": round(float(ig), 3) if np.isfinite(ig) else None,
                "p_değeri": float(p) if np.isfinite(p) else None,
                "Brier_ETAS": round(brier(pe, ya), 5), "Brier_model": round(brier(pm, ya), 5),
            })
            rel[f"M{m:g}"] = {
                "etas": reliability_curve(pe, ya, n_bins=8).to_dict("records"),
                "model": reliability_curve(pm, ya, n_bins=8).to_dict("records")}
        by_year = []
        for yr in self.test_years:
            t = test[test["year"] == yr]
            yv = t[f"y_count_m{self.train_mag:g}"].to_numpy()
            by_year.append({"yıl": yr, "gözlenen": int(yv.sum()),
                            "ETAS": round(t["lam_etas"].sum(), 1), "model": round(t["lam_model"].sum(), 1),
                            "logL_ETAS": round(poisson_loglik(t["lam_etas"].to_numpy(), yv), 1),
                            "logL_model": round(poisson_loglik(t["lam_model"].to_numpy(), yv), 1)})
        out = {"by_threshold": rows, "by_year": by_year, "reliability": rel}
        (self.dir / "evaluation.json").write_text(json.dumps(out, indent=2, default=float))
        tab[["cell_id", "year", "lat_c", "lon_c", "lam_etas", "lam_model"] +
            [f"y_count_m{m:g}" for m in self.mags]].to_csv(self.dir / "predictions.csv.gz", index=False)
        return {"thresholds": len(rows)}

    # ------------------------------------------------------------ gözlemlerin katkısı
    def _fit_eval_cols(self, tab, cols, prior_scale: float = 1.0):
        """Verilen sütunlarla eğit (μ' sabit), sınamada eşik eşik ETAS'a karşı bilgi kazancı."""
        e = self.etas_model()
        b = e.b
        train = tab[tab["year"].isin(self.train_years)]
        test = tab[tab["year"].isin(self.test_years)]
        Zt, st = self._design(train, cols)
        if prior_scale != 1.0:
            pm, ps = self._priors(Zt.columns, self.cfg["model"]["prior_sd"])
            beta, _, _ = pfit.fit(Zt.to_numpy(), np.log(train["etas_train"]).to_numpy(),
                                  train["y_fit"].to_numpy(), pm * prior_scale, ps)
        else:
            beta, _ = self._fit_poisson(Zt, np.log(train["etas_train"]), train["y_fit"].to_numpy(),
                                        self.cfg["model"]["prior_sd"])
        Zs, _ = self._design(test, cols, st)
        lam = pfit.predict(Zs.to_numpy(), np.log(test["etas_train"]).to_numpy(), beta)
        out = {}
        for m in self.mags:
            s = 10 ** (-b * (m - self.train_mag))
            le, lm = test["etas_train"].to_numpy() * s, lam * s
            y = test[f"y_count_m{m:g}"].to_numpy()
            n = y.sum()
            if n == 0:
                out[f"M≥{m:g}"] = (None, None)
                continue
            idx = np.repeat(np.arange(len(y)), y.astype(int))
            per = np.log(lm[idx]) - np.log(le[idx]) - (lm.sum() - le.sum()) / n
            out[f"M≥{m:g}"] = (round(float(per.mean()), 4),
                               float(stats.ttest_1samp(per, 0).pvalue) if n > 1 else None)
        return out, per

    def stage_obs_eval(self):
        res = {}
        if self.obs_enabled():
            tab = pd.read_csv(self.dir / "feature_table.csv.gz", parse_dates=["window_start"])
            m = json.loads((self.dir / "model.json").read_text())
            e = self.etas_model()
            scale = 10 ** (-e.b * (self.train_mag - self.grid.mag_edges[0]))
            tab["etas_train"] = np.exp(tab["etas_log_rate"]) * scale
            raw_col = f"y_count_m{self.train_mag:g}"
            tr = tab["year"].isin(self.train_years)
            k = tab.loc[tr, raw_col].sum() / max(tab.loc[tr, "y_wcount_train"].sum(), 1e-9)
            tab["y_fit"] = tab["y_wcount_train"] * k if self.cfg["model"]["use_uncertainty_weights"] else tab[raw_col]
            cols = m["columns"]
            obs_cols = [c for c in cols if c in OBS_FEATURES]
            groups = {}
            for c in obs_cols:
                groups.setdefault(OBS_FEATURES[c].group, []).append(c)
            variants = {"gözlemsiz (önceki model)": [c for c in cols if c not in obs_cols],
                        "tüm gözlemlerle": cols}
            for g, gc in groups.items():
                variants[f"yalnız {g} eklenmiş"] = [c for c in cols if c not in obs_cols or c in gc]
            rows = []
            for name, vc in variants.items():
                r, _ = self._fit_eval_cols(tab, vc)
                rows.append({"model": name, **{f"{k_} kazanç": v[0] for k_, v in r.items()},
                             **{f"{k_} p": v[1] for k_, v in r.items()}})
            # öncülün tek başına etkisi: öncül sıfır vs literatür (aynı veriyle)
            r0, _ = self._fit_eval_cols(tab, cols, prior_scale=0.0)
            rows.append({"model": "tüm gözlemler, öncül = 0 (literatürsüz)",
                         **{f"{k_} kazanç": v[0] for k_, v in r0.items()},
                         **{f"{k_} p": v[1] for k_, v in r0.items()}})
            res["yearly_ablation"] = rows
            res["coverage"] = json.loads((self.dir / "obs_coverage.json").read_text())
        st = self.cfg.get("observations", {}).get("shortterm", {})
        if st.get("use", False):
            res["shortterm"] = self._shortterm(st)
        (self.dir / "obs_evaluation.json").write_text(json.dumps(res, indent=2, default=float))
        return {"variants": len(res.get("yearly_ablation", []))}

    def _shortterm(self, cfg):
        from .observations import shortterm as stm
        obs = self.root / "data" / "obs"
        boxes = stm.make_boxes(TURKEY, cfg.get("box_deg", 2.0))
        start = pd.Timestamp(f"{self.train_years[0]}-01-01")
        end = pd.Timestamp(f"{self.test_years[-1] + 1}-01-01")
        train_end = pd.Timestamp(f"{self.test_years[0]}-01-01")
        sources = {
            "tec_daily.csv.gz": (["tec_z_max", "tec_z_min", "tec_frac_pos"], 300.0),
            "ulf_daily.csv.gz": (["ulf_zh_rel_z", "ulf_z_rel_z"], cfg.get("ulf_radius_km", 200.0)),
            "thermal_daily.csv.gz": (["thermal_dts_z", "thermal_dsa_z"], 80.0),
        }
        out = {}
        for mag in cfg.get("target_mags", [4.0, 4.5]):
            panel = stm.build_panel(self.catalog, TURKEY, boxes, start, end, mag)
            wins = pd.DatetimeIndex(sorted(panel["win"].unique()))
            feats = {}
            for fn, (cols, rad) in sources.items():
                if not (obs / fn).exists():
                    continue
                d = pd.read_csv(obs / fn, parse_dates=["date"])
                if "disturbed" in d.columns:
                    d = d[~d["disturbed"].astype(bool)]
                for c in cols:
                    if c not in d.columns:
                        continue
                    db = stm.daily_to_box(d[["date", "lat", "lon", c]].dropna(), c, boxes, rad)
                    agg = "min" if c == "tec_z_min" else ("mean" if c.endswith("frac_pos") else "max")
                    arr = stm.lead_feature(db, wins, agg)
                    # panel sırası: pencere × kutu
                    feats[c] = arr.reshape(-1)
            if not feats:
                continue
            res, coefs = stm.evaluate(panel, feats, train_end, end)
            out[f"M≥{mag:g}"] = {"results": res.to_dict("records"), "coefs": coefs.to_dict("records"),
                                 "n_boxes": len(boxes), "window_days": stm.WIN_DAYS, "lead_days": stm.LEAD_DAYS}
        return out

    def stage_report(self):
        from .report import write_report, write_layers
        write_report(self)
        if self.cfg.get("map", {}).get("write_layers", False):
            write_layers(self)

    def run(self):
        t_all = time.time()
        for st in STAGES:
            if self.done(st):
                log.info("[%s] tamamlanmış, atlanıyor", st)
                continue
            log.info("[%s] başlıyor", st)
            t = time.time()
            info = getattr(self, f"stage_{st}")() or {}
            self.mark(st, seconds=round(time.time() - t, 1), **info)
            log.info("[%s] bitti (%.1f s)", st, time.time() - t)
        log.info("Tüm aşamalar bitti (%.1f s). Çıktılar: %s", time.time() - t_all, self.dir)
