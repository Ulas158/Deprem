"""Parametre kombinasyonu araması: arama 2015-2019, dokunulmamış sınama 2020-2025.

    python scripts/train.py --config configs/train_search.toml          # ortak veri (bir kez)
    python scripts/search_params.py --config configs/train_search.toml  # arama + sınama + rapor

Adımlar (ayrıntı: src/deprem/search.py):
  1. gerçek arama (Bayesçi optimizasyon, TPE) yalnız arama yıllarında
  2. aynı arama, özellikleri karıştırılmış verilerde (şans sınırı)
  3. ilk K kombinasyon + karşılaştırma modelleri sınama yıllarında (Holm düzeltmesi)
  4. en iyi kombinasyon için deprem deprem (zaman+konum) puanı
  5. katalog belirsizliği: katalog hatalarıyla yeniden örneklenip sınama tekrarlanır
  6. rapor: reports/parametre_aramasi.md

HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # paralel süreçlerde iş parçacığı taşmasını önle

import argparse
import json
import logging
import pickle
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from deprem import search as sr
from deprem.pipeline import Run

ROOT = Path(__file__).resolve().parents[1]
log = logging.getLogger("search")


def _null_worker(args):
    data_path, seed, years, n_trials = args
    data = pickle.load(open(data_path, "rb"))
    tab = sr.shuffled_table(data, seed)
    df = sr.run_search(data, years, n_trials, seed=1000 + seed, tab=tab)
    return seed, df


def _replicate_worker(args):
    from deprem.features.table import FEATURE_REGISTRY, assign_cell, build_table, make_cells
    from deprem.longwalk import perturb
    from deprem.region import TURKEY
    cfg_path, data_path, rep, cfg, years, mags_ew = args
    run = Run(Path(cfg_path), ROOT)
    data = pickle.load(open(data_path, "rb"))
    rng = np.random.default_rng(10_000 + rep)
    cat = perturb(run.catalog, rng)
    from deprem.catalog.weights import compute_weights
    cat["weight"] = compute_weights(cat)
    cells = make_cells(TURKEY, run.cfg["data"]["cell_deg"])
    starts = [run.window(y)[0] for y in data.years]
    feats = run.cfg["model"]["features"]
    feats = list(FEATURE_REGISTRY) if feats == "all" else feats
    cat_groups = [g for g in ["counts", "moment_release", "time_since", "b_value", "rate_change", "depth",
                              "clustering", "swarms_fluid"] if g in feats]
    fns = [FEATURE_REGISTRY[g] for g in cat_groups if g in FEATURE_REGISTRY]
    t = build_table(cat, TURKEY, run.cfg["data"]["cell_deg"], starts, pd.DateOffset(years=1),
                    target_mags=run.mags, mc=run.cfg["data"]["mc"], faults=run.faults, features=fns)
    t["year"] = t["window_start"].dt.year
    t = t.sort_values(["year", "cell_id"]).reset_index(drop=True)
    tab = data.tab.copy()
    for c in t.columns:
        if c in tab.columns and c not in ("cell_id", "year", "window_start"):
            tab[c] = t[c].to_numpy()
    cat["cell_id"] = assign_cell(cat, TURKEY, cells)
    wl = []
    for t0 in starts:
        sel = cat[(cat["time"] >= t0) & (cat["time"] < t0 + pd.DateOffset(years=1)) & (cat["mag"] >= run.train_mag)
                  & (cat["cell_id"] >= 0)]
        wl.append(sel.groupby("cell_id")["weight"].sum().reindex(cells["cell_id"]).fillna(0).to_numpy())
    tab["y_wcount_train"] = np.concatenate(wl)
    k = tab[f"y_count_m{run.train_mag:g}"].sum() / max(tab["y_wcount_train"].sum(), 1e-9)
    tab["y_w"] = tab["y_wcount_train"] * k
    tab = sr.add_friction_proxy(tab)
    res = sr.evaluate(data, cfg, years, tab=tab)
    out = {"rep": rep, **{f"M≥{m:g}": res["pooled"][m] for m in data.mags}, "score": res["score"]}
    # deprem deprem (bozulmuş katalogla)
    ew = eventwalk(run, data, cfg, years, mags_ew, catalog=cat, tab=tab)
    for tg, g in ew.groupby("target"):
        out[f"olay {tg} toplam"] = float((g["time_ll_model"] + g["space_ll_model"] - g["time_ll_etas"]
                                         - g["space_ll_etas"]).mean())
    return out


def eventwalk(run, data, cfg, years, mags, catalog=None, tab=None):
    from deprem.forecast.etas import ETAS
    from deprem.forecast.eventwalk import score_year
    from deprem.pipeline import ts
    pj = json.loads((run.dir / "walk_etas_params.json").read_text())
    c = run.cfg["etas"]
    cat = run.catalog if catalog is None else catalog
    tab = data.tab if tab is None else tab
    n_cells = (run.grid.lat_edges.size - 1) * (run.grid.lon_edges.size - 1)
    frames = []
    for Y in years:
        lam = sr.fit_predict(data, cfg, Y, tab=tab)
        eta = np.zeros(n_cells)
        eta[tab.loc[tab["year"] == Y, "cell_id"].to_numpy()] = np.log(lam) - data.off_test[Y]
        t0, t1 = run.window(Y)
        e = ETAS(m0=c["m0"], bandwidth_km=c["bandwidth_km"], r_max_km=c["r_max_km"], dt_max_days=c["dt_max_days"])
        train = cat[(cat["mag"] >= e.m0) & (cat["time"] >= ts(c["fit_start"])) & (cat["time"] < t0)]
        e._background(train)
        e.params = pj[str(Y)]["params"]
        e.b = pj[str(Y)]["b"]
        for mg in mags:
            d = score_year(e, cat, t0, t1, run.grid, eta, mg)
            d["year"], d["target"] = Y, f"M≥{mg:g}"
            frames.append(d)
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "configs/train_search.toml")
    ap.add_argument("--search-years", default=None, help="varsayılan: yapılandırmadaki [search] search_years")
    ap.add_argument("--test-years", default="2020-2025")
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--null-runs", type=int, default=4)
    ap.add_argument("--top-k", type=int, default=5)
    ap.add_argument("--replicates", type=int, default=8)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--skip", nargs="*", default=[], help="atlanacak adımlar: null, replicates, eventwalk")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    span = lambda s: list(range(int(s.split("-")[0]), int(s.split("-")[1]) + 1))
    run = Run(args.config, ROOT)
    sy = span(args.search_years or run.cfg.get("search", {}).get("search_years", "2015-2019"))
    ty = span(args.test_years)
    out = run.dir / "search"
    out.mkdir(exist_ok=True)
    dp = out / "data.pkl"
    if not dp.exists():
        pickle.dump(sr.load(run), open(dp, "wb"))
    data = pickle.load(open(dp, "rb"))
    assert not set(sy) & set(ty), "arama ve sınama yılları çakışamaz"
    groups = sr.groups_of(data)
    summary = {"search_years": sy, "test_years": ty, "groups": {g: c for g, c in groups.items()},
               "n_trials": args.trials}

    # 1) gerçek arama --------------------------------------------------------
    t = time.time()
    default = sr.default_config(groups)
    etas_only = dict(default, groups=[])
    trials_p = out / "trials.csv"
    if trials_p.exists():
        trials = pd.read_csv(trials_p)
    else:
        trials = sr.run_search(data, sy, args.trials, seed=0, extra_configs=[default])
        trials.to_csv(trials_p, index=False)
    log.info("arama bitti: %d deneme, en iyi %.4f (%.0f s)", len(trials), trials["score"].iloc[0], time.time() - t)

    # 2) şans sınırı ---------------------------------------------------------
    null_best = []
    if "null" not in args.skip:
        jobs = [(str(dp), k, sy, args.trials) for k in range(args.null_runs)]
        with ProcessPoolExecutor(args.workers) as ex:
            for k, df in ex.map(_null_worker, jobs):
                df.to_csv(out / f"null_trials_{k}.csv", index=False)
                null_best.append(float(df["score"].iloc[0]))
        log.info("karıştırılmış aramaların en iyi puanları: %s", null_best)
    summary["null_best_scores"] = null_best
    summary["real_best_score"] = float(trials["score"].iloc[0])

    # 3) sınama --------------------------------------------------------------
    top = []
    seen = set()
    for _, r in trials.iterrows():
        cfg = json.loads(r["config"])
        key = json.dumps({k: cfg[k] for k in ("groups", "ard", "energy_prior", "window", "uncertainty_weights")},
                         sort_keys=True)
        if key in seen:
            continue
        seen.add(key)
        top.append((f"aramada {len(top) + 1}. sıradaki", cfg, float(r["score"])))
        if len(top) >= args.top_k:
            break
    best = top[0][1]
    compare = [("yalnız ETAS'ın kendi düzeltmesi (özelliksiz)", etas_only, None),
               ("önceki model: tüm parametreler, varsayılan ayar", default, None)]
    # en iyi kombinasyonun çeşitlemeleri: enerji öncülü aç/kapa, ARD aç/kapa
    compare.append(("en iyi + enerji öncülü " + ("kapalı" if best["energy_prior"] else "açık"),
                    dict(best, energy_prior=not best["energy_prior"]), None))
    compare.append(("en iyi + veriyle önem ağırlığı (ARD) " + ("kapalı" if best["ard"] else "açık"),
                    dict(best, ard=not best["ard"]), None))
    # doluluk özetleri: ortalama / moment ağırlıklı / en dolu parça / zincirleme (temel: kapasite grupları dışı)
    cap_groups = {"doluluk oranı", "doluluk + Coulomb", "zincirleme kırılma", "kapasite ve kayaç",
                  "GPS kilitlenme / krip", "sürtünme vekili"}
    base_groups = [g for g in default["groups"] if g not in cap_groups]
    variants = {"doluluk: basit ortalama": ["fill_mean"], "doluluk: moment ağırlıklı ortalama": ["fill_m0w"],
                "doluluk: en dolu parça": ["fill_max", "fill_p90"],
                "doluluk + Coulomb aktarımı": ["fillc_mean", "fillc_max"],
                "zincirleme kırılma (olası uzunluk/büyüklük)": ["cascade_len_km", "cascade_mw"],
                "kapasite + kayaç sertliği + kabuk": ["cap_log_m0", "rigidity_gpa", "moho_km"],
                "GPS kilitlenme / krip": ["locking_mean", "creep_mm_yr"],
                "sürtünme vekili": ["fric_velocity_strengthening"]}
    var_rows = []
    for name, cols in variants.items():
        cols = [c for c in cols if c in data.feature_cols]
        if not cols:
            continue
        d2 = sr.SearchData(**{**data.__dict__, "feature_cols": [c for c in data.feature_cols
                                                               if sr.column_group(c) not in cap_groups] + cols})
        gname = sr.column_group(cols[0])
        cfgv = dict(default, groups=base_groups + [gname])
        cfg0 = dict(default, groups=base_groups)
        rs, rt = sr.evaluate(d2, cfgv, sy), sr.evaluate(d2, cfgv, ty)
        r0s, r0t = sr.evaluate(d2, cfg0, sy), sr.evaluate(d2, cfg0, ty)
        var_rows.append({"özet": name, "arama dönemi katkısı": rs["score"] - r0s["score"],
                         "sınama dönemi katkısı": rt["score"] - r0t["score"],
                         **{f"sınama M≥{m:g} katkı": rt["pooled"][m] - r0t["pooled"][m] for m in data.score_mags}})
    summary["capacity_variants"] = var_rows

    rows, pvals = [], []
    for name, cfg, sscore in top + compare:
        rs = sr.evaluate(data, cfg, sy) if sscore is None else None
        rt = sr.evaluate(data, cfg, ty)
        row = {"model": name, "arama puanı": sscore if sscore is not None else rs["score"],
               "sınama puanı": rt["score"]}
        for m in data.mags:
            lo, hi, p = sr.year_block_ci(rt["per_year"][m])
            row[f"M≥{m:g} kazanç"] = rt["pooled"][m]
            row[f"M≥{m:g} %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
            row[f"M≥{m:g} geçtiği yıl"] = f"{sum(1 for v in rt['per_year'][m] if len(v) and v.mean() > 0)}/" \
                                         f"{sum(1 for v in rt['per_year'][m] if len(v))}"
            row[f"M≥{m:g} p"] = p
        # birleşik p: M≥3/4/5 ortalaması için yıl blok bootstrap
        comb = [np.concatenate([rt["per_year"][m][i] for m in data.score_mags]) for i in range(len(ty))]
        lo, hi, p = sr.year_block_ci(comb)
        row["birleşik %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
        row["birleşik p (tek yönlü)"] = p
        row["config"] = json.dumps(cfg, ensure_ascii=False)
        pvals.append(p)
        rows.append(row)
    adj = sr.holm(pvals)
    for r, a in zip(rows, adj):
        r["Holm düzeltilmiş p"] = float(a)
    summary["test_table"] = rows
    _, coef = sr.fit_predict(data, best, ty[-1], return_coef=True)
    summary["best_coef_last_year"] = coef
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float))

    # 4) deprem deprem ----------------------------------------------------------
    if "eventwalk" not in args.skip:
        m0 = run.cfg["etas"]["m0"]
        ew_mags = sorted({max(float(m), m0) for m in data.score_mags})
        ew = eventwalk(run, data, best, ty, ew_mags)
        ew.to_csv(out / "eventwalk_best.csv.gz", index=False)
        from deprem.forecast.eventwalk import summarize
        summary["eventwalk"] = [{"hedef": tg, **summarize(g)} for tg, g in ew.groupby("target", sort=False)]
        (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float))

    # 5) katalog belirsizliği ---------------------------------------------------
    if "replicates" not in args.skip and args.replicates > 0:
        m0 = run.cfg["etas"]["m0"]
        rep_mags = sorted({max(float(m), m0) for m in data.score_mags})[-2:]
        jobs = [(str(args.config), str(dp), k, best, ty, rep_mags) for k in range(args.replicates)]
        with ProcessPoolExecutor(min(args.workers, 3)) as ex:
            reps = list(ex.map(_replicate_worker, jobs))
        pd.DataFrame(reps).to_csv(out / "replicates.csv", index=False)
        summary["replicates"] = reps
        (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float))
    log.info("bitti: %s", out)


if __name__ == "__main__":
    main()
