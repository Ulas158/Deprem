"""Bölge bazında parametre ayarı (kısmi havuzlama): arama yılları + dokunulmamış 2020-2025.

    python scripts/train.py --config configs/train_search.toml              # ortak veri (bir kez)
    python scripts/search_regional.py --config configs/train_search.toml    # kısa: 2012+, M≥3
    python scripts/search_regional.py --config configs/train_search_long.toml --trials 150  # uzun: 1900+, M≥5,5

Adımlar (yöntem: src/deprem/regional.py):
  1. ortak (tek Türkiye) ayar araması — önceki PR'daki arama (varsa sonucu yeniden kullanılır)
  2. bölgesel arama: en iyi ortak ayarın üstüne bölge sapmaları (τ, τ0, τ_r, hangi gruplar)
  3. şans kontrolü: aynı bölgesel arama, bölge etiketleri karıştırılmış (sahte) bölgelerle
  4. sınama 2020-2025: ETAS'a ve ortak ayara karşı, yıl blok bootstrap, Holm düzeltmesi
  5. rapor: reports/bolgesel_ayar.md (scripts/regional_report.py)

HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import logging
import pickle
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from deprem import regional as rg
from deprem import search as sr
from deprem.features.table import make_cells
from deprem.pipeline import Run
from deprem.region import TURKEY

ROOT = Path(__file__).resolve().parents[1]
log = logging.getLogger("regional")


def _null_worker(args):
    data_path, cfg, seed, years, n_trials, reg = args
    data = pickle.load(open(data_path, "rb"))
    fake = rg.fake_regions(np.asarray(reg), seed)
    df = rg.run_search(data, cfg, years, fake, n_trials, seed=2000 + seed)
    return seed, df


def _global_worker(args):
    data_path, years, n_trials, seed = args
    data = pickle.load(open(data_path, "rb"))
    groups = sr.groups_of(data)
    return sr.run_search(data, years, n_trials, seed=seed, extra_configs=[sr.default_config(groups)])


def _regional_worker(args):
    data_path, cfg, years, reg, n_trials, extra = args
    data = pickle.load(open(data_path, "rb"))
    return rg.run_search(data, cfg, years, np.asarray(reg), n_trials, seed=0, extra=extra)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "configs/train_search.toml")
    ap.add_argument("--test-years", default="2020-2025")
    ap.add_argument("--global-trials", type=int, default=300)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--null-runs", type=int, default=4)
    ap.add_argument("--top-k", type=int, default=3)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    span = lambda s: list(range(int(s.split("-")[0]), int(s.split("-")[1]) + 1))
    run = Run(args.config, ROOT)
    sy = span(run.cfg.get("search", {}).get("search_years", "2015-2019"))
    ty = span(args.test_years)
    assert not set(sy) & set(ty), "arama ve sınama yılları çakışamaz"
    sdir = run.dir / "search"
    sdir.mkdir(exist_ok=True)
    out = run.dir / "regional"
    out.mkdir(exist_ok=True)
    dp = sdir / "data.pkl"
    if not dp.exists():
        pickle.dump(sr.load(run), open(dp, "wb"))
    data = pickle.load(open(dp, "rb"))
    groups = sr.groups_of(data)
    cells = make_cells(TURKEY, run.cfg["data"]["cell_deg"])
    reg = rg.cell_regions(cells)
    summary = {"search_years": sy, "test_years": ty, "regions": rg.REGIONS, "other": rg.OTHER,
               "n_cells_by_region": {n: int((reg == k).sum()) for k, n in enumerate(rg.REGION_NAMES)},
               "n_trials": args.trials}

    # 1) ortak ayar araması (önceki PR) + 2) bölgesel arama ---------------------
    t = time.time()
    gp = sdir / "trials.csv"
    with ProcessPoolExecutor(args.workers) as ex:
        if not gp.exists():
            pd.DataFrame(ex.submit(_global_worker, (str(dp), sy, args.global_trials, 0)).result()).to_csv(
                gp, index=False)
        gtrials = pd.read_csv(gp)
        best_global = json.loads(gtrials["config"].iloc[0])
        log.info("ortak ayar: arama puanı %.4f (%.0f s)", gtrials["score"].iloc[0], time.time() - t)
        near_global = {"tau": 0.01, "tau0": 0.01, "tau_r": {}, "vary": []}
        start = [rg.default_rcfg(best_global["groups"]), near_global,
                 {"tau": 0.01, "tau0": 0.3, "tau_r": {}, "vary": []}]
        rp = out / "trials.csv"
        fut_real = None if rp.exists() else ex.submit(
            _regional_worker, (str(dp), best_global, sy, reg.tolist(), args.trials, start))
        # 3) şans kontrolü: sahte bölgeler ----------------------------------------
        null_paths = [out / f"null_trials_{k}.csv" for k in range(args.null_runs)]
        jobs = [(str(dp), best_global, k, sy, args.trials, reg.tolist())
                for k, p in enumerate(null_paths) if not p.exists()]
        null_futs = [ex.submit(_null_worker, j) for j in jobs]
        if fut_real is not None:
            fut_real.result().to_csv(rp, index=False)
        for f in null_futs:
            k, df = f.result()
            df.to_csv(out / f"null_trials_{k}.csv", index=False)
    rtrials = pd.read_csv(rp)
    log.info("bölgesel arama: en iyi %.4f (%.0f s)", rtrials["score"].iloc[0], time.time() - t)
    g_search = sr.evaluate(data, best_global, sy)["score"]
    summary["global_search_score"] = g_search
    summary["global_config"] = best_global
    summary["regional_best_search_score"] = float(rtrials["score"].iloc[0])
    summary["null_best_scores"] = [float(pd.read_csv(p)["score"].max()) for p in null_paths]
    summary["null_rcfgs"] = [pd.read_csv(p)["rcfg"].iloc[0] for p in null_paths]

    # 4) sınama ----------------------------------------------------------------
    top, seen = [], set()
    for _, r in rtrials.iterrows():
        rc = json.loads(r["rcfg"])
        key = json.dumps(sorted(rc["vary"])) + f"{np.log10(rc['tau']):.1f}"
        if key in seen:
            continue
        seen.add(key)
        top.append((f"bölgesel: aramada {len(top) + 1}. sıradaki", best_global, rc, float(r["score"])))
        if len(top) >= args.top_k:
            break
    best_rc = top[0][2]
    etas_only = dict(sr.default_config(groups), groups=[])
    compare = [
        ("yalnız ETAS'ın kendi düzeltmesi (özelliksiz)", etas_only, None, None),
        ("ortak ayar (tek Türkiye, önceki PR'daki gibi aranmış)", best_global, None, None),
        ("bölgesel: yalnız bölge sabiti (katsayılar ortak)", best_global,
         dict(best_rc, vary=[], tau=0.01), None),
        ("havuzlamasız: her bölge tamamen ayrı", best_global,
         {"tau": 10.0, "tau0": 3.0, "tau_r": {}, "vary": list(best_global["groups"])}, None),
        ("önceki varsayılan model + bölgesel en iyi ayar", sr.default_config(groups),
         dict(best_rc, vary=[g for g in best_rc["vary"] if g in groups]), None),
    ]
    results, rows = {}, []
    for name, cfg, rc, s in top + compare:
        rs = rg.evaluate(data, cfg, rc, sy, reg)["score"] if s is None else s
        rt = rg.evaluate(data, cfg, rc, ty, reg)
        results[name] = rt
        row = {"model": name, "arama puanı": rs, "sınama puanı": rt["score"]}
        for m in data.mags:
            lo, hi, p = sr.year_block_ci(rt["per_year"][m])
            row[f"M≥{m:g} kazanç"] = rt["pooled"][m]
            row[f"M≥{m:g} %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
        comb = [np.concatenate([rt["per_year"][m][i] for m in data.score_mags]) for i in range(len(ty))]
        lo, hi, p = sr.year_block_ci(comb)
        row["ETAS'a göre %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
        row["ETAS'a göre p"] = p
        row["rcfg"] = json.dumps(rc, ensure_ascii=False) if rc else ""
        rows.append(row)
    gname = compare[1][0]
    g = results[gname]
    for row in rows:
        rt = results[row["model"]]
        comb_a = [np.concatenate([rt["per_year"][m][i] for m in data.score_mags]) for i in range(len(ty))]
        comb_b = [np.concatenate([g["per_year"][m][i] for m in data.score_mags]) for i in range(len(ty))]
        lo, hi, p = rg.paired_ci(comb_a, comb_b)
        row["ortak ayara göre fark"] = row["sınama puanı"] - g["score"]
        row["ortak ayara göre %95 GA"] = f"[{lo:.3f}, {hi:.3f}]"
        row["ortak ayara göre p"] = p
    for key in ("ETAS'a göre p", "ortak ayara göre p"):
        idx = [i for i, r in enumerate(rows) if not (key == "ortak ayara göre p" and r["model"] == gname)]
        adj = sr.holm([rows[i][key] for i in idx])
        for i, a in zip(idx, adj):
            rows[i][key.replace(" p", " Holm p")] = float(a)
    summary["test_table"] = rows
    summary["region_breakdown"] = {
        "ortak ayar": rg.region_breakdown(g, data.mags),
        "bölgesel en iyi": rg.region_breakdown(results[top[0][0]], data.mags),
        "havuzlamasız": rg.region_breakdown(results[compare[3][0]], data.mags)}
    # sahte bölgelerin sınama yıllarındaki kazancı (şans kontrolünün sınama tarafı)
    fake_test = []
    for k, rc_s in enumerate(summary["null_rcfgs"]):
        fr = rg.fake_regions(reg, k)
        fake_test.append(rg.evaluate(data, best_global, json.loads(rc_s), ty, fr)["score"] - g["score"])
    summary["null_test_gain_vs_global"] = fake_test
    _, coef = rg.fit_predict(data, best_global, best_rc, ty[-1], reg, return_coef=True)
    summary["best_coef_last_year"] = coef
    summary["best_rcfg"] = best_rc
    (out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float))
    log.info("bitti: %s", out)


if __name__ == "__main__":
    main()
