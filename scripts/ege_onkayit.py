"""Ön kayıtlı Ege sınaması (docs/onkayit/2026-11-25_ege_ustunluk.md).

İki adım, sırayla ve ayrı commit'lerle:

    python scripts/train.py --config configs/train_ege_2026.toml --until features
    python scripts/ege_onkayit.py secim     # ortak ayarı seç + dondur (2026'ya bakmaz)
    git commit docs/onkayit/2026-11-25_ege_dondurulmus_ayar.json ...
    python scripts/ege_onkayit.py sinama    # 2026 puanı (bir kez)

``secim`` yalnız 2015-2019 (arama) ve 2020-2025 (PR #13 sonucunun yeniden üretimi)
yıllarını kullanır. ``sinama`` dondurulmuş ayarı okur, 2012-2025 ile bir kez eğitir ve
2026-01-01 → 2026-11-25 Ege olaylarında ETAS'la karşılaştırır.

HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

from deprem import regional as rg
from deprem import search as sr
from deprem.features.table import assign_cell, make_cells
from deprem.pipeline import Run
from deprem.region import TURKEY

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/train_ege_2026.toml"
FROZEN = ROOT / "docs/onkayit/2026-11-25_ege_dondurulmus_ayar.json"
RESULT = ROOT / "docs/onkayit/2026-11-25_ege_sonuc.json"
LIVE = ROOT / "forecasts/ege"
TEST_YEAR = 2026
EGE = rg.REGION_NAMES.index("Ege / Batı Anadolu")
N_BOOT, SEED = 10_000, 2026
log = logging.getLogger("ege")


# ------------------------------------------------------------------ ölçüt (ön kayıt bölüm 5)
def ig(n, lm, le):
    """Bölgeye kısıtlı olay başına bilgi kazancı (Poisson log-olabilirlik farkı / olay)."""
    n, lm, le = (np.asarray(a, float) for a in (n, lm, le))
    N = n.sum()
    if N == 0:
        return np.nan
    return float(((n * np.log(lm / le)).sum() - (lm.sum() - le.sum())) / N)


def cell_bootstrap(n, lm, le, n_boot=N_BOOT, seed=SEED):
    """Hücre blokları bootstrap'ı: %95 GA ve tek yönlü p = P(IG ≤ 0)."""
    n, lm, le = (np.asarray(a, float) for a in (n, lm, le))
    rng = np.random.default_rng(seed)
    a = n * np.log(lm / le) - (lm - le)
    pick = rng.integers(0, len(n), (n_boot, len(n)))
    num, den = a[pick].sum(1), n[pick].sum(1)
    ok = den > 0
    bs = num[ok] / den[ok]
    return float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975)), float((bs <= 0).mean())


def verdict(g, p, n):
    if n < 30:
        return "kararsız (30'dan az olay)"
    if g <= 0:
        return "kaldı"
    if p >= 0.05:
        return "kararsız"
    return "geçti (güçlü)" if g >= 0.10 else "geçti (zayıf)"


# ------------------------------------------------------------------ veri
def load_data(run):
    dp = run.dir / "search" / "data.pkl"
    dp.parent.mkdir(exist_ok=True)
    if not dp.exists():
        pickle.dump(sr.load(run), open(dp, "wb"))
    data = pickle.load(open(dp, "rb"))
    # y_w ölçeği yalnız 2026 öncesi yıllardan (sınama yılının sayıları eğitime sızmasın)
    tab = data.tab
    tm = data.train_mag
    old = tab["year"] < TEST_YEAR
    k = tab.loc[old, f"y_count_m{tm:g}"].sum() / max(tab.loc[old, "y_wcount_train"].sum(), 1e-9)
    tab["y_w"] = tab["y_wcount_train"] * k
    return data


def region_of(run):
    cells = make_cells(TURKEY, run.cfg["data"]["cell_deg"])
    return rg.cell_regions(cells)


def band_counts(run, year, lo, hi):
    """[lo, hi) büyüklük bandındaki olay sayısı, hücre başına (yıl penceresinde)."""
    cat = run.catalog
    cells = make_cells(TURKEY, run.cfg["data"]["cell_deg"])
    t0, t1 = run.window(year)
    sel = cat[(cat["time"] >= t0) & (cat["time"] < t1) & (cat["mag"] >= lo) & (cat["mag"] < hi)].copy()
    cid = assign_cell(sel, TURKEY, cells)
    cid = cid[cid >= 0]
    return np.bincount(cid, minlength=len(cells)).astype(float)


def year_arrays(data, cfg, Y, reg, run=None, cache=None):
    """Y yılı için hücre başına (n, λ_model, λ_ETAS) M≥3, M≥4 ve 2,5-3 bandı."""
    tab = data.tab[data.tab["year"] == Y]
    cid = tab["cell_id"].to_numpy()
    lam = sr.fit_predict(data, cfg, Y, cache=cache)
    base = np.exp(data.off_test[Y])
    b = data.b[Y]
    out = {}
    for mg in (3.0, 4.0):
        s = 10 ** (-b * (mg - data.train_mag))
        out[f"M≥{mg:g}"] = (tab[f"y_count_m{mg:g}"].to_numpy(float), lam * s, base * s)
    if run is not None:
        s = 10 ** (-b * (2.5 - data.train_mag)) - 1.0
        out["2,5≤M<3"] = (band_counts(run, Y, 2.5, 3.0)[cid], lam * s, base * s)
    return out, reg[cid], lam / base


def summarize(arrs, mask):
    res = {}
    for k, (n, lm, le) in arrs.items():
        n, lm, le = n[mask], lm[mask], le[mask]
        g = ig(n, lm, le)
        lo, hi, p = cell_bootstrap(n, lm, le) if n.sum() else (np.nan, np.nan, np.nan)
        res[k] = {"olay": int(n.sum()), "IG": g, "GA95": [lo, hi], "p_tek_yonlu": p,
                  "model_beklenen": float(lm.sum()), "etas_beklenen": float(le.sum())}
    return res


def old_metric(data, cfg, years, reg, cache=None):
    """PR #13'teki ölçüt: olay başına kazanç, sayı farkı terimi tüm Türkiye'den."""
    out = {}
    for mg in (3.0, 4.0):
        g_all = []
        for Y in years:
            tab = data.tab[data.tab["year"] == Y]
            lam = sr.fit_predict(data, cfg, Y, cache=cache)
            g = sr.per_event_gains(data, lam, Y)[mg]
            r = np.repeat(reg[tab["cell_id"].to_numpy()], tab[f"y_count_m{mg:g}"].to_numpy().astype(int))
            g_all.append(g[r == EGE])
        g_all = np.concatenate(g_all)
        out[f"M≥{mg:g}"] = {"olay": int(len(g_all)), "kazanç": float(g_all.mean()) if len(g_all) else np.nan}
    return out


# ------------------------------------------------------------------ adımlar
def secim(args):
    run = Run(CONFIG, ROOT)
    data = load_data(run)
    reg = region_of(run)
    sy = list(range(2015, 2020))
    groups = sr.groups_of(data)
    tp = run.dir / "search" / "trials.csv"
    t = time.time()
    if not tp.exists():
        sr.run_search(data, sy, args.trials, seed=0,
                      extra_configs=[sr.default_config(groups)]).to_csv(tp, index=False)
    trials = pd.read_csv(tp)
    best = json.loads(trials["config"].iloc[0])
    log.info("ortak ayar seçildi: arama puanı %.4f (%.0f s)", trials["score"].iloc[0], time.time() - t)
    # PR #13 Ege sonucunun yeniden üretimi (2020-2025; yeni veri değil)
    ty = list(range(2020, 2026))
    cache = {}
    rep_old = old_metric(data, best, ty, reg, cache={})
    per = {}
    for Y in ty:
        arrs, rr, _ = year_arrays(data, best, Y, reg, run=run, cache=cache)
        for k, v in arrs.items():
            per.setdefault(k, []).append(tuple(a[rr == EGE] for a in v))
    rep_new = {}
    for k, lst in per.items():
        n, lm, le = (np.concatenate([x[i] for x in lst]) for i in range(3))
        lo, hi, p = cell_bootstrap(n, lm, le)
        rep_new[k] = {"olay": int(n.sum()), "IG": ig(n, lm, le), "GA95": [lo, hi], "p_tek_yonlu": p}
    overall = sr.evaluate(data, best, ty)
    doc = {
        "on_kayit": "docs/onkayit/2026-11-25_ege_ustunluk.md",
        "olusturma": pd.Timestamp.now(tz="UTC").isoformat(),
        "not": "Bu ayar 2015-2019 arama yıllarıyla seçildi; 2026 verisi hiçbir hesapta kullanılmadı.",
        "arama": {"yillar": sy, "deneme": int(len(trials)), "tohum": 0, "puan": float(trials["score"].iloc[0])},
        "ayar": best,
        "yeniden_uretim_2020_2025": {
            "PR13_ege_M3_eski_olcut": 0.315,
            "eski_olcut_ege": rep_old,
            "bolgeye_kisitli_ege": rep_new,
            "tum_turkiye_sinama_puani": overall["score"],
        },
    }
    FROZEN.write_text(json.dumps(doc, indent=2, ensure_ascii=False, default=float))
    log.info("yazıldı: %s", FROZEN)
    print(json.dumps(doc["yeniden_uretim_2020_2025"], indent=2, ensure_ascii=False, default=float))


def sinama(args):
    frozen = json.loads(FROZEN.read_text())
    cfg = frozen["ayar"]
    run = Run(CONFIG, ROOT)
    data = load_data(run)
    reg = region_of(run)
    assert TEST_YEAR in data.off_test, "2026 ETAS'ı yok"
    arrs, rr, mult = year_arrays(data, cfg, TEST_YEAR, reg, run=run, cache={})
    t0, t1 = run.window(TEST_YEAR)
    prim = summarize(arrs, rr == EGE)
    rest = summarize({k: v for k, v in arrs.items() if k != "2,5≤M<3"}, rr != EGE)
    old = old_metric(data, cfg, [TEST_YEAR], reg, cache={})
    m3 = prim["M≥3"]
    res = {
        "on_kayit": frozen["on_kayit"],
        "hesaplama": pd.Timestamp.now(tz="UTC").isoformat(),
        "pencere": [t0.isoformat(), t1.isoformat()],
        "katalog_kaynaklari_2026": run.catalog[(run.catalog["time"] >= t0)]["source"].value_counts().to_dict()
        if "source" in run.catalog else {},
        "birincil": {"bolge": "Ege / Batı Anadolu", **m3,
                     "karar": verdict(m3["IG"], m3["p_tek_yonlu"], m3["olay"])},
        "ikincil": {"ege": {k: v for k, v in prim.items() if k != "M≥3"},
                    "ege_eski_olcut": old, "ege_disi_turkiye": rest},
        "etas_2026_b": data.b[TEST_YEAR],
    }
    RESULT.write_text(json.dumps(res, indent=2, ensure_ascii=False, default=float))
    log.info("yazıldı: %s", RESULT)
    # canlı test için dondurulmuş hücre çarpanları (model / ETAS, 2026)
    LIVE.mkdir(parents=True, exist_ok=True)
    tab = data.tab[data.tab["year"] == TEST_YEAR]
    pd.DataFrame({"cell_id": tab["cell_id"].to_numpy(), "lat_c": tab["lat_c"].to_numpy(),
                  "lon_c": tab["lon_c"].to_numpy(), "bolge": [rg.REGION_NAMES[i] for i in rr],
                  "model_carpani": mult}).to_csv(LIVE / "model_carpani_2026.csv", index=False, float_format="%.6g")
    pj = json.loads((run.dir / "walk_etas_params.json").read_text())[str(TEST_YEAR)]
    from deprem.longwalk import background_share
    from deprem.forecast.etas import ETAS
    c = run.cfg["etas"]
    e = ETAS(m0=c["m0"], bandwidth_km=c["bandwidth_km"], r_max_km=c["r_max_km"], dt_max_days=c["dt_max_days"])
    e.params, e.b = pj["params"], pj["b"]
    cat = run.catalog
    train = cat[(cat["mag"] >= e.m0) & (cat["time"] >= pd.Timestamp(c["fit_start"], tz="UTC"))
                & (cat["time"] < t0)]
    e._background(train)
    U = background_share(e, run.grid)
    (LIVE / "etas_2026.json").write_text(json.dumps(
        {"params": e.params, "b": e.b, "m0": e.m0, "dt_max_days": e.dt_max_days,
         "fit": [c["fit_start"], t0.isoformat()], "not": "2026 ileri yürüyen ETAS (ön kayıt bölüm 4)"},
        indent=2))
    pd.DataFrame({"cell_id": np.arange(len(U)), "arka_plan_payi": U}).to_csv(
        LIVE / "etas_arka_plan_2026.csv", index=False, float_format="%.6g")
    print(json.dumps(res, indent=2, ensure_ascii=False, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("adim", choices=["secim", "sinama"])
    ap.add_argument("--trials", type=int, default=300)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    {"secim": secim, "sinama": sinama}[args.adim](args)


if __name__ == "__main__":
    main()
