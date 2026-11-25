"""Parametre ve model ayarı kombinasyonu araması (aşırı uyuma karşı korumalı).

Soru: hangi parametre grupları, hangi önem ağırlıklarıyla ve hangi model
ayarlarıyla bir sonraki yılın depremlerini (yer + sayı) ETAS'tan daha iyi tahmin eder?

Kurallar (ön kayıt):

1. **Arama yalnız arama döneminde** yapılır (``search_years``, varsayılan 2015-2019).
   Her arama yılı Y için model yalnız Y öncesi yıllarla eğitilir ve Y tahmin edilir
   (ileriye kayan). Puan: olay başına bilgi kazancı (ETAS'a göre, M≥3/4/5 ortalaması).
2. **Son sınama dönemi** (``test_years``, 2020-2025) aramada hiç kullanılmaz. Arama
   bittikten sonra en iyi tek kombinasyon (önceden seçilmiş) burada bir kez sınanır.
3. **Şans eseri iyi görünme**: aynı arama, özelliklerin hücreler arasında rastgele
   karıştırıldığı (bilgisi yok edilmiş) verilerde de yapılır. Gerçek aramanın en iyi
   puanı, karıştırılmış aramaların en iyi puanlarını geçmiyorsa arama dönemi puanı
   anlamsızdır.
4. **Çoklu test**: sınamaya giden ilk K kombinasyon için Holm düzeltmesi.

Model: λ(hücre, Y) = ETAS(hücre, Y) · exp(β0 + Σ β_k z_k), β_k ~ N(m_k, (w_g · s)²).
w_g, grubun "önem ağırlığı"dır (0'a yakın = katkısı kısılmış). ``ard`` seçeneği
ağırlıkları veriden belirler (iki aşamalı uyarlamalı ridge: ilk uydurmada etkisi küçük
çıkan parametrelerin katkısı ikinci uydurmada kısılır).
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

from .observations import fit as pfit
from .observations.registry import OBS_FEATURES, prior_for

log = logging.getLogger(__name__)

# ------------------------------------------------------------------ sütun → grup
CATALOG_GROUPS = [
    ("sayılar", r"^n_mc_|^n_m5_all$"),
    ("moment boşalımı", r"^log_m0_|^max_mag_all$"),
    ("son depremden geçen süre", r"^yrs_since_"),
    ("b-değeri", r"^b_10y$"),
    ("fay geometrisi", r"^fault_"),
    ("hız değişimi", r"^rate_|_rate_z$|^z_rate"),
    ("derinlik", r"depth"),
    ("kümelenme", r"^frac_dep|^aft|cluster|^n_aftershock|dependent"),
    ("gelgit ve periyodik yükler", r"tidal|seasonal|hour_of_day"),
    ("sürüler ve akışkan", r"^swarm|hydraulic"),
    ("uzak tetikleme", r"^dyn_"),
    ("insan kaynaklı sahalar", r"^site_|reservoir"),
    ("akışkan (statik)", r"heat_flow|thermal_spring|radon|groundwater"),
    ("Coulomb gerilme aktarımı", r"^dcff_"),
    ("kayma açığı (enerji birikimi)", r"^slip_deficit"),
    ("kapasite ve kayaç", r"^cap_log_m0$|^rigidity_gpa$|^moho_km$"),
    ("doluluk oranı", r"^fill_"),
    ("doluluk + Coulomb", r"^fillc_"),
    ("zincirleme kırılma", r"^cascade_"),
    ("GPS kilitlenme / krip", r"^locking_mean$|^creep_mm_yr$"),
    ("sürtünme vekili", r"^fric_"),
]
ENERGY_GROUPS = {"kayma açığı (enerji birikimi)", "doluluk oranı", "doluluk + Coulomb", "zincirleme kırılma",
                 "kapasite ve kayaç", "Coulomb gerilme aktarımı"}
SMALL_EFFECT_GROUPS = {"gelgit ve periyodik yükler", "uzak tetikleme"}


def column_group(c: str) -> str:
    base = re.sub(r"__(var|eksik)$", "", c)
    if base in OBS_FEATURES:
        return OBS_FEATURES[base].group
    for g, pat in CATALOG_GROUPS:
        if re.search(pat, base):
            return g
    return "diğer"


# ------------------------------------------------------------------ veri
@dataclass
class SearchData:
    tab: pd.DataFrame                 # yıl, hücre sıralı
    off_test: dict                    # Y -> log ETAS (train_mag) dizisi (o yılın hücreleri)
    off_train: dict                   # (Y, yp) -> log ETAS dizisi
    b: dict                           # Y -> b
    mags: list
    train_mag: float
    feature_cols: list                # μ'den bağımsız sütunlar
    dcff_by_mu: dict                  # μ' -> sütun listesi
    years: list = field(default_factory=list)
    score_mags: tuple = (3.0, 4.0, 5.0)
    refit_every: int = 1              # katsayılar kaç yılda bir yeniden eğitilir (uzun dönem için 5)


def add_friction_proxy(tab: pd.DataFrame) -> pd.DataFrame:
    """Sürtünme davranışı vekili: yüksek mikrodeprem üretkenliği (biriken momente göre),
    yüksek b-değeri ve tekrarlayan deprem aileleri hız-güçlenen (sürünmeye yatkın, daha az
    enerji biriktiren) fay davranışını gösterir (Scholz 1998; Schorlemmer vd. 2005)."""
    parts = []
    if "n_mc_10y" in tab and "fault_log_m0_rate" in tab:
        r = np.log10(tab["n_mc_10y"] + 1) - tab["fault_log_m0_rate"].fillna(tab["fault_log_m0_rate"].min())
        parts.append(r)
    if "b_10y" in tab:
        parts.append(tab["b_10y"])
    rep = [c for c in tab.columns if c.startswith("rep_") and "n_fam" in c]
    if rep:
        parts.append(np.log1p(tab[rep[0]].fillna(0)))
    if not parts:
        return tab
    z = [(p - p.groupby(tab["year"]).transform("mean")) / p.groupby(tab["year"]).transform("std") for p in parts]
    tab["fric_velocity_strengthening"] = pd.concat(z, axis=1).mean(axis=1, skipna=True)
    return tab


def load(run) -> SearchData:
    """Koşu klasöründen özellik tablosu ve ileriye kayan ETAS ofsetlerini yükle."""
    from .features.table import make_cells
    from .region import TURKEY
    wf = run._wf_etas()
    comp = run._train_offset_mode() == "compensator"
    wtr = np.load(run.dir / ("walk_etas_train_comp.npz" if comp else "walk_etas_train.npz"))
    tab = pd.read_csv(run.dir / "feature_table.csv.gz", parse_dates=["window_start"])
    tab = tab[tab["year"].isin(list(wf))].sort_values(["year", "cell_id"]).reset_index(drop=True)
    tab = add_friction_proxy(tab)
    n_cells = len(make_cells(TURKEY, run.cfg["data"]["cell_deg"]))
    mmin = run.grid.mag_edges[0]
    tm = run.train_mag
    off_test, off_train, bs = {}, {}, {}
    for Y in sorted(wf):
        v = wf[Y]
        b = float(v[-1])
        bs[Y] = b
        cid = tab.loc[tab["year"] == Y, "cell_id"].to_numpy()
        off_test[Y] = np.log(v[:-1][:n_cells][cid] * 10 ** (-b * (tm - mmin)) + 1e-9)
        for yp in sorted(wf):
            k = f"{Y}_{yp}"
            if yp < Y and k in wtr.files:
                cp_ = tab.loc[tab["year"] == yp, "cell_id"].to_numpy()
                off_train[(Y, yp)] = np.log(wtr[k][:n_cells][cp_] * 10 ** (-b * (tm - mmin)) + 1e-9)
    mu_any = run.cfg["coulomb"]["mu_eff_grid"][0]
    cols = run.feature_columns(tab, mu_any)
    base_cols = [c for c in cols if not c.startswith("dcff_")]
    if "fric_velocity_strengthening" in tab and "fric_velocity_strengthening" not in base_cols:
        base_cols.append("fric_velocity_strengthening")
    dcff = {mu: [f"dcff_max_mu{mu:g}", f"dcff_mean_mu{mu:g}", f"dcff_frac_pos_mu{mu:g}"]
            for mu in run.cfg["coulomb"]["mu_eff_grid"]}
    # sabit / tamamen boş sütunları at
    base_cols = [c for c in base_cols if tab[c].notna().any() and tab[c].nunique(dropna=True) > 1]
    k = tab[f"y_count_m{tm:g}"].sum() / max(tab["y_wcount_train"].sum(), 1e-9)
    tab["y_w"] = tab["y_wcount_train"] * k
    sc = run.cfg.get("search", {})
    return SearchData(tab=tab, off_test=off_test, off_train=off_train, b=bs, mags=list(run.mags),
                      train_mag=tm, feature_cols=base_cols, dcff_by_mu=dcff, years=sorted(wf),
                      score_mags=tuple(float(m) for m in sc.get("score_mags", run.mags[:3])),
                      refit_every=int(sc.get("refit_every", 1)))


def groups_of(data: SearchData) -> dict:
    g = {}
    for c in data.feature_cols + data.dcff_by_mu[next(iter(data.dcff_by_mu))]:
        g.setdefault(column_group(c), []).append(c)
    return g


# ------------------------------------------------------------------ tek kombinasyonun değerlendirilmesi
def design(tab, cols, st=None):
    X = tab[cols].astype(float)
    if st is None:
        miss = [c for c in cols if X[c].isna().any()]
        st = {"mean": X.mean().fillna(0), "std": X.std().replace(0, 1).fillna(1), "miss": miss,
              "sparse": {c: bool(X[c].isna().mean() > 0.5) for c in miss}}
    Z = ((X - st["mean"]) / st["std"]).clip(-5, 5).fillna(0.0)
    extra = {}
    for c in st["miss"]:
        extra[f"{c}__var" if st["sparse"][c] else f"{c}__eksik"] = (
            X[c].notna() if st["sparse"][c] else X[c].isna()).astype(float)
    if extra:
        Z = pd.concat([Z, pd.DataFrame(extra, index=Z.index)], axis=1)
    return Z, st


def priors(columns, cfg) -> tuple[np.ndarray, np.ndarray]:
    s = cfg["prior_sd"]
    w = cfg.get("weights", {})
    m, sd = [], []
    for c in columns:
        g = column_group(c)
        base = re.sub(r"__(var|eksik)$", "", c)
        p = prior_for(base) if base == c else None
        mean, s0 = (p if p else (0.0, s))
        if p:
            s0 = max(s0, 1e-3) * (s / 0.5)  # literatür sd'si de genel ölçekle birlikte
        if cfg.get("energy_prior"):
            if g in ENERGY_GROUPS and base == c:
                mean, s0 = mean + 0.15, s0 * 2.0
            elif g in SMALL_EFFECT_GROUPS:
                s0 = s0 * 0.3
        s0 = s0 * w.get(g, 1.0)
        if cfg.get("literature_prior", True) is False:
            mean = 0.0
        m.append(mean)
        sd.append(max(s0, 1e-4))
    return np.array(m), np.array(sd)


def selected_columns(data: SearchData, cfg) -> list:
    on = set(cfg["groups"])
    cols = [c for c in data.feature_cols if column_group(c) in on]
    if "Coulomb gerilme aktarımı" in on:
        cols += data.dcff_by_mu[cfg.get("mu_eff", 0.4)]
    return cols


def fit_predict(data: SearchData, cfg, Y: int, return_coef: bool = False, tab=None, fit_year=None,
                cache=None):
    """Y yılını tahmin et. Katsayılar ``fit_year`` (varsayılan Y) öncesi yıllarla, o yılın
    ETAS'ıyla eğitilir; ``cache`` verilirse aynı fit_year için tekrar kullanılır."""
    tab = data.tab if tab is None else tab
    cols = selected_columns(data, cfg)
    FY = Y if fit_year is None else fit_year
    win = cfg.get("window")
    te = tab[tab["year"] == Y]
    if cache is not None and FY in cache:
        beta, cols_z, st = cache[FY]
        if beta is None:
            lam = np.exp(data.off_test[Y] + cols_z)
            return (lam, {"(sabit)": float(cols_z)}) if return_coef else lam
        Zs, _ = design(te, cols, st)
        Zs = Zs.reindex(columns=cols_z, fill_value=0.0)
        lam = pfit.predict(Zs.to_numpy(), data.off_test[Y], beta)
        return (lam, dict(zip(["(sabit)"] + list(cols_z), beta.tolist()))) if return_coef else lam
    tr_mask = (tab["year"] < FY) & tab["year"].isin([yp for (yy, yp) in data.off_train if yy == FY])
    if win:
        tr_mask &= tab["year"] >= FY - win
    tr = tab[tr_mask]
    off_tr = np.concatenate([data.off_train[(FY, int(yp))] for yp in sorted(tr["year"].unique())])
    ycol = "y_w" if cfg.get("uncertainty_weights", True) else f"y_count_m{data.train_mag:g}"
    y = tr[ycol].to_numpy()
    if not cols:
        beta = np.array([np.log(max(y.sum(), 1e-9) / np.exp(off_tr).sum())])
        if cache is not None:
            cache[FY] = (None, float(beta[0]), None)
        lam = np.exp(data.off_test[Y] + beta[0])
        return (lam, {"(sabit)": float(beta[0])}) if return_coef else lam
    Z, st = design(tr, cols)
    pm, ps = priors(Z.columns, cfg)
    beta, sd, _ = pfit.fit(Z.to_numpy(), off_tr, y, pm, ps)
    if cfg.get("ard"):
        # uyarlamalı ridge: ilk uydurmada etkisi küçük parametrelerin katkısını kıs
        mag = np.abs(beta[1:] - pm)
        rel = np.clip(mag / max(np.median(mag[mag > 0]) if (mag > 0).any() else 1.0, 1e-6), 0.05, 3.0)
        beta, sd, _ = pfit.fit(Z.to_numpy(), off_tr, y, pm, ps * rel)
    if cache is not None:
        cache[FY] = (beta, list(Z.columns), st)
    Zs, _ = design(te, cols, st)
    Zs = Zs.reindex(columns=Z.columns, fill_value=0.0)
    lam = pfit.predict(Zs.to_numpy(), data.off_test[Y], beta)
    if return_coef:
        return lam, dict(zip(["(sabit)"] + list(Z.columns), beta.tolist()))
    return lam


def per_event_gains(data: SearchData, lam, Y, tab=None):
    tab = data.tab if tab is None else tab
    te = tab[tab["year"] == Y]
    base = np.exp(data.off_test[Y])
    out = {}
    for mg in data.mags:
        s = 10 ** (-data.b[Y] * (mg - data.train_mag))
        le, lm = base * s, lam * s
        y = te[f"y_count_m{mg:g}"].to_numpy()
        n = int(y.sum())
        if n == 0:
            out[mg] = np.array([])
            continue
        idx = np.repeat(np.arange(len(y)), y.astype(int))
        out[mg] = np.log(lm[idx]) - np.log(le[idx]) - (lm.sum() - le.sum()) / n
    return out


def evaluate(data: SearchData, cfg, years, tab=None, score_mags=None) -> dict:
    score_mags = data.score_mags if score_mags is None else score_mags
    per = {mg: [] for mg in data.mags}
    cache = {}
    k = max(int(getattr(data, "refit_every", 1)), 1)
    y0 = min(years)
    for Y in years:
        fy = Y - ((Y - y0) % k)  # blok başı: katsayılar k yılda bir, yalnız blok öncesiyle eğitilir
        lam = fit_predict(data, cfg, Y, tab=tab, fit_year=fy, cache=cache)
        for mg, v in per_event_gains(data, lam, Y, tab=tab).items():
            per[mg].append(v)
    pooled = {mg: (np.concatenate(v).mean() if sum(len(x) for x in v) else np.nan) for mg, v in per.items()}
    score = float(np.nanmean([pooled[m] for m in score_mags if m in pooled]))
    return {"score": score, "pooled": pooled, "per_year": per}


# ------------------------------------------------------------------ arama uzayı
def default_config(groups) -> dict:
    return {"groups": sorted(groups), "weights": {}, "prior_sd": 0.5, "mu_eff": 0.4, "ard": False,
            "energy_prior": False, "uncertainty_weights": True, "window": None}


def suggest(trial, groups, mu_grid) -> dict:
    on, w = [], {}
    for g in sorted(groups):
        if trial.suggest_categorical(f"g::{g}", [True, False]):
            on.append(g)
            w[g] = trial.suggest_float(f"w::{g}", 0.05, 3.0, log=True)
    return {"groups": on, "weights": w,
            "prior_sd": trial.suggest_float("prior_sd", 0.02, 1.5, log=True),
            "mu_eff": trial.suggest_categorical("mu_eff", list(mu_grid)),
            "ard": trial.suggest_categorical("ard", [False, True]),
            "energy_prior": trial.suggest_categorical("energy_prior", [False, True]),
            "uncertainty_weights": trial.suggest_categorical("uncertainty_weights", [True, False]),
            "window": trial.suggest_categorical("window", [0, 3, 5]) or None}


def run_search(data: SearchData, years, n_trials: int, seed: int = 0, tab=None, extra_configs=()):
    import optuna
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    groups = groups_of(data)
    sampler = optuna.samplers.TPESampler(seed=seed, n_startup_trials=min(60, n_trials // 3), multivariate=True)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    rows = []

    def obj(trial):
        cfg = suggest(trial, groups, list(data.dcff_by_mu))
        r = evaluate(data, cfg, years, tab=tab)
        rows.append({"trial": trial.number, "score": r["score"],
                     **{f"M≥{m:g}": r["pooled"][m] for m in data.mags}, "config": json.dumps(cfg, ensure_ascii=False)})
        return r["score"] if np.isfinite(r["score"]) else -9.0

    for c in extra_configs:  # bilinen başlangıç noktaları (ör. tüm parametreler, varsayılan ayar)
        params = {}
        for g in groups:
            params[f"g::{g}"] = g in c["groups"]
            if g in c["groups"]:
                params[f"w::{g}"] = float(np.clip(c["weights"].get(g, 1.0), 0.05, 3.0))
        params.update(prior_sd=c["prior_sd"], mu_eff=c["mu_eff"], ard=c["ard"], energy_prior=c["energy_prior"],
                      uncertainty_weights=c["uncertainty_weights"], window=c["window"] or 0)
        study.enqueue_trial(params)
    study.optimize(obj, n_trials=n_trials)
    return pd.DataFrame(rows).sort_values("score", ascending=False).reset_index(drop=True)


def shuffled_table(data: SearchData, seed: int) -> pd.DataFrame:
    """Özellikleri her yıl içinde hücreler arasında aynı permütasyonla karıştır (etiket,
    ETAS ofseti yerinde kalır). Özelliklerin konum bilgisi yok edilir, kendi aralarındaki
    ilişki korunur: bu verideki "en iyi" puan saf şans sınırıdır."""
    rng = np.random.default_rng(seed)
    tab = data.tab.copy()
    feat = sorted(set(data.feature_cols) | {c for v in data.dcff_by_mu.values() for c in v})
    for y, idx in tab.groupby("year").groups.items():
        idx = np.asarray(idx)
        perm = rng.permutation(idx)
        tab.loc[idx, feat] = tab.loc[perm, feat].to_numpy()
    return tab


def year_block_ci(per_year: list[np.ndarray], n_boot=4000, seed=1):
    rng = np.random.default_rng(seed)
    s = np.array([p.sum() for p in per_year])
    n = np.array([len(p) for p in per_year])
    ok = n > 0
    s, n = s[ok], n[ok]
    if len(s) < 2:
        return (np.nan, np.nan, np.nan)
    pick = rng.integers(0, len(s), (n_boot, len(s)))
    bs = s[pick].sum(1) / n[pick].sum(1)
    p_one = float((bs <= 0).mean())  # tek yönlü: kazanç ≤ 0 olasılığı (bootstrap)
    return float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975)), p_one


def holm(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    m = len(p)
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (m - k) * p[i]))
        adj[i] = run
    return adj
