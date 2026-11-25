"""Kısa vadeli öncül sinyal hipotezlerinin testi (iyonosfer, ULF, termal).

Soru: "Bir bölgede son 10 günde anomali görüldüyse, sonraki 5 günde o
bölgede deprem olasılığı, depremlerin kendi kümelenmesinden beklenenden
yüksek mi?"

Birimler: 2° × 2° kutular × 5 günlük pencereler. Hedef: kümelenmeden
ayıklanmış (Gardner-Knopoff) bağımsız M ≥ ``target_mag`` depremler. Böylece
artçı dizileri sinyali şişiremez.

Sıfır (null) model: kutunun eğitim yıllarındaki arka plan oranı (ofset) +
son 30 ve 5 gündeki M≥3 sayıları (kısa vadeli kümelenme). Alternatif:
null + anomali özelliği, öncül = ``registry`` başlangıç etkisi.
Eğitim 2012-2019, sınama 2020-2025 (modele hiç gösterilmez).

Ölçütler:
* sınamada log-olabilirlik farkı ve olay başına bilgi kazancı (p-değeri),
* alarm olasılık kazancı: eğitimde %90'lık eşiği aşan anomali sonrası
  pencerelerdeki olay payı / alarm altındaki uzay-zaman payı.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import stats

from . import fit as pfit
from .registry import OBS_FEATURES, apply_transform

log = logging.getLogger(__name__)

LEAD_DAYS = 10
WIN_DAYS = 5


def make_boxes(region, step=2.0):
    lats = np.arange(region.min_lat, region.max_lat - 1e-9, step)
    lons = np.arange(region.min_lon, region.max_lon - 1e-9, step)
    b = pd.DataFrame([(a, o) for a in lats for o in lons], columns=["lat0", "lon0"])
    b["lat_c"], b["lon_c"] = b["lat0"] + step / 2, b["lon0"] + step / 2
    b["box"] = np.arange(len(b))
    b.attrs["step"] = step
    return b


def _box_of(lat, lon, region, boxes):
    step = boxes.attrs["step"]
    nlon = int(np.ceil((region.max_lon - region.min_lon) / step - 1e-9))
    i = np.floor((np.asarray(lat) - region.min_lat) / step).astype(int)
    j = np.floor((np.asarray(lon) - region.min_lon) / step).astype(int)
    ok = (i >= 0) & (j >= 0) & (j < nlon) & (i * nlon + j < len(boxes))
    return np.where(ok, i * nlon + j, -1)


def _windows(start, end):
    return pd.date_range(start, end, freq=f"{WIN_DAYS}D")


def daily_to_box(df: pd.DataFrame, col: str, boxes: pd.DataFrame, max_km: float) -> pd.DataFrame:
    """Gözlem düğümlerinden kutu merkezine: ``max_km`` içindeki düğümlerin en büyüğü (günlük)."""
    nodes = df[["lat", "lon"]].drop_duplicates().to_numpy()
    dx = (boxes["lon_c"].to_numpy()[:, None] - nodes[:, 1][None]) * 85.0
    dy = (boxes["lat_c"].to_numpy()[:, None] - nodes[:, 0][None]) * 111.0
    near = np.hypot(dx, dy) <= max_km
    piv = df.pivot_table(index="date", columns=["lat", "lon"], values=col, aggfunc="max")
    piv = piv.reindex(columns=pd.MultiIndex.from_arrays([nodes[:, 0], nodes[:, 1]]))
    A = piv.to_numpy()
    out = np.full((len(piv), len(boxes)), np.nan)
    for b in range(len(boxes)):
        sel = near[b]
        if sel.any():
            with np.errstate(all="ignore"):
                out[:, b] = np.nanmax(A[:, sel], axis=1) if col != "tec_z_min" else np.nanmin(A[:, sel], axis=1)
    return pd.DataFrame(out, index=pd.to_datetime(piv.index), columns=boxes["box"])


def lead_feature(daily_box: pd.DataFrame, windows, agg="max") -> np.ndarray:
    """Her pencere başlangıcından önceki ``LEAD_DAYS`` günün özeti (pencere günleri hariç)."""
    out = np.full((len(windows), daily_box.shape[1]), np.nan)
    D = daily_box.sort_index()
    for k, w in enumerate(windows):
        s = D.loc[w - pd.Timedelta(days=LEAD_DAYS): w - pd.Timedelta(days=1)]
        if len(s) >= LEAD_DAYS // 2:
            with np.errstate(all="ignore"):
                out[k] = (s.max() if agg == "max" else s.min() if agg == "min" else s.mean()).to_numpy()
    return out


def build_panel(cat: pd.DataFrame, region, boxes, start, end, target_mag: float) -> pd.DataFrame:
    win = _windows(start, end)
    cat = cat.copy()
    cat["box"] = _box_of(cat["latitude"], cat["longitude"], region, boxes)
    cat = cat[cat["box"] >= 0]
    rows = []
    t = cat["time"].dt.tz_localize(None) if cat["time"].dt.tz is not None else cat["time"]
    cat = cat.assign(tn=t)
    m3 = cat[cat["mag"] >= 3.0]
    tgt = cat[(cat["mag"] >= target_mag) & cat["independent"]]
    for w in win:
        w1 = w + pd.Timedelta(days=WIN_DAYS)
        y = tgt[(tgt["tn"] >= w) & (tgt["tn"] < w1)].groupby("box").size()
        n30 = m3[(m3["tn"] >= w - pd.Timedelta(days=30)) & (m3["tn"] < w)].groupby("box").size()
        n5 = m3[(m3["tn"] >= w - pd.Timedelta(days=5)) & (m3["tn"] < w)].groupby("box").size()
        d = pd.DataFrame({"box": boxes["box"], "win": w})
        d["y"] = d["box"].map(y).fillna(0).to_numpy()
        d["n30"] = d["box"].map(n30).fillna(0).to_numpy()
        d["n5"] = d["box"].map(n5).fillna(0).to_numpy()
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def _ig(lam_a, lam_b, y):
    n = y.sum()
    if n == 0:
        return np.nan, np.nan
    idx = np.repeat(np.arange(len(y)), y.astype(int))
    per = np.log(lam_a[idx]) - np.log(lam_b[idx]) - (lam_a.sum() - lam_b.sum()) / n
    p = stats.ttest_1samp(per, 0).pvalue if n > 1 else np.nan
    return float(per.mean()), float(p)


def evaluate(panel: pd.DataFrame, feats: dict[str, np.ndarray], train_end: pd.Timestamp,
             test_end: pd.Timestamp, alarm_q: float = 0.9) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Her hipotez özelliği tek başına ve birlikte: null'a karşı sınama sonuçları."""
    tr = (panel["win"] < train_end).to_numpy()
    te = ((panel["win"] >= train_end) & (panel["win"] < test_end)).to_numpy()
    # arka plan: eğitimdeki kutu oranı (Laplace düzeltmeli)
    nwin_tr = panel.loc[tr, "win"].nunique()
    bg = panel[tr].groupby("box")["y"].sum().reindex(panel["box"].unique()).fillna(0)
    rate = (bg + 0.1) / max(nwin_tr, 1)
    off = np.log(panel["box"].map(rate).to_numpy())
    base = np.column_stack([np.log1p(panel["n30"]), np.log1p(panel["n5"])])
    y = panel["y"].to_numpy()
    results, coefs = [], []
    sets = {k: [k] for k in feats}
    groups = {}
    for k in feats:
        groups.setdefault(OBS_FEATURES[k].group, []).append(k)
    sets.update({f"[{g}] hepsi": v for g, v in groups.items() if len(v) > 1})
    tec_thermal = [k for k in feats if OBS_FEATURES[k].group in ("İyonosfer TEC", "Termal")]
    if tec_thermal:
        sets["[TEC + termal] hepsi"] = tec_thermal
    for name, cols in sets.items():
        F = np.column_stack([apply_transform(c, feats[c]) for c in cols])
        avail = np.isfinite(F).all(1)
        rtr, rte = tr & avail, te & avail
        if rtr.sum() < 100 or rte.sum() < 50 or y[rtr].sum() < 5:
            results.append({"hipotez": name, "durum": "yetersiz veri"})
            continue
        Bz_tr, Bz_te = pfit.standardize(base[rtr], base[rte])
        Fz_tr, Fz_te = pfit.standardize(F[rtr], F[rte])
        b0, _, _ = pfit.fit(Bz_tr, off[rtr], y[rtr], np.zeros(2), np.ones(2))
        pm = np.array([OBS_FEATURES[c].prior_mean for c in cols])
        ps = np.array([OBS_FEATURES[c].prior_sd for c in cols])
        b1, s1, _ = pfit.fit(np.column_stack([Bz_tr, Fz_tr]), off[rtr], y[rtr],
                             np.r_[0.0, 0.0, pm], np.r_[1.0, 1.0, ps])
        lam0 = pfit.predict(Bz_te, off[rte], b0)
        lam1 = pfit.predict(np.column_stack([Bz_te, Fz_te]), off[rte], b1)
        # yalnız öncülle (eğitimsiz) tahmin: başlangıç etkisinin kendisi ne kadar iyi?
        lamp = pfit.predict(np.column_stack([Bz_te, Fz_te]), off[rte], np.r_[b0, pm])
        yt = y[rte]
        ig, p = _ig(lam1, lam0, yt)
        igp, pp = _ig(lamp, lam0, yt)
        row = {"hipotez": name, "durum": "test edildi", "eğitim_olay": int(y[rtr].sum()),
               "sınama_olay": int(yt.sum()),
               "logL_null": round(float((yt * np.log(lam0) - lam0).sum()), 2),
               "logL_hipotez": round(float((yt * np.log(lam1) - lam1).sum()), 2),
               "bilgi_kazancı_olay_başına": round(ig, 4), "p_değeri": p,
               "yalnız_öncül_bilgi_kazancı": round(igp, 4)}
        if len(cols) == 1:
            thr = np.nanquantile(F[rtr, 0], alarm_q) if OBS_FEATURES[cols[0]].prior_mean >= 0 \
                else np.nanquantile(F[rtr, 0], 1 - alarm_q)
            alarm = F[rte, 0] >= thr if OBS_FEATURES[cols[0]].prior_mean >= 0 else F[rte, 0] <= thr
            frac_t = alarm.mean()
            frac_e = (yt * alarm).sum() / max(yt.sum(), 1)
            row["alarm_olasılık_kazancı"] = round(frac_e / frac_t, 2) if frac_t > 0 else None
            # binom testi: alarm altındaki olay payı şans eseri mi?
            row["alarm_p"] = float(stats.binomtest(int((yt * alarm).sum()), int(yt.sum()), frac_t,
                                                   alternative="greater").pvalue) if yt.sum() else None
        results.append(row)
        for c, mp, sp, bp, sd in zip(cols, pm, ps, b1[3:], s1[3:]):
            coefs.append({"set": name, "özellik": c, "öncül_ort": mp, "öncül_sd": sp,
                          "öğrenilen": round(float(bp), 3), "sonsal_sd": round(float(sd), 3)})
    return pd.DataFrame(results), pd.DataFrame(coefs)
