"""Öğrenme için (hücre × zaman penceresi) özellik tablosu.

Her satır bir mekân hücresi ve bir tahmin penceresidir (ör. 0,5° × 1 yıl).
Özellikler **yalnızca pencere başlangıcından önceki** veriden hesaplanır;
etiket ise pencere içinde olan biter. Bu, gelecekten sızıntıyı (leakage)
yapısal olarak engeller. Yeni değişkenler ``FEATURES`` listesine bir
fonksiyon eklenerek tanımlanır; her biri aynı imzayı kullanır, böylece
ileride yüzlerce değişken eklenip tek tek geçmiş veriyle sınanabilir.

Değerlendirme kuralı: model, ``time_split`` ile ayrılan eğitim yıllarında
öğrenir; sınama yılları ona hiç gösterilmez. Sonuçlar ``deprem.forecast``
içindeki referans modellerle (düzleştirilmiş sismisite, ileride ETAS)
aynı olabilirlik ölçütleriyle karşılaştırılmalıdır.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd

from ..faults.moment import mw_to_m0
from ..geo import point_to_polyline_km
from ..region import Region
from ..stats.gr import b_value_mle


def make_cells(region: Region, cell_deg: float) -> pd.DataFrame:
    lats = np.arange(region.min_lat, region.max_lat - 1e-9, cell_deg)
    lons = np.arange(region.min_lon, region.max_lon - 1e-9, cell_deg)
    LA, LO = np.meshgrid(lats, lons, indexing="ij")
    cells = pd.DataFrame({"lat0": LA.ravel(), "lon0": LO.ravel()})
    cells["lat_c"] = cells["lat0"] + cell_deg / 2
    cells["lon_c"] = cells["lon0"] + cell_deg / 2
    cells["cell_id"] = np.arange(len(cells))
    cells.attrs["cell_deg"] = cell_deg
    return cells


def assign_cell(events: pd.DataFrame, region: Region, cells: pd.DataFrame) -> np.ndarray:
    d = cells.attrs["cell_deg"]
    n_lon = int(round((cells["lon0"].max() - region.min_lon) / d)) + 1
    i = np.floor((events["latitude"].to_numpy() - region.min_lat) / d).astype(int)
    j = np.floor((events["longitude"].to_numpy() - region.min_lon) / d).astype(int)
    n_lat = len(cells) // n_lon
    ok = (i >= 0) & (i < n_lat) & (j >= 0) & (j < n_lon)
    out = np.full(len(events), -1)
    out[ok] = i[ok] * n_lon + j[ok]
    return out


# ---------------------------------------------------------------- özellik kaydı
# İmza: f(past_events, cells, t0, ctx) -> DataFrame (index = cell_id)
FeatureFn = Callable[[pd.DataFrame, pd.DataFrame, pd.Timestamp, dict], pd.DataFrame]
FEATURE_REGISTRY: dict[str, FeatureFn] = {}


def register_feature(name: str):
    """Yeni özellik grubu eklemek için dekoratör::

        @register_feature("benim_ozelligim")
        def f(past, cells, t0, ctx):
            return pd.DataFrame({...}).reindex(cells["cell_id"])
    """
    def deco(fn):
        FEATURE_REGISTRY[name] = fn
        return fn
    return deco


@register_feature("counts")
def f_counts(past, cells, t0, ctx):
    out = {}
    mc = ctx["mc"]
    p = past[past["mag"] >= mc]
    for yrs in (1, 5, 10, 30):
        w = p[p["time"] >= t0 - pd.DateOffset(years=yrs)]
        out[f"n_mc_{yrs}y"] = w.groupby("cell_id").size()
    w5 = past[past["mag"] >= 5.0]
    out["n_m5_all"] = w5.groupby("cell_id").size()
    return pd.DataFrame(out).reindex(cells["cell_id"]).fillna(0)


@register_feature("moment_release")
def f_moment_release(past, cells, t0, ctx):
    out = {}
    for yrs in (5, 30):
        w = past[past["time"] >= t0 - pd.DateOffset(years=yrs)]
        out[f"log_m0_{yrs}y"] = np.log10(
            w.assign(m0=mw_to_m0(w["mw_best"])).groupby("cell_id")["m0"].sum())
    out["max_mag_all"] = past.groupby("cell_id")["mag"].max()
    return pd.DataFrame(out).reindex(cells["cell_id"])


@register_feature("time_since")
def f_time_since(past, cells, t0, ctx):
    out = {}
    for m in (4.0, 5.0, 6.0):
        last = past[past["mag"] >= m].groupby("cell_id")["time"].max()
        out[f"yrs_since_m{m:g}"] = (t0 - last).dt.days / 365.25
    return pd.DataFrame(out).reindex(cells["cell_id"])


@register_feature("b_value")
def f_b_value(past, cells, t0, ctx, min_n: int = 30):
    mc = ctx["mc"]
    w = past[(past["time"] >= t0 - pd.DateOffset(years=10)) & (past["mag"] >= mc)]
    b = w.groupby("cell_id")["mag"].apply(
        lambda m: b_value_mle(m, mc)[0] if len(m) >= min_n else np.nan)
    return pd.DataFrame({"b_10y": b}).reindex(cells["cell_id"])


@register_feature("faults")
def f_faults(past, cells, t0, ctx):
    """Zamandan bağımsız fay özellikleri (hücre merkezine 25 km içindeki faylar)."""
    if "fault_static" not in ctx:
        faults = ctx.get("faults")
        if faults is None or faults.empty:
            ctx["fault_static"] = pd.DataFrame(index=cells["cell_id"])
        else:
            rows = []
            lat, lon = cells["lat_c"].to_numpy(), cells["lon_c"].to_numpy()
            dist = np.column_stack([point_to_polyline_km(lat, lon, c) for c in faults["coords"]])
            near = dist <= 25.0
            rate = faults["m0_rate_Nm_yr"].to_numpy()
            sr = faults["slip_rate_mm_yr"].to_numpy()
            for k in range(len(cells)):
                sel = near[k]
                rows.append({
                    "fault_dist_km": dist[k].min(),
                    "fault_log_m0_rate": np.log10(rate[sel].sum()) if sel.any() and rate[sel].sum() > 0 else np.nan,
                    "fault_max_slip_mm_yr": sr[sel].max() if sel.any() else 0.0,
                    "fault_n": int(sel.sum()),
                })
            ctx["fault_static"] = pd.DataFrame(rows, index=cells["cell_id"])
    return ctx["fault_static"]


@register_feature("rate_change")
def f_rate_change(past, cells, t0, ctx):
    """Sismik sessizlik / hızlanma: son 1 yıl hızı vs önceki 9 yıl (Poisson z)."""
    p = past[past["mag"] >= ctx["mc"]]
    n1 = p[p["time"] >= t0 - pd.DateOffset(years=1)].groupby("cell_id").size()
    n10 = p[p["time"] >= t0 - pd.DateOffset(years=10)].groupby("cell_id").size()
    df = pd.DataFrame({"n1": n1, "n10": n10}).reindex(cells["cell_id"]).fillna(0)
    prev = df["n10"] - df["n1"]
    z = (df["n1"] - prev / 9) / np.sqrt(df["n1"] + prev / 81 + 1e-9)
    return pd.DataFrame({"rate_z_1y": z.where(df["n10"] > 0)})


@register_feature("depth")
def f_depth(past, cells, t0, ctx):
    w = past[past["time"] >= t0 - pd.DateOffset(years=10)]
    g = w.groupby("cell_id")["depth_km"]
    return pd.DataFrame({"depth_mean_10y": g.mean(), "depth_p90_10y": g.quantile(0.9)}
                        ).reindex(cells["cell_id"])


@register_feature("clustering")
def f_clustering(past, cells, t0, ctx):
    """Bağımlı (artçı/öncü/sürü) olay oranı ve ana şok başına artçı sayısı."""
    if "independent" not in past.columns:
        return pd.DataFrame(index=cells["cell_id"])
    w = past[past["time"] >= t0 - pd.DateOffset(years=10)]
    g = w.groupby("cell_id")
    frac = 1 - g["independent"].mean()
    big = w[w["independent"] & (w["mag"] >= 4)].groupby("cell_id").size()
    dep = w[~w["independent"]].groupby("cell_id").size()
    prod = (dep / big).where(big > 0)
    return pd.DataFrame({"frac_dependent_10y": frac, "aftershocks_per_m4_10y": prod}
                        ).reindex(cells["cell_id"])


FEATURES: list[FeatureFn] = list(FEATURE_REGISTRY.values())


def build_table(catalog: pd.DataFrame, region: Region, cell_deg: float,
                window_starts: list[pd.Timestamp], window: pd.DateOffset,
                target_mags, mc: float, faults: pd.DataFrame | None = None,
                features: list[FeatureFn] | None = None,
                ctx_extra: dict | None = None) -> pd.DataFrame:
    """Hücre × pencere tablosu: özellikler + etiketler.

    ``target_mags`` bir veya birkaç büyüklük eşiğidir (ör. [3, 4, 5]). Her
    eşik m için etiketler: ``y_count_m{m}`` (pencere içi M>=m olay sayısı) ve
    ``y_any_m{m}`` (en az bir olay oldu mu). Ayrıca ``y_max_mag``.
    Küçük eşiklerde kalibre edilen bir modelin büyük eşiklere ne kadar
    taşındığı, aynı tabloda eşik eşik ölçülebilir
    (bkz. ``deprem.features.evaluate``).
    """
    target_mags = np.atleast_1d(target_mags).astype(float)
    features = features or list(FEATURE_REGISTRY.values())
    cells = make_cells(region, cell_deg)
    cat = catalog.copy()
    cat["cell_id"] = assign_cell(cat, region, cells)
    cat = cat[cat["cell_id"] >= 0]
    ctx = {"mc": mc, "faults": faults, **(ctx_extra or {})}
    parts = []
    for t0 in window_starts:
        t1 = t0 + window
        past = cat[cat["time"] < t0]
        fut = cat[(cat["time"] >= t0) & (cat["time"] < t1)]
        feats = pd.concat([fn(past, cells, t0, ctx) for fn in features], axis=1)
        for m in target_mags:
            cnt = fut[fut["mag"] >= m].groupby("cell_id").size()
            cnt = cnt.reindex(cells["cell_id"]).fillna(0).to_numpy()
            feats[f"y_count_m{m:g}"] = cnt
            feats[f"y_any_m{m:g}"] = (cnt > 0).astype(int)
        feats["y_max_mag"] = fut.groupby("cell_id")["mag"].max().reindex(cells["cell_id"]).to_numpy()
        feats.insert(0, "window_start", t0)
        feats.insert(1, "lat_c", cells["lat_c"].to_numpy())
        feats.insert(2, "lon_c", cells["lon_c"].to_numpy())
        parts.append(feats.reset_index().rename(columns={"index": "cell_id"}))
    return pd.concat(parts, ignore_index=True)


def time_split(table: pd.DataFrame, test_start: pd.Timestamp):
    """Zamana göre ayır: eğitim < test_start <= sınama. Rastgele bölme yasak."""
    train = table[table["window_start"] < test_start]
    test = table[table["window_start"] >= test_start]
    return train, test


# Gaz/sıvı ve dış etki özellikleri kendilerini yukarıdaki kayda ekler
# (bkz. deprem.external). ``ctx_extra`` ile uzak depremler, saha listesi vb. verilir.
from .. import external as _external  # noqa: E402,F401
