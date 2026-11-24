"""Tarihsel kataloğu ham dosyalardan üretir.

Kullanım (depo kökünden):
    python -m tarihsel.deprem_tarihsel.build
veya
    cd tarihsel && python -m deprem_tarihsel.build

Çıktılar ``tarihsel/data`` altında:
    tarihsel_tum_kayitlar.csv  Bütün kaynaklardaki bütün kayıtlar (yedek; hiçbir şey atılmaz)
    tarihsel_katalog.csv       Her olay için tek (tercih edilen) kayıt
    tamlik.json                Stepp (1972) tamlık analizi sonucu
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from .config import DEFAULT, Config
from .sources import READERS
from .timeutil import to_epoch_s

DATA = Path(__file__).resolve().parent.parent / "data"

SCHEMA_COLUMNS = [
    # aletsel katalogla ortak çekirdek
    "event_id", "time", "latitude", "longitude", "depth_km", "mag", "mag_type", "source", "agency",
    "time_err_s", "horiz_err_km", "depth_err_km", "mag_err",
    # tarihsel ekler
    "time_min", "time_max", "time_precision", "decimal_year",
    "mag_is_estimated", "mag_method", "mag_orig", "mag_orig_type",
    "intensity_io", "area", "country", "deaths", "reference", "source_event_id",
    "suspect", "notes",
    # yineleme
    "duplicate_group", "is_preferred", "n_sources", "group_agencies", "group_mag_min", "group_mag_max",
    # ağırlıklar
    "mag_err_used", "mag_err_assigned", "horiz_err_used", "horiz_err_assigned",
    "w_mag", "w_loc", "w_time", "w_unc", "in_completeness_window", "w_completeness",
    "w_suspect", "training_weight",
]


# ---------------------------------------------------------------------------
def load_all(cfg: Config = DEFAULT) -> pd.DataFrame:
    recs = []
    for name in cfg.source_priority:
        recs += READERS[name](cfg)
    df = pd.DataFrame(recs)
    # Tarihsel büyüklüklerin hepsi tahmindir (makrosismik / dönüşüm / uzman yargısı).
    df["mag_is_estimated"] = df["mag"].notna()
    return df


def _haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dphi = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def find_duplicates(df: pd.DataFrame, cfg: Config = DEFAULT) -> pd.DataFrame:
    """Farklı kataloglardaki aynı olayı gruplar. Bir grupta aynı katalogdan iki kayıt olmaz."""
    df = df.reset_index(drop=True).copy()
    tmin = to_epoch_s(df["time_min"])
    tmax = to_epoch_s(df["time_max"])
    tmid = (tmin + tmax) / 2
    lat, lon = df["latitude"].to_numpy(float), df["longitude"].to_numpy(float)
    mag = df["mag"].to_numpy(float)
    ag = df["agency"].to_numpy()
    slack = cfg.dup_time_slack_days * 86400

    pairs = []
    for i in range(len(df)):
        j = np.arange(i + 1, len(df))
        j = j[ag[j] != ag[i]]
        if not len(j):
            continue
        gap = np.maximum(0, np.maximum(tmin[j] - tmax[i], tmin[i] - tmax[j]))
        ok = gap <= slack
        d = _haversine_km(lat[i], lon[i], lat[j], lon[j])
        ok &= d <= cfg.dup_max_km
        dm = np.abs(mag[j] - mag[i])
        ok &= ~(dm > cfg.dup_max_dmag)  # NaN -> geçer
        for k in np.nonzero(ok)[0]:
            score = (d[k] / cfg.dup_max_km + abs(tmid[j[k]] - tmid[i]) / (365.25 * 86400)
                     + (0 if np.isnan(dm[k]) else dm[k] / cfg.dup_max_dmag))
            pairs.append((score, i, int(j[k])))
    pairs.sort()

    parent = list(range(len(df)))
    members = {i: {ag[i]} for i in range(len(df))}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for _, i, j in pairs:
        ri, rj = find(i), find(j)
        if ri == rj or members[ri] & members[rj]:
            continue
        parent[rj] = ri
        members[ri] |= members.pop(rj)

    roots = np.array([find(i) for i in range(len(df))])
    prio = {a: k for k, a in enumerate(cfg.source_priority)}
    df["_root"] = roots
    df["_prio"] = df["agency"].map(prio)
    df["_hasmag"] = df["mag"].notna()
    df["_bad"] = df["suspect"].astype(bool)
    order = df.sort_values(["_root", "_bad", "_hasmag", "_prio"], ascending=[True, True, False, True])
    preferred = set(order.groupby("_root").head(1).index)
    df["is_preferred"] = df.index.isin(preferred)
    # Grup kimliği: tercih edilen kaydın event_id'si
    gid = order.groupby("_root").head(1).set_index("_root")["event_id"]
    df["duplicate_group"] = df["_root"].map(gid)
    df["n_sources"] = df.groupby("_root")["agency"].transform("count")
    df["group_agencies"] = df.groupby("_root")["agency"].transform(lambda s: ",".join(sorted(s)))
    # Kaynaklar arası büyüklük tutarsızlığı: tercih edilen kaydın ağırlığına yansıtılır.
    df["group_mag_min"] = df.groupby("_root")["mag"].transform("min")
    df["group_mag_max"] = df.groupby("_root")["mag"].transform("max")
    return df.drop(columns=["_root", "_prio", "_hasmag", "_bad"])


# ---------------------------------------------------------------------------
def stepp_completeness(df: pd.DataFrame, cfg: Config = DEFAULT) -> dict:
    """Stepp (1972): her Mw alt sınırı için oranın kararlı kaldığı en uzun dönem.

    Bitiş yılından geriye doğru T uzunluğunda pencereler alınır; λ_T = N/T.
    Referans oran, en az ``stepp_min_events`` olay içeren en kısa penceredir.
    λ_T, λ_ref - k·sqrt(λ_ref/T) altına düşene kadar katalog tam kabul edilir.
    """
    ev = df[df["is_preferred"] & ~df["suspect"].astype(bool) & df["mag"].notna()]
    years = ev["decimal_year"].to_numpy()
    mags = ev["mag"].to_numpy()
    end = cfg.max_year + 1
    first = int(np.floor(years.min())) if len(years) else end
    result = {}
    for mc in cfg.completeness_mag_bins:
        y = years[mags >= mc]
        Ts = list(range(cfg.stepp_window_years, end - first + 1, cfg.stepp_window_years))
        rates = [(T, int(np.sum(y >= end - T)) / T) for T in Ts]
        ref = next(((T, r) for T, r in rates if r * T >= cfg.stepp_min_events), None)
        entry = {"rates": [{"T_years": T, "start_year": end - T, "n": int(round(r * T)), "rate_per_yr": round(r, 5)}
                           for T, r in rates]}
        if ref is None:
            entry.update(start_year=cfg.default_completeness[mc], method="varsayılan (yetersiz olay)")
        else:
            lam = ref[1]
            start = end - ref[0]
            for T, r in rates:
                if T < ref[0]:
                    continue
                if r < lam - cfg.stepp_tolerance_sigma * math.sqrt(lam / T):
                    break
                start = end - T
            entry.update(start_year=int(start), method="stepp", ref_rate_per_yr=round(lam, 5))
        result[mc] = entry
    # Büyük büyüklükler için başlangıç, küçükler için olandan geç olamaz (tutarlılık).
    prev = None
    for mc in sorted(result):
        s = result[mc]["start_year"]
        if prev is not None and s > prev:
            result[mc]["start_year"] = prev
            result[mc]["method"] += " (monotonluk düzeltmesi)"
        prev = result[mc]["start_year"]
    return result


def in_window(mag: float, year: float, table: dict) -> bool:
    if mag is None or (isinstance(mag, float) and math.isnan(mag)):
        return False
    bins = [mc for mc in sorted(table) if mag >= mc]
    if not bins:
        return False
    return year >= table[bins[-1]]["start_year"]


# ---------------------------------------------------------------------------
def add_weights(df: pd.DataFrame, completeness: dict, cfg: Config = DEFAULT) -> pd.DataFrame:
    df = df.copy()
    dm = df["agency"].map(cfg.default_mag_err)
    dh = df["agency"].map(cfg.default_horiz_err_km)
    df["mag_err_assigned"] = df["mag"].notna() & df["mag_err"].isna()
    # σM: kaynak değeri (yoksa atanan), kaynaklar arası yarım-aralıkla karesel birleştirilir.
    spread = ((df["group_mag_max"] - df["group_mag_min"]) / 2).fillna(0.0)
    base = df["mag_err"].fillna(dm)
    df["mag_err_used"] = np.where(df["mag"].isna(), cfg.missing_mag_err, np.sqrt(base**2 + spread**2))
    df["horiz_err_assigned"] = df["horiz_err_km"].isna()
    df["horiz_err_used"] = df["horiz_err_km"].fillna(dh)
    t_days = df["time_err_s"] / 86400.0

    qm = (df["mag_err_used"] / cfg.scale_mag) ** 2
    ql = (df["horiz_err_used"] / cfg.scale_horiz_km) ** 2
    qt = (t_days / cfg.scale_time_days) ** 2
    qz = (df["depth_err_km"].fillna(0.0) / cfg.scale_depth_km) ** 2
    df["w_mag"] = 1 / (1 + qm)
    df["w_loc"] = 1 / (1 + ql + qz)
    df["w_time"] = 1 / (1 + qt)
    df["w_unc"] = np.maximum(cfg.w_floor, 1 / (1 + qm + ql + qz + qt))
    df["in_completeness_window"] = [in_window(m, y, completeness) for m, y in zip(df["mag"], df["decimal_year"])]
    df["w_completeness"] = np.where(df["in_completeness_window"], 1.0, cfg.incomplete_factor)
    df["w_suspect"] = np.where(df["suspect"].astype(bool), cfg.suspect_factor, 1.0)
    df["training_weight"] = df["w_unc"] * df["w_completeness"] * df["w_suspect"] * df["is_preferred"].astype(float)
    for c in ["w_mag", "w_loc", "w_time", "w_unc", "training_weight", "mag_err_used", "horiz_err_used"]:
        df[c] = df[c].round(4)
    return df


def build(cfg: Config = DEFAULT, out_dir: Path = DATA) -> dict:
    df = load_all(cfg)
    df = find_duplicates(df, cfg)
    completeness = stepp_completeness(df, cfg)
    df = add_weights(df, completeness, cfg)
    df = df.sort_values(["decimal_year", "agency"]).reset_index(drop=True)[SCHEMA_COLUMNS]

    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "tarihsel_tum_kayitlar.csv", index=False)
    pref = df[df["is_preferred"]].reset_index(drop=True)
    pref.to_csv(out_dir / "tarihsel_katalog.csv", index=False)
    meta = {
        "config": {k: (v if not isinstance(v, tuple) else list(v)) for k, v in asdict(cfg).items()},
        "completeness": {str(k): v for k, v in completeness.items()},
        "counts": {
            "all_records": int(len(df)),
            "events": int(len(pref)),
            "by_agency_all": df["agency"].value_counts().to_dict(),
            "by_agency_preferred": pref["agency"].value_counts().to_dict(),
            "events_mag_ge_6": int((pref["mag"] >= 6).sum()),
            "events_mag_ge_7": int((pref["mag"] >= 7).sum()),
        },
    }
    (out_dir / "tamlik.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2, default=str))
    return meta


if __name__ == "__main__":
    m = build()
    print(json.dumps(m["counts"], ensure_ascii=False, indent=2))
    for k, v in m["completeness"].items():
        print(f"Mw>={k}: {v['start_year']} ({v['method']})")
