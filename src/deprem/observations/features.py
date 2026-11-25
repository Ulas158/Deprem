"""Yıllık (hücre × yıl) gözlem özelliklerini tek tabloda topla.

Her yıl ``t0`` için yalnız ``t0`` öncesi veriler kullanılır. Hiç verisi
olmayan özellik sütun olarak eklenmez (model onu görmez, raporda
"veri bekliyor" yazar). Kısmen eksik olanlar NaN kalır; eğitim hattı eksik
bayrağı ekler.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from . import dvv as dvv_mod
from . import gnss, gonaf, insar, repeaters
from .registry import OBS_FEATURES, apply_transform

log = logging.getLogger(__name__)


def yearly_table(root: Path, cells: pd.DataFrame, years: list[int], catalog: pd.DataFrame,
                 include: list[str] | None = None) -> pd.DataFrame:
    obs = Path(root) / "data" / "obs"
    raw = Path(root) / "data" / "raw"
    sty = pd.read_csv(obs / "gnss_station_year.csv") if (obs / "gnss_station_year.csv").exists() else None
    dv = None
    if (obs / "dvv_series.csv.gz").exists():
        dv = pd.read_csv(obs / "dvv_series.csv.gz", parse_dates=["time"])
    ins = insar.cell_features(raw / "insar", cells)
    frames = []
    for y in years:
        t0 = pd.Timestamp(f"{y}-01-01", tz="UTC")
        parts = [pd.DataFrame(index=cells["cell_id"])]
        if sty is not None:
            parts.append(gnss.cell_features(sty, cells, y))
        past = catalog[catalog["time"] < t0]
        fam = repeaters.find_families(past) if len(past) else pd.DataFrame()
        parts.append(repeaters.cell_features(fam, cells, t0))
        if dv is not None:
            parts.append(dvv_mod.cell_features(dv, cells, t0.tz_localize(None)))
        parts.append(gonaf.cell_features(raw / "gonaf", cells, t0))
        parts.append(ins)
        d = pd.concat(parts, axis=1)
        d["year"] = y
        frames.append(d.reset_index().rename(columns={"index": "cell_id"}))
    tab = pd.concat(frames, ignore_index=True)
    feat_cols = [c for c in tab.columns if c in OBS_FEATURES]
    if include is not None:
        feat_cols = [c for c in feat_cols if c in include or OBS_FEATURES[c].group in include]
    # tamamen boş veya sabit sütunları at
    keep = [c for c in feat_cols if tab[c].notna().sum() > 0 and tab[c].nunique(dropna=True) > 1]
    for c in keep:
        tab[c] = apply_transform(c, tab[c])
    log.info("Gözlem özellikleri: %s", keep)
    return tab[["cell_id", "year"] + keep]
