"""GEM Küresel Aktif Fay veritabanı (Styron & Pagani 2020).

Türkiye segmentleri büyük ölçüde EMME (Earthquake Model of the Middle East,
Danciu vd. 2018) fay modelinden gelir ve her segment için kayma hızı,
eğim, kayma açısı (rake) ve sismojenik derinlik aralığı içerir. Değerler
``(en_iyi, min, maks)`` biçiminde dizgedir.

Kaynak: https://github.com/GEMScienceTools/gem-global-active-faults
Lisans: CC BY-SA 4.0
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from ..geo import polyline_length_km
from ..region import Region

GEM_URL = ("https://raw.githubusercontent.com/GEMScienceTools/gem-global-active-faults/"
           "master/geojson/gem_active_faults_harmonized.geojson")

_TRIPLE = re.compile(r"\(([^,]*),([^,]*),([^)]*)\)")


def parse_triple(s):
    """'(0.2,0.1,0.4)' -> (0.2, 0.1, 0.4); eksik alanlar NaN."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return (np.nan, np.nan, np.nan)
    if isinstance(s, (int, float)):
        return (float(s), np.nan, np.nan)
    mt = _TRIPLE.match(str(s).replace(" ", ""))
    if not mt:
        try:
            return (float(s), np.nan, np.nan)
        except ValueError:
            return (np.nan, np.nan, np.nan)
    return tuple(float(v) if v not in ("", "None") else np.nan for v in mt.groups())


def download(dest: Path) -> Path:
    dest = Path(dest)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(GEM_URL, timeout=300)
        r.raise_for_status()
        dest.write_bytes(r.content)
    return dest


def _lines(geom):
    if geom is None:
        return []
    if geom["type"] == "LineString":
        return [geom["coordinates"]]
    if geom["type"] == "MultiLineString":
        return geom["coordinates"]
    return []


def load_faults(path: Path, region: Region) -> pd.DataFrame:
    """Bölgeye değen fay segmentlerini tablo olarak yükle.

    Her satır bir segmenttir; ``coords`` sütunu [(lon, lat), ...] listesidir.
    """
    data = json.loads(Path(path).read_text())
    rows = []
    for f in data["features"]:
        p = f["properties"]
        for line in _lines(f["geometry"]):
            c = np.asarray(line, dtype=float)[:, :2]
            inside = ((c[:, 1] >= region.min_lat) & (c[:, 1] <= region.max_lat)
                      & (c[:, 0] >= region.min_lon) & (c[:, 0] <= region.max_lon))
            if not inside.any():
                continue
            sr, sr_min, sr_max = parse_triple(p.get("net_slip_rate"))
            dip = parse_triple(p.get("average_dip"))[0]
            rake = parse_triple(p.get("average_rake"))[0]
            up = parse_triple(p.get("upper_seis_depth"))[0]
            lo = parse_triple(p.get("lower_seis_depth"))[0]
            rows.append({
                "fault_id": p.get("catalog_id"), "name": p.get("name"),
                "catalog": p.get("catalog_name"), "slip_type": p.get("slip_type"),
                "slip_rate_mm_yr": sr, "slip_rate_min": sr_min, "slip_rate_max": sr_max,
                "dip": dip, "rake": rake, "upper_depth_km": up, "lower_depth_km": lo,
                "coords": [tuple(x) for x in c.tolist()],
            })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["upper_depth_km"] = df["upper_depth_km"].fillna(0.0)
    df["lower_depth_km"] = df["lower_depth_km"].fillna(15.0)
    df["dip"] = df["dip"].fillna(90.0)
    df["length_km"] = df["coords"].map(polyline_length_km)
    return df.reset_index(drop=True)
