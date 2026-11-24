#!/usr/bin/env python3
"""Harita için statik veri katmanlarını üretir (web/data/*.geojson).

Kaynaklar:
  - GEM Global Active Faults (CC BY-SA 4.0)
    https://github.com/GEMScienceTools/gem-global-active-faults
  - Natural Earth 10m ülke sınırları ve yerleşimler (kamu malı)
    https://github.com/nvkelso/natural-earth-vector

Kullanım:
  python3 web/tools/build_data.py <gem.geojson> <ne_countries.geojson> <ne_places.geojson>
Ek bağımlılık gerektirmez (yalnızca standart kütüphane).
"""
import json
import math
import re
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "data"

# Türkiye ve çevresi
FAULT_BBOX = (24.0, 34.0, 46.0, 43.5)      # lon_min, lat_min, lon_max, lat_max
BASE_BBOX = (14.0, 28.0, 56.0, 49.0)


def in_bbox(pt, b):
    return b[0] <= pt[0] <= b[2] and b[1] <= pt[1] <= b[3]


def dp(points, eps):
    """Douglas-Peucker sadeleştirme."""
    if len(points) < 3:
        return points
    (x1, y1), (x2, y2) = points[0], points[-1]
    dx, dy = x2 - x1, y2 - y1
    norm = math.hypot(dx, dy) or 1e-12
    idx, dmax = 0, 0.0
    for i in range(1, len(points) - 1):
        px, py = points[i]
        d = abs(dy * px - dx * py + x2 * y1 - y2 * x1) / norm
        if d > dmax:
            idx, dmax = i, d
    if dmax > eps:
        return dp(points[: idx + 1], eps)[:-1] + dp(points[idx:], eps)
    return [points[0], points[-1]]


def rnd(points, nd):
    out = []
    for x, y in points:
        p = [round(x, nd), round(y, nd)]
        if not out or out[-1] != p:
            out.append(p)
    return out


def triple(s):
    """'(1.55,0.8,2.22)' -> (1.55, 0.8, 2.22); boş alanlar None."""
    if not s:
        return (None, None, None)
    vals = [v.strip() for v in re.sub(r"[()]", "", str(s)).split(",")]
    vals += [""] * (3 - len(vals))
    out = []
    for v in vals[:3]:
        try:
            out.append(float(v))
        except ValueError:
            out.append(None)
    return tuple(out)


def build_faults(path):
    src = json.load(open(path))
    feats = []
    for f in src["features"]:
        g = f.get("geometry")
        if not g:
            continue
        lines = [g["coordinates"]] if g["type"] == "LineString" else g["coordinates"]
        if not any(in_bbox(p, FAULT_BBOX) for l in lines for p in l):
            continue
        p = f["properties"]
        rate, rmin, rmax = triple(p.get("net_slip_rate"))
        dip = triple(p.get("average_dip"))[0]
        rake = triple(p.get("average_rake"))[0]
        for l in lines:
            coords = rnd(l, 4)
            if len(coords) < 2:
                continue
            feats.append({
                "type": "Feature",
                "geometry": {"type": "LineString", "coordinates": coords},
                "properties": {
                    "name": p.get("name") or "",
                    "slip_type": p.get("slip_type") or "",
                    "slip_rate": rate, "slip_rate_min": rmin, "slip_rate_max": rmax,
                    "dip": dip, "rake": rake,
                    "catalog": p.get("catalog_name") or "",
                },
            })
    return {"type": "FeatureCollection", "features": feats}


def build_countries(path):
    src = json.load(open(path))
    feats = []
    for f in src["features"]:
        g = f["geometry"]
        polys = [g["coordinates"]] if g["type"] == "Polygon" else g["coordinates"]
        keep = []
        for poly in polys:
            if not any(in_bbox(p, BASE_BBOX) for p in poly[0]):
                continue
            rings = []
            for ring in poly:
                # kapalı halkayı ikiye bölerek sadeleştir (uçlar çakışık olduğundan)
                m = len(ring) // 2
                r = rnd(dp(ring[: m + 1], 0.004)[:-1] + dp(ring[m:], 0.004), 3)
                if len(r) >= 4:
                    rings.append(r)
            if rings:
                keep.append(rings)
        if keep:
            p = f["properties"]
            feats.append({
                "type": "Feature",
                "geometry": {"type": "MultiPolygon", "coordinates": keep},
                "properties": {"name": p.get("NAME_TR") or p.get("NAME"), "iso": p.get("ADM0_A3")},
            })
    return {"type": "FeatureCollection", "features": feats}


def build_cities(path):
    src = json.load(open(path))
    feats = []
    for f in src["features"]:
        p = f["properties"]
        tur = p.get("adm0_a3") == "TUR"
        if not ((tur and p["pop_max"] >= 150000) or (not tur and p.get("adm0cap") and in_bbox([p["longitude"], p["latitude"]], FAULT_BBOX))):
            continue
        feats.append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [round(p["longitude"], 3), round(p["latitude"], 3)]},
            "properties": {"name": p["name"], "pop": p["pop_max"], "tr": tur},
        })
    feats.sort(key=lambda f: -f["properties"]["pop"])
    return {"type": "FeatureCollection", "features": feats}


def dump(obj, name):
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    path.write_text(json.dumps(obj, ensure_ascii=False, separators=(",", ":")))
    print(f"{name}: {len(obj['features'])} öğe, {path.stat().st_size/1024:.0f} KB")


if __name__ == "__main__":
    gem, countries, places = sys.argv[1:4]
    dump(build_faults(gem), "faults.geojson")
    dump(build_countries(countries), "countries.geojson")
    dump(build_cities(places), "cities.geojson")
