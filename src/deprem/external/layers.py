"""Gaz/sıvı ve dış etki katmanlarını web haritası için yazma.

Biçim: ``web/data/layers/<id>.layer.json`` + ``<id>.geojson`` (bkz. web/README.md,
"Veri katmanı ekleme"). Koordinatlar 4 ondalığa yuvarlanır.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

GROUP_HUMAN = "İnsan kaynaklı etkiler"
GROUP_FLUID = "Gaz ve sıvılar"
GROUP_LOADS = "Dış yükler (gelgit, mevsim, uzak depremler)"

STATUS_COLORS = {"aktif": "#35e08a", "pasif": "#ff9f43", "kapali": "#8796ab",
                 "bilinmiyor": "#b18cff", "dogal": "#4fc3f7"}
STATUS_LABELS = {"aktif": "Aktif", "pasif": "Pasif / durduruldu", "kapali": "Kapalı",
                 "bilinmiyor": "Durumu bilinmiyor", "dogal": "Doğal"}


def _clean(v):
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else round(float(v), 4)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, pd.Timestamp):
        return v.strftime("%Y-%m-%d")
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return None
    return v


def points_geojson(df: pd.DataFrame, props: list[str], lat="lat", lon="lon", note="") -> dict:
    feats = []
    for _, r in df.iterrows():
        feats.append({"type": "Feature",
                      "geometry": {"type": "Point", "coordinates": [round(float(r[lon]), 4), round(float(r[lat]), 4)]},
                      "properties": {p: _clean(r[p]) for p in props if p in df.columns}})
    return {"type": "FeatureCollection", "metadata": {"note": note}, "features": feats}


def cells_geojson(df: pd.DataFrame, cell_deg: float, props: list[str], note="") -> dict:
    feats = []
    h = cell_deg / 2
    for _, r in df.iterrows():
        x, y = float(r["lon_c"]), float(r["lat_c"])
        ring = [[round(x - h, 4), round(y - h, 4)], [round(x + h, 4), round(y - h, 4)],
                [round(x + h, 4), round(y + h, 4)], [round(x - h, 4), round(y + h, 4)],
                [round(x - h, 4), round(y - h, 4)]]
        feats.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [ring]},
                      "properties": {p: _clean(r[p]) for p in props if p in df.columns}})
    return {"type": "FeatureCollection", "metadata": {"note": note}, "features": feats}


def write_layer(out_dir: Path, spec: dict, geojson: dict | None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    lid = spec["id"]
    spec = {"data": f"{lid}.geojson", "visible": False, **spec}
    (out_dir / f"{lid}.layer.json").write_text(json.dumps(spec, ensure_ascii=False, indent=2) + "\n")
    if geojson is not None:
        (out_dir / f"{lid}.geojson").write_text(json.dumps(geojson, ensure_ascii=False, separators=(",", ":")))


def site_layer_specs() -> list[tuple[dict, list[str]]]:
    """(tanım, kategori listesi) çiftleri: her saha türü ayrı aç/kapa katmanı."""
    tip = [{"property": "category", "label": "Tür"},
           {"property": "status", "label": "Durum"},
           {"property": "start_year", "label": "Başlangıç yılı", "digits": 0},
           {"property": "end_year", "label": "Bitiş yılı", "digits": 0},
           {"property": "fluid_process", "label": "Olası etki"},
           {"property": "notes", "label": "Not"},
           {"property": "loc_accuracy_km", "label": "Konum belirsizliği", "unit": "km", "digits": 0},
           {"property": "verified", "label": "Doğrulandı mı"},
           {"property": "source", "label": "Kaynak"}]
    base = {"group": GROUP_HUMAN, "style": {"type": "circle", "radius": 6, "outline": "#0b111a"},
            "color": {"property": "status", "type": "categorical", "categories": STATUS_COLORS, "default": "#8aa0bb"},
            "legend": {"title": "Durum", "categoryLabels": STATUS_LABELS},
            "tooltip": {"title": "name", "fields": tip},
            "note": "Konumların çoğu yaklaşık ve doğrulanmamıştır; tahmin değildir."}
    def leg(title):
        return {"legend": {"title": f"{title} · durum", "categoryLabels": STATUS_LABELS}}

    return [
        ({**base, **leg("Barajlar"), "id": "ext_barajlar", "label": "Barajlar ve rezervuarlar", "order": 300,
          "description": "Su yükü ve gözenek basıncı; dolum yılıyla"}, ["baraj"]),
        ({**base, **leg("Madenler"), "id": "ext_madenler", "label": "Madenler", "order": 310,
          "description": "Kömür, metal, bor, çözelti madenciliği; aktif/pasif ve tarihçe"},
         ["maden_komur", "maden_metal", "maden_endustriyel", "maden_cozelti"]),
        ({**base, **leg("Jeotermal"), "id": "ext_jeotermal", "label": "Jeotermal sahalar", "order": 320,
          "description": "Akışkan çekimi ve geri basma"}, ["jeotermal"]),
        ({**base, **leg("Petrol/gaz"), "id": "ext_petrol_gaz", "label": "Petrol, gaz sahaları ve gaz depoları", "order": 330,
          "description": "Akışkan çekimi, CO2 ve su enjeksiyonu, döngüsel gaz depolama"},
         ["petrol", "dogalgaz", "gaz_depolama"]),
        ({**base, **leg("Gaz çıkışları"), "id": "ext_gaz_cikislari", "label": "Doğal gaz çıkışları", "order": 400, "group": GROUP_FLUID,
          "description": "Deniz tabanı / fay boyu gaz çıkışları (literatür)"}, ["dogal_gaz_cikisi"]),
    ]
