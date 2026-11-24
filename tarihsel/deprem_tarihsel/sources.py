"""Kaynak kataloğu okuyucuları. Her biri ortak şemadaki sözlüklerin listesini döndürür.

Ortak alanlar (aletsel katalog şemasıyla uyumlu):
    event_id, time, latitude, longitude, depth_km, mag, mag_type, source, agency,
    time_err_s, horiz_err_km, depth_err_km, mag_err
Tarihsel ekler: time_min, time_max, time_precision, decimal_year, mag_is_estimated,
    mag_method, mag_orig, mag_orig_type, intensity_io, area, country, deaths,
    reference, source_event_id, suspect, notes
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import pandas as pd

from .config import Config
from .timeutil import astronomical_year, summarize

RAW = Path(__file__).resolve().parent.parent / "data" / "raw"


def _num(v):
    if v is None:
        return None
    if isinstance(v, str):
        v = v.strip()
        if not v:
            return None
        try:
            return float(v)
        except ValueError:
            return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def parse_intensity(v):
    """'8-9' -> 8.5, '10' -> 10.0, 'F' / boş -> None."""
    if v is None:
        return None
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return None
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)", s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2
    return _num(s)


def ms_to_mw(ms: float) -> float:
    """Scordilis (2006) küresel Ms->Mw bağıntısı."""
    if ms <= 6.1:
        return 0.67 * ms + 2.07
    return 0.99 * ms + 0.08


def _horiz(lat_unc, lon_unc):
    a, b = _num(lat_unc), _num(lon_unc)
    if a is None and b is None:
        return None
    if a is None or b is None:
        return a if b is None else b
    return round(math.sqrt(a * b), 1)  # elipsin geometrik ortalama yarıçapı


def _in_region(lat, lon, cfg: Config) -> bool:
    return (
        lat is not None
        and lon is not None
        and cfg.lat_min <= lat <= cfg.lat_max
        and cfg.lon_min <= lon <= cfg.lon_max
    )


def _base(agency, sid, year, mo, da, ho, mi, se, lat, lon):
    rec = {
        "event_id": f"{agency}:{sid}",
        "source_event_id": str(sid),
        "source": "historical",
        "agency": agency,
        "latitude": lat,
        "longitude": lon,
    }
    rec.update(summarize(year, mo, da, ho, mi, se))
    return rec


# ---------------------------------------------------------------------------
def read_epica(cfg: Config, path: Path = RAW / "EPICA_v1.1.xlsx") -> list[dict]:
    df = pd.read_excel(path, sheet_name="EPICA")
    out = []
    for r in df.itertuples(index=False):
        lat, lon = _num(r.Lat), _num(r.Lon)
        if not _in_region(lat, lon, cfg) or int(r.Year) > cfg.max_year:
            continue
        rec = _base("EPICA", r.EqID, int(r.Year), _num(r.Mo), _num(r.Da), _num(r.Ho), _num(r.Mi), None, lat, lon)
        mw, mwu = _num(r.Mw), _num(r.MwUnc)
        tmw = None if pd.isna(r.TMw) else str(r.TMw)
        rec.update(
            depth_km=_num(r.H),
            depth_err_km=_num(r.HUnc),
            horiz_err_km=_horiz(r.LatUnc, r.LonUnc),
            mag=mw,
            mag_type="mw" if mw is not None else None,
            mag_err=mwu,
            mag_orig=mw,
            mag_orig_type="mw" if mw is not None else None,
            mag_method={"wm": "makrosismik+katalog ağırlıklı ort.", "MMw": "makrosismik veri noktaları",
                        "CMw": "kaynak katalog", "nd": None}.get(tmw, tmw),
            intensity_io=parse_intensity(r.Io),
            area=None if pd.isna(r.Ax) else str(r.Ax),
            country=None,
            deaths=None,
            reference="; ".join(str(x) for x in (r.MDPsSource, r.CatSource) if not pd.isna(x)) or None,
            suspect=False,
            notes=f"TEpi={r.TEpi}; TEpiUnc={r.TEpiUnc}; TMw={tmw}",
        )
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
def read_share_cet(cfg: Config, path: Path = RAW / "SHARE_CET.xls") -> list[dict]:
    df = pd.read_excel(path, sheet_name="catalogue")
    out = []
    for r in df.itertuples(index=False):
        lat, lon = _num(r.Lat), _num(r.Lon)
        if not _in_region(lat, lon, cfg) or int(r.Year) > cfg.max_year:
            continue
        sid = int(r.En) if _num(r.En) is not None else f"{r.Year}-{r.Lat}-{r.Lon}"
        rec = _base("SHARE_CET", sid, int(r.Year), _num(r.Mo), _num(r.Da), _num(r.Ho), _num(r.Mi), _num(r.Se), lat, lon)
        mw, mwu = _num(r.Mw), _num(r.MwUnc)
        tmw = None if pd.isna(r.TMw) else str(r.TMw)
        rec.update(
            depth_km=_num(r.H),
            depth_err_km=None,
            horiz_err_km=_horiz(r.LatUnc, r.LonUnc),
            mag=mw,
            mag_type="mw" if mw is not None else None,
            mag_err=mwu,
            mag_orig=mw,
            mag_orig_type="mw" if mw is not None else None,
            mag_method={"co": "bölgesel katalog", "cc": "bölgesel katalog, Io veya başka M türünden dönüştürülmüş",
                        "nd": None}.get(tmw, tmw),
            intensity_io=parse_intensity(r.Io),
            area=None if pd.isna(r.Ax) else str(r.Ax),
            country=None,
            deaths=None,
            reference="SHARE Central & Eastern Turkey catalogue",
            suspect=False,
            notes=f"TMw={tmw}",
        )
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
def read_ghec(cfg: Config, path: Path = RAW / "GEM-GHEC-v1.txt") -> list[dict]:
    lines = [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    out = []
    for ln in lines[1:]:
        f = dict(zip(header, ln.split("\t")))
        lat, lon = _num(f.get("Lat")), _num(f.get("Lon"))
        year = int(f["Year"])
        if not _in_region(lat, lon, cfg) or year > cfg.max_year:
            continue
        rec = _base("GEM_GHEC", f["En"], year, _num(f.get("Mo")), _num(f.get("Da")), _num(f.get("Ho")),
                    _num(f.get("Mi")), _num(f.get("Se")), lat, lon)
        m, mu, mt = _num(f.get("M")), _num(f.get("MUnc")), (f.get("MType") or "").strip()
        mag, err, mtype, note = m, mu, None, ""
        if m is not None:
            if mt == "w":
                mtype = "mw"
            elif mt == "s":
                mag = round(ms_to_mw(m), 2)
                err = math.sqrt((mu or cfg.default_mag_err["GEM_GHEC"]) ** 2 + cfg.ms_to_mw_err**2)
                err = round(err, 2)
                mtype = "mw"
                note = "Ms->Mw (Scordilis 2006); "
            else:
                mtype = "mw"
                err = round(math.sqrt((mu or cfg.default_mag_err["GEM_GHEC"]) ** 2 + cfg.unknown_mag_type_err**2), 2)
                note = "büyüklük türü belirsiz, Mw sayıldı; "
        rec.update(
            depth_km=_num(f.get("Dep")),
            depth_err_km=None,
            horiz_err_km=_horiz(f.get("LatUnc"), f.get("LonUnc")),
            mag=mag,
            mag_type=mtype,
            mag_err=err,
            mag_orig=m,
            mag_orig_type={"w": "mw", "s": "ms"}.get(mt, mt or None),
            mag_method=(f.get("MDet") or None),
            intensity_io=parse_intensity(f.get("Io")),
            area=(f.get("Area") or None),
            country=None,
            deaths=None,
            reference=f"{f.get('Source')}; M: {f.get('Msource')}",
            suspect=False,
            notes=(note + (f.get("Remarks") or "").strip('"')).strip("; ") or None,
        )
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
NOAA_COUNTRIES = {
    "TURKEY", "GREECE", "SYRIA", "IRAQ", "IRAN", "ARMENIA", "GEORGIA", "AZERBAIJAN",
    "CYPRUS", "LEBANON", "BULGARIA", "RUSSIA", "ISRAEL", "JORDAN", "EGYPT",
}


def read_noaa(cfg: Config, paths=None) -> list[dict]:
    paths = paths or sorted(RAW.glob("ncei_tr_pre1900_p*.json"))
    items = []
    for p in paths:
        items += json.loads(Path(p).read_text())["items"]
    out = []
    for it in items:
        lat, lon = _num(it.get("latitude")), _num(it.get("longitude"))
        year = astronomical_year(int(it["year"]), historical_bc=True)
        if not _in_region(lat, lon, cfg) or year > cfg.max_year:
            continue
        rec = _base("NOAA_NCEI", it["id"], year, it.get("month"), it.get("day"), it.get("hour"),
                    it.get("minute"), it.get("second"), lat, lon)
        notes = []
        if int(it["year"]) < 0:
            notes.append(f"MÖ {abs(int(it['year']))} (astronomik yıl {year})")
        mag, err, orig, otype = None, None, None, None
        base = cfg.default_mag_err["NOAA_NCEI"]
        if _num(it.get("eqMagMw")) is not None:
            mag = orig = _num(it["eqMagMw"])
            otype = "mw"
            err = base
        elif _num(it.get("eqMagMs")) is not None:
            orig, otype = _num(it["eqMagMs"]), "ms"
            mag = round(ms_to_mw(orig), 2)
            err = round(math.sqrt(base**2 + cfg.ms_to_mw_err**2), 2)
            notes.append("Ms->Mw (Scordilis 2006)")
        elif _num(it.get("eqMagnitude")) is not None:
            orig = mag = _num(it.get("eqMagUnk") or it["eqMagnitude"])
            otype = "unk"
            err = round(math.sqrt(base**2 + cfg.unknown_mag_type_err**2), 2)
            notes.append("büyüklük türü belirsiz, Mw sayıldı")
        country = it.get("country")
        suspect = country not in NOAA_COUNTRIES
        if suspect:
            notes.append("ülke adı koordinatla uyumsuz; konum kuşkulu")
        rec.update(
            depth_km=_num(it.get("eqDepth")),
            depth_err_km=None,
            horiz_err_km=None,
            mag=mag,
            mag_type="mw" if mag is not None else None,
            mag_err=err,
            mag_orig=orig,
            mag_orig_type=otype,
            mag_method=None,
            intensity_io=_num(it.get("intensity")),
            area=it.get("locationName"),
            country=country,
            deaths=_num(it.get("deaths")),
            reference="NOAA NCEI Global Significant Earthquake Database",
            suspect=suspect,
            notes="; ".join(notes) or None,
        )
        out.append(rec)
    return out


READERS = {
    "EPICA": read_epica,
    "SHARE_CET": read_share_cet,
    "GEM_GHEC": read_ghec,
    "NOAA_NCEI": read_noaa,
}
