"""Açık deprem kataloglarından veri çekme.

Desteklenen kaynaklar:

* ``usgs``  USGS ComCat FDSN olay servisi. Türkiye için küresel ağa dayanır,
  tamamlılık yaklaşık M 4 civarındadır; uzun geçmiş (1900+) sağlar.
* ``emsc``  EMSC FDSN olay servisi. Türkiye olaylarının büyük kısmı KOERI
  ve AFAD çözümlerinden gelir; küçük büyüklüklerde çok daha zengindir.
* ``afad``  AFAD Deprem Dairesi açık API'si (servisnet.afad.gov.tr). Ulusal
  ağın resmi kataloğu, M~1'e kadar iner. **Sınırlama:** Kasım 2026 itibarıyla
  açık API yalnızca son ~2 haftayı döndürüyor; daha eski istekler boş liste
  veriyor. Bu yüzden (a) AFAD verisi düzenli aralıklarla çekilip önbellekte
  biriktirilmeli, (b) geçmiş katalog AFAD web kataloğundan
  (deprem.afad.gov.tr/event-catalog) elle dışa aktarılıp
  :func:`load_afad_export` ile okunmalı.

Her istek zaman dilimlerine bölünür ve ham yanıt diske önbelleğe alınır.
Böylece analiz tekrar çalıştırıldığında servisler yeniden sorgulanmaz ve
sonuçlar o anki veri anlık görüntüsüne göre tekrarlanabilir kalır.
"""

from __future__ import annotations

import io
import json
import logging
import time as _time
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import requests

from ..region import Region
from .schema import empty_catalog, normalize

log = logging.getLogger(__name__)

USER_AGENT = "deprem-research-prototype/0.1"


@dataclass
class FetchSpec:
    region: Region
    start: pd.Timestamp
    end: pd.Timestamp
    min_mag: float
    chunk: str = "90D"  # pandas frekans dizgesi


def _time_chunks(start: pd.Timestamp, end: pd.Timestamp, freq: str):
    edges = list(pd.date_range(start, end, freq=freq))
    if not edges or edges[0] != start:
        edges.insert(0, start)
    if edges[-1] != end:
        edges.append(end)
    return list(zip(edges[:-1], edges[1:]))


def _get(url: str, params: dict, retries: int = 4, timeout: int = 120) -> requests.Response:
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout,
                             headers={"User-Agent": USER_AGENT})
            if r.status_code in (200, 204, 404):
                return r
            last = RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except requests.RequestException as e:  # ağ hatası
            last = e
        _time.sleep(2 ** attempt)
    raise RuntimeError(f"İstek başarısız: {url} {params}") from last


# --------------------------------------------------------------------- USGS
USGS_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"
USGS_LIMIT = 20000


def _usgs_chunk(region: Region, t0, t1, min_mag) -> str:
    params = {
        "format": "csv", "starttime": t0.isoformat(), "endtime": t1.isoformat(),
        "minlatitude": region.min_lat, "maxlatitude": region.max_lat,
        "minlongitude": region.min_lon, "maxlongitude": region.max_lon,
        "minmagnitude": min_mag, "orderby": "time-asc", "limit": USGS_LIMIT,
    }
    return _get(USGS_URL, params).text


def parse_usgs(text: str) -> pd.DataFrame:
    if not text.strip():
        return empty_catalog()
    raw = pd.read_csv(io.StringIO(text))
    if raw.empty:
        return empty_catalog()
    raw = raw[raw["type"] == "earthquake"]
    df = pd.DataFrame({
        "time": raw["time"], "latitude": raw["latitude"], "longitude": raw["longitude"],
        "depth_km": raw["depth"], "mag": raw["mag"], "mag_type": raw["magType"],
        "source": "usgs", "agency": raw["magSource"], "event_id": raw["id"],
        "horiz_err_km": raw.get("horizontalError"), "depth_err_km": raw.get("depthError"),
        "mag_err": raw.get("magError"),
    })
    return normalize(df)


# --------------------------------------------------------------------- EMSC
EMSC_URL = "https://www.seismicportal.eu/fdsnws/event/1/query"


def _emsc_chunk(region: Region, t0, t1, min_mag) -> str:
    params = {
        "format": "text", "starttime": t0.strftime("%Y-%m-%dT%H:%M:%S"),
        "endtime": t1.strftime("%Y-%m-%dT%H:%M:%S"),
        "minlat": region.min_lat, "maxlat": region.max_lat,
        "minlon": region.min_lon, "maxlon": region.max_lon,
        "minmag": min_mag, "orderby": "time-asc", "limit": 20000,
    }
    return _get(EMSC_URL, params).text


def parse_emsc(text: str) -> pd.DataFrame:
    if not text.strip():
        return empty_catalog()
    raw = pd.read_csv(io.StringIO(text), sep="|")
    if raw.empty:
        return empty_catalog()
    raw.columns = [c.lstrip("#").strip() for c in raw.columns]
    df = pd.DataFrame({
        "time": raw["Time"], "latitude": raw["Latitude"], "longitude": raw["Longitude"],
        "depth_km": raw["Depth/km"], "mag": raw["Magnitude"], "mag_type": raw["MagType"],
        "source": "emsc", "agency": raw["Author"], "event_id": raw["EventID"].astype(str),
    })
    return normalize(df)


# --------------------------------------------------------------------- AFAD
AFAD_URL = "https://servisnet.afad.gov.tr/apigateway/deprem/apiv2/event/filter"


def _afad_chunk(region: Region, t0, t1, min_mag) -> str:
    params = {
        "start": t0.strftime("%Y-%m-%dT%H:%M:%S"), "end": t1.strftime("%Y-%m-%dT%H:%M:%S"),
        "minlat": region.min_lat, "maxlat": region.max_lat,
        "minlon": region.min_lon, "maxlon": region.max_lon,
        "minmag": min_mag, "format": "json", "orderby": "time",
    }
    return _get(AFAD_URL, params).text


def parse_afad(text: str) -> pd.DataFrame:
    if not text.strip():
        return empty_catalog()
    rows = json.loads(text)
    if not rows:
        return empty_catalog()
    raw = pd.DataFrame(rows)
    df = pd.DataFrame({
        # AFAD apiv2 zamanları UTC olarak döner (yerel saat değil)
        "time": pd.to_datetime(raw["date"]).dt.tz_localize("UTC"),
        "latitude": raw["latitude"], "longitude": raw["longitude"],
        "depth_km": raw["depth"], "mag": raw["magnitude"], "mag_type": raw["type"],
        "source": "afad", "agency": "AFAD", "event_id": raw["eventID"].astype(str),
    })
    return normalize(df)


_AFAD_EXPORT_ALIASES = {
    "time": ["date", "Date", "Tarih", "Olus Zamani", "eventDate", "Origin Time"],
    "latitude": ["latitude", "Latitude", "Enlem"],
    "longitude": ["longitude", "Longitude", "Boylam"],
    "depth_km": ["depth", "Depth", "Derinlik", "Depth(km)", "Derinlik(km)"],
    "mag": ["magnitude", "Magnitude", "Büyüklük", "Buyukluk", "Mag"],
    "mag_type": ["type", "Type", "Tip", "Magnitude Type", "Büyüklük Tipi"],
    "event_id": ["eventID", "EventID", "Event ID", "Deprem Kodu"],
}


def load_afad_export(path: Path) -> pd.DataFrame:
    """AFAD web kataloğundan elle indirilen CSV/XLSX dosyasını oku.

    Sütun adları dışa aktarma biçimine göre değişebildiği için yaygın
    adlar eşlenir; eşlenemeyen zorunlu sütun hata verir. Zamanların UTC
    olduğu varsayılır; dosyada yerel saat varsa önce dönüştürülmeli.
    """
    path = Path(path)
    raw = pd.read_excel(path) if path.suffix.lower() in (".xls", ".xlsx") else pd.read_csv(path)
    cols = {}
    for key, names in _AFAD_EXPORT_ALIASES.items():
        hit = next((n for n in names if n in raw.columns), None)
        if hit is None and key in ("time", "latitude", "longitude", "mag"):
            raise ValueError(f"AFAD dosyasında '{key}' sütunu bulunamadı: {list(raw.columns)}")
        cols[key] = raw[hit] if hit else pd.NA
    df = pd.DataFrame(cols)
    df["time"] = pd.to_datetime(df["time"]).dt.tz_localize("UTC")
    df["source"] = "afad"
    df["agency"] = "AFAD"
    if df["event_id"].isna().all():
        df["event_id"] = [f"afadx{i}" for i in range(len(df))]
    df["event_id"] = df["event_id"].astype(str)
    return normalize(df)


# ------------------------------------------------------------- ISC (AFAD arşivi)
ISC_URL = "https://isc.ac.uk/fdsnws/event/1/query"


def _isc_afad_chunk(region: Region, t0, t1, min_mag) -> str:
    """ISC Bülteni'nden yalnızca AFAD'ın kendi çözümleri (contributor=AFAD).

    AFAD'ın ISC'ye gönderdiği çözümler 2010'dan itibaren bulunur; ISC'nin
    gözden geçirme gecikmesi nedeniyle son ~1 yıl eksik olabilir.
    """
    params = {
        "format": "text", "contributor": "AFAD",
        "starttime": t0.strftime("%Y-%m-%dT%H:%M:%S"), "endtime": t1.strftime("%Y-%m-%dT%H:%M:%S"),
        "minlatitude": region.min_lat, "maxlatitude": region.max_lat,
        "minlongitude": region.min_lon, "maxlongitude": region.max_lon,
        "minmagnitude": min_mag, "nodata": 404,
    }
    r = _get(ISC_URL, params, timeout=300)
    return r.text if r.status_code == 200 else ""


def parse_isc(text: str, source: str = "afad") -> pd.DataFrame:
    lines = [ln for ln in text.splitlines()
             if ln.strip() and (ln.startswith("#EventID") or not ln.startswith("#"))]
    if len(lines) < 2:
        return empty_catalog()
    raw = pd.read_csv(io.StringIO("\n".join(lines)), sep="|")
    raw.columns = [c.lstrip("#").strip() for c in raw.columns]
    if "EventType" in raw.columns:
        # patlama, maden patlaması vb. at
        raw = raw[raw["EventType"].fillna("earthquake").str.contains("earthquake")]
    df = pd.DataFrame({
        "time": raw["Time"].astype(str).str.rstrip("Z") + "Z",
        "latitude": raw["Latitude"], "longitude": raw["Longitude"],
        "depth_km": raw["Depth/km"], "mag": raw["Magnitude"], "mag_type": raw["MagType"],
        "source": source, "agency": raw["Author"],
        "event_id": "isc" + raw["EventID"].astype(str),
    })
    return normalize(df)


SOURCES = {
    "usgs": (_usgs_chunk, parse_usgs, "csv"),
    "emsc": (_emsc_chunk, parse_emsc, "txt"),
    "afad": (_afad_chunk, parse_afad, "json"),
    "isc_afad": (_isc_afad_chunk, parse_isc, "txt"),
}


def _fetch_range(source, spec, t0, t1, cache_dir, refresh, now, depth=0):
    chunk_fn, parse_fn, ext = SOURCES[source]
    cache = None
    if cache_dir is not None:
        d = Path(cache_dir) / source / f"{spec.region.name}_M{spec.min_mag:g}"
        d.mkdir(parents=True, exist_ok=True)
        cache = d / f"{t0:%Y%m%dT%H%M%S}_{t1:%Y%m%dT%H%M%S}.{ext}"
    if cache is not None and cache.exists() and not refresh and t1 < now:
        text = cache.read_text()
    else:
        log.info("%s: %s - %s çekiliyor", source, t0, t1)
        text = chunk_fn(spec.region, t0, t1, spec.min_mag)
    part = parse_fn(text)
    if len(part) >= SERVICE_LIMIT and depth < 8:
        # servis sınırına ulaşıldı: dilimi ikiye böl
        mid = t0 + (t1 - t0) / 2
        return pd.concat([
            _fetch_range(source, spec, t0, mid, cache_dir, refresh, now, depth + 1),
            _fetch_range(source, spec, mid, t1, cache_dir, refresh, now, depth + 1)])
    if cache is not None and t1 < now:
        cache.write_text(text)
    return part


SERVICE_LIMIT = 20000


def fetch(source: str, spec: FetchSpec, cache_dir: Path | None = None,
          refresh: bool = False) -> pd.DataFrame:
    """Bir kaynaktan katalog çek; zaman dilimlerini önbellekle.

    Henüz bitmemiş (geleceğe uzanan) dilim önbelleğe yazılmaz. Servis
    sınırına ulaşan dilim otomatik olarak ikiye bölünür.
    """
    now = pd.Timestamp.now(tz="UTC")
    frames = [_fetch_range(source, spec, t0, t1, cache_dir, refresh, now)
              for t0, t1 in _time_chunks(spec.start, spec.end, spec.chunk)]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return empty_catalog()
    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["source", "event_id"])
    return normalize(out)
