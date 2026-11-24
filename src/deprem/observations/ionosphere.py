"""İyonosfer toplam elektron içeriği (TEC) anomalileri — TEST EDİLECEK HİPOTEZ.

Kaynak: CODE (Bern Üniversitesi) küresel iyonosfer haritaları (GIM, IONEX),
2,5° enlem × 5° boylam, 1-2 saatlik. Jeomanyetik etkinlik: GFZ Kp indeksi.

Hipotez (Pulinets & Ouzounov 2011; Liu vd. 2006, 2009): büyük depremlerden
1-10 gün önce episantr üzerinde TEC'de olağan dışı artış/azalış olur.
Eleştiri: Dautermann vd. 2007 (JGR), Masci vd. 2017 (zamansal kümelenme ve
güneş/jeomanyetik etkinlikle karışma); ileriye dönük testlerde tutarlı
sonuç yok. Bu yüzden kanıt düzeyi **zayıf** sayılır.

Bastırılan karışıklıklar:
* Jeomanyetik fırtına günleri (günlük en büyük Kp ≥ 4) dışarıda bırakılır.
* Güneş etkinliği ve mevsim etkisini azaltmak için her düğümün anomalisinden
  aynı enlemdeki tüm boylamların (kuşak) ortalama anomalisi çıkarılır:
  yerel sinyal = yerel − küresel kuşak.

Anomali tanımı (Liu vd. 2006): her epok için önceki 15 günün aynı saatindeki
ortanca M ve çeyrekler arası açıklık IQR; z = (TEC − M) / (IQR/1.349).
Yalnız ``t`` gününden önceki veriler kullanılır (sızıntı yok).
"""

from __future__ import annotations

import gzip
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ..region import Region
from .net import download

log = logging.getLogger(__name__)

CODE = "https://download.aiub.unibe.ch/CODE"
KP_URL = "https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt"


def _urls(day: pd.Timestamp) -> list[tuple[str, str]]:
    y, doy = day.year, day.dayofyear
    new = f"COD0OPSFIN_{y}{doy:03d}0000_01D_01H_GIM.INX.gz"
    old = f"CODG{doy:03d}0.{y % 100:02d}I.Z"
    if day >= pd.Timestamp("2022-11-27"):
        return [(f"{CODE}/{y}/{new}", new), (f"{CODE}/{y}/{old}", old)]
    return [(f"{CODE}/{y}/{old}", old), (f"{CODE}/{y}/{new}", new)]


def _read_bytes(path: Path) -> str:
    if path.suffix == ".gz":
        return gzip.decompress(path.read_bytes()).decode("ascii", "ignore")
    # .Z (LZW): gzip -dc çözer
    return subprocess.run(["gzip", "-dc", str(path)], capture_output=True, check=True).stdout.decode(
        "ascii", "ignore")


def parse_ionex(text: str):
    """IONEX'ten TEC haritaları: (epok listesi, enlemler, boylamlar, dizi[epok, enlem, boylam]) TECU."""
    lines = text.splitlines()
    exp = -1
    lat1 = lat2 = dlat = lon1 = lon2 = dlon = None
    i = 0
    for i, ln in enumerate(lines):
        lab = ln[60:].strip()
        if lab == "EXPONENT":
            exp = int(ln[:6])
        elif lab.startswith("LAT1 / LAT2 / DLAT"):
            lat1, lat2, dlat = map(float, ln[:60].split())
        elif lab.startswith("LON1 / LON2 / DLON"):
            lon1, lon2, dlon = map(float, ln[:60].split())
        elif lab == "END OF HEADER":
            break
    lats = np.arange(lat1, lat2 + dlat / 2, dlat) if dlat < 0 else np.arange(lat1, lat2 + dlat / 2, dlat)
    lons = np.arange(lon1, lon2 + dlon / 2, dlon)
    nlon = len(lons)
    maps, epochs = [], []
    cur = None
    row = -1
    buf: list[float] = []
    k = i + 1
    while k < len(lines):
        ln = lines[k]
        lab = ln[60:].strip()
        if cur is not None and row >= 0:
            buf += [float(ln[j:j + 5]) for j in range(0, len(ln.rstrip()), 5) if ln[j:j + 5].strip()]
            if len(buf) >= nlon:
                v = np.array(buf[:nlon])
                v[v == 9999] = np.nan
                cur[row] = v * 10.0 ** exp
                row = -1
        elif lab == "START OF TEC MAP":
            cur = np.full((len(lats), nlon), np.nan)
        elif lab == "EPOCH OF CURRENT MAP" and cur is not None:
            f = list(map(int, ln[:36].split()))
            epochs.append(pd.Timestamp(*f[:6]))
        elif lab.startswith("LAT/LON1/LON2/DLON/H") and cur is not None:
            lat = float(ln[2:8])
            row = int(round((lat - lat1) / dlat))
            buf = []
        elif lab == "END OF TEC MAP" and cur is not None:
            maps.append(cur)
            cur = None
        elif lab.startswith("START OF RMS MAP"):
            break
        k += 1
    return epochs, lats, lons, np.array(maps)


def fetch_kp(raw: Path) -> pd.DataFrame:
    p = download(KP_URL, raw / "Kp_since_1932.txt")
    rows = []
    for ln in p.read_text().splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.split()
        kp = [float(x) for x in f[7:15]]
        rows.append({"date": pd.Timestamp(int(f[0]), int(f[1]), int(f[2])), "kp_max": max(kp),
                     "ap": float(f[23]), "f107": float(f[25])})
    return pd.DataFrame(rows)


def _day_nodes(raw: Path, day: pd.Timestamp, region: Region):
    """Bir gün için bölge düğümleri ve enlem kuşağı ortalamaları (saatlik, 0-23 UT)."""
    cache = raw / "parsed" / f"{day:%Y%m%d}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["local"], z["zonal"], z["lats"], z["lons"]
    path = None
    for url, name in _urls(day):
        path = download(url, raw / "ionex" / str(day.year) / name, min_bytes=1000)
        if path is not None:
            break
    if path is None:
        return None
    try:
        epochs, lats, lons, maps = parse_ionex(_read_bytes(path))
    except Exception as e:
        log.warning("IONEX %s okunamadı: %s", day.date(), e)
        return None
    if maps.size == 0:
        return None
    hours = np.array([(e - day).total_seconds() / 3600 for e in epochs])
    target = np.arange(24)
    li = np.where((lats >= region.min_lat - 2.5) & (lats <= region.max_lat + 2.5))[0]
    lo = np.where((lons >= region.min_lon - 5) & (lons <= region.max_lon + 5))[0]
    # saatliğe enterpole et (eski dosyalar 2 saatlik)
    def interp(series):
        return np.array([np.interp(target, hours, s) for s in series])
    sub = maps[:, li][:, :, lo]  # epok, enlem, boylam
    local = interp(sub.reshape(len(hours), -1).T).T.reshape(24, len(li), len(lo))
    zon = np.nanmean(maps[:, li, :], axis=2)  # epok, enlem
    zonal = interp(zon.T).T
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, local=local, zonal=zonal, lats=lats[li], lons=lons[lo])
    return local, zonal, lats[li], lons[lo]


def build_daily(raw: Path, region: Region, start: pd.Timestamp, end: pd.Timestamp,
                window: int = 15, workers: int = 8) -> pd.DataFrame:
    """Gün × düğüm TEC anomali özellikleri.

    Sütunlar: date, lat, lon, tec_z_max (günün en büyük pozitif z'si, yerel −
    kuşak), tec_z_min (en büyük negatif), tec_frac_pos / tec_frac_neg (|z|>2
    saat oranı), kp_max, disturbed.
    """
    days = pd.date_range(start, end - pd.Timedelta(days=1), freq="D")
    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(lambda d: _day_nodes(raw, d, region), days))
    ok = [r for r in res if r is not None]
    if not ok:
        return pd.DataFrame()
    lats, lons = ok[0][2], ok[0][3]
    nl, nlo = len(lats), len(lons)
    L = np.full((len(days), 24, nl, nlo), np.nan)
    Z = np.full((len(days), 24, nl), np.nan)
    for k, r in enumerate(res):
        if r is not None and r[0].shape == (24, nl, nlo):
            L[k], Z[k] = r[0], r[1]
    kp = fetch_kp(raw.parent / "kp").set_index("date")["kp_max"]
    kpd = kp.reindex(days).to_numpy()
    rows = []
    for k in range(window, len(days)):
        past = L[k - window:k]
        pz = Z[k - window:k]
        # fırtına günlerini referanstan çıkar
        good = ~(kpd[k - window:k] >= 4)
        if good.sum() < 8:
            continue
        med = np.nanmedian(past[good], axis=0)
        iqr = np.nanpercentile(past[good], 75, axis=0) - np.nanpercentile(past[good], 25, axis=0)
        zmed = np.nanmedian(pz[good], axis=0)
        ziqr = np.nanpercentile(pz[good], 75, axis=0) - np.nanpercentile(pz[good], 25, axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            zl = (L[k] - med) / (iqr / 1.349 + 0.5)
            zz = (Z[k] - zmed) / (ziqr / 1.349 + 0.5)
        z = zl - zz[:, :, None]  # yerel − kuşak
        for a in range(nl):
            for b in range(nlo):
                s = z[:, a, b]
                if np.isnan(s).all():
                    continue
                rows.append((days[k], lats[a], lons[b], np.nanmax(s), np.nanmin(s),
                             np.nanmean(s > 2), np.nanmean(s < -2), kpd[k]))
    df = pd.DataFrame(rows, columns=["date", "lat", "lon", "tec_z_max", "tec_z_min",
                                     "tec_frac_pos", "tec_frac_neg", "kp_max"])
    df["disturbed"] = df["kp_max"] >= 4
    return df
