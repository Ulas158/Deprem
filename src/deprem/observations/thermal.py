"""Yüzey sıcaklığı (termal) anomalileri — TEST EDİLECEK HİPOTEZ.

Kaynak: NASA POWER günlük ``TS_MIN`` (yer yüzeyi "deri" sıcaklığı günlük
en düşük, yaklaşık gece) ve ``T2M_MIN`` (2 m hava). POWER, MERRA-2
yeniden analizine dayanır (≈0,5° × 0,625°); doğrudan uydu kızılötesi
ölçümü değildir. Gerçek uydu LST (MODIS) NASA Earthdata hesabı ister; o
bağlanınca aynı arayüzle değiştirilebilir.

Hipotez (Tronin 2006; Tramutoli vd. 2005 RST/RETIRA; Ouzounov vd.): büyük
depremlerden günler-haftalar önce fay üzerinde gece yüzey sıcaklığı artar
(gaz salımı, gizli ısı). Eleştiri: bulut, mevsim ve hava olaylarıyla
karışma; körlemesine ileriye dönük testlerde tutarlı sonuç yok. Kanıt
düzeyi **zayıf**.

Yöntem (RETIRA benzeri): ΔT(x,t) = TS_MIN(x,t) − bölge ortalaması(t)
(sinoptik hava etkisini çıkarır), ayrıca ΔT − ΔT2M (yüzey-hava farkı).
Anomali: aynı takvim penceresinin (±15 gün) **önceki yıllardaki**
ortalaması ve standart sapmasına göre z. Yalnız geçmiş yıllar → sızıntı yok.
"""

from __future__ import annotations

import json
import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ..region import Region
from .net import download

log = logging.getLogger(__name__)
POWER = "https://power.larc.nasa.gov/api/temporal/daily/point"


def grid_points(region: Region, step: float = 1.0):
    lats = np.arange(region.min_lat + step / 2, region.max_lat, step)
    lons = np.arange(region.min_lon + step / 2, region.max_lon, step)
    return [(round(a, 2), round(b, 2)) for a in lats for b in lons]


def fetch_point(raw: Path, lat: float, lon: float, start: pd.Timestamp, end: pd.Timestamp):
    params = {"parameters": "TS_MIN,T2M_MIN", "community": "AG", "latitude": lat, "longitude": lon,
              "start": f"{start:%Y%m%d}", "end": f"{end:%Y%m%d}", "format": "JSON"}
    p = download(POWER, raw / f"{lat:.2f}_{lon:.2f}_{start:%Y%m%d}_{end:%Y%m%d}.json", params=params,
                 min_bytes=500, timeout=300)
    if p is None:
        return None
    d = json.loads(p.read_text())["properties"]["parameter"]
    df = pd.DataFrame({k: pd.Series(v) for k, v in d.items()})
    df.index = pd.to_datetime(df.index, format="%Y%m%d")
    df = df.replace(-999.0, np.nan)
    df["lat"], df["lon"] = lat, lon
    return df


def build_daily(raw: Path, region: Region, start: pd.Timestamp, end: pd.Timestamp,
                half_window: int = 15, min_years: int = 3) -> pd.DataFrame:
    pts = grid_points(region)
    # Önceki yıllardan iklim referansı için başlangıcı 5 yıl geri al
    s0 = start - pd.DateOffset(years=5)
    with ThreadPoolExecutor(4) as ex:
        frames = [f for f in ex.map(lambda p: fetch_point(raw, p[0], p[1], s0, end), pts) if f is not None]
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames).rename_axis("date").reset_index()
    ts = df.pivot(index="date", columns=["lat", "lon"], values="TS_MIN")
    ta = df.pivot(index="date", columns=["lat", "lon"], values="T2M_MIN")
    dts = ts.sub(ts.mean(axis=1), axis=0)         # bölgesel sinoptik etkiyi çıkar
    dsa = (ts - ta).sub((ts - ta).mean(axis=1), axis=0)
    out = []
    doy = dts.index.dayofyear.to_numpy()
    yr = dts.index.year.to_numpy()
    for name, M in (("thermal_dts", dts), ("thermal_dsa", dsa)):
        A = M.to_numpy()
        Z = np.full_like(A, np.nan)
        for k in range(len(A)):
            if dts.index[k] < start:
                continue
            dd = np.abs(((doy - doy[k]) + 182) % 365 - 182)
            ref = (yr < yr[k]) & (dd <= half_window)
            if len(np.unique(yr[ref])) < min_years:
                continue
            mu, sd = np.nanmean(A[ref], axis=0), np.nanstd(A[ref], axis=0)
            Z[k] = (A[k] - mu) / np.where(sd > 0, sd, np.nan)
        z = pd.DataFrame(Z, index=M.index, columns=M.columns).stack(["lat", "lon"], future_stack=True)
        out.append(z.rename(name + "_z"))
    res = pd.concat(out, axis=1).reset_index()
    return res[res["date"] >= start].dropna(how="all", subset=["thermal_dts_z", "thermal_dsa_z"])
