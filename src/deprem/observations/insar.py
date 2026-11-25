"""Uydu radarı (InSAR) yüzey deformasyonu.

Durum: **veri bekliyor.** COMET LiCSAR açık arşivi (gws-access.jasmin.ac.uk)
Türkiye çerçeveleri için yalnız ham interferogramları ve plaka hareketi
modelini (``vlos_eur``, gözlem değil) sunuyor; gözlenmiş hız haritası yok.
Ham interferogramlardan hız üretmek (LiCSBAS) çerçeve başına onlarca GB
indirme ve saatlerce hesap ister.

Bu modül iki hazır kaynağı okur, ikisi de ``data/raw/insar/`` altına konur:

1. **LiCSBAS çıktısı** veya herhangi bir görüş doğrultusu (LOS) hız GeoTIFF'i
   (``*.vel*.tif``, mm/yıl).
2. **Weiss vd. 2020 (GRL) Anadolu hız alanı** Doğu/Kuzey GeoTIFF'leri
   (``*E*.tif`` ve ``*N*.tif``) — gerinim hızı hesaplanır.

Özellikler (hücre başına, zamandan bağımsız):
* ``insar_los_rate_mm_yr``: hücredeki ortanca LOS hızı
* ``insar_strain_rate_ns``: Doğu/Kuzey hız alanının gradyanından ikinci
  değişmez gerinim hızı (nanogerinim/yıl)

Sızıntı uyarısı: hız haritaları genellikle 2015-2020+ verisinden üretilir;
sınama yıllarıyla çakışır. Raporda ayrıca işaretlenir.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def _sample_cells(path: Path, cells: pd.DataFrame, cell_deg: float, stat="median"):
    import rasterio
    out = np.full(len(cells), np.nan)
    with rasterio.open(path) as r:
        for k, c in enumerate(cells.itertuples()):
            win = rasterio.windows.from_bounds(c.lon0, c.lat0, c.lon0 + cell_deg, c.lat0 + cell_deg, r.transform)
            try:
                a = r.read(1, window=win, boundless=True, fill_value=np.nan).astype(float)
            except Exception:
                continue
            if r.nodata is not None:
                a[a == r.nodata] = np.nan
            a[a == 0] = np.nan
            if np.isfinite(a).sum() > 20:
                out[k] = np.nanmedian(a) if stat == "median" else np.nanmax(a)
    return out


def _strain_map(e_path: Path, n_path: Path, cells, cell_deg):
    import rasterio
    with rasterio.open(e_path) as re_, rasterio.open(n_path) as rn:
        ve, vn = re_.read(1).astype(float), rn.read(1).astype(float)
        tr = re_.transform
    for a in (ve, vn):
        a[(a == 0) | ~np.isfinite(a)] = np.nan
    lat_c = tr.f + tr.e * (np.arange(ve.shape[0]) + 0.5)
    dx_km = abs(tr.a) * 111.195 * np.cos(np.radians(lat_c))[:, None]
    dy_km = abs(tr.e) * 111.195
    dve_dy, dve_dx = np.gradient(ve)
    dvn_dy, dvn_dx = np.gradient(vn)
    exx, eyy = dve_dx / dx_km, -dvn_dy / dy_km
    exy = 0.5 * (-dve_dy / dy_km + dvn_dx / dx_km)
    inv2 = np.sqrt(exx ** 2 + eyy ** 2 + 2 * exy ** 2) * 1e3  # mm/yıl/km → nanogerinim/yıl
    out = np.full(len(cells), np.nan)
    for k, c in enumerate(cells.itertuples()):
        r0, c0 = ~tr * (c.lon0, c.lat0 + cell_deg)
        r1, c1 = ~tr * (c.lon0 + cell_deg, c.lat0)
        a = inv2[int(max(c0, 0)):int(max(c1, 0)), int(max(r0, 0)):int(max(r1, 0))]
        if a.size and np.isfinite(a).sum() > 20:
            out[k] = np.nanmedian(a)
    return out


def cell_features(raw: Path, cells: pd.DataFrame) -> pd.DataFrame:
    raw = Path(raw)
    cell_deg = cells.attrs.get("cell_deg", 0.5)
    df = pd.DataFrame(index=cells["cell_id"])
    if not raw.exists():
        return df
    vel = sorted(raw.glob("*vel*.tif"))
    if vel:
        vals = np.vstack([_sample_cells(p, cells, cell_deg) for p in vel])
        with np.errstate(all="ignore"):
            df["insar_los_rate_mm_yr"] = np.nanmean(vals, axis=0)
    e = sorted(p for p in raw.glob("*.tif") if "_E" in p.stem or p.stem.endswith("E"))
    n = sorted(p for p in raw.glob("*.tif") if "_N" in p.stem or p.stem.endswith("N"))
    if e and n:
        df["insar_strain_rate_ns"] = _strain_map(e[0], n[0], cells, cell_deg)
    log.info("InSAR: %d hız, %d D/K dosyası", len(vel), min(len(e), len(n)))
    return df
