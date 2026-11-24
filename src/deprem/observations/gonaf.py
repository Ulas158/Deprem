"""GONAF kuyu gerinimölçerleri (Marmara).

Durum: **veri bekliyor.** GONAF (Geophysical Observatory at the North
Anatolian Fault; GFZ Potsdam + AFAD, ICDP) Marmara çevresinde 5 kuyuda
nanogerinim duyarlılıkta gerinim ve kuyu içi sismometre ölçer
(Bohnhoff vd. 2017, Sci. Drilling). Veriler açık FDSN'de değil; GFZ'ye
veri talebiyle alınır.

Veri gelince ``data/raw/gonaf/`` altına:
* ``stations.csv``: ``station,lat,lon``
* her istasyon için ``<station>.csv``: ``time,<kanal1>,<kanal2>,...``
  (gerinim, nanogerinim; gelgit ve basınç düzeltmesi yapılmış olması
  tercih edilir).

Özellik ``gonaf_strain_anom``: ``t0`` öncesi son 1 yılın gerinim hızının,
önceki 3 yılın doğrusal eğiliminden sapmasının gürültüye oranı (kanalların
en büyüğü), 50 km içindeki hücrelere atanır.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load(raw: Path):
    raw = Path(raw)
    if not (raw / "stations.csv").exists():
        return None
    st = pd.read_csv(raw / "stations.csv")
    series = {}
    for s in st["station"]:
        p = raw / f"{s}.csv"
        if p.exists():
            d = pd.read_csv(p, parse_dates=["time"]).set_index("time").resample("1D").mean()
            series[s] = d
    return st, series


def cell_features(raw: Path, cells: pd.DataFrame, t0: pd.Timestamp, radius_km: float = 50.0) -> pd.DataFrame:
    idx = cells["cell_id"]
    data = load(raw)
    if data is None:
        return pd.DataFrame(index=idx)
    st, series = data
    t0n = t0.tz_localize(None) if t0.tzinfo else t0
    vals = {}
    for s, d in series.items():
        w = d[(d.index < t0n) & (d.index >= t0n - pd.DateOffset(years=4))]
        if len(w) < 600:
            continue
        tt = (w.index - t0n).days.to_numpy() / 365.25
        early = tt < -1
        best = 0.0
        for c in w.columns:
            y = w[c].to_numpy(float)
            ok = np.isfinite(y)
            if (ok & early).sum() < 300 or (ok & ~early).sum() < 150:
                continue
            p = np.polyfit(tt[ok & early], y[ok & early], 1)
            res = y - np.polyval(p, tt)
            sl = np.polyfit(tt[ok & ~early], res[ok & ~early], 1)[0]
            noise = np.nanstd(np.diff(res[ok & early])) * np.sqrt(365)
            best = max(best, abs(sl) / max(noise, 1e-9))
        vals[s] = best
    out = np.full(len(cells), np.nan)
    for s, v in vals.items():
        r = st[st["station"] == s].iloc[0]
        d = np.hypot((cells["lon_c"] - r["lon"]) * 85.0, (cells["lat_c"] - r["lat"]) * 111.0).to_numpy()
        sel = d <= radius_km
        out[sel] = np.fmax(out[sel], v)
    return pd.DataFrame({"gonaf_strain_anom": out}, index=idx)
