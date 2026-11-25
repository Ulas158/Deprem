"""ULF manyetik alan anomalileri — TEST EDİLECEK HİPOTEZ.

Kaynak: INTERMAGNET dakikalık verisi (BGS GIN servisi). Yerel istasyon İznik
(IZN, Boğaziçi Üniversitesi KRDAE), uzak referans Panagyurishte (PAG,
Bulgaristan). İstanbul (ISK) verisi serviste boş olduğundan kullanılmıyor.

Hipotez (Hayakawa vd. 1996; Hattori 2004): deprem hazırlık sürecinde
litosferden yayılan ULF (0,01-0,1 Hz) sinyaller, gece saatlerinde düşey
bileşenin (Z) yatay bileşene (H) oranını yükseltir. Eleştiri: Masci 2011
(jeomanyetik etkinlik ve aletsel kaymalarla açıklanabilir), Thomas vd. 2009
(Loma Prieta sinyali aletsel). Kanıt düzeyi **zayıf**.

Sınırlama: dakikalık veride en yüksek frekans 8,3 mHz'dir; klasik ULF bandı
(10-100 mHz) için saniyelik veri gerekir. Burada 2-8 mHz (Pc5/Pi2 alt bandı)
kullanılır. Bu, hipotezin zayıflatılmış bir sınamasıdır ve raporda belirtilir.

Yöntem: her gün yerel gece (Türkiye UT+2/3 → 20-02 UT) için Z ve H
spektral gücü (Welch), polarizasyon oranı R = Z/H; referans istasyonla
oranlanarak küresel jeomanyetik etkinlik bastırılır (gradyan yöntemi);
anomali = önceki 30 günün ortalamasına göre z-değeri (sızıntı yok).
"""

from __future__ import annotations

import io
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import welch

from .net import download

log = logging.getLogger(__name__)

GIN = "https://imag-data.bgs.ac.uk/GIN_V1/GINServices"
LOCAL, REFERENCE = "IZN", "PAG"
# İznik konumu (anomali bu noktaya atanır; etki yarıçapı özellik tarafında)
STATION_LATLON = {"IZN": (40.50, 29.72), "PAG": (42.51, 24.18)}


def fetch_month(raw: Path, obs: str, month: pd.Timestamp) -> Path | None:
    days = month.days_in_month
    params = {"Request": "GetData", "format": "iaga2002", "testObsys": "0",
              "observatoryIagaCode": obs, "samplesPerDay": "minute",
              "publicationState": "best-avail", "dataStartDate": f"{month:%Y-%m-%d}",
              "dataDuration": str(days)}
    return download(GIN, raw / obs / f"{obs}_{month:%Y%m}.iaga", params=params, min_bytes=1000)


def read_iaga(path: Path) -> pd.DataFrame:
    lines = [ln for ln in Path(path).read_text(errors="ignore").splitlines()
             if ln[:2] in ("19", "20") and ln[4] == "-"]
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=r"\s+", header=None,
                     names=["date", "time", "doy", "X", "Y", "Z", "F"])
    df["t"] = pd.to_datetime(df["date"] + " " + df["time"])
    for c in "XYZF":
        df.loc[df[c] >= 88888, c] = np.nan
    df["H"] = np.hypot(df["X"], df["Y"])
    return df[["t", "H", "Z"]].set_index("t")


def _night_power(df: pd.DataFrame, day: pd.Timestamp, band=(0.002, 0.008)):
    s = df.loc[day - pd.Timedelta(hours=4): day + pd.Timedelta(hours=2) - pd.Timedelta(minutes=1)]
    if len(s) < 300 or s.isna().mean().max() > 0.2:
        return np.nan, np.nan
    out = []
    for c in ("Z", "H"):
        x = s[c].interpolate(limit_direction="both").to_numpy()
        x = x - np.polyval(np.polyfit(np.arange(len(x)), x, 1), np.arange(len(x)))
        f, p = welch(x, fs=1 / 60.0, nperseg=128)
        sel = (f >= band[0]) & (f <= band[1])
        out.append(np.trapezoid(p[sel], f[sel]))
    return out[0], out[1]


def build_daily(raw: Path, start: pd.Timestamp, end: pd.Timestamp, window: int = 30) -> pd.DataFrame:
    months = pd.date_range(start.to_period("M").to_timestamp(), end, freq="MS")
    series = {}
    for obs in (LOCAL, REFERENCE):
        frames = []
        for m in months:
            p = fetch_month(raw, obs, m)
            if p is not None:
                try:
                    frames.append(read_iaga(p))
                except Exception as e:
                    log.warning("%s %s okunamadı: %s", obs, m.date(), e)
        series[obs] = pd.concat(frames).sort_index() if frames else pd.DataFrame(columns=["H", "Z"])
    days = pd.date_range(start, end - pd.Timedelta(days=1), freq="D")
    rows = []
    for d in days:
        zl, hl = _night_power(series[LOCAL], d)
        zr, hr = _night_power(series[REFERENCE], d)
        rows.append({"date": d, "zh_local": zl / hl if hl else np.nan,
                     "zh_ref": zr / hr if hr else np.nan, "z_pow_local": zl, "z_pow_ref": zr})
    df = pd.DataFrame(rows)
    # gradyan: yerel oran / referans oran (log), ve yerel Z gücü / referans Z gücü
    df["log_zh_rel"] = np.log10(df["zh_local"]) - np.log10(df["zh_ref"])
    df["log_z_rel"] = np.log10(df["z_pow_local"]) - np.log10(df["z_pow_ref"])
    for c in ("log_zh_rel", "log_z_rel"):
        past = df[c].shift(1).rolling(window, min_periods=15)
        df[c.replace("log_", "ulf_") + "_z"] = (df[c] - past.mean()) / past.std()
    lat, lon = STATION_LATLON[LOCAL]
    df["lat"], df["lon"] = lat, lon
    return df
