"""Ortam gürültüsünden sismik hız değişimi (dv/v).

Kabuktaki sismik dalga hızı gerilmeye ve gözenek basıncına duyarlıdır:
gerilme artınca çatlaklar kapanır ve hız artar; hasar ve akışkan basıncı
artınca hız düşer. Sürekli yer gürültüsünün kendi kendisiyle (tek istasyon,
bileşenler arası) korelasyonundaki kod dalgalarının zamanla "esnemesi"
bu değişimi 10⁻⁴ duyarlılıkla ölçer (Brenguier vd. 2008, Science, Parkfield;
Wegler & Sens-Schönfelder 2007; Hobiger vd. 2012 tek istasyon yöntemi).
Deprem öncesi hız düşüşü raporları vardır (ör. Chen vd. 2010), ancak
genellemesi sınırlıdır; büyük depremlerin eş-sismik hız düşüşü ise iyi
bilinir. Bu özellik **güçlü fiziksel temelli** sayılır, öncül etkisi temkinli
tutulur.

Uygulama (hafif sürüm, bu ortamda çalışabilir boyutta):
* KOERI (KO ağı) FDSN'den her ``step_days`` günde bir, yerel gece 1 saatlik
  3 bileşen kayıt.
* 0,5-2 Hz bant, tek bit normalizasyon, spektral beyazlatma.
* ZN, ZE, NE çapraz korelasyonları (±20 s), 5 örneklik kayan yığın.
* Referans: istasyonun **ilk iki yılının** yığını (eğitim dönemi).
* Esneme (stretching) yöntemiyle 2-15 s kod penceresinde dv/v.
Ham kayıt işlendikten sonra silinir; yalnız korelasyonlar saklanır.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from ..region import Region
from .net import UA

log = logging.getLogger(__name__)

FDSN = "https://eida.koeri.boun.edu.tr/fdsnws"
FS = 20.0
MAXLAG = 20.0


def stations(region: Region, start_before: str = "2011-06-01") -> pd.DataFrame:
    r = requests.get(f"{FDSN}/station/1/query", params={"network": "KO", "level": "channel",
                     "channel": "HH?,BH?", "format": "text"}, headers=UA, timeout=120)
    r.raise_for_status()
    rows = [ln.split("|") for ln in r.text.splitlines() if not ln.startswith("#") and ln.strip()]
    df = pd.DataFrame(rows).iloc[:, [0, 1, 2, 3, 4, 5, 6, 7, 8, 15, 16]]
    df.columns = ["net", "sta", "loc", "cha", "lat", "lon", "elev", "depth", "azimuth", "start", "end"]
    df["lat"], df["lon"] = df["lat"].astype(float), df["lon"].astype(float)
    df["start"] = pd.to_datetime(df["start"], errors="coerce")
    g = df.groupby("sta").agg(lat=("lat", "first"), lon=("lon", "first"), start=("start", "min"),
                              chans=("cha", lambda c: ",".join(sorted(set(c)))))
    g = g[(g["start"] <= pd.Timestamp(start_before)) & g["chans"].str.contains("HHZ|BHZ")]
    g = g[(g["lat"].between(region.min_lat, region.max_lat)) & (g["lon"].between(region.min_lon, region.max_lon))]
    return g.reset_index()


def pick_spread(st: pd.DataFrame, n: int, prefer: list[str] | None = None) -> pd.DataFrame:
    """Uzaklığı en büyükleyen açgözlü seçim (en uzak nokta örneklemesi)."""
    prefer = [p for p in (prefer or []) if p in set(st["sta"])]
    chosen = list(prefer) or [st.iloc[0]["sta"]]
    xy = np.column_stack([st["lon"] * 0.77, st["lat"]])
    while len(chosen) < min(n, len(st)):
        cidx = st.index[st["sta"].isin(chosen)]
        d = np.min(np.linalg.norm(xy[:, None, :] - xy[cidx][None, :, :], axis=2), axis=1)
        d[cidx] = -1
        chosen.append(st.iloc[int(np.argmax(d))]["sta"])
    return st[st["sta"].isin(chosen)].reset_index(drop=True)


def _get_hour(sta: str, t: pd.Timestamp):
    from obspy import read
    params = {"network": "KO", "station": sta, "channel": "HH?,BH?",
              "starttime": t.strftime("%Y-%m-%dT%H:%M:%S"),
              "endtime": (t + pd.Timedelta(hours=1)).strftime("%Y-%m-%dT%H:%M:%S")}
    for k in range(3):
        try:
            r = requests.get(f"{FDSN}/dataselect/1/query", params=params, headers=UA, timeout=120)
            if r.status_code == 204 or not r.content:
                return None
            r.raise_for_status()
            return read(BytesIO(r.content))
        except Exception:
            if k == 2:
                return None
    return None


def _prep(st):
    """HH tercih et, 0,5-2 Hz, FS'ye indir, tek bit + beyazlatma; Z,N,E dizileri."""
    if st is None:
        return None
    hh = st.select(channel="HH?")
    st = hh if len(hh) >= 3 else st.select(channel="BH?")
    st.merge(fill_value=0)
    out = {}
    for c in "ZNE12":
        tr = st.select(component=c)
        if len(tr):
            out[c if c in "ZNE" else ("N" if c == "1" else "E")] = tr[0]
    if not all(k in out for k in "ZNE"):
        return None
    n = int(3600 * FS) - 1
    arr = {}
    for k, tr in out.items():
        tr.detrend("demean")
        tr.taper(0.02)
        tr.filter("bandpass", freqmin=0.5, freqmax=2.0, zerophase=True)
        tr.resample(FS) if tr.stats.sampling_rate != FS else None
        x = tr.data[:n].astype(float)
        if len(x) < n * 0.9:
            return None
        x = np.pad(x, (0, n - len(x)))
        x = np.sign(x)  # tek bit
        X = np.fft.rfft(x, 2 * n)
        f = np.fft.rfftfreq(2 * n, 1 / FS)
        band = (f >= 0.5) & (f <= 2.0)
        X = np.where(band, X / (np.abs(X) + 1e-9), 0)  # beyazlatma
        arr[k] = X
    return arr, n


def _ccf(arr, n):
    lag = int(MAXLAG * FS)
    res = {}
    for a, b in (("Z", "N"), ("Z", "E"), ("N", "E")):
        c = np.fft.irfft(arr[a] * np.conj(arr[b]), 2 * n)
        c = np.concatenate([c[-lag:], c[:lag + 1]])
        res[a + b] = c / (np.abs(c).max() + 1e-12)
    return res


def station_ccfs(sta: str, out_dir: Path, start: pd.Timestamp, end: pd.Timestamp,
                 step_days: int = 5, hour_ut: int = 0) -> Path:
    path = out_dir / f"ccf_{sta}.npz"
    done = {}
    if path.exists():
        z = np.load(path, allow_pickle=True)
        done = {pd.Timestamp(t): {k: z[k][i] for k in ("ZN", "ZE", "NE")} for i, t in enumerate(z["t"])}
    times = [t for t in pd.date_range(start, end, freq=f"{step_days}D") + pd.Timedelta(hours=hour_ut)
             if t not in done]
    for i, t in enumerate(times):
        p = _prep(_get_hour(sta, t))
        if p is not None:
            done[t] = _ccf(*p)
        if (i + 1) % 100 == 0 or i == len(times) - 1:
            ts = sorted(done)
            np.savez_compressed(path, t=np.array([str(x) for x in ts]),
                                **{k: np.array([done[x][k] for x in ts]) for k in ("ZN", "ZE", "NE")})
            log.info("dv/v %s: %d/%d saat işlendi", sta, i + 1, len(times))
    return path


def stretching(ref, cur, t_lag, tmin=2.0, tmax=15.0, eps=0.02, n=201):
    """Esneme: cur(t) ≈ ref(t·(1+ε)). Hız artınca varışlar erken gelir → ε = dv/v (pozitif = hız artışı)."""
    w = (np.abs(t_lag) >= tmin) & (np.abs(t_lag) <= tmax)
    best, bcc = 0.0, -1
    for e in np.linspace(-eps, eps, n):
        r = np.interp(t_lag * (1 + e), t_lag, ref)
        cc = np.corrcoef(r[w], cur[w])[0, 1]
        if cc > bcc:
            best, bcc = e, cc
    return best, bcc


def station_dvv(path: Path, ref_years: float = 2.0, stack: int = 5) -> pd.DataFrame:
    z = np.load(path, allow_pickle=True)
    t = pd.to_datetime(z["t"])
    lag = np.arange(-int(MAXLAG * FS), int(MAXLAG * FS) + 1) / FS
    rows = []
    comps = {}
    for k in ("ZN", "ZE", "NE"):
        C = z[k]
        ref_sel = t < t[0] + pd.DateOffset(years=int(ref_years))
        ref = C[ref_sel].mean(0)
        S = pd.DataFrame(C).rolling(stack, min_periods=3).mean().to_numpy()
        dv = np.full(len(t), np.nan)
        cc = np.full(len(t), np.nan)
        for i in range(len(t)):
            if np.isfinite(S[i]).all():
                dv[i], cc[i] = stretching(ref, S[i], lag)
        comps[k] = (dv, cc)
    dv = np.nanmean([np.where(comps[k][1] > 0.5, comps[k][0], np.nan) for k in comps], axis=0)
    return pd.DataFrame({"time": t, "dvv": dv})


def run(raw: Path, obs: Path, region: Region, max_stations: int, start: pd.Timestamp,
        end: pd.Timestamp, step_days: int = 5, workers: int = 8):
    st = stations(region)
    # Marmara ve Doğu Anadolu öncelikli
    st = pick_spread(st, max_stations, prefer=["ADVT", "GEML", "MRMX", "SVRH", "KMRS", "ERZN"])
    out = obs / "dvv"
    out.mkdir(parents=True, exist_ok=True)
    st.to_csv(out / "stations.csv", index=False)
    with ThreadPoolExecutor(workers) as ex:
        paths = list(ex.map(lambda s: station_ccfs(s, out, start, end, step_days), st["sta"]))
    frames = []
    for s, p in zip(st.itertuples(), paths):
        d = station_dvv(p)
        d["sta"], d["lat"], d["lon"] = s.sta, s.lat, s.lon
        frames.append(d)
    res = pd.concat(frames)
    res.to_csv(obs / "dvv_series.csv.gz", index=False)
    return obs / "dvv_series.csv.gz"


def cell_features(dvv: pd.DataFrame, cells: pd.DataFrame, t0: pd.Timestamp,
                  radius_km: float = 60.0) -> pd.DataFrame:
    """Son 6 ayın dv/v ortalaması − önceki 2 yılın ortalaması (yüzde), yalnız t0 öncesi."""
    idx = cells["cell_id"]
    if dvv is None or dvv.empty:
        return pd.DataFrame(index=idx)
    rows = []
    for s, g in dvv.groupby("sta"):
        g = g[g["time"] < t0]
        recent = g[g["time"] >= t0 - pd.DateOffset(months=6)]["dvv"]
        before = g[(g["time"] < t0 - pd.DateOffset(months=6)) & (g["time"] >= t0 - pd.DateOffset(months=30))]["dvv"]
        if recent.notna().sum() >= 10 and before.notna().sum() >= 30:
            rows.append((g["lat"].iloc[0], g["lon"].iloc[0], 100 * (recent.mean() - before.mean())))
    out = np.full(len(cells), np.nan)
    if rows:
        a = np.array(rows)
        dx = (cells["lon_c"].to_numpy()[:, None] - a[:, 1][None]) * 85.0
        dy = (cells["lat_c"].to_numpy()[:, None] - a[:, 0][None]) * 111.0
        d = np.hypot(dx, dy)
        w = np.where(d <= radius_km, np.exp(-(d / (radius_km / 2)) ** 2), 0)
        ok = w.sum(1) > 0
        out[ok] = (w[ok] * a[:, 2][None]).sum(1) / w[ok].sum(1)
    return pd.DataFrame({"dvv_change_pct": out}, index=idx)
