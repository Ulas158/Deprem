"""İstasyon çiftleri arası ortam gürültüsü korelasyonu ile sismik hız değişimi (dv/v).

PR 9'daki tek istasyon yöntemi (gecede 1 saat, bileşenler arası korelasyon)
ölçümlerin yalnız %11'inde güvenilir sonuç verdi; ölçüm gürültüsü (%0,2-0,7)
beklenen öncül etkiden (~%0,1) büyüktü. Bu modül standart yolu izler
(Brenguier vd. 2008; Bensen vd. 2007; Weaver vd. 2011):

* Aynı bölgedeki istasyonların düşey bileşenleri (Z) her örnek günde birlikte
  indirilir (tek FDSN isteği; KOERI KO ağı).
* 0,1-1 Hz bant (ikincil mikrosismik + yerel gürültü), 5 Hz'e seyreltme,
  30 dakikalık pencerelerde tek bit + spektral beyazlatma.
* Her istasyon çifti için ZZ çapraz korelasyonu (±180 s), günlük yığın.
* Geriye dönük 60 günlük kayan yığın (yalnız geçmiş: tahmin anında kullanılabilir),
  nedensel ve nedensel olmayan yanlar ortalanır (saat hatasına dayanıklılık).
* Yüzey dalgası varışından başlayan 40 s'lik erken kod penceresinde esneme ile dv/v,
  Weaver vd. (2011) formülüyle her ölçümün belirsizliği.
* Referans: çiftin ilk iki yılının yığını (yalnız geçmiş veri).
* Bölgesel seri: çiftlerin ters varyans ağırlıklı ortalaması.

Ham kayıt işlendikten sonra silinir; yalnız günlük korelasyonlar saklanır.
Hesap yükü: veri indirme baskındır (yaklaşık 20 MB / örnek gün / 12 istasyon).
"""

from __future__ import annotations

import itertools
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from .net import UA

log = logging.getLogger(__name__)

FDSN = "https://eida.koeri.boun.edu.tr/fdsnws"
FS = 5.0
FMIN, FMAX = 0.1, 1.0
WIN_S = 1800
MAXLAG = 180.0
NLAG = int(MAXLAG * FS)

# Doğu Marmara pilotu: 2009'dan beri kesintisiz geniş bantlı KOERI istasyonları
# (İzmit-Çınarcık-Adalar kesimi, 1999 kırığının batı ucu ve Marmara sismik boşluğu).
MARMARA_EAST = ["ADVT", "ARMT", "YLV", "HRTX", "MDNY", "CAVI", "SILT", "BGKT",
                "KLYT", "ISK", "GULT", "EDC"]


def lags() -> np.ndarray:
    return np.arange(-NLAG, NLAG + 1) / FS


def station_table(stas: list[str]) -> pd.DataFrame:
    r = requests.get(f"{FDSN}/station/1/query", params={"network": "KO", "station": ",".join(stas),
                     "level": "station", "format": "text"}, headers=UA, timeout=120)
    r.raise_for_status()
    rows = [ln.split("|") for ln in r.text.splitlines() if ln.strip() and not ln.startswith("#")]
    df = pd.DataFrame(rows).iloc[:, [1, 2, 3]]
    df.columns = ["sta", "lat", "lon"]
    df[["lat", "lon"]] = df[["lat", "lon"]].astype(float)
    return df.drop_duplicates("sta").set_index("sta").loc[stas].reset_index()


def dist_km(lat1, lon1, lat2, lon2) -> float:
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return float(6371 * 2 * np.arcsin(np.sqrt(a)))


def pairs_of(stas: list[str]) -> list[tuple[str, str]]:
    return list(itertools.combinations(stas, 2))


# ------------------------------------------------------------------ indirme + ön işleme
def fetch(stas: list[str], t0: pd.Timestamp, hours: float, retries: int = 3):
    from obspy import read
    params = {"network": "KO", "station": ",".join(stas), "channel": "HHZ,BHZ",
              "starttime": t0.strftime("%Y-%m-%dT%H:%M:%S"),
              "endtime": (t0 + pd.Timedelta(hours=hours)).strftime("%Y-%m-%dT%H:%M:%S")}
    for k in range(retries):
        try:
            r = requests.get(f"{FDSN}/dataselect/1/query", params=params, headers=UA, timeout=900)
            if r.status_code == 204 or not r.content:
                return None
            r.raise_for_status()
            return read(BytesIO(r.content))
        except Exception as e:  # ağ hatası: tekrar dene
            if k == retries - 1:
                log.warning("indirilemedi %s: %s", t0, e)
    return None


def whiten_windows(x: np.ndarray, valid: np.ndarray) -> np.ndarray | None:
    """x: 5 Hz, 0,1-1 Hz süzülmüş dizi. Pencere başına tek bit + beyazlatılmış spektrum."""
    n = int(WIN_S * FS)
    nw = len(x) // n
    if nw == 0:
        return None
    nfft = 2 * n
    f = np.fft.rfftfreq(nfft, 1 / FS)
    band = (f >= FMIN) & (f <= FMAX)
    taper = np.ones(band.sum())
    k = max(1, int(0.1 * len(taper)))
    taper[:k] = np.hanning(2 * k)[:k]
    taper[-k:] = np.hanning(2 * k)[k:]
    out = np.zeros((nw, len(f)), complex)
    ok = np.zeros(nw, bool)
    for i in range(nw):
        seg = x[i * n:(i + 1) * n]
        if valid[i * n:(i + 1) * n].mean() < 0.98 or np.std(seg) == 0:
            continue
        X = np.fft.rfft(np.sign(seg), nfft)
        Y = np.zeros_like(X)
        Y[band] = X[band] / (np.abs(X[band]) + 1e-12) * taper
        out[i], ok[i] = Y, True
    return np.where(ok[:, None], out, np.nan)


def preprocess(st, stas: list[str], t0: pd.Timestamp, hours: float) -> dict[str, np.ndarray]:
    from obspy import UTCDateTime
    n_total = int(hours * 3600 * FS)
    res = {}
    for s in stas:
        sel = st.select(station=s, channel="HHZ") or st.select(station=s, channel="BHZ")
        if not len(sel):
            continue
        sel = sel.copy()
        try:
            sel.merge(method=1, fill_value=None)
            tr = sel[0]
            mask = np.ma.getmaskarray(tr.data)
            tr.data = np.ma.filled(tr.data.astype(float), 0.0)
            tr.detrend("demean")
            tr.taper(0.01)
            tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
            fac = int(round(tr.stats.sampling_rate / FS))
            tr.data = tr.data[::fac]  # 1 Hz üst sınırlı süzgeçten sonra güvenli
            m = mask[::fac]
            off = int(round((tr.stats.starttime - UTCDateTime(t0.to_pydatetime())) * FS))
        except Exception as e:
            log.debug("%s %s: %s", s, t0, e)
            continue
        x = np.zeros(n_total)
        v = np.zeros(n_total, bool)
        a, b = max(off, 0), min(off + len(tr.data), n_total)
        if b - a <= 0:
            continue
        x[a:b] = tr.data[a - off:b - off]
        v[a:b] = ~m[a - off:b - off]
        # deprem/ani sinyal: pencere içi genliği medyanın 20 katını geçen yerleri geçersiz say
        env = np.abs(x)
        med = np.median(env[v]) if v.any() else 0
        if med > 0:
            v &= env < 20 * med
        w = whiten_windows(x, v)
        if w is not None:
            res[s] = w
    return res


def day_ccfs(spec: dict[str, np.ndarray], pairs: list[tuple[str, str]]):
    n = int(WIN_S * FS)
    C = np.full((len(pairs), 2 * NLAG + 1), np.nan, np.float32)
    nwin = np.zeros(len(pairs), np.int16)
    for j, (a, b) in enumerate(pairs):
        if a not in spec or b not in spec:
            continue
        A, B = spec[a], spec[b]
        ok = np.isfinite(A[:, 0].real) & np.isfinite(B[:, 0].real)
        if ok.sum() < 2:
            continue
        c = np.fft.irfft((A[ok] * np.conj(B[ok])).mean(0), 2 * n)
        c = np.concatenate([c[-NLAG:], c[:NLAG + 1]])
        C[j] = c / (np.abs(c).max() + 1e-12)
        nwin[j] = ok.sum()
    return C, nwin


def process_day(stas, pairs, day: pd.Timestamp, hour_ut: int, hours: float, out: Path) -> int:
    path = out / f"{day:%Y-%m-%d}.npz"
    if path.exists():
        return -1
    st = fetch(stas, day + pd.Timedelta(hours=hour_ut), hours)
    if st is None:
        C, nw = np.full((len(pairs), 2 * NLAG + 1), np.nan, np.float32), np.zeros(len(pairs), np.int16)
    else:
        spec = preprocess(st, stas, day + pd.Timedelta(hours=hour_ut), hours)
        C, nw = day_ccfs(spec, pairs)
    tmp = path.with_name(path.stem + ".tmp.npz")
    np.savez_compressed(tmp, C=C.astype(np.float16), nwin=nw)
    tmp.replace(path)
    return int((nw > 0).sum())


def run_ccfs(out: Path, stas: list[str], start: str, end: str, step_days: int = 2,
             hour_ut: int = 20, hours: float = 4.0, workers: int = 12) -> Path:
    """Örnek günlerin çift korelasyonlarını ``out/days/`` altına yaz (kaldığı yerden sürer)."""
    out = Path(out)
    (out / "days").mkdir(parents=True, exist_ok=True)
    pairs = pairs_of(stas)
    sta = station_table(stas)
    sta.to_csv(out / "stations.csv", index=False)
    pd.DataFrame(pairs, columns=["a", "b"]).to_csv(out / "pairs.csv", index=False)
    days = [d for d in pd.date_range(start, end, freq=f"{step_days}D")
            if not (out / "days" / f"{d:%Y-%m-%d}.npz").exists()]
    log.info("%d örnek gün işlenecek", len(days))
    done = 0
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(process_day, stas, pairs, d, hour_ut, hours, out / "days"): d for d in days}
        for f in as_completed(futs):
            done += 1
            try:
                k = f.result()
            except Exception as e:
                log.warning("%s: %s", futs[f], e)
                continue
            if done % 25 == 0 or done == len(days):
                log.info("%d/%d gün (%s: %d çift)", done, len(days), futs[f].date(), k)
    return out


def load_ccfs(out: Path) -> tuple[pd.DatetimeIndex, np.ndarray, pd.DataFrame, pd.DataFrame]:
    out = Path(out)
    files = sorted((out / "days").glob("????-??-??.npz"))
    t = pd.to_datetime([f.stem for f in files])
    C = np.stack([np.load(f)["C"].astype(np.float32) for f in files])  # (gün, çift, gecikme)
    return t, C, pd.read_csv(out / "pairs.csv"), pd.read_csv(out / "stations.csv")


# ------------------------------------------------------------------ dv/v
def coda_window(d_km: float) -> tuple[float, float]:
    """Yüzey dalgasının (~3,5 km/s) varışından başlayan 40 s'lik pencere (erken kod).

    Seçim, yalnız 2009-2011 verisinde (değerlendirme depremlerine bakmadan) çiftleri
    rastgele iki yarıya bölüp bölgesel serinin gürültüsünü ölçerek yapıldı:
    geç kod (d/2,5+10 s'den 60 s) %0,083, erken pencere %0,053; 60 günlük yığınla %0,037.
    Erken pencere doğrudan dalgayı da içerdiği için gürültü kaynağının mevsimsel
    değişimine daha duyarlıdır; bu, mevsim sınamasında görünür.
    """
    t1 = max(8.0, d_km / 3.5)
    return t1, min(t1 + 40.0, MAXLAG - 20.0)


def stretching(ref, cur, t, t1, t2, eps=0.015, n=151):
    """cur(t) ≈ ref(t·(1+ε)) → ε = dv/v (pozitif = hız artışı). Kaba ızgara + parabol."""
    w = (t >= t1) & (t <= t2)
    grid = np.linspace(-eps, eps, n)
    cur_w = cur[w] - cur[w].mean()
    cc = np.empty(n)
    for i, e in enumerate(grid):
        r = np.interp(t[w] * (1 + e), t, ref)
        r = r - r.mean()
        cc[i] = (r @ cur_w) / (np.linalg.norm(r) * np.linalg.norm(cur_w) + 1e-12)
    i = int(np.argmax(cc))
    best = grid[i]
    if 0 < i < n - 1:
        y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
        den = y0 - 2 * y1 + y2
        if den < 0:
            best = grid[i] + 0.5 * (y0 - y2) / den * (grid[1] - grid[0])
    return float(best), float(cc[i])


def weaver_sigma(cc: float, t1: float, t2: float) -> float:
    """Weaver vd. (2011) esneme belirsizliği."""
    cc = min(max(cc, 1e-3), 0.999)
    T = 1.0 / (FMAX - FMIN)
    wc = 2 * np.pi * np.sqrt(FMIN * FMAX)
    return float(np.sqrt(1 - cc ** 2) / (2 * cc) * np.sqrt(6 * np.sqrt(np.pi / 2) * T / (wc ** 2 * (t2 ** 3 - t1 ** 3))))


def symmetric(c: np.ndarray) -> np.ndarray:
    return 0.5 * (c[NLAG:] + c[NLAG::-1])


def pair_dvv(t: pd.DatetimeIndex, C: np.ndarray, d_km: float, stack_days: int = 60,
             ref_years: float = 2.0, min_cc: float = 0.5) -> pd.DataFrame:
    """C: (gün, gecikme). Geriye dönük kayan yığın; referans ilk ``ref_years`` yıl."""
    good = np.isfinite(C).all(1)
    if good.sum() < 50:
        return pd.DataFrame()
    tg, Cg = t[good], C[good]
    ref_end = tg[0] + pd.DateOffset(years=int(ref_years))
    ref = symmetric(Cg[tg < ref_end].mean(0))
    tl = np.arange(NLAG + 1) / FS
    t1, t2 = coda_window(d_km)
    rows = []
    cs = np.vstack([np.zeros(Cg.shape[1]), np.cumsum(Cg, 0)])
    tv = tg.values
    lo = np.searchsorted(tv, (tg - pd.Timedelta(days=stack_days)).values, side="right")
    for i in range(len(tg)):
        k = i + 1 - lo[i]
        if k < max(3, stack_days // 4):  # yığında yeterli örnek gün yok
            continue
        cur = symmetric((cs[i + 1] - cs[lo[i]]) / k)
        e, cc = stretching(ref, cur, tl, t1, t2)
        rows.append((tg[i], e, cc, weaver_sigma(cc, t1, t2), k))
    df = pd.DataFrame(rows, columns=["time", "dvv", "cc", "sigma", "n_days"])
    df.loc[df["cc"] < min_cc, ["dvv", "sigma"]] = np.nan
    return df


def all_pairs_dvv(out: Path, stack_days: int = 60, min_snr: float = 3.0) -> pd.DataFrame:
    t, C, pairs, sta = load_ccfs(out)
    pos = sta.set_index("sta")
    frames = []
    for j, p in enumerate(pairs.itertuples()):
        a, b = pos.loc[p.a], pos.loc[p.b]
        d = dist_km(a.lat, a.lon, b.lat, b.lon)
        Cj = C[:, j, :]
        # çift kalitesi: tüm yığında kod penceresi SNR'si
        good = np.isfinite(Cj).all(1)
        if good.sum() < 50:
            continue
        full = symmetric(Cj[good].mean(0))
        t1, t2 = coda_window(d)
        tl = np.arange(NLAG + 1) / FS
        noise = np.std(full[tl > MAXLAG - 15]) + 1e-12
        snr = np.abs(full[(tl >= t1) & (tl <= t2)]).max() / noise
        df = pair_dvv(t, Cj, d, stack_days)
        if df.empty:
            continue
        df["a"], df["b"], df["dist_km"], df["snr"] = p.a, p.b, d, snr
        df["lat"], df["lon"] = (a.lat + b.lat) / 2, (a.lon + b.lon) / 2
        frames.append(df)
    res = pd.concat(frames, ignore_index=True)
    res = res[res["snr"] >= min_snr]
    return res


def regional_series(pairs: pd.DataFrame) -> pd.DataFrame:
    """Çiftlerin ters varyans ağırlıklı ortalaması ve belirsizliği."""
    p = pairs.dropna(subset=["dvv", "sigma"])
    p = p[p["sigma"] > 0]
    w = 1 / p["sigma"] ** 2
    g = p.assign(w=w, wx=w * p["dvv"]).groupby("time")
    out = pd.DataFrame({"dvv": g["wx"].sum() / g["w"].sum(), "sigma": 1 / np.sqrt(g["w"].sum()),
                        "n_pairs": g.size()})
    return out.reset_index()


def cell_features(pairs: pd.DataFrame, cells: pd.DataFrame, t0: pd.Timestamp,
                  radius_km: float = 40.0) -> pd.DataFrame:
    """Son 6 ay dv/v ortalaması − önceki 2 yıl (yüzde), çift orta noktalarından; yalnız t0 öncesi."""
    idx = cells["cell_id"]
    if pairs is None or pairs.empty:
        return pd.DataFrame(index=idx)
    rows = []
    for (a, b), g in pairs[pairs["time"] < t0].groupby(["a", "b"]):
        recent = g[g["time"] >= t0 - pd.DateOffset(months=6)]["dvv"]
        before = g[(g["time"] < t0 - pd.DateOffset(months=6)) & (g["time"] >= t0 - pd.DateOffset(months=30))]["dvv"]
        if recent.notna().sum() >= 20 and before.notna().sum() >= 80:
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
    return pd.DataFrame({"dvv_pair_change_pct": out}, index=idx)
