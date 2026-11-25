"""GPS hızlarından fay kilitlenme oranı ve yavaş kayma (krip) hızı.

Her fay için fayın iki yanındaki GPS istasyonlarının fayla paralel hızı
(doğrultu atımlı faylar) ya da faya dik hızı (eğim atımlı faylar) şu elastik
modelle ters çözülür (Savage & Burford 1973; yüzeyde krip eklenmiş hali):

    v(x) = c0 + c1·x + (s − C)/π · atan(x / D) + (C/2) · sign(x)

* x: faya dik işaretli uzaklık (km), D: kilitli derinlik (fay veritabanından)
* s: fayın uzun dönem kayma hızı, C: yüzeye kadar sürünen (krip) kısım
* c0, c1: blok ötelemesi ve bölgesel gerinim (fayla ilgisiz eğilim)

Bilinmeyenler doğrusal olduğu için Gauss öncüllü en küçük kareler kullanılır:
s ~ N(fay veritabanı hızı, aralıktan sd), C ~ N(0, s/2). Kilitlenme oranı
L = (s − C)/s; belirsizlik sonsal dağılımdan örneklenerek hesaplanır.
L = 1 tam kilitli (enerji biriktiriyor), L = 0 tamamen sürünüyor (enerjiyi
yavaş yavaş boşaltıyor).

Sınır: istasyonlar seyrek (çoğu yerde 20-50 km aralık), bu yüzden çözüm çoğu
fayda öncüle yakın kalır; yalnız iki yanında yeterli istasyon olan faylar
için hesaplanır, diğerleri boş bırakılır. Yalnız ``year`` öncesi GPS verisi
kullanılır.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..geo import project_to_polyline

SIGMA_V = 1.5      # mm/yıl, istasyon hızı belirsizliği (renkli gürültü ve model hatası dahil)
MAX_DIST = 60.0    # km
ALONG_PAD = 25.0   # km, segment uçlarının ötesindeki istasyonlar da (fay devam ediyor)
MIN_SLIP = 2.0     # mm/yıl altındaki faylarda GPS sinyali gürültüden ayrılamaz


def station_velocities(sty: pd.DataFrame, year: int) -> pd.DataFrame:
    """Her istasyonun ``year`` öncesi yıllık hız tahminlerinin ortancası."""
    d = sty[sty["year"] <= year]  # satır yılı Y: [Y−5, Y) verisi
    if d.empty:
        return pd.DataFrame(columns=["site", "lat", "lon", "ve", "vn", "n"])
    g = d.groupby("site")
    return pd.DataFrame({"lat": g["lat"].first(), "lon": g["lon"].first(),
                         "ve": g["ve"].median(), "vn": g["vn"].median(),
                         "n": g.size()}).reset_index()


def _signed(lat, lon, coords):
    """Faya dik işaretli uzaklık (km; doğrultunun sağı +) ve yerel doğrultu (derece)."""
    d, along, L = project_to_polyline(lat, lon, coords)
    c = np.asarray(coords, float)
    lat0 = np.radians(c[:, 1].mean())
    kx, ky = 111.195 * np.cos(lat0), 111.195
    x, y = c[:, 0] * kx, c[:, 1] * ky
    cum = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    k = np.clip(np.searchsorted(cum, along) - 1, 0, len(c) - 2)
    dx, dy = x[k + 1] - x[k], y[k + 1] - y[k]
    px, py = np.asarray(lon) * kx - x[k], np.asarray(lat) * ky - y[k]
    side = np.sign(dx * py - dy * px)  # + sol
    strike = (np.degrees(np.arctan2(dx, dy)) + 360) % 360
    return -side * d, along, L, strike


def fit_fault(row: pd.Series, vel: pd.DataFrame, n_samp: int = 400, rng=None) -> dict | None:
    sr = float(row["slip_rate_mm_yr"])
    if not np.isfinite(sr) or sr < MIN_SLIP or vel.empty:
        return None
    rng = rng or np.random.default_rng(0)
    x, along, L, strike = _signed(vel["lat"].to_numpy(), vel["lon"].to_numpy(), row["coords"])
    inside = (np.abs(x) <= MAX_DIST) & (along > -ALONG_PAD) & (along < L + ALONG_PAD)
    if inside.sum() < 5 or (x[inside] > 0).sum() < 2 or (x[inside] < 0).sum() < 2:
        return None
    x, strike = x[inside], strike[inside]
    ve, vn = vel["ve"].to_numpy()[inside], vel["vn"].to_numpy()[inside]
    phi = np.radians(strike)
    rake = float(row["rake"]) if pd.notna(row["rake"]) else 180.0
    strike_slip = abs(np.sin(np.radians(rake))) < 0.7
    if strike_slip:
        v = ve * np.sin(phi) + vn * np.cos(phi)   # doğrultu yönündeki bileşen
        sign = -1.0 if np.cos(np.radians(rake)) < 0 else 1.0  # sağ yanal: sağdaki blok doğrultunun tersine
    else:
        v = ve * np.cos(phi) - vn * np.sin(phi)   # doğrultuya dik (sağa) bileşen
        sign = 1.0 if np.sin(np.radians(rake)) < 0 else -1.0  # normal: açılma
    D = max(float(row["lower_depth_km"]) - float(row["upper_depth_km"]), 3.0)
    # bilinmeyenler: c0, c1, a = (s − C), C   (işaret ``sign`` ile fay yönüne çevrilir)
    G = np.column_stack([np.ones_like(x), x / 100.0, sign * np.arctan(x / D) / np.pi,
                         sign * 0.5 * np.sign(x)])
    lo, hi = row.get("slip_rate_min", np.nan), row.get("slip_rate_max", np.nan)
    s_sd = max((hi - lo) / 4 if np.isfinite(lo) and np.isfinite(hi) and hi > lo else 0.3 * sr, 0.15 * sr, 0.5)
    # öncül: a + C = s ~ N(sr, s_sd); C ~ N(0, sr/2) (krip çoğu fayda küçük)
    # (a, C) için ortak Gauss öncül: a = s − C
    m0 = np.array([0.0, 0.0, sr, 0.0])
    P = np.zeros((4, 4))
    P[0, 0] = 1 / 50.0 ** 2
    P[1, 1] = 1 / 20.0 ** 2
    # s = a + C bileşeni
    u = np.array([0, 0, 1.0, 1.0])
    P += np.outer(u, u) / s_sd ** 2
    w = np.array([0, 0, 0, 1.0])
    P += np.outer(w, w) / (sr / 2) ** 2
    A = G.T @ G / SIGMA_V ** 2 + P
    b = G.T @ v / SIGMA_V ** 2 + P @ m0
    cov = np.linalg.inv(A)
    mu = cov @ b
    smp = rng.multivariate_normal(mu, cov, n_samp)
    a, C = smp[:, 2], smp[:, 3]
    C = np.clip(C, 0, None)
    a = np.clip(a, 0, None)
    s = a + C
    Lk = np.where(s > 0, a / np.maximum(s, 1e-6), 1.0)
    # öncüle göre bilgi: verinin C belirsizliğini ne kadar daralttığı
    shrink = 1 - np.sqrt(cov[3, 3]) / (sr / 2)
    return {"fault_id": row["fault_id"], "n_sta": int(inside.sum()),
            "locking": float(np.median(Lk)), "locking_sd": float(np.std(Lk)),
            "creep_mm_yr": float(np.median(C)), "creep_sd": float(np.std(C)),
            "slip_geod_mm_yr": float(np.median(s)), "data_info": float(shrink),
            "rms_mm_yr": float(np.sqrt(np.mean((v - G @ mu) ** 2)))}


def locking_table(faults: pd.DataFrame, sty: pd.DataFrame, year: int) -> pd.DataFrame:
    vel = station_velocities(sty, year)
    rng = np.random.default_rng(year)
    rows = [r for r in (fit_fault(f, vel, rng=rng) for _, f in faults.iterrows()) if r]
    return pd.DataFrame(rows)


# Literatürdeki krip ölçümleri (karşılaştırma için; modele girmez)
KNOWN_CREEP = pd.DataFrame([
    {"yer": "KAF İsmetpaşa", "lat": 40.87, "lon": 32.60, "krip_mm_yil": 8.0,
     "kaynak": "Çakır vd. 2005; Kaneko vd. 2013 (InSAR, kripmetre ~6-9 mm/yıl)"},
    {"yer": "KAF Destek", "lat": 40.95, "lon": 36.45, "krip_mm_yil": 6.0,
     "kaynak": "Kaneko vd. 2013 (InSAR, ~5-8 mm/yıl)"},
    {"yer": "Marmara İzmit (1999 sonrası)", "lat": 40.72, "lon": 30.0, "krip_mm_yil": 5.0,
     "kaynak": "Aslan vd. 2019 (deprem sonrası sığ krip, azalan)"},
    {"yer": "DAF Hazar-Palu", "lat": 38.5, "lon": 39.5, "krip_mm_yil": 4.0,
     "kaynak": "Cavalié & Jónsson 2014 (sığ kısmi krip)"},
])
