"""Ay ve Güneş için düşük doğruluklu konum hesabı ve gelgit fazları.

Ek bağımlılık olmadan, olay zamanlarında gelgit fazını hesaplamak için
yeterli doğrulukta (Ay ~0,3°, Güneş ~0,01°) analitik seriler kullanılır
(Montenbruck & Pfleger, *Astronomy on the Personal Computer*, MiniMoon;
Meeus, *Astronomical Algorithms*, bölüm 25 ve 47 ana terimleri).

Gelgit büyüklüğü olarak derece-2 gelgit potansiyelinin düşey ivmesi
kullanılır: g ∝ Σ (M/r³)·(3cos²ψ − 1). Bu, katı-yer gelgitinin fay
üzerindeki Coulomb gerilmesi değildir (o, fay yönelimi ve okyanus
yüklemesi gerektirir); fazın kaba bir vekilidir. Tam gelgit gerilmesi
ileride SPOTL/ETERNA gibi araçlarla eklenebilir.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ARC = 206264.806
J2000 = 2451545.0
MOON_EARTH_MASS = 0.0123000371
SUN_EARTH_MASS = 332946.0487
AU_KM = 149597870.7


def julian_day(times) -> np.ndarray:
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    secs = (t - pd.Timestamp("1970-01-01", tz="UTC")).total_seconds().to_numpy()
    return secs / 86400.0 + 2440587.5


def _frac(x):
    return x - np.floor(x)


def _ecl_to_eq(lam, beta, eps):
    x = np.cos(beta) * np.cos(lam)
    y = np.cos(beta) * np.sin(lam) * np.cos(eps) - np.sin(beta) * np.sin(eps)
    z = np.cos(beta) * np.sin(lam) * np.sin(eps) + np.sin(beta) * np.cos(eps)
    return np.arctan2(y, x), np.arcsin(np.clip(z, -1, 1))


def moon_position(jd):
    """(ekliptik boylam, enlem [rad], uzaklık [km])."""
    T = (np.asarray(jd) - J2000) / 36525.0
    L0 = _frac(0.606433 + 1336.855225 * T)
    l = 2 * np.pi * _frac(0.374897 + 1325.552410 * T)
    ls = 2 * np.pi * _frac(0.993133 + 99.997361 * T)
    D = 2 * np.pi * _frac(0.827361 + 1236.853086 * T)
    F = 2 * np.pi * _frac(0.259086 + 1342.227825 * T)
    dL = (22640 * np.sin(l) - 4586 * np.sin(l - 2 * D) + 2370 * np.sin(2 * D)
          + 769 * np.sin(2 * l) - 668 * np.sin(ls) - 412 * np.sin(2 * F)
          - 212 * np.sin(2 * l - 2 * D) - 206 * np.sin(l + ls - 2 * D)
          + 192 * np.sin(l + 2 * D) - 165 * np.sin(ls - 2 * D) - 125 * np.sin(D)
          - 110 * np.sin(l + ls) + 148 * np.sin(l - ls) - 55 * np.sin(2 * F - 2 * D))
    S = F + (dL + 412 * np.sin(2 * F) + 541 * np.sin(ls)) / ARC
    h = F - 2 * D
    N = (-526 * np.sin(h) + 44 * np.sin(l + h) - 31 * np.sin(-l + h) - 23 * np.sin(ls + h)
         + 11 * np.sin(-ls + h) - 25 * np.sin(-2 * l + F) + 21 * np.sin(-l + F))
    lam = 2 * np.pi * _frac(L0 + dL / 1296000.0)
    beta = (18520.0 * np.sin(S) + N) / ARC
    r = (385000.56 - 20905.36 * np.cos(l) - 3699.11 * np.cos(2 * D - l)
         - 2955.97 * np.cos(2 * D) - 569.93 * np.cos(2 * l) + 246.16 * np.cos(2 * l - 2 * D))
    return lam, beta, r


def sun_position(jd):
    d = np.asarray(jd) - J2000
    g = np.radians(357.529 + 0.98560028 * d)
    q = 280.459 + 0.98564736 * d
    lam = np.radians(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    r = (1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)) * AU_KM
    return lam, np.zeros_like(lam), r


def _obliquity(jd):
    return np.radians(23.439 - 0.00000036 * (np.asarray(jd) - J2000))


def _gmst_rad(jd):
    d = np.asarray(jd) - J2000
    return np.radians(np.mod(280.46061837 + 360.98564736629 * d, 360.0))


def tidal_gravity(jd, lat, lon) -> np.ndarray:
    """Ay + Güneş derece-2 düşey gelgit ivmesi (göreli birim, Ay ortalaması ≈ 1)."""
    jd = np.asarray(jd, dtype=float)
    phi = np.radians(np.asarray(lat, dtype=float))
    lam_obs = np.radians(np.asarray(lon, dtype=float))
    eps = _obliquity(jd)
    th = _gmst_rad(jd) + lam_obs
    out = 0.0
    for (lam, beta, r), mass in ((moon_position(jd), MOON_EARTH_MASS),
                                 (sun_position(jd), SUN_EARTH_MASS)):
        ra, dec = _ecl_to_eq(lam, beta, eps)
        H = th - ra
        cpsi = np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.cos(H)
        out = out + mass / (r / 384400.0) ** 3 * (3 * cpsi ** 2 - 1) / MOON_EARTH_MASS
    return out


def elongation(jd) -> np.ndarray:
    """Ay − Güneş ekliptik boylam farkı [rad] (0 yeni ay, π dolunay)."""
    return np.mod(moon_position(jd)[0] - sun_position(jd)[0], 2 * np.pi)


def semidiurnal_phase(times, lat, lon, half_window_h: float = 18.0,
                      step_min: float = 6.0, chunk: int = 4000) -> np.ndarray:
    """Tanaka (2002) tarzı gelgit fazı [derece, −180..180].

    Olayın etrafındaki gelgit eğrisinde en yakın tepe 0°, önceki çukur −180°,
    sonraki çukur +180° kabul edilir; aradaki zaman doğrusal eşlenir.
    """
    jd0 = julian_day(times)
    lat = np.broadcast_to(np.asarray(lat, dtype=float), jd0.shape)
    lon = np.broadcast_to(np.asarray(lon, dtype=float), jd0.shape)
    offs = np.arange(-half_window_h * 60, half_window_h * 60 + 1e-9, step_min) / 1440.0
    mid = len(offs) // 2
    phase = np.full(len(jd0), np.nan)
    for s in range(0, len(jd0), chunk):
        sl = slice(s, s + chunk)
        jd = jd0[sl, None] + offs[None, :]
        g = tidal_gravity(jd, lat[sl, None], lon[sl, None])
        dg = np.diff(g, axis=1)
        is_max = np.zeros_like(g, dtype=bool)
        is_min = np.zeros_like(g, dtype=bool)
        is_max[:, 1:-1] = (dg[:, :-1] > 0) & (dg[:, 1:] <= 0)
        is_min[:, 1:-1] = (dg[:, :-1] < 0) & (dg[:, 1:] >= 0)
        idx = np.arange(len(offs))
        for k in range(g.shape[0]):
            mx = idx[is_max[k]]
            mn = idx[is_min[k]]
            prev_min = mn[mn <= mid]
            next_min = mn[mn > mid]
            if len(prev_min) == 0 or len(next_min) == 0:
                continue
            a, b = prev_min[-1], next_min[0]
            peaks = mx[(mx > a) & (mx < b)]
            if len(peaks) == 0:
                continue
            p = peaks[np.argmax(g[k, peaks])]
            if mid <= p:
                phase[s + k] = -180.0 * (p - mid) / max(p - a, 1)
            else:
                phase[s + k] = 180.0 * (mid - p) / max(b - p, 1)
    return phase


def fortnightly_phase(times) -> np.ndarray:
    """Dolunay/yeni ay (kuvvetli gelgit) 0°, dördün (zayıf gelgit) ±180° [derece]."""
    e = elongation(julian_day(times))
    return np.degrees(np.angle(np.exp(2j * e)))


def annual_phase(times) -> np.ndarray:
    """Yıl içi faz [derece]; 1 Ocak 0°."""
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    doy = (t.dayofyear - 1 + t.hour / 24.0) / np.where(t.is_leap_year, 366.0, 365.0)
    return np.degrees(np.angle(np.exp(2j * np.pi * np.asarray(doy))))


def lunar_semidiurnal_phase(times, lon) -> np.ndarray:
    """Yalnızca Ay'ın yarı günlük (M2 benzeri) fazı [derece]: 2 × yerel Ay saat açısı.

    0° Ay'ın meridyen geçişi (veya karşı geçişi). Güneş bileşenini (S2, tam 12 sa)
    içermez; bu yüzden gece/gündüz algılama farkından etkilenmez.
    """
    jd = julian_day(times)
    lam, beta, _ = moon_position(jd)
    ra, _ = _ecl_to_eq(lam, beta, _obliquity(jd))
    H = _gmst_rad(jd) + np.radians(np.asarray(lon, dtype=float)) - ra
    return np.degrees(np.angle(np.exp(2j * H)))


def local_hour_phase(times, lon) -> np.ndarray:
    """Yerel güneş saati fazı [derece] (gece yarısı 0°). Algılama/gürültü kontrolü için."""
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    h = (t.hour + t.minute / 60.0 + np.asarray(lon, dtype=float) / 15.0) % 24
    return np.degrees(np.angle(np.exp(1j * np.radians(np.asarray(h) * 15.0))))
