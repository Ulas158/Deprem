"""Küçük coğrafi yardımcılar (ek bağımlılık olmadan)."""

from __future__ import annotations

import numpy as np

EARTH_RADIUS_KM = 6371.0


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = (np.sin((lat2 - lat1) / 2) ** 2
         + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2)
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def polyline_length_km(coords) -> float:
    """[(lon, lat), ...] çizgisinin uzunluğu."""
    c = np.asarray(coords, dtype=float)
    if len(c) < 2:
        return 0.0
    return float(haversine_km(c[:-1, 1], c[:-1, 0], c[1:, 1], c[1:, 0]).sum())


def point_to_polyline_km(lat, lon, coords) -> np.ndarray:
    """Noktaların bir çoklu çizgiye en yakın uzaklığı (yerel düzlem yaklaşımı)."""
    lat = np.atleast_1d(np.asarray(lat, dtype=float))
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    c = np.asarray(coords, dtype=float)
    lat0 = np.radians(c[:, 1].mean())
    kx = 111.195 * np.cos(lat0)
    ky = 111.195
    px, py = lon * kx, lat * ky
    ax, ay = c[:-1, 0] * kx, c[:-1, 1] * ky
    bx, by = c[1:, 0] * kx, c[1:, 1] * ky
    if len(c) == 1:
        return np.hypot(px - c[0, 0] * kx, py - c[0, 1] * ky)
    dx, dy = bx - ax, by - ay
    seg2 = np.where(dx * dx + dy * dy == 0, 1e-12, dx * dx + dy * dy)
    tt = ((px[:, None] - ax) * dx + (py[:, None] - ay) * dy) / seg2
    tt = np.clip(tt, 0, 1)
    qx, qy = ax + tt * dx, ay + tt * dy
    return np.hypot(px[:, None] - qx, py[:, None] - qy).min(axis=1)


def project_to_polyline(lat, lon, coords):
    """Noktaları çoklu çizgiye izdüşür.

    Dönen: (faya_dik_uzaklık_km, çizgi_boyunca_konum_km, çizgi_uzunluğu_km).
    Yerel düzlem yaklaşımı (segment uzunlukları ≲ 500 km için yeterli).
    """
    lat = np.atleast_1d(np.asarray(lat, dtype=float))
    lon = np.atleast_1d(np.asarray(lon, dtype=float))
    c = np.asarray(coords, dtype=float)
    lat0 = np.radians(c[:, 1].mean())
    kx, ky = 111.195 * np.cos(lat0), 111.195
    px, py = lon * kx, lat * ky
    x, y = c[:, 0] * kx, c[:, 1] * ky
    if len(c) == 1:
        return np.hypot(px - x[0], py - y[0]), np.zeros_like(px), 0.0
    dx, dy = np.diff(x), np.diff(y)
    seg_len = np.hypot(dx, dy)
    cum = np.concatenate([[0], np.cumsum(seg_len)])
    seg2 = np.where(seg_len == 0, 1e-12, seg_len ** 2)
    tt = np.clip(((px[:, None] - x[:-1]) * dx + (py[:, None] - y[:-1]) * dy) / seg2, 0, 1)
    qx, qy = x[:-1] + tt * dx, y[:-1] + tt * dy
    d = np.hypot(px[:, None] - qx, py[:, None] - qy)
    k = d.argmin(axis=1)
    rows = np.arange(len(px))
    along = cum[k] + tt[rows, k] * seg_len[k]
    return d[rows, k], along, float(cum[-1])
