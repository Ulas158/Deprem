"""Kısmi tarihlerden (yalnız yıl, yıl+ay, ...) zaman aralığı ve 1σ hesabı.

Takvim notu: numpy proleptik Gregoryen takvim kullanır. Kaynak kataloglar
tarihleri kaynağın verdiği biçimde (1582 öncesi çoğunlukla Jülyen) aktarır;
burada dönüşüm yapılmaz. Bu fark (en çok ~10 gün) ``time_err_s`` yanında
küçüktür ama gün düzeyinde eşleştirmelerde akılda tutulmalıdır.

Yıl notu: NOAA MÖ yılları negatif verir ve 0 yılı yoktur (MÖ 1 = -1).
``astronomical_year`` bunu astronomik yıla (MÖ 1 = 0) çevirir.
"""

from __future__ import annotations

import math

import numpy as np

SECONDS_PER_DAY = 86400.0


def astronomical_year(year: int, historical_bc: bool) -> int:
    if historical_bc and year < 0:
        return year + 1
    return year


def _ystr(y: int) -> str:
    return f"{y:04d}" if y >= 0 else f"-{abs(y):04d}"


def _valid(v) -> bool:
    return v is not None and not (isinstance(v, float) and math.isnan(v))


def interval(year, month=None, day=None, hour=None, minute=None, second=None):
    """Verilen kesinliğe göre [başlangıç, bitiş) aralığını ve kesinlik etiketini döndürür."""
    y = int(year)
    if not _valid(month) or int(month) == 0:
        start = np.datetime64(_ystr(y), "Y").astype("datetime64[s]")
        end = (np.datetime64(_ystr(y), "Y") + np.timedelta64(1, "Y")).astype("datetime64[s]")
        return start, end, "year"
    m = int(month)
    if not _valid(day) or int(day) == 0:
        mm = np.datetime64(f"{_ystr(y)}-{m:02d}", "M")
        return mm.astype("datetime64[s]"), (mm + np.timedelta64(1, "M")).astype("datetime64[s]"), "month"
    d = np.datetime64(f"{_ystr(y)}-{m:02d}-{int(day):02d}", "D")
    if not _valid(hour):
        return d.astype("datetime64[s]"), (d + np.timedelta64(1, "D")).astype("datetime64[s]"), "day"
    h = d.astype("datetime64[h]") + np.timedelta64(int(hour), "h")
    if not _valid(minute):
        return h.astype("datetime64[s]"), (h + np.timedelta64(1, "h")).astype("datetime64[s]"), "hour"
    mi = h.astype("datetime64[m]") + np.timedelta64(int(minute), "m")
    if not _valid(second):
        return mi.astype("datetime64[s]"), (mi + np.timedelta64(1, "m")).astype("datetime64[s]"), "minute"
    s = mi.astype("datetime64[s]") + np.timedelta64(int(float(second)), "s")
    return s, s + np.timedelta64(1, "s"), "second"


def iso(t: np.datetime64) -> str:
    """ISO 8601 (UTC, 'Z'); MÖ yıllar dört haneli negatif yazılır: -0363-01-01T00:00:00Z."""
    s = str(t.astype("datetime64[s]"))
    if s.startswith("-"):
        y, rest = s[1:].split("-", 1)
        s = f"-{int(y):04d}-{rest}"
    return s + "Z"


def to_epoch_s(values) -> np.ndarray:
    """iso() çıktılarını (pandas'a ihtiyaç duymadan) epoch saniyesine çevirir."""
    return np.array([np.datetime64(v.rstrip("Z"), "s") for v in values]).astype("int64")


def summarize(year, month=None, day=None, hour=None, minute=None, second=None) -> dict:
    """Orta nokta zamanı, aralık sınırları, kesinlik ve 1σ (tekdüze dağılım: genişlik/√12)."""
    start, end, precision = interval(year, month, day, hour, minute, second)
    width_s = float((end - start) / np.timedelta64(1, "s"))
    mid = start + np.timedelta64(int(width_s // 2), "s")
    sigma_s = width_s / math.sqrt(12.0)
    # Ondalık yıl (modelleme için sürekli zaman ekseni).
    y0 = np.datetime64("1970-01-01T00:00:00", "s")
    dec_year = 1970.0 + float((mid - y0) / np.timedelta64(1, "s")) / (365.2425 * SECONDS_PER_DAY)
    return {
        "time": iso(mid),
        "time_min": iso(start),
        "time_max": iso(end),
        "time_precision": precision,
        "time_err_s": round(sigma_s, 1),
        "decimal_year": round(dec_year, 5),
    }
