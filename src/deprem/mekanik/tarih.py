"""Marmara'da Kuzey Anadolu Fayı kuzey kolunun tarihsel kırılmaları.

Kırık uçları boylam olarak verilir ve iz boyunca x_km'ye çevrilir. Tarihsel
depremlerin kırık uzunlukları belirsizdir; her uç için ``belirsiz_km``
verilir ve topluluk üyelerinde bu aralıkta rastgele kaydırılır.

Kaynaklar: Ambraseys & Jackson 2000 (GJI), Ambraseys 2002 (BSSA),
Parsons 2004 (JGR), Pondard vd. 2007 (GJI), Aksoy vd. 2010 (1912),
Barka vd. 2002 (1999 İzmit), Akyüz vd. 2002 (1999 Düzce),
EPICA/AHEAD tarihsel katalogu (büyüklükler). 2025 Silivri: AFAD/USGS.
Not: 1957 Abant ve 1967 Mudurnu güney koldadır; model ağında değildir.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# tarih, ad, batı ucu boylam, doğu ucu boylam, Mw, Mw belirsizliği, uç belirsizliği km, güven
KIRILMALAR = [
    ("1509-09-10", "1509 İstanbul (Kıyamet-i Suğra)", 28.00, 29.10, 7.2, 0.3, 25, "düşük"),
    ("1659-02-17", "1659 Tekirdağ", 26.95, 27.55, 7.0, 0.3, 25, "düşük"),
    ("1719-05-25", "1719 İzmit", 29.35, 30.60, 7.4, 0.3, 25, "orta"),
    ("1754-09-02", "1754 İzmit Körfezi", 29.10, 29.80, 6.9, 0.3, 20, "düşük"),
    ("1766-05-22", "1766 Mayıs Doğu Marmara", 28.70, 29.40, 7.1, 0.3, 20, "orta"),
    ("1766-08-05", "1766 Ağustos Ganos", 26.60, 27.80, 7.4, 0.3, 20, "orta"),
    ("1894-07-10", "1894 Adalar-İzmit Körfezi", 29.00, 29.75, 7.0, 0.3, 20, "orta"),
    ("1912-08-09", "1912 Mürefte-Şarköy", 26.40, 27.85, 7.3, 0.15, 10, "iyi"),
    ("1943-06-20", "1943 Hendek-Adapazarı", 30.45, 30.85, 6.5, 0.2, 10, "orta"),
    ("1963-09-18", "1963 Çınarcık", 29.05, 29.25, 6.3, 0.2, 8, "orta"),
    ("1999-08-17", "1999 İzmit (Kocaeli)", 29.25, 30.95, 7.6, 0.05, 3, "iyi"),
    ("1999-11-12", "1999 Düzce", 30.96, 31.40, 7.2, 0.05, 3, "iyi"),
    ("2022-11-23", "2022 Düzce-Gölyaka", 30.92, 31.06, 6.1, 0.05, 5, "iyi"),
    ("2025-04-23", "2025 Silivri açığı", 28.05, 28.35, 6.2, 0.05, 5, "iyi"),
]


def ondalik_yil(t: str) -> float:
    ts = pd.Timestamp(t)
    return ts.year + (ts.dayofyear - 0.5) / (366 if ts.is_leap_year else 365)


def boylam_to_x(lon, cells: pd.DataFrame):
    """Kabuk hücrelerinin boylam-x eşlemesinden doğrusal enterpolasyon."""
    k = cells[cells["z_ust"] == 0].sort_values("x_km")
    return np.interp(lon, k["lon"].to_numpy(), k["x_km"].to_numpy())


def tablo(cells: pd.DataFrame) -> pd.DataFrame:
    df = pd.DataFrame(KIRILMALAR, columns=["tarih", "ad", "lon0", "lon1", "mw", "mw_sd", "uc_km", "guven"])
    df["t"] = df["tarih"].map(ondalik_yil)
    df["x0"] = boylam_to_x(df["lon0"], cells)
    df["x1"] = boylam_to_x(df["lon1"], cells)
    df["uzunluk_km"] = df["x1"] - df["x0"]
    return df
