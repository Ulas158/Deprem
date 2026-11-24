"""Çalışma bölgeleri.

Ana çalışma alanı Türkiye ve yakın çevresidir (TURKEY). Alt bölgeler,
sonuçları tektonik olarak farklı kuşaklar arasında karşılaştırmak için
tanımlanmıştır; sınırlar kasıtlı olarak kaba kutulardır ve ileride
poligon veya faya-uzaklık tabanlı seçimle değiştirilebilir.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Region:
    name: str
    min_lat: float
    max_lat: float
    min_lon: float
    max_lon: float
    max_depth_km: float = 70.0

    def contains(self, df: pd.DataFrame) -> pd.Series:
        return (
            df["latitude"].between(self.min_lat, self.max_lat)
            & df["longitude"].between(self.min_lon, self.max_lon)
        )

    def area_km2(self) -> float:
        lat0 = np.deg2rad((self.min_lat + self.max_lat) / 2)
        dy = (self.max_lat - self.min_lat) * 111.195
        dx = (self.max_lon - self.min_lon) * 111.195 * np.cos(lat0)
        return float(dx * dy)


# Türkiye ve komşu deniz alanları (Ege, Doğu Akdeniz kıyısı, Karadeniz kıyısı)
TURKEY = Region("Türkiye", min_lat=35.5, max_lat=42.5, min_lon=25.5, max_lon=45.0)

# Karşılaştırma için kaba alt bölgeler (çakışabilirler)
SUBREGIONS = {
    "marmara": Region("Marmara (KAF batı)", 40.0, 41.3, 26.0, 30.5),
    "kaf_orta": Region("KAF orta (Bolu-Erzincan)", 40.0, 41.2, 30.5, 39.5),
    "dogu_anadolu": Region("Doğu Anadolu Fay Zonu", 36.0, 39.0, 35.8, 41.0),
    "bati_anadolu": Region("Batı Anadolu / Ege graben", 36.5, 40.0, 26.0, 30.5),
    "dogu": Region("Doğu Anadolu platosu (Van-Kars)", 38.0, 41.5, 41.0, 45.0),
}
