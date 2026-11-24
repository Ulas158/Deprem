"""İleriye dönük test edilebilir ızgara tabanlı olasılık modelleri.

Bütün modeller (bugünkü basit referanslar, ileride ETAS, fay-tabanlı ve
Coulomb-tabanlı modeller) aynı arayüzü uygular: bir test dönemi için
her (enlem, boylam, büyüklük) hücresinde **beklenen olay sayısı** üretir.
Böylece CSEP tarzı testlerle (``deprem.forecast.tests``) aynı veri
üzerinde karşılaştırılabilirler.

Kural: Model yalnızca ``fit_end`` öncesi veriyle eğitilir; test dönemi
verisi modele asla gösterilmez.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..region import Region
from ..analysis.seismicity import grid_edges


@dataclass
class ForecastGrid:
    region: Region
    cell_deg: float
    mag_edges: np.ndarray  # büyüklük kutu kenarları, son kutu açık uçlu sayılır

    def __post_init__(self):
        self.lat_edges, self.lon_edges = grid_edges(self.region, self.cell_deg)

    @property
    def shape(self):
        return (self.lat_edges.size - 1, self.lon_edges.size - 1, self.mag_edges.size - 1)

    def count(self, df: pd.DataFrame) -> np.ndarray:
        """Gözlenen olayları aynı hücrelere say."""
        mags = df["mag"].clip(upper=self.mag_edges[-1] - 1e-6)
        h, _ = np.histogramdd(
            np.column_stack([df["latitude"], df["longitude"], mags]),
            bins=[self.lat_edges, self.lon_edges, self.mag_edges])
        return h

    def cell_areas_km2(self) -> np.ndarray:
        lat_c = (self.lat_edges[:-1] + self.lat_edges[1:]) / 2
        dy = self.cell_deg * 111.195
        dx = self.cell_deg * 111.195 * np.cos(np.radians(lat_c))
        return np.repeat((dx * dy)[:, None], self.lon_edges.size - 1, axis=1)


class Forecast(ABC):
    name: str = "base"

    @abstractmethod
    def fit(self, catalog: pd.DataFrame, fit_start: pd.Timestamp,
            fit_end: pd.Timestamp) -> "Forecast":
        ...

    @abstractmethod
    def expected(self, grid: ForecastGrid, start: pd.Timestamp,
                 end: pd.Timestamp) -> np.ndarray:
        """Şekli grid.shape olan beklenen olay sayıları."""


def _gr_mag_fractions(mag_edges, b, mmin):
    """GR'ye göre büyüklük kutularına düşen olay oranları (son kutu açık)."""
    lo = mag_edges[:-1]
    hi = np.append(mag_edges[1:-1], np.inf)
    lo = np.maximum(lo, mmin)
    f = 10 ** (-b * (lo - mmin)) - np.where(np.isinf(hi), 0, 10 ** (-b * (hi - mmin)))
    return np.clip(f, 0, None)


class UniformPoisson(Forecast):
    """Referans 0: bölge içinde mekânda düzgün, zamanda sabit oranlı Poisson."""

    name = "uniform_poisson"

    def __init__(self, mmin: float, b: float = 1.0):
        self.mmin, self.b = mmin, b

    def fit(self, catalog, fit_start, fit_end):
        d = catalog[(catalog["time"] >= fit_start) & (catalog["time"] < fit_end)
                    & (catalog["mag"] >= self.mmin)]
        self.rate_per_day = len(d) / max((fit_end - fit_start).days, 1)
        return self

    def expected(self, grid, start, end):
        days = (end - start).days
        area = grid.cell_areas_km2()
        spatial = area / area.sum()
        mag = _gr_mag_fractions(grid.mag_edges, self.b, self.mmin)
        return self.rate_per_day * days * spatial[:, :, None] * mag[None, None, :]


class SmoothedSeismicity(Forecast):
    """Referans 1: zamanda sabit, mekânda Gauss çekirdekle düzleştirilmiş Poisson.

    Frankel (1995) / Helmstetter vd. (2007) türü düzleştirilmiş sismisite;
    ETAS ve diğer modellerin geçmesi gereken asgari çıta budur.
    """

    name = "smoothed_seismicity"

    def __init__(self, mmin: float, b: float = 1.0, bandwidth_km: float = 25.0,
                 floor: float = 1e-3):
        self.mmin, self.b, self.bw, self.floor = mmin, b, bandwidth_km, floor

    def fit(self, catalog, fit_start, fit_end):
        d = catalog[(catalog["time"] >= fit_start) & (catalog["time"] < fit_end)
                    & (catalog["mag"] >= self.mmin)]
        self._train = d[["latitude", "longitude"]].to_numpy()
        self.rate_per_day = len(d) / max((fit_end - fit_start).days, 1)
        return self

    def expected(self, grid, start, end):
        lat_c = (grid.lat_edges[:-1] + grid.lat_edges[1:]) / 2
        lon_c = (grid.lon_edges[:-1] + grid.lon_edges[1:]) / 2
        LA, LO = np.meshgrid(lat_c, lon_c, indexing="ij")
        dens = np.zeros_like(LA)
        ky = 111.195
        kx = 111.195 * np.cos(np.radians(LA))
        for la, lo in self._train:
            d2 = ((LA - la) * ky) ** 2 + ((LO - lo) * kx) ** 2
            dens += np.exp(-d2 / (2 * self.bw ** 2))
        area = grid.cell_areas_km2()
        dens = dens / dens.sum()
        uniform = area / area.sum()
        spatial = (1 - self.floor) * dens + self.floor * uniform
        days = (end - start).days
        mag = _gr_mag_fractions(grid.mag_edges, self.b, self.mmin)
        return self.rate_per_day * days * spatial[:, :, None] * mag[None, None, :]
