"""CRUST1.0 kabuk modeli (Laske vd. 2013): kayaç sertliği (rijitlik) ve kabuk kalınlığı.

1°×1° hücrelerde 9 katman (su, buz, 3 tortul, üst/orta/alt kabuk, manto) için
üst sınır yüksekliği (km, deniz seviyesine göre), Vp, Vs (km/s) ve yoğunluk (g/cm³).
Dosyalar ``data/ext/crust1/`` altına konur (crust1.bnds, .vs, .rho; ``download``).

Rijitlik (makaslama modülü) μ = ρ·Vs². Bir derinlik aralığında kalınlık ağırlıklı
ortalama alınır. Yüksek μ = daha sert kayaç: aynı gerilme düşümü için daha az kayma,
yani aynı gerilmeyi taşımak için daha az kayma açığı yeter.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import requests

FILES = ("crust1.bnds", "crust1.vs", "crust1.rho", "crust1.vp")
URL = "https://raw.githubusercontent.com/junlysky/crust1.0/master/"  # özgün dağıtımın aynası


def download(dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    for f in FILES:
        p = dest / f
        if not p.exists():
            r = requests.get(URL + f, timeout=300)
            r.raise_for_status()
            p.write_bytes(r.content)
    return dest


class Crust1:
    def __init__(self, folder: Path):
        folder = Path(folder)
        self.bnds = np.loadtxt(folder / "crust1.bnds").reshape(180, 360, 9)
        self.vs = np.loadtxt(folder / "crust1.vs").reshape(180, 360, 9)
        self.rho = np.loadtxt(folder / "crust1.rho").reshape(180, 360, 9)

    @staticmethod
    def _idx(lat, lon):
        i = np.clip((90 - np.asarray(lat, float)).astype(int), 0, 179)
        j = np.clip((np.asarray(lon, float) + 180).astype(int), 0, 359)
        return i, j

    def moho_depth_km(self, lat, lon):
        i, j = self._idx(lat, lon)
        return -self.bnds[i, j, 8]

    def rigidity_pa(self, lat, lon, z_top_km, z_bot_km):
        """[z_top, z_bot] (km, deniz seviyesinden aşağı) aralığında ortalama μ (Pa)."""
        i, j = self._idx(lat, lon)
        top = -self.bnds[i, j, :]                       # katman üstleri, derinlik (km)
        bot = np.concatenate([top[..., 1:], np.full(top.shape[:-1] + (1,), 1e4)], axis=-1)
        mu = (self.rho[i, j, :] * 1e3) * (self.vs[i, j, :] * 1e3) ** 2
        zt = np.asarray(z_top_km, float)[..., None]
        zb = np.asarray(z_bot_km, float)[..., None]
        zb = np.maximum(zb, zt + 0.1)
        ov = np.clip(np.minimum(bot, zb) - np.maximum(top, zt), 0, None)
        ov = np.where(self.vs[i, j, :] > 0.5, ov, 0.0)  # su/buz katmanı hariç
        return (ov * mu).sum(-1) / np.maximum(ov.sum(-1), 1e-9)
