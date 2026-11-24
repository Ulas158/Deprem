"""Belirsizliğe dayalı eğitim ağırlıkları.

İlke (proje kararı, 24 Kasım 2026):

* Ağırlık bir kaydın **ölçüm belirsizliğinden** türetilir (ters varyans
  tarzı). Küçük depremler sırf küçük oldukları için düşürülmez; modern
  ağlarda iyi ölçülürler ve sayıca bilginin çoğunu taşırlar.
* Eski/tarihsel kayıtlar belirsizlikleri büyük olduğu için kendiliğinden
  düşük ağırlık alır.
* Büyük depremler ağırlıkla şişirilmez (az oldukları için model onları
  ezberler); başarıları ayrı ölçütlerle izlenir
  (``deprem.features.evaluate.transfer_report``, eşik başına).

Ağırlık, her belirsizliğin bir **ölçek** değerine göre normalize edilmiş
karelerinin toplamından hesaplanır:

    w = 1 / (1 + Σ_k (σ_k / s_k)²)

``s_k`` "bu kadar belirsizlik ağırlığı yarıya indirir" ölçeğidir ve
parametriktir. Belirsizliği verilmeyen alanlar kaynak/dönem varsayılanlarıyla
doldurulur; hangi değerin varsayılan olduğu ``*_imputed`` sütunlarında tutulur.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .schema import UNCERTAINTY_COLUMNS


@dataclass
class WeightConfig:
    # Ağırlığı yarıya indiren belirsizlik ölçekleri
    scales: dict = field(default_factory=lambda: {
        # Ölçekler tarihsel katalog thread'iyle hizalı (24 Kasım 2026)
        "time_err_s": 30 * 86400.0,  # 30 gün: yıllık pencerelerde bundan küçük hata önemsiz
        "horiz_err_km": 20.0,
        "depth_err_km": 30.0,
        "mag_err": 0.3,
    })
    # Kaynak vermediğinde kullanılan varsayılan 1σ değerleri:
    # (başlangıç yılı, {alan: σ}) listesi; o yıldan sonraki olaylara uygulanır.
    defaults_by_era: list = field(default_factory=lambda: [
        (1900, {"time_err_s": 10.0, "horiz_err_km": 30.0, "depth_err_km": 20.0, "mag_err": 0.4}),
        (1964, {"time_err_s": 3.0, "horiz_err_km": 15.0, "depth_err_km": 15.0, "mag_err": 0.3}),
        (1998, {"time_err_s": 1.0, "horiz_err_km": 7.0, "depth_err_km": 10.0, "mag_err": 0.2}),
        (2008, {"time_err_s": 0.5, "horiz_err_km": 3.0, "depth_err_km": 5.0, "mag_err": 0.15}),
    ])
    # İsteğe bağlı kaynak başına çarpanlar (ör. tek istasyonlu çözümler için)
    source_factor: dict = field(default_factory=dict)
    # Ağırlıkların bu değerden küçüğü bu değere çekilir (sıfır ağırlık yok)
    floor: float = 1e-3


def fill_uncertainties(df: pd.DataFrame, cfg: WeightConfig | None = None) -> pd.DataFrame:
    cfg = cfg or WeightConfig()
    out = df.copy()
    years = out["time"].dt.year.to_numpy()
    eras = sorted(cfg.defaults_by_era, key=lambda e: e[0])
    for col in UNCERTAINTY_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
        default = np.full(len(out), np.nan)
        for start_year, vals in eras:
            if col in vals:
                default[years >= start_year] = vals[col]
        if eras and col in eras[0][1]:
            default[years < eras[0][0]] = eras[0][1][col] * 2  # 1900 öncesi: iki kat
        miss = out[col].isna().to_numpy()
        out[f"{col}_imputed"] = miss
        out.loc[miss, col] = default[miss]
    return out


def compute_weights(df: pd.DataFrame, cfg: WeightConfig | None = None) -> pd.Series:
    """Her olay için (0, 1] aralığında eğitim ağırlığı."""
    cfg = cfg or WeightConfig()
    filled = fill_uncertainties(df, cfg)
    s = np.zeros(len(filled))
    for col, scale in cfg.scales.items():
        s += (filled[col].to_numpy(dtype=float) / scale) ** 2
    w = 1.0 / (1.0 + s)
    if cfg.source_factor:
        w *= filled["source"].map(cfg.source_factor).fillna(1.0).to_numpy()
    return pd.Series(np.clip(w, cfg.floor, 1.0), index=df.index, name="weight")
