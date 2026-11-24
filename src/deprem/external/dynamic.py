"""Uzak büyük depremlerin yüzey dalgalarıyla dinamik tetikleme.

Büyük bir depremin yüzey dalgaları binlerce km uzakta kısa süreli
(saniye–dakika) gerilme salınımları (~0,01–0,1 MPa) yaratır. Bu, özellikle
jeotermal/volkanik ve akışkan bakımından zengin bölgelerde dakikalar–günler
içinde küçük deprem artışlarına yol açabilir (Hill vd. 1993 Landers;
Gomberg vd. 2001 Denali; Velasco vd. 2008; Brodsky & van der Elst 2014).
Büyük (M≥6) uzak tetikleme örnekleri tartışmalı ve çok nadirdir.

Ölçüt: her uzak olay için, yerel olayların ``after_days`` içindeki sayısı
ile önceki ``before_days`` içindeki sayısı karşılaştırılır; tüm uzak
olaylar üzerinde toplanıp β istatistiği (Matthews & Reasenberg 1988)
hesaplanır. β > 2 kabaca anlamlı artış sayılır.

Genlik vekili: Ms formülünden A* = M − 1,66·log10(Δ°); yalnızca
A* ≥ ``min_amp`` olan uzak olaylar kullanılır (dalga genliği vekili).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..catalog.clean import epoch_seconds
from ..geo import haversine_km


def amplitude_proxy(mag, dist_km):
    deg = np.maximum(np.asarray(dist_km) / 111.195, 1.0)
    return np.asarray(mag) - 1.66 * np.log10(deg)


def beta_statistic(n_after, n_total, p_after):
    n_total = np.asarray(n_total, dtype=float)
    exp = n_total * p_after
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(n_total > 0, (np.asarray(n_after) - exp) / np.sqrt(exp * (1 - p_after)), np.nan)


def remote_triggering(local: pd.DataFrame, remote: pd.DataFrame, lat_c, lon_c, groups,
                      min_dist_km: float = 1000.0, min_amp: float = 4.0,
                      before_days: float = 30.0, after_days: float = 5.0) -> pd.DataFrame:
    """Grup (hücre) başına uzak tetikleme β değeri.

    ``local``: yerel olaylar (``time`` ve grup sütunu ``groups`` adıyla).
    ``lat_c``, ``lon_c``: grup kimliğiyle indekslenmiş merkez koordinatları.
    """
    lt = epoch_seconds(local["time"])
    rt = epoch_seconds(remote["time"])
    lg = local[groups].to_numpy()
    ids = lat_c.index.to_numpy()
    n_after = pd.Series(0.0, index=ids)
    n_tot = pd.Series(0.0, index=ids)
    n_used = pd.Series(0, index=ids)
    bd, ad = before_days * 86400.0, after_days * 86400.0
    for k, (_, q) in enumerate(remote.iterrows()):
        d = haversine_km(lat_c.to_numpy(), lon_c.to_numpy(), q["latitude"], q["longitude"])
        use = (d >= min_dist_km) & (amplitude_proxy(q["mag"], d) >= min_amp)
        if not use.any():
            continue
        t0 = rt[k]
        sel_b = (lt >= t0 - bd) & (lt < t0)
        sel_a = (lt >= t0) & (lt < t0 + ad)
        ca = pd.Series(lg[sel_a]).value_counts().reindex(ids, fill_value=0)
        cb = pd.Series(lg[sel_b]).value_counts().reindex(ids, fill_value=0)
        u = pd.Series(use, index=ids)
        n_after += ca.where(u, 0)
        n_tot += (ca + cb).where(u, 0)
        n_used += u.astype(int)
    p = after_days / (after_days + before_days)
    return pd.DataFrame({"dyn_beta": beta_statistic(n_after, n_tot, p),
                         "dyn_n_after": n_after, "dyn_n_total": n_tot, "dyn_n_remote": n_used},
                        index=ids)
