"""Hücre özellikleri: ısı akısı (deprem derinliğinden), levha hızı (GPS), sığ krip önceliği."""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import termal

# GPS profil uyumundan (bu çalışma, 2008-2026 NGL hızları) ve literatürden
# (Reilinger vd. 2006; Ergintav vd. 2014) kuzey kol kayma hızı, mm/yıl
V_PL_SEGMENT = {0: 18.5, 1: 18.5, 2: 22.0, 3: 18.0}
V_PL_SD = 2.0

# Sığ krip öncülü: tekrarlayan depremler ve deniz tabanı jeodezisi Orta
# Marmara'da (Batı Yükselimi - Orta Çukur) kısmi sığ krip olduğunu
# düşündürüyor (Schmittbuhl vd. 2016; Bohnhoff vd. 2017; Yamamoto vd. 2019).
# Kumburgaz ve Adalar için deniz tabanı ölçümleri kilitli olduğunu gösteriyor
# (Lange vd. 2019; Ergintav vd. 2014). Boylam aralığı, ek krip derinliği (km) aralığı
SIG_KRIP = [((27.75, 28.05), (0.0, 6.0))]


def isi_akisi_profili(d90: pd.DataFrame, cells: pd.DataFrame, d_min=8.0, d_max=20.0, rng=None):
    """Her hücre için ısı akısı. D90: AFAD ve EMSC ortalaması; toplulukta aralarından örneklenir."""
    a = d90["d90_afad"].to_numpy()
    e = d90["d90_emsc"].to_numpy()
    h = d90["d90_hepsi"].to_numpy()
    lo = np.where(np.isfinite(a) & np.isfinite(e), np.minimum(a, e), h)
    hi = np.where(np.isfinite(a) & np.isfinite(e), np.maximum(a, e), h)
    if rng is None:
        d = 0.5 * (lo + hi)
    else:
        u = rng.uniform()  # tek bir kaynak tercihi tüm fay için (sistematik fark)
        d = lo + u * (hi - lo)
    x = d90["x_km"].to_numpy()
    ok = np.isfinite(d)
    d = np.interp(cells["x_km"].to_numpy(), x[ok], d[ok])
    # 20 km'lik pürüzsüzleştirme zaten bin genişliğinde; aşırı uçları kırp
    d = np.clip(d, d_min, d_max)
    return termal.isi_akisi_derinlikten(d), d


def levha_hizi(cells: pd.DataFrame, rng=None):
    seg = cells["seg"].to_numpy()
    v = np.array([V_PL_SEGMENT[s] for s in seg], float)
    if rng is not None:
        dv = {s: rng.normal(0, V_PL_SD) for s in V_PL_SEGMENT}
        v = v + np.array([dv[s] for s in seg])
    return v / 1000.0


def sig_krip_maskesi(cells: pd.DataFrame, rng=None):
    """Ek sığ krip: bu hücreler hıza göre güçlenen kabul edilir."""
    m = np.zeros(len(cells), bool)
    for (l0, l1), (z0, z1) in SIG_KRIP:
        zk = 0.5 * (z0 + z1) if rng is None else rng.uniform(z0, z1)
        m |= cells["lon"].between(l0, l1).to_numpy() & (cells["z_km"].to_numpy() < zk)
    return m
