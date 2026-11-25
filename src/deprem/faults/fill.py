"""Eksik eleman parametrelerini doldurma (ortak yardımcı).

Kural (proje kararı): ölçümü olmayan bir fay elemanının değeri boş
bırakılmaz; ancak **doldurulduğu açıkça işaretlenir** ve belirsizliği
gerçek veriye uzaklıkla büyür. Eğitimde doldurulan değerler öğrenilebilir
parametre olur: doldurma değeri öncül ortalama, belirsizlik öncül genişliği
(``trainable_vector``, ``gaussian_log_prior``).

Doldurma sırası (her parametre için ayrı):

1. ``FILL_ALONG`` (1): aynı fay üzerinde, doğrultu boyunca sağdaki ve
   soldaki en yakın gerçek değerin uzaklık ağırlıklı ortalaması (bir taraf
   yoksa diğer taraf). Aynı ``along_km`` konumundaki farklı derinlikler
   birlikte ele alınır.
2. ``FILL_NEIGHBOR`` (2): fay üzerinde hiç veri yoksa çevredeki faylardaki
   en yakın ``k`` gerçek değerin ters uzaklık ağırlıklı ortalaması
   (``neighbor_km`` içinde).
3. ``FILL_REGIONAL`` (3): o da yoksa tüm gerçek değerlerin ortancası.
4. ``FILL_DEFAULT`` (4): hiç ölçüm yoksa ve ``defaults`` içinde literatür
   değeri varsa o değer (geniş belirsizlikle).
5. Hiçbiri yoksa NaN kalır ve bayrak ``FILL_NONE`` (9) olur: hiçbir
   ölçümü olmayan bir büyüklük (ör. radon) uydurulmaz.

Her ``col`` için iki sütun eklenir: ``<col>__fill`` (0 = gerçek veri) ve
``<col>__sd`` (belirsizlik; gerçek veride ``measurement_sd``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

FILL_REAL, FILL_ALONG, FILL_NEIGHBOR, FILL_REGIONAL, FILL_DEFAULT, FILL_NONE = 0, 1, 2, 3, 4, 9
FILL_LABELS = {0: "gerçek veri", 1: "fay boyunca komşulardan", 2: "komşu faylardan",
               3: "bölgesel ortanca", 4: "literatür varsayılanı", 9: "veri yok"}

# Hiç ölçüm yokken kullanılabilecek kaba literatür değerleri: (ortalama, sd)
LITERATURE_DEFAULTS = {
    "pore_pressure_lambda": (0.4, 0.2),        # hidrostatik ≈ 0,37
    "vp_vs_ratio": (1.75, 0.06),
    "heat_flow_mw_m2": (80.0, 30.0),           # Anadolu ortalaması mertebesi; batıda yüksek
    "curie_depth_km": (15.0, 6.0),
    "crustal_thickness_km": (36.0, 7.0),
    "friction_a_minus_b": (-0.003, 0.005),
    "locking_ratio": (0.7, 0.3),
}


def _xy(lat, lon, lat0=39.0, lon0=35.0):
    return np.column_stack([(np.asarray(lon) - lon0) * 111.195 * np.cos(np.radians(lat0)),
                            (np.asarray(lat) - lat0) * 111.195])


def fill_missing(elements: pd.DataFrame, cols, length_scale_km: float = 20.0,
                 neighbor_km: float = 100.0, k: int = 8, measurement_sd: dict | None = None,
                 defaults: dict | None = LITERATURE_DEFAULTS,
                 fault_col: str = "fault_id", along_col: str = "along_km") -> pd.DataFrame:
    """``cols`` sütunlarını doldurup bayrak ve belirsizlik sütunlarıyla döndür.

    Belirsizlik: sd = s_bölge · √(1 − exp(−d / L)), d en yakın gerçek veriye
    uzaklık, L = ``length_scale_km``, s_bölge gerçek değerlerin standart
    sapması (gerçek veri az ise varsayılan sd). Komşu fay ve bölgesel
    doldurmada sd en az s_bölge'dir.
    """
    out = elements.copy()
    measurement_sd = measurement_sd or {}
    defaults = defaults or {}
    xy = _xy(out["lat"], out["lon"])
    fid = out[fault_col].to_numpy() if fault_col in out else np.zeros(len(out))
    along = out[along_col].to_numpy(dtype=float) if along_col in out else np.zeros(len(out))
    for col in [cols] if isinstance(cols, str) else cols:
        v = out[col].to_numpy(dtype=float) if col in out else np.full(len(out), np.nan)
        real = np.isfinite(v)
        flag = np.where(real, FILL_REAL, FILL_NONE).astype(int)
        dist = np.where(real, 0.0, np.nan)
        filled = v.copy()
        s_reg = np.nanstd(v[real]) if real.sum() >= 2 else np.nan
        if not (np.isfinite(s_reg) and s_reg > 0):
            s_reg = defaults.get(col, (np.nan, np.nan))[1]
        if not np.isfinite(s_reg) and real.any():
            s_reg = max(0.5 * abs(float(np.nanmedian(v[real]))), 1e-6)
        if real.any():
            # 1. fay boyunca sol/sağ
            for f in np.unique(fid):
                idx = np.flatnonzero(fid == f)
                r = idx[real[idx]]
                m = idx[~real[idx]]
                if len(r) == 0 or len(m) == 0:
                    continue
                # aynı doğrultu konumundaki derinlikleri birleştir
                a_r = pd.Series(v[r]).groupby(along[r]).mean()
                pos, val = a_r.index.to_numpy(), a_r.to_numpy()
                am = along[m]
                j = np.searchsorted(pos, am)
                left = np.clip(j - 1, 0, len(pos) - 1)
                right = np.clip(j, 0, len(pos) - 1)
                dl = np.where(j > 0, am - pos[left], np.inf)
                dr = np.where(j < len(pos), pos[right] - am, np.inf)
                wl = np.where(np.isfinite(dl), 1 / np.maximum(dl, 1e-3), 0)
                wr = np.where(np.isfinite(dr), 1 / np.maximum(dr, 1e-3), 0)
                filled[m] = (wl * val[left] + wr * val[right]) / (wl + wr)
                dist[m] = np.minimum(dl, dr)
                flag[m] = FILL_ALONG
            # 2. komşu faylar
            rest = np.flatnonzero(flag == FILL_NONE)
            if len(rest):
                ri = np.flatnonzero(real)
                tree = cKDTree(xy[ri])
                kk = min(k, len(ri))
                d, nn = tree.query(xy[rest], k=kk)
                d, nn = d.reshape(len(rest), kk), nn.reshape(len(rest), kk)
                w = np.where(d <= neighbor_km, 1 / np.maximum(d, 0.5) ** 2, 0)
                ok = w.sum(axis=1) > 0
                filled[rest[ok]] = (w[ok] * v[ri][nn[ok]]).sum(axis=1) / w[ok].sum(axis=1)
                dist[rest[ok]] = d[ok, 0]
                flag[rest[ok]] = FILL_NEIGHBOR
                # 3. bölgesel ortanca
                far = rest[~ok]
                filled[far] = np.nanmedian(v[real])
                dist[far] = d[~ok, 0]
                flag[far] = FILL_REGIONAL
        elif col in defaults:
            mean, sd = defaults[col]
            filled[:] = mean
            flag[:] = FILL_DEFAULT
        sd = np.full(len(out), np.nan)
        sd[flag == FILL_REAL] = measurement_sd.get(col, 0.0)
        a = flag == FILL_ALONG
        sd[a] = s_reg * np.sqrt(1 - np.exp(-dist[a] / length_scale_km))
        nb = flag == FILL_NEIGHBOR
        sd[nb] = s_reg * np.maximum(1.0, np.sqrt(dist[nb] / length_scale_km))
        sd[flag == FILL_REGIONAL] = 2 * s_reg
        if col in defaults:
            sd[flag == FILL_DEFAULT] = defaults[col][1]
        out[col] = filled
        out[f"{col}__fill"] = flag
        out[f"{col}__sd"] = sd
    return out


def fill_summary(filled: pd.DataFrame, cols) -> pd.DataFrame:
    rows = []
    for c in [cols] if isinstance(cols, str) else cols:
        vc = filled[f"{c}__fill"].value_counts(normalize=True)
        rows.append({"parametre": c, **{FILL_LABELS[k]: round(float(vc.get(k, 0)), 3) for k in FILL_LABELS}})
    return pd.DataFrame(rows)


def trainable_vector(filled: pd.DataFrame, col: str, min_sd: float = 1e-6):
    """Eğitimde öğrenilecek doldurulmuş değerler.

    Döndürür: (index, başlangıç = öncül ortalama, öncül sd). Gerçek veriler
    ve hiç doldurulamayanlar dışarıda kalır.
    """
    f = filled[f"{col}__fill"].to_numpy()
    sel = (f != FILL_REAL) & (f != FILL_NONE)
    idx = filled.index[sel]
    mu = filled.loc[idx, col].to_numpy(dtype=float)
    sd = np.maximum(filled.loc[idx, f"{col}__sd"].to_numpy(dtype=float), min_sd)
    return idx, mu, sd


def gaussian_log_prior(x, mu, sd) -> float:
    """Doldurulmuş değerler için öncül log-yoğunluğu (sabit terim hariç): −½Σ((x−μ)/σ)²."""
    z = (np.asarray(x) - np.asarray(mu)) / np.asarray(sd)
    return float(-0.5 * np.sum(z * z))


def training_weight(filled: pd.DataFrame, col: str) -> np.ndarray:
    """Basit modeller için ağırlık: gerçek veri 1, doldurulan 1/(1+(sd/s)²), veri yok 0."""
    f = filled[f"{col}__fill"].to_numpy()
    sd = filled[f"{col}__sd"].to_numpy(dtype=float)
    s = np.nanmax(sd[f != FILL_REAL]) if (f != FILL_REAL).any() else 1.0
    s = s if np.isfinite(s) and s > 0 else 1.0
    w = 1 / (1 + (np.nan_to_num(sd) / s) ** 2)
    w[f == FILL_REAL] = 1.0
    w[f == FILL_NONE] = 0.0
    return w
