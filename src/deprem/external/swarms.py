"""Deprem sürüleri ve akışkan göçü (gözenek basıncı difüzyonu) göstergeleri.

Fiziksel arka plan: akışkan basıncı p arttığında fay üzerindeki etkin
normal gerilme (σn − p) azalır ve Coulomb kırılma ölçütü
τ ≥ μ(σn − p) + C daha kolay sağlanır. Basınç cephesi kayaç içinde
difüzyonla ilerlediğinde tetiklenen olaylar r ≈ √(4πDt) zarfı içinde
kalır (Shapiro vd. 1997). Burada D hidrolik yayınımdır (tipik olarak
0,01–10 m²/s). Ana şok–artçı dizilerinden farklı olarak sürülerde en büyük
olay dizinin başında değildir ve büyüklükler birbirine yakındır.

Örnek çalışmalar: Matsushiro 1965–67 sürüsü (yeraltından su çıkışı),
Kuzey Bohemya/Vogtland sürüleri (CO2 kaynaklı), Hainzl 2004, Shelly vd.
2016 (Long Valley), Ross vd. 2020 (Cahuilla sürüsü), L'Aquila 2009 öncesi
akışkan rolü (Terakawa vd. 2010), Marmara'da gaz çıkışı–sismisite
ilişkisi (Géli vd. 2018).

Bu modülün çıktıları **hipotez göstergeleridir**; sürü sınıflaması basit
bir kurala dayanır ve konum hataları (birkaç km) küçük dizilerde D
tahminini büyük ölçüde etkiler.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ..catalog.clean import epoch_seconds


def _xyz(lat, lon, depth, lat0=39.0, lon0=35.0):
    x = (np.asarray(lon) - lon0) * 111.195 * np.cos(np.radians(lat0))
    y = (np.asarray(lat) - lat0) * 111.195
    return np.column_stack([x, y, np.asarray(depth, dtype=float)])


def link_clusters(events: pd.DataFrame, r_km: float = 10.0, dt_days: float = 5.0) -> np.ndarray:
    """Tek bağlantılı uzay-zaman kümeleme: iki olay ≤ r_km ve ≤ dt_days ise aynı küme.

    Döndürür: her olay için küme kimliği (events sırasıyla).
    """
    n = len(events)
    if n == 0:
        return np.zeros(0, dtype=int)
    t = epoch_seconds(events["time"]) / 86400.0
    order = np.argsort(t, kind="stable")
    xyz = _xyz(events["latitude"], events["longitude"], events["depth_km"].fillna(10.0))
    xy = xyz[order, :2]
    ts = t[order]
    parent = np.arange(n)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    # zamanı dt_days ölçeğine göre 3. boyut yaparak tek KD-ağacı:
    # |Δt|·(r/dt) ≤ r ve |Δxy| ≤ r koşulu, Chebyshev normunda yaklaşık eşdeğer.
    pts = np.column_stack([xy, ts * (r_km / dt_days)])
    tree = cKDTree(pts)
    pairs = tree.query_pairs(r_km, p=np.inf, output_type="ndarray")
    if len(pairs):
        d = np.hypot(*(xy[pairs[:, 0]] - xy[pairs[:, 1]]).T)
        pairs = pairs[d <= r_km]
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    roots = np.array([find(i) for i in range(n)])
    out = np.empty(n, dtype=int)
    out[order] = roots
    return out


def cluster_stats(events: pd.DataFrame, cluster: np.ndarray, min_n: int = 5,
                  swarm_gap: float = 0.5, swarm_pos: float = 0.1) -> pd.DataFrame:
    """Küme başına özet: büyüklük, süre, sürü bayrağı, difüzyon katsayısı.

    Sürü kuralı: (en büyük − ikinci en büyük) ≤ ``swarm_gap`` VE en büyük olay
    dizinin ilk ``swarm_pos`` kesiminde değil. ``D_m2_s``: Shapiro zarfı
    r² = 4πDt ile her olay için D_i = r_i²/(4π t_i) (t_i > 1 sa) değerlerinin
    %90'lık dilimi; r ilk olaydan 3B uzaklıktır. En az 10 olay gerekir.
    """
    ev = events.assign(_c=cluster)
    t_all = epoch_seconds(ev["time"])
    ev["_t"] = t_all
    xyz = _xyz(ev["latitude"], ev["longitude"], ev["depth_km"].fillna(10.0))
    ev["_x"], ev["_y"], ev["_z"] = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    rows = []
    for cid, g in ev.groupby("_c"):
        if len(g) < min_n:
            continue
        g = g.sort_values("_t")
        m = g["mag"].to_numpy()
        top = np.sort(m)[::-1]
        imax = int(np.argmax(m))
        is_swarm = (top[0] - top[1] <= swarm_gap) and (imax / (len(g) - 1) >= swarm_pos)
        dt = g["_t"].to_numpy() - g["_t"].iloc[0]
        r = np.sqrt((g["_x"] - g["_x"].iloc[0]) ** 2 + (g["_y"] - g["_y"].iloc[0]) ** 2
                    + (g["_z"] - g["_z"].iloc[0]) ** 2).to_numpy() * 1000.0
        ok = dt > 3600
        D = np.quantile(r[ok] ** 2 / (4 * np.pi * dt[ok]), 0.9) if ok.sum() >= 10 else np.nan
        rows.append({
            "cluster": cid, "n": len(g), "t_start": g["time"].iloc[0], "t_end": g["time"].iloc[-1],
            "duration_days": dt[-1] / 86400.0, "lat": g["latitude"].mean(), "lon": g["longitude"].mean(),
            "depth_km": g["depth_km"].mean(), "mag_max": top[0], "mag_gap": top[0] - top[1],
            "max_pos": imax / (len(g) - 1), "is_swarm": bool(is_swarm), "D_m2_s": D,
            "radius_km": float(np.quantile(r, 0.9) / 1000.0),
        })
    return pd.DataFrame(rows)
