"""Yapay zekâ ile zenginleştirilmiş katalog (PİLOT).

Ham sürekli sismogramlarda derin öğrenme faz seçicisi PhaseNet (Zhu &
Beroza 2019; SeisBench, Woollam vd. 2022) P ve S varışlarını işaretler;
basit bir ızgara-arama ilişkilendiricisi bunları depremlere dönüştürür.
Güney Kaliforniya'da bu yaklaşım katalogdaki olay sayısını ~10 kat artırdı
(Ross vd. 2019, Science) ve öncü/sürü süreçlerini görünür kıldı.

Bu ortamda tam 15 yıllık Türkiye kataloğu üretmek günler süren hesap ve
terabaytlarca ham veri ister. Bu modül:
* seçilen bir bölgede (varsayılan Marmara) birkaç günlük pilot çalıştırır,
* AFAD/EMSC kataloğundaki olay sayısıyla karşılaştırıp zenginleşme oranını
  (``mlcat_rate_ratio``) raporlar,
* tam çalıştırma için aynı kodu daha uzun süreye ölçeklemeye hazırdır.
Sonuçlar eğitime, tüm yıllar ve bölgeler için üretilene kadar **girmez**.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .dvv import FDSN, _get_hour, stations as ko_stations

log = logging.getLogger(__name__)

MARMARA = dict(min_lat=40.2, max_lat=41.2, min_lon=27.0, max_lon=30.0)
VP, VS = 6.0, 3.5


def _picks_for_hour(model, sta, t):
    st = _get_hour(sta, t)
    if st is None:
        return []
    hh = st.select(channel="HH?")
    st = hh if len(hh) >= 3 else st.select(channel="BH?")
    st.merge(fill_value=0)
    try:
        out = model.classify(st, P_threshold=0.3, S_threshold=0.3)
    except Exception as e:  # kanal/örnekleme sorunları
        log.warning("%s %s: %s", sta, t, e)
        return []
    return [(sta, p.phase, pd.Timestamp(p.peak_time.datetime), p.peak_value) for p in out.picks]


def associate(picks: pd.DataFrame, sta_xy: dict, box: dict, min_sta: int = 4, tol: float = 1.5,
              depth_km: float = 10.0):
    """Izgara-arama: her aday kaynak için P'den geri hesaplanan köken zamanlarını kümele."""
    lats = np.arange(box["min_lat"], box["max_lat"] + 1e-9, 0.1)
    lons = np.arange(box["min_lon"], box["max_lon"] + 1e-9, 0.1)
    P = picks[picks["phase"] == "P"].copy()
    if P.empty:
        return pd.DataFrame()
    P["ts"] = (P["time"] - P["time"].min()).dt.total_seconds()
    t_ref = P["time"].min()
    best = {}
    for la in lats:
        for lo in lons:
            d = np.array([np.hypot((sta_xy[s][1] - lo) * 85.0, (sta_xy[s][0] - la) * 111.0) for s in P["sta"]])
            tt = np.hypot(d, depth_km) / VP
            order = np.argsort(P["ts"].to_numpy() - tt)
            to = (P["ts"].to_numpy() - tt)[order]
            sts = P["sta"].to_numpy()[order]
            i = 0
            while i < len(to):
                j = np.searchsorted(to, to[i] + tol)
                n = len(set(sts[i:j]))
                if n >= min_sta:
                    key = int(round(to[i] / 5.0))
                    if key not in best or best[key][0] < n:
                        best[key] = (n, la, lo, to[i])
                    i = j
                else:
                    i += 1
    ev = pd.DataFrame([{"time": t_ref + pd.Timedelta(seconds=v[3]), "latitude": v[1], "longitude": v[2],
                        "n_sta": v[0]} for v in best.values()])
    return ev.sort_values("time").reset_index(drop=True) if not ev.empty else ev


def pilot(raw: Path, obs: Path, days: int = 1, start: str = "2024-03-01", box: dict | None = None,
          max_sta: int = 10):
    import seisbench.models as sbm
    from deprem.region import Region
    box = box or MARMARA
    model = sbm.PhaseNet.from_pretrained("instance")  # İtalya verisiyle eğitilmiş: bölgesel ağlara uygun
    reg = Region("pilot", box["min_lat"] - 0.3, box["max_lat"] + 0.3, box["min_lon"] - 0.3, box["max_lon"] + 0.3)
    st = ko_stations(reg, start_before=start)
    st = st.head(max_sta)
    sta_xy = {r.sta: (r.lat, r.lon) for r in st.itertuples()}
    hours = pd.date_range(start, periods=24 * days, freq="h")
    picks = []
    for t in hours:
        for s in st["sta"]:
            picks += _picks_for_hour(model, s, t)
        log.info("PhaseNet %s: toplam %d seçim", t, len(picks))
    pk = pd.DataFrame(picks, columns=["sta", "phase", "time", "prob"])
    ev = associate(pk, sta_xy, box)
    out = obs / "mlcat"
    out.mkdir(parents=True, exist_ok=True)
    pk.to_csv(out / "picks.csv.gz", index=False)
    ev.to_csv(out / "events.csv", index=False)
    # karşılaştırma: birleşik katalog
    cat_p = Path(__file__).resolve().parents[3] / "data/processed/turkey_catalog.csv.gz"
    summary = {"start": start, "days": days, "stations": list(sta_xy), "box": box,
               "n_ml_events": int(len(ev)), "n_picks": int(len(pk))}
    if cat_p.exists():
        cat = pd.read_csv(cat_p)
        cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601").dt.tz_localize(None)
        t1 = pd.Timestamp(start) + pd.Timedelta(days=days)
        c = cat[(cat["time"] >= pd.Timestamp(start)) & (cat["time"] < t1)
                & cat["latitude"].between(box["min_lat"], box["max_lat"])
                & cat["longitude"].between(box["min_lon"], box["max_lon"])]
        summary["n_catalog_events"] = int(len(c))
        summary["mlcat_rate_ratio"] = round(len(ev) / max(len(c), 1), 2)
        if len(ev):
            matched = sum(((ev["time"] - t).abs() < pd.Timedelta(seconds=10)).any() for t in c["time"])
            summary["catalog_events_recovered"] = int(matched)
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    log.info("ML katalog pilotu: %s", summary)
    return summary
