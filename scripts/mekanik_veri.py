"""Marmara mekanik kırılma modeli için veri: katalog (AFAD, USGS, EMSC) ve GPS.

Kullanım: python scripts/mekanik_veri.py
Çıktı: data/mekanik/katalog_marmara.parquet, data/raw/gnss/tenv3/*.tenv3
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from deprem.catalog.clean import basic_clean, merge_sources, moment_magnitude
from deprem.catalog.sources import FetchSpec, fetch
from deprem.observations import gnss
from deprem.region import Region

ROOT = Path(__file__).resolve().parents[1]
MARMARA = Region("Marmara_mekanik", 39.9, 41.3, 25.8, 31.9, max_depth_km=40.0)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    out = ROOT / "data" / "mekanik"
    out.mkdir(parents=True, exist_ok=True)
    raw = ROOT / "data" / "raw"
    end = pd.Timestamp.now(tz="UTC").floor("D")
    frames = []
    for src, start, mmin, chunk in [("usgs", "1900-01-01", 2.5, "3650D"),
                                    ("emsc", "1998-01-01", 1.5, "365D"),
                                    ("isc_afad", "2005-01-01", 1.5, "365D"),
                                    ("afad", str((end - pd.Timedelta(days=13)).date()), 1.0, "14D")]:
        try:
            spec = FetchSpec(MARMARA, pd.Timestamp(start, tz="UTC"), end, mmin, chunk)
            d = fetch(src, spec, raw)
            logging.info("%s: %d kayıt", src, len(d))
            frames.append(d)
        except Exception as e:  # bir kaynak çökerse diğerleri devam etsin
            logging.warning("%s alınamadı: %s", src, e)
    cat = merge_sources([basic_clean(f, MARMARA) for f in frames if len(f)])
    cat["mw"] = moment_magnitude(cat)
    cat.to_parquet(out / "katalog_marmara.parquet")
    logging.info("birleşik katalog: %d", len(cat))

    st = gnss.station_list(raw / "gnss", MARMARA, pad_deg=0.6, min_years=2.5)
    st.to_csv(out / "gnss_istasyonlar.csv", index=False)
    gnss.fetch_all(raw / "gnss", st["site"].tolist())
    logging.info("GNSS istasyon: %d", len(st))


if __name__ == "__main__":
    main()
