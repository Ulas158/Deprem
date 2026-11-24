"""Türkiye kataloglarını çek, birleştir, temizle ve kaydet.

Kullanım:
    python scripts/fetch_catalogs.py                    # varsayılanlar
    python scripts/fetch_catalogs.py --afad-export data/raw/afad_export.csv

Çıktı: data/processed/turkey_catalog.parquet (veya .csv)
Ham yanıtlar data/raw/<kaynak>/ altında önbelleğe alınır.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from deprem.catalog.clean import basic_clean, merge_sources, moment_magnitude
from deprem.catalog.sources import FetchSpec, fetch, load_afad_export
from deprem.region import TURKEY

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--end", default=pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%d"))
    ap.add_argument("--emsc-start", default="1998-01-01")
    ap.add_argument("--emsc-min-mag", type=float, default=2.0)
    ap.add_argument("--usgs-start", default="1900-01-01")
    ap.add_argument("--usgs-min-mag", type=float, default=4.0)
    ap.add_argument("--afad-days", type=int, default=14,
                    help="AFAD açık API'sinden çekilecek son gün sayısı")
    ap.add_argument("--afad-export", type=Path, nargs="*", default=[],
                    help="AFAD web kataloğundan elle indirilmiş dosyalar")
    ap.add_argument("--isc-afad-start", default="2010-01-01",
                    help="AFAD'ın ISC'deki çözümlerinin başlangıcı")
    ap.add_argument("--no-afad", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    raw = ROOT / "data" / "raw"
    end = pd.Timestamp(args.end, tz="UTC")
    frames = []

    usgs = fetch("usgs", FetchSpec(TURKEY, pd.Timestamp(args.usgs_start, tz="UTC"), end,
                                   args.usgs_min_mag, chunk="3650D"), raw)
    logging.info("USGS: %d olay", len(usgs))
    frames.append(usgs)

    emsc = fetch("emsc", FetchSpec(TURKEY, pd.Timestamp(args.emsc_start, tz="UTC"), end,
                                   args.emsc_min_mag, chunk="30D"), raw)
    logging.info("EMSC: %d olay", len(emsc))
    frames.append(emsc)

    if not args.no_afad:
        try:
            afad = fetch("afad", FetchSpec(TURKEY, end - pd.Timedelta(days=args.afad_days), end,
                                           0.0, chunk="1D"), raw)
            logging.info("AFAD (açık API, son %d gün): %d olay", args.afad_days, len(afad))
            frames.append(afad)
        except Exception as e:  # ağ izni yoksa diğer kaynaklarla devam et
            logging.warning("AFAD çekilemedi: %s", e)
    if not args.no_afad:
        isc = fetch("isc_afad", FetchSpec(TURKEY, pd.Timestamp(args.isc_afad_start, tz="UTC"), end,
                                          0.0, chunk="10D"), raw)
        logging.info("AFAD arşivi (ISC üzerinden): %d olay", len(isc))
        frames.append(isc)
    for p in args.afad_export:
        frames.append(load_afad_export(p))

    frames = [basic_clean(f, TURKEY) for f in frames]
    merged = merge_sources(frames, priority=("afad", "emsc", "usgs"))
    merged["mw_best"] = moment_magnitude(merged)
    out = ROOT / "data" / "processed"
    out.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out / "turkey_catalog.csv.gz", index=False)
    logging.info("Birleşik katalog: %d olay -> %s", len(merged), out / "turkey_catalog.csv.gz")


if __name__ == "__main__":
    main()
