"""Katalog dışı gözlemleri indir ve özellik tablolarına çevir.

    python scripts/fetch_observations.py gnss insar tec geomag thermal repeaters
    python scripts/fetch_observations.py all
    python scripts/fetch_observations.py dvv --dvv-max-stations 12      # dalga biçimi, yavaş
    python scripts/fetch_observations.py mlcat --mlcat-days 2           # PhaseNet pilotu, yavaş

Ham veri data/raw/<kaynak>/ altına önbelleğe alınır; işlenmiş tablolar
data/obs/ altına yazılır ve eğitim hattı (scripts/train.py) bunları okur.
Hiçbir dosya yoksa ilgili özellik "eksik" sayılır, model yine çalışır.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd

from deprem.region import TURKEY

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OBS = ROOT / "data" / "obs"
YEARS = list(range(2010, 2027))


def do_gnss(args):
    from deprem.observations import gnss
    raw = RAW / "gnss"
    sites = gnss.station_list(raw, TURKEY)
    logging.info("GNSS: bölgede %d istasyon", len(sites))
    gnss.fetch_all(raw, sites["site"].tolist())
    sty = gnss.station_year_table(raw, sites, YEARS)
    OBS.mkdir(parents=True, exist_ok=True)
    sty.to_csv(OBS / "gnss_station_year.csv", index=False)
    logging.info("GNSS: %d istasyon-yıl satırı", len(sty))


def do_insar(args):
    raw = RAW / "insar"
    tifs = sorted(raw.glob("*.tif")) if raw.exists() else []
    if not tifs:
        logging.warning("InSAR: data/raw/insar/ altında hız GeoTIFF'i yok. Açık arşivde hazır hız "
                        "haritası bulunmuyor; LiCSBAS çıktısı veya Weiss vd. 2020 D/K hız dosyaları "
                        "buraya konunca eğitim hattı otomatik okur (bkz. deprem/observations/insar.py).")
    else:
        logging.info("InSAR: %d GeoTIFF bulundu; eğitim hattı doğrudan okuyacak", len(tifs))


def do_tec(args):
    from deprem.observations import ionosphere
    raw = RAW / "tec"
    ionosphere.fetch_kp(RAW / "kp")
    daily = ionosphere.build_daily(raw, TURKEY, pd.Timestamp(args.start), pd.Timestamp(args.end))
    OBS.mkdir(parents=True, exist_ok=True)
    daily.to_csv(OBS / "tec_daily.csv.gz", index=False)


def do_geomag(args):
    from deprem.observations import geomag
    d = geomag.build_daily(RAW / "geomag", pd.Timestamp(args.start), pd.Timestamp(args.end))
    OBS.mkdir(parents=True, exist_ok=True)
    d.to_csv(OBS / "ulf_daily.csv.gz", index=False)


def do_thermal(args):
    from deprem.observations import thermal
    d = thermal.build_daily(RAW / "thermal", TURKEY, pd.Timestamp(args.start), pd.Timestamp(args.end))
    OBS.mkdir(parents=True, exist_ok=True)
    d.to_csv(OBS / "thermal_daily.csv.gz", index=False)


def do_repeaters(args):
    from deprem.observations import repeaters
    cat = pd.read_csv(ROOT / "data/processed/turkey_catalog.csv.gz")
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    # Not: eğitim hattı aileleri her yıl için yalnız geçmiş veriyle yeniden bulur;
    # bu dosya inceleme ve harita içindir.
    fam = repeaters.find_families(cat)
    OBS.mkdir(parents=True, exist_ok=True)
    fam.to_csv(OBS / "repeater_events.csv.gz", index=False)


def do_dvv(args):
    from deprem.observations import dvv
    out = dvv.run(RAW / "waveforms", OBS, TURKEY, max_stations=args.dvv_max_stations,
                  start=pd.Timestamp(args.dvv_start), end=pd.Timestamp(args.end),
                  step_days=args.dvv_step_days)
    logging.info("dv/v: %s", out)


def do_mlcat(args):
    from deprem.observations import mlcatalog
    mlcatalog.pilot(RAW / "waveforms", OBS, days=args.mlcat_days)


TASKS = {"gnss": do_gnss, "insar": do_insar, "tec": do_tec, "geomag": do_geomag,
         "thermal": do_thermal, "repeaters": do_repeaters, "dvv": do_dvv, "mlcat": do_mlcat}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tasks", nargs="+", choices=list(TASKS) + ["all"])
    ap.add_argument("--start", default="2010-01-01")
    ap.add_argument("--end", default=pd.Timestamp.now().strftime("%Y-%m-%d"))
    ap.add_argument("--dvv-start", default="2012-01-01")
    ap.add_argument("--dvv-step-days", type=int, default=5)
    ap.add_argument("--dvv-max-stations", type=int, default=12)
    ap.add_argument("--mlcat-days", type=int, default=1)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    tasks = [t for t in TASKS if t not in ("dvv", "mlcat")] if "all" in args.tasks else args.tasks
    for t in tasks:
        logging.info("== %s", t)
        TASKS[t](args)


if __name__ == "__main__":
    main()
