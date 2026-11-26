"""Model değişkenlerinin etki payı: her parçayı çıkar (veya ekle), sıralı modeli baştan yürüt.

    python scripts/etki_payi.py --jobs 4            # bütün sürümler (runs/etki/)
    python scripts/etki_payi.py partial no_cfs      # yalnız bunlar
    python scripts/etki_payi_rapor.py               # rapor: reports/etki_payi.md

Her sürüm, fay fay sıralı modelin (src/deprem/faultwalk.py) aynı yürüyüşüdür: her deprem
yalnız ondan önceki depremlerle eğitilmiş modelle puanlanır. Ek olarak her aralık için
Marmara kutusunun kendi beklenen sayısı (bölge kompansatörü) yazılır; böylece Marmara puanı
"Marmara'da ne zaman ve nerede" sorusunun kendi olabilirliğidir.
HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
OUT = ROOT / "runs" / "etki"
FW = ROOT / "runs" / "faultwalk"
CACHE = FW / "walkdata.pkl"
ENERGY = FW / "energy.pkl"
DVV = ROOT / "data" / "dvv" / "east_marmara_dvv.csv"
FREEZE_YEAR = 1900
# Marmara kutusu (projenin Marmara canlı kaydıyla aynı: hücre merkezi enlem 40,0-41,3; boylam 25,5-30,5)
MARMARA = dict(lat=(39.99, 41.31), lon=(25.49, 30.51))

ALL = ["partial", "uniform", "etas", "no_fast", "no_slow", "no_trig", "no_far", "no_cross", "no_cfs",
       "no_load", "load_time", "no_fmap", "no_smap", "no_bmap", "no_slip", "pooled", "fixed",
       "partial_bval", "partial_chain", "partial_deficit", "partial_dvv"]


def region_mask(W):
    import numpy as np
    la, lo = MARMARA["lat"], MARMARA["lon"]
    return (W.lat_c >= la[0]) & (W.lat_c <= la[1]) & (W.lon_c >= lo[0]) & (W.lon_c <= lo[1]) & np.ones_like(W.lat_c, bool)


def _init(counter):
    with counter.get_lock():
        core = counter.value % os.cpu_count()
        counter.value += 1
    try:
        os.sched_setaffinity(0, {core})
    except (AttributeError, OSError):
        pass


def _job(spec):
    name, every = spec
    from deprem.faultwalk import VARIANTS, walk
    from deprem.faultwalk_data import time_matrix, year_start_days
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    out = OUT / f"{name}.csv"
    if out.exists():
        return name
    with open(CACHE, "rb") as fh:
        W = pickle.load(fh)
    kw = dict(refit_every=every, progress=name, maxiter=30, cell_probs=False, region=region_mask(W))
    v = VARIANTS["partial" if name == "fixed" else name]
    if name == "fixed":
        kw["freeze_at"] = year_start_days(FREEZE_YEAR)
    if v.load_mode == "time":
        W.extra["Ltime"] = time_matrix(ROOT, W)
    covs = [c for c in v.cov.split(",") if c]
    if "dvv" in covs:
        import pandas as pd
        from deprem.faultwalk_data import dvv_covariate
        W.extra.setdefault("cov", {})["dvv"] = dvv_covariate(W, pd.read_csv(DVV))
    elif covs:
        from deprem.fault_energy import build_energy, covariates
        E = build_energy(ROOT, W, cache=ENERGY)
        W.extra.setdefault("cov", {}).update(covariates(W, E))
    df, snaps = walk(W, v, **kw)
    df.to_csv(out, index=False)
    with open(OUT / f"{name}.params.json", "w") as fh:
        json.dump(snaps, fh, ensure_ascii=False)
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--every", type=int, default=1, help="kaç depremde bir yeniden eğitilsin (hepsinde aynı)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    if not CACHE.exists():
        from deprem.faultwalk_data import build
        build(ROOT, cache=CACHE)
    names = args.variants or ALL
    if any(n.startswith("partial_") and n != "partial_dvv" for n in names) and not ENERGY.exists():
        from deprem.fault_energy import build_energy
        with open(CACHE, "rb") as fh:
            build_energy(ROOT, pickle.load(fh), cache=ENERGY)
    import multiprocessing as mp
    counter = mp.Value("i", 0)
    with ProcessPoolExecutor(args.jobs, initializer=_init, initargs=(counter,)) as ex:
        for n in ex.map(_job, [(n, args.every) for n in names]):
            logging.info("bitti: %s", n)


if __name__ == "__main__":
    main()
