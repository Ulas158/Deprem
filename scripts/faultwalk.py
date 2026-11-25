"""Depremden depreme, fay fay öğrenen modelin sıralı yürüyüşü.

    python scripts/faultwalk.py --prepare                 # veri (katalog, fay birimleri, Coulomb)
    python scripts/faultwalk.py etas partial pooled ...   # sürümleri çalıştır (paralel: --jobs 4)
    python scripts/faultwalk.py --nulls 10                # şans kontrolü (karıştırılmış fay/Coulomb)
    python scripts/faultwalk_report.py                    # rapor

Çıktılar runs/faultwalk/ altında. Yöntem: src/deprem/faultwalk.py
HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "faultwalk"
CACHE = OUT / "walkdata.pkl"
FREEZE_YEAR = 1900


_CORE = None


def _init(counter):
    """Her işçi süreci tek bir çekirdeğe sabitlenir (iş parçacığı çekişmesini önler)."""
    global _CORE
    with counter.get_lock():
        _CORE = counter.value % os.cpu_count()
        counter.value += 1
    try:
        os.sched_setaffinity(0, {_CORE})
    except (AttributeError, OSError):
        pass


def _job(spec):
    import numpy as np
    from deprem.faultwalk import VARIANTS, walk
    from deprem.faultwalk_data import year_start_days
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    name, kind, seed, every = spec
    with open(CACHE, "rb") as fh:
        W = pickle.load(fh)
    out = OUT / f"{name}.csv"
    if out.exists():
        return name
    kw = dict(refit_every=every, progress=name, maxiter=30)
    if kind == "fixed":
        v = VARIANTS["partial"]
        kw["freeze_at"] = year_start_days(FREEZE_YEAR)
    elif kind in VARIANTS:
        v = VARIANTS[kind]
    elif kind == "null_fault":            # fay birimlerini hücreler arasında karıştır
        v = VARIANTS["partial"]
        rng = np.random.default_rng(seed)
        uoc = W.unit_of_cell.copy()
        fay = np.where(W.units["kind"].to_numpy()[uoc] == "fay")[0]
        uoc[fay] = rng.permutation(uoc[fay])
        kw["unit_override"] = uoc
        kw["cell_probs"] = False
    elif kind == "null_cfs":              # Coulomb işaretlerini kaynaklar arasında karıştır
        v = VARIANTS["partial"]
        rng = np.random.default_rng(seed)
        kw["cfs_override"] = W.cfs[rng.permutation(W.cfs.shape[0])]
        kw["cell_probs"] = False
    else:
        raise ValueError(kind)
    df, snaps = walk(W, v, **kw)
    df.to_csv(out, index=False)
    with open(OUT / f"{name}.params.json", "w") as fh:
        json.dump(snaps, fh, ensure_ascii=False)
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*")
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--every", type=int, default=1)
    ap.add_argument("--nulls", type=int, default=0)
    ap.add_argument("--null-every", type=int, default=20)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    OUT.mkdir(parents=True, exist_ok=True)
    if args.prepare or not CACHE.exists():
        from deprem.faultwalk_data import build
        if CACHE.exists():
            CACHE.unlink()
        build(ROOT, cache=CACHE)
    specs = [(v, v, 0, args.every) for v in args.variants]
    if args.nulls:
        e = args.null_every
        specs += [(f"partial_e{e}", "partial", 0, e), (f"pooled_e{e}", "pooled", 0, e),
                  (f"no_cfs_e{e}", "no_cfs", 0, e)]
        for s in range(args.nulls):
            specs += [(f"null_fault_{s}", "null_fault", s, e), (f"null_cfs_{s}", "null_cfs", s, e)]
    if not specs:
        return
    import multiprocessing as mp
    counter = mp.Value("i", 0)
    with ProcessPoolExecutor(args.jobs, initializer=_init, initargs=(counter,)) as ex:
        for name in ex.map(_job, specs):
            logging.info("bitti: %s", name)


if __name__ == "__main__":
    main()
