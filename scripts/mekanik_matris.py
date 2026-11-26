"""Hücre ağını ve elastik etkileşim matrisini hesapla (bir kez; ~30 dk)."""
from pathlib import Path
import time
import numpy as np
import os
from deprem.mekanik import geometri as g

ROOT = Path(__file__).resolve().parents[1]
SURUM = os.environ.get("DEPREM_MEKANIK_SURUM", "")
out = ROOT / "data" / ("mekanik" + (f"_{SURUM}" if SURUM else ""))
out.mkdir(parents=True, exist_ok=True)
iz = g.fay_izleri(ROOT / "data/raw/gem_active_faults.geojson", yan_kollar="v2" in SURUM)
# "kaba": 2 km kutu (ikinci tur arama ağı)
c = g.hucreler(iz, ds=2.0, dz=2.0) if "kaba" in SURUM else g.hucreler(iz)
c.to_parquet(out / "hucreler.parquet")
t = time.time()
Kt, Kn = g.etkilesim_matrisi(c)
np.save(out / "K_tau.npy", Kt)
np.save(out / "K_sn.npy", Kn)
print("bitti", len(c), time.time() - t)
