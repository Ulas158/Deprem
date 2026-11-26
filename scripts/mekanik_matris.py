"""Hücre ağını ve elastik etkileşim matrisini hesapla (bir kez; ~30 dk)."""
from pathlib import Path
import time
import numpy as np
from deprem.mekanik import geometri as g

ROOT = Path(__file__).resolve().parents[1]
out = ROOT / "data" / "mekanik"
out.mkdir(parents=True, exist_ok=True)
iz = g.fay_izleri(ROOT / "data/raw/gem_active_faults.geojson")
c = g.hucreler(iz)
c.to_parquet(out / "hucreler.parquet")
t = time.time()
Kt, Kn = g.etkilesim_matrisi(c)
np.save(out / "K_tau.npy", Kt)
np.save(out / "K_sn.npy", Kn)
print("bitti", len(c), time.time() - t)
