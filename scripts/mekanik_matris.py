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
if "mta" in SURUM:
    # MTA Türkiye Diri Fay Haritası 2026 (Elmacı vd. 2026); ham veri data/mta/tdfh_marmara.geojson
    from deprem.mekanik import mta
    import pandas as pd
    iz, sinir = mta.fay_izleri_mta(ROOT / "data/mta/tdfh_marmara.geojson", ROOT / "data/raw/gem_active_faults.geojson")
    pd.DataFrame(sinir, columns=["kol", "kesim", "lon0", "lon1"]).to_csv(out / "kesimler.csv", index=False)
else:
    iz = g.fay_izleri(ROOT / "data/raw/gem_active_faults.geojson", yan_kollar="v2" in SURUM)
# "kaba": 2 km kutu (ikinci tur arama ağı)
c = g.hucreler(iz, ds=2.0, dz=2.0) if "kaba" in SURUM else g.hucreler(iz)
if "mta" in SURUM:
    # hücreye MTA kesim adı (kendi kolunda, boylamına göre)
    c["kesim"] = ""
    for r in sinir:
        m = (c.kol == r[0]) & (c.lon >= r[2]) & (c.lon <= r[3])
        c.loc[m, "kesim"] = r[1]
    # derinlik (D90) profili: eski ağın x'inden boylama, oradan yeni x'e taşınır
    eski = pd.read_parquet(ROOT / "data/mekanik_v2_kaba/hucreler.parquet")
    e0 = eski[(eski.kol == 0) & (eski.z_ust == 0)].sort_values("x_km")
    y0 = c[(c.kol == 0) & (c.z_ust == 0)].sort_values("lon")
    d90 = pd.read_csv(ROOT / "data/mekanik_v2_kaba/d90_profil.csv")
    lon_d = np.interp(d90.x_km, e0.x_km, e0.lon)
    d90["x_km"] = np.interp(lon_d, y0.lon, y0.x_km)
    d90.to_csv(out / "d90_profil.csv", index=False)
c.to_parquet(out / "hucreler.parquet")
t = time.time()
Kt, Kn = g.etkilesim_matrisi(c)
np.save(out / "K_tau.npy", Kt)
np.save(out / "K_sn.npy", Kn)
print("bitti", len(c), time.time() - t)
