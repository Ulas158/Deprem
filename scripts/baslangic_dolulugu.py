"""Pencerelerin 1964 başlangıç durumu büyük modelden (Ulaş: "1964'ten önceki depremleri ne yapacaksın").

Büyük model (2 km, MTA) 1500'den 1964'e kadar koşulur: 1509'dan itibaren tarihsel büyük depremler ve
1905-1963 USGS M>=5 depremleri işlenir. 1964'teki doluluk (gerilme / kırılma eşiği) her hücre için
kaydedilir; pencereler 1964'te bu doluluktan başlar (ozumseme.PENCERE_BASLANGIC). Böylece 1964 öncesi
depremlerin faya bıraktığı gerilme (ör. 1912 Mürefte'nin Ganos'u boşaltması) pencerelere aktarılır.
Parametreler: --takim literatur (ilk tur); büyük model eğitildikten sonra en iyi takımla yeniden üretilir.
Çıktı: data/mekanik_mta_kaba/doluluk_1964.parquet
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
import mekanik_ikinci_tur as M  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--takim", default="literatur")
    ap.add_argument("--t", type=float, default=1964.0)
    a = ap.parse_args()
    th = dict(oz.baslangic_theta(), zayif=0.0)
    if a.takim != "literatur":
        th.update(json.loads(Path(a.takim).read_text()))
    ag, ol, _ = M._hazirla("buyuk", "mta_kaba", True)
    ol = [o for o in ol if o.t < a.t]
    A = M.AYAR["buyuk"]
    m, _ = oz.Kosu(ag, th, ol, 1, A["isinma"], A["t_bas"], a.t, a.t).calistir()
    c = ag.cells
    d = pd.DataFrame({"kol": c.kol.to_numpy(), "x_km": c.x_km.to_numpy(), "z_km": c.z_km.to_numpy(),
                      "lon": c.lon.to_numpy(), "lat": c.lat.to_numpy(), "doluluk": oz.doluluk(m)})
    d = d[np.isfinite(d.doluluk)]
    f = ROOT / "data" / "mekanik_mta_kaba" / f"doluluk_{int(a.t)}.parquet"
    d.to_parquet(f)
    (ROOT / "reports" / "tables" / f"doluluk_{int(a.t)}_ozet.json").write_text(json.dumps(
        {"takim": a.takim, "n": int(len(d)), "ortalama": float(d.doluluk.mean()),
         "kol0_boylam_0.2": {f"{l:.1f}": round(float(d[(d.kol == 0) & (np.abs(d.lon - l) < 0.1)].doluluk.mean()), 3)
                             for l in np.arange(26.3, 31.5, 0.2)}}, indent=1))
    print(f, len(d), d.doluluk.describe().round(3).to_dict())


if __name__ == "__main__":
    main()
