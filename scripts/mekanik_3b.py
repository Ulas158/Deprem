"""Modelin 3B görünümü: Marmara kıyıları üstünde fay duvarı, hücre hücre doluluk.

Kullanım: python scripts/mekanik_3b.py [--durum runs/mekanik/isinma_kalin.pkl] [--senaryo kalin]
Çıktı: reports/figures/mekanik_3b.png ve reports/mekanik_3b/veri.json (etkileşimli sayfa için)
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import deney, geometri  # noqa: E402

KUTU = (26.0, 31.6, 40.15, 41.45)
ABARTI = 2.0  # derinlik abartması (görsellik için)
ETIKET = [("Saros-Ganos", 26.55), ("Tekirdağ", 27.5), ("Orta Marmara", 27.95), ("Kumburgaz", 28.35),
          ("Adalar", 28.85), ("Çınarcık", 29.12), ("İzmit", 29.85), ("Sapanca-Akyazı", 30.4),
          ("Düzce", 31.15)]
SEHIR = [("İstanbul", 28.98, 41.01), ("Tekirdağ", 27.51, 40.98), ("İzmit", 29.92, 40.77),
         ("Bursa", 29.06, 40.19), ("Yalova", 29.27, 40.65), ("Düzce", 31.16, 40.84)]


def xy(lon, lat):
    p = geometri.to_local_m(np.asarray(lat), np.asarray(lon), 0 * np.asarray(lat, float))
    return p[:, 0] / 1000.0, p[:, 1] / 1000.0


def kiyilar():
    import json as js
    f = js.loads((ROOT / "data/raw/ne/ne_10m_coastline.geojson").read_text())
    out = []
    for ft in f["features"]:
        g = ft["geometry"]
        parts = [g["coordinates"]] if g["type"] == "LineString" else g["coordinates"]
        for p in parts:
            c = np.asarray(p)
            ok = (c[:, 0] > KUTU[0]) & (c[:, 0] < KUTU[1]) & (c[:, 1] > KUTU[2]) & (c[:, 1] < KUTU[3])
            if ok.sum() < 2:
                continue
            # kutu dışına çıkan kısımları böl
            idx = np.flatnonzero(ok)
            for seg in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
                if len(seg) >= 2:
                    out.append(c[seg])
    return out


def hucre_kose(cells):
    """Her hücre için 4 köşe (x, y km; z km aşağı negatif, abartılı)."""
    x, y = xy(cells.lon.to_numpy(), cells.lat.to_numpy())
    phi = np.radians(cells.strike.to_numpy())
    sx, sy = np.sin(phi), np.cos(phi)
    h = cells.L_km.to_numpy() / 2
    z0 = -cells.z_ust.to_numpy() * ABARTI
    z1 = -(cells.z_ust + cells.W_km).to_numpy() * ABARTI
    k = np.stack([np.stack([x - sx * h, y - sy * h, z0], 1), np.stack([x + sx * h, y + sy * h, z0], 1),
                  np.stack([x + sx * h, y + sy * h, z1], 1), np.stack([x - sx * h, y - sy * h, z1], 1)], 1)
    return k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--durum", default=str(ROOT / "runs/mekanik/isinma_kalin.pkl"))
    ap.add_argument("--senaryo", default="kalin")
    ap.add_argument("--etiket", default="800 yıllık ısınma sonu, örnek model durumu")
    a = ap.parse_args()
    import mekanik_model as mm
    v = deney.Veri(ROOT)
    pj = json.loads((mm.RUN / "secilen_parametreler.json").read_text())
    p = mm.Parametreler(**{k: val for k, val in pj.items() if k in mm.Parametreler.__dataclass_fields__
                           and k not in ("isi_akisi", "v_pl", "extra")})
    m, _ = mm.model_senaryo(v, p, a.senaryo)
    st = pickle.loads(Path(a.durum).read_bytes())
    st = st[-1] if isinstance(st, list) else st
    m.yukle(st) if "slip" in st and hasattr(m, "yukle") else None
    c = v.cells
    sinif = np.where(m.viskoz, 2, np.where(m.vs, 1, 0))  # 0 kilitli, 1 yavaş kayan, 2 sünek
    dol = np.where(m.vw, np.clip(m.S / np.maximum(m.S_fail, 1), -0.1, 1.0), np.nan)
    K = hucre_kose(c)
    kiyi = kiyilar()

    # ---- etkileşimli sayfa için veri
    out = ROOT / "reports" / "mekanik_3b"
    out.mkdir(parents=True, exist_ok=True)
    kx = [[round(float(t), 2) for t in xy(k[:, 0], k[:, 1])[0]] for k in kiyi]
    ky = [[round(float(t), 2) for t in xy(k[:, 0], k[:, 1])[1]] for k in kiyi]
    ana = c[(c.z_ust == 0) & (c.katman == "kabuk")].sort_values("x_km")
    et = []
    for ad, lo in ETIKET:
        r = ana.iloc[np.argmin(abs(ana.lon - lo))]
        ex, ey = xy([r.lon], [r.lat])
        et.append({"ad": ad, "x": round(float(ex[0]), 1), "y": round(float(ey[0]), 1)})
    sh = []
    for ad, lo, la in SEHIR:
        ex, ey = xy([lo], [la])
        sh.append({"ad": ad, "x": round(float(ex[0]), 1), "y": round(float(ey[0]), 1)})
    veri = {"kose": np.round(K, 2).reshape(len(c), -1).tolist(), "sinif": sinif.tolist(),
            "dol": [None if not np.isfinite(d) else round(float(d), 3) for d in dol],
            "x_km": np.round(c.x_km.to_numpy(), 1).tolist(), "z": c.z_ust.tolist(),
            "seg": c.seg_ad.tolist(), "kiyi_x": kx, "kiyi_y": ky, "etiket": et, "sehir": sh,
            "abarti": ABARTI, "durum": a.etiket}
    (out / "veri.json").write_text(json.dumps(veri, ensure_ascii=False, separators=(",", ":")))

    # ---- PNG
    cmap = matplotlib.colormaps["YlOrRd"]
    renk = np.where(sinif[:, None] == 2, np.array([[0.72, 0.72, 0.75, 0.35]]),
                    np.where(sinif[:, None] == 1, np.array([[0.35, 0.65, 0.75, 0.9]]),
                             cmap(np.nan_to_num(np.clip(dol, 0, 1)))))
    fig = plt.figure(figsize=(13, 7.5), dpi=150)
    ax = fig.add_subplot(111, projection="3d")
    ax.add_collection3d(Poly3DCollection(K, facecolors=renk, edgecolors=(0, 0, 0, 0.04), linewidths=0.1))
    seg = []
    for k in kiyi:
        kx_, ky_ = xy(k[:, 0], k[:, 1])
        pts = np.stack([kx_, ky_, np.zeros_like(kx_)], 1)
        seg += [pts[i:i + 2] for i in range(len(pts) - 1)]
    ax.add_collection3d(Line3DCollection(seg, colors="#3b5b7a", linewidths=0.6))
    for e in et:
        ax.text(e["x"], e["y"], 4, e["ad"], fontsize=7.5, ha="center", color="#222")
    for s in sh:
        ax.scatter([s["x"]], [s["y"]], [0], s=8, color="#111")
        ax.text(s["x"], s["y"], 1.5, s["ad"], fontsize=7, color="#444")
    x0, y0 = xy([KUTU[0]], [KUTU[2]])
    x1, y1 = xy([KUTU[1]], [KUTU[3]])
    ax.set_xlim(x0[0], x1[0])
    ax.set_ylim(y0[0], y1[0])
    ax.set_zlim(-40 * ABARTI, 8)
    ax.set_box_aspect((x1[0] - x0[0], y1[0] - y0[0], 48 * ABARTI))
    ax.view_init(elev=32, azim=-60)
    ax.set_zticks([0, -10 * ABARTI, -20 * ABARTI, -30 * ABARTI, -40 * ABARTI])
    ax.set_zticklabels(["0", "10", "20", "30", "40 km"], fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"Marmara fay modeli, 3B: 1 km × 1 km hücreler, renk = doluluk ({a.etiket})", fontsize=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.02)
    cb.set_label("doluluk: birikmiş gerilme / kırılma sınırı (kilitli hücreler)")
    fig.text(0.02, 0.04, "Mavi: yavaşça kayan (krip) hücreler   Gri: sünek akan kuşak (20-40 km)   "
             f"Derinlik {ABARTI:g} kat abartılı", fontsize=8, color="#333")
    fig.tight_layout()
    fig.savefig(ROOT / "reports/figures/mekanik_3b.png", bbox_inches="tight", pad_inches=0.1)
    print("ok", len(c), np.nanmean(dol))


if __name__ == "__main__":
    main()
