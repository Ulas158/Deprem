"""MTA Türkiye Diri Fay Haritası 2026'dan Kuzey Anadolu Fayı kuzey kolu izi.

Kaynak: Elmacı vd. 2026, Türkiye Diri Fay Haritası, MTA Özel Yayın Serisi-60 (ISBN 978-625-8107-74-6).
Ham veri depoya girmez (data/mta/ altında, .gitignore); burada yalnız türetilmiş iz kurulur.

Her MTA kesimi birçok küçük parçadan (paralel kollar dahil) oluşur. Ana kolun tek izi, bütün
kesimlerin parçaları 0,5 km aralıkla sıklaştırılıp batıdan doğuya boylam kutularında (0,01°)
süreklilik izlenerek çıkarılır (ana_kol_izleri). Kademeli geçişler tek çizgiye iner; 2 km hücre
ağında bu ihmal edilebilir.

Kesimler, eski GEM tabanlı ağın 4 ana kol parçasıyla (seg 0-3) aynı gruplara toplanır ki
levha hızı ve sığ krip gibi seg/boylam tabanlı ayarlar değişmeden çalışsın. Yan kollar (seg 4-8)
şimdilik GEM'den kalır.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

KM_PER_DEG = 111.195

# (MTA SEGMENT_ADI, bu boylamdan batıdaki noktaları at, bu boylamdan doğudakileri at, grup)
ANA_KOL = [
    ("Saros Fayı", 26.20, None, 0),
    ("Ganos Fayı", None, None, 0),
    ("Tekirdağ Fayı", None, None, 0),
    ("Kumburgaz Fayı", None, None, 1),
    ("Avcılar Fayı", None, None, 1),
    ("Adalar Fayı", None, None, 1),
    ("Darıca Fayı", None, None, 2),
    ("Karamürsel Fayı", None, None, 2),
    ("Gölcük Fayı", None, None, 2),
    ("Tepetarla Fayı", None, None, 2),
    ("Arifiye Fayı", None, None, 2),
    ("Karadere Fayı", None, 31.00, 2),  # KD uzantısı Düzce birleşiminin doğusunda atılır
    ("Düzce Fayı", 31.00, None, 3),
]
# Orta (güney) kol: GEM'deki Biga-Bandırma (seg 4), Bandırma-Gemlik (5), İznik-Geyve (6) gruplarıyla
ORTA_KOL = [
    ("Sinekçi Fayı", None, None, 4),
    ("Edincik Fayı", None, None, 4),
    ("Bandırma Fayı", None, None, 5),
    ("Zeytinbağı Fayı", None, None, 5),
    ("Gençali Fayı", None, None, 5),
    ("Gemlik Fayı", None, None, 6),
    ("İznik - Mekece Fayı", None, None, 6),
    ("Geyve Fayı", None, None, 6),
]
GRUP_AD = {0: "Ganos-Saros", 1: "Orta Marmara-Çınarcık", 2: "İzmit-Sapanca-Akyazı", 3: "Düzce",
           4: "Biga-Bandırma", 5: "Bandırma-Gemlik", 6: "İznik-Geyve"}


def _parcalar(geojson: Path, segment_adi: str):
    g = json.loads(Path(geojson).read_text())
    out = []
    for f in g["features"]:
        p, gm = f["properties"], f["geometry"]
        if gm is None or p.get("SEGMENT_ADI") != segment_adi or "KUZEY ANADOLU" not in (p.get("FAY_ADI") or ""):
            continue
        ls = [gm["coordinates"]] if gm["type"] == "LineString" else gm["coordinates"]
        out += [np.asarray(c, float)[:, :2] for c in ls]
    return out


def _sik(c, ds_km=0.5):
    kx = KM_PER_DEG * np.cos(np.radians(40.7))
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(c[:, 0]) * kx, np.diff(c[:, 1]) * KM_PER_DEG))])
    if s[-1] == 0:
        return c
    t = np.linspace(0, s[-1], max(int(s[-1] / ds_km), 1) + 1)
    return np.column_stack([np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])])


def kesim_izi(geojson: Path, segment_adi: str, dlon=0.01):
    P = np.vstack([_sik(c) for c in _parcalar(geojson, segment_adi)])
    b = np.floor(P[:, 0] / dlon).astype(int)
    kut = np.unique(b)
    lon = (kut + 0.5) * dlon
    lat = np.array([np.median(P[b == k, 1]) for k in kut])
    if len(lat) >= 3:
        lat = np.array([np.median(lat[max(i - 1, 0):i + 2]) for i in range(len(lat))])
    # uçlar: gerçek en batı/doğu nokta
    return np.column_stack([np.r_[P[:, 0].min(), lon, P[:, 0].max()],
                            np.r_[P[np.argmin(P[:, 0]), 1], lat, P[np.argmax(P[:, 0]), 1]]])


def ana_kol_izleri(geojson: Path, dlon=0.01, lon_bas=26.20, liste=None):
    """Batıdan doğuya tek sürekli iz. Bütün ana kol kesimlerinin noktaları birlikte izlenir:
    her boylam kutusunda bir önceki kutulardan doğrusal uzatmaya en yakın nokta seçilir (paralel
    kollar arasında atlamamak için), sonra 5 kutuluk medyan + 3 kutuluk ortalamayla pürüzsüzleştirilir.
    Her kutunun kesim adı, seçilen noktanın kesimidir (7 kutuluk çoğunlukla düzeltilir)."""
    liste = ANA_KOL if liste is None else liste
    P, E = [], []
    for k, (ad, l0, l1, _) in enumerate(liste):
        for c in _parcalar(geojson, ad):
            q = _sik(c)
            q = q[(q[:, 0] >= (l0 if l0 is not None else -999)) & (q[:, 0] <= (l1 if l1 is not None else 999))]
            if len(q) == 0:
                continue
            P.append(q)
            E.append(np.full(len(q), k))
    P, E = np.vstack(P), np.concatenate(E)
    keep = P[:, 0] >= lon_bas
    P, E = P[keep], E[keep]
    b = np.floor(P[:, 0] / dlon).astype(int)
    lon, lat, etk = [], [], []
    for kb in np.unique(b):
        ii = np.flatnonzero(b == kb)
        if len(lat) >= 3:
            egim = (lat[-1] - lat[-3]) / (lon[-1] - lon[-3])
            tah = lat[-1] + egim * ((kb + 0.5) * dlon - lon[-1])
        elif lat:
            tah = lat[-1]
        else:
            # başlangıç: en batıdaki kesimin (Saros) medyanı
            tah = np.median(P[ii, 1])
        j = ii[np.argmin(np.abs(P[ii, 1] - tah))]
        # kutu içindeki aynı koldan noktaların medyanı (seçilen noktaya 1 km'den yakın)
        yak = ii[np.abs(P[ii, 1] - P[j, 1]) < 1.0 / KM_PER_DEG]
        lon.append((kb + 0.5) * dlon)
        lat.append(float(np.median(P[yak, 1])))
        etk.append(int(np.bincount(E[yak]).argmax()))
    lon, lat, etk = np.array(lon), np.array(lat), np.array(etk)
    n = len(lat)
    lat = np.array([np.median(lat[max(i - 2, 0):i + 3]) for i in range(n)])
    lat = np.array([np.mean(lat[max(i - 1, 0):i + 2]) for i in range(n)])
    etk = np.array([np.bincount(etk[max(i - 3, 0):i + 4]).argmax() for i in range(n)])
    etk = np.maximum.accumulate(etk)  # batıdan doğuya sıra korunur
    out = []
    for k in np.unique(etk):
        ii = np.flatnonzero(etk == k)
        i0 = max(ii[0] - 1, 0)  # bir önceki kesimin son noktasıyla bağlan
        ad, _, _, gr = liste[k]
        out.append({"ad": ad, "coords": np.column_stack([lon[i0:ii[-1] + 1], lat[i0:ii[-1] + 1]]), "grup": gr})
    return out


def fay_izleri_mta(mta_geojson: Path, gem_path: Path) -> tuple[list[dict], list[tuple]]:
    """geometri.hucreler için iz listesi, eski GEM ağıyla aynı seg numaraları:
    ana kol seg 0-3 (MTA), orta kol seg 4-6 (MTA), İmralı seg 7 ve Armutlu-Çınarcık güney seg 8 (GEM;
    MTA haritasında deniz içindeki bu iki kol yok). İkinci dönen: MTA kesim sınırları
    [(kol, ad, lon0, lon1)]."""
    from . import geometri as g
    gem = {iz["ad"]: iz for iz in g.fay_izleri(gem_path, yan_kollar=True)}
    izler, sinir = [], []
    for kol, liste, gruplar, lon_bas, x_bas in [(0, ANA_KOL, range(4), 26.20, None),
                                                (1, ORTA_KOL, range(4, 7), -999.0, 600.0)]:
        kes = ana_kol_izleri(mta_geojson, lon_bas=lon_bas, liste=liste)
        for j, gr in enumerate(gruplar):
            parca = [k for k in kes if k["grup"] == gr]
            c = np.vstack([parca[0]["coords"]] + [p["coords"][1:] for p in parca[1:]])
            ref = gem.get(GRUP_AD[gr], {})
            izler.append({"fault_id": f"MTA_{gr}", "ad": GRUP_AD[gr], "coords": c, "kol": kol,
                          "x_bas": x_bas if j == 0 else None, "gem_slip_rate": ref.get("gem_slip_rate", np.nan)})
            sinir += [(kol, p["ad"], float(p["coords"][0, 0]), float(p["coords"][-1, 0])) for p in parca]
    izler += [gem["İmralı"], gem["Armutlu-Çınarcık güney"]]
    return izler, sinir
