"""Fay fay, depremden depreme öğrenen model için veri hazırlığı.

Bütün diziler bir kez hesaplanır ve ``faultwalk.py`` içindeki model tarafından
kullanılır. Hiçbir dizi, ait olduğu andan sonraki bir veriyi kullanmaz:

* yıl haritası S[y]: yalnız y yılının 1 Ocak'ından önceki depremlerden
* fay doluluğu L[s, u]: yalnız s adımının başından önceki büyük depremlerden
* Coulomb işareti s[j, c]: yalnız j depreminin kendi kırılmasından

Birimler
--------
* Uzay: 0,5° hücreler (önceki çalışmalarla aynı). Her olay bir hücreye düşer.
* Fay birimi: her GEM/EMME fay segmenti ayrı bir birimdir ("her fay ayrı").
  Bir hücre, 25 km içindeki en yüksek moment birikim hızlı segmente bağlanır;
  yakınında fay olmayan hücreler bölgesel "fay dışı" birimlere bağlanır.
* Büyük fay sistemleri: uçları 12 km içinde birleşen segmentlerden oluşan ve
  toplam moment birikim hızı yüksek bağlı bileşenler (Kuzey Anadolu, Doğu
  Anadolu, ...). Diğer segmentler coğrafi bölge sistemlerine düşer.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .faults import gem
from .faults.moment import add_geometry, mw_to_m0
from .geo import point_to_polyline_km
from .region import TURKEY

log = logging.getLogger(__name__)

T0 = pd.Timestamp("1700-01-01", tz="UTC")
M0 = 5.5
CELL = 0.5
ERAS = [(1700, 1800), (1800, 1850), (1850, 1900)]   # tespit oranı öğrenilen eksik dönemler; 1900+ = 1
SUB = 3                                             # hücre içi alt örnekleme (3x3)


def days(t) -> np.ndarray:
    return ((pd.to_datetime(t, utc=True) - T0) / pd.Timedelta(days=1)).to_numpy(float) \
        if not isinstance(t, pd.Timestamp) else (t - T0) / pd.Timedelta(days=1)


def year_start_days(y: int) -> float:
    return (pd.Timestamp(f"{y}-01-01", tz="UTC") - T0) / pd.Timedelta(days=1)


def xy_km(lat, lon):
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    return (lon - 35.0) * 111.195 * np.cos(np.radians(39.0)), (lat - 39.0) * 111.195


# ------------------------------------------------------------------ büyük fay sistemleri
# yaklaşık fay izleri (boylam, enlem); büyük sistemlerin koridorlarını tanımlamak için
NAF_TRACE = [(26.0, 40.6), (27.5, 40.8), (29.0, 40.75), (30.3, 40.7), (31.2, 40.8), (32.5, 40.9), (34.0, 41.0),
             (35.5, 40.9), (36.5, 40.6), (38.0, 40.1), (39.5, 39.75), (40.7, 39.4)]
EAF_TRACE = [(40.7, 39.4), (39.5, 38.6), (38.6, 38.3), (37.5, 37.9), (36.9, 37.6), (36.5, 37.0), (36.2, 36.4)]


def fault_systems(f: pd.DataFrame) -> pd.Series:
    """Segmentleri büyük fay sistemlerine ayır.

    Kuzey Anadolu ve Doğu Anadolu: segment orta noktası yaklaşık ana ize 25 km içinde ve
    kayma hızı ≥ 2 mm/yıl. Ege graben sistemi: Batı Anadolu (boylam < 30,5; enlem 36,8-39,6)
    normal faylar. Geri kalanlar coğrafi "diğer" sistemlere düşer.
    """
    out = []
    for lat, lon, sr, st in zip(f["mid_lat"], f["mid_lon"], f["slip_rate_mm_yr"].fillna(0), f["slip_type"].fillna("")):
        dn = point_to_polyline_km([lat], [lon], NAF_TRACE)[0]
        de = point_to_polyline_km([lat], [lon], EAF_TRACE)[0]
        if dn <= 25 and sr >= 2 and dn <= de:
            out.append("Kuzey Anadolu Fayı")
        elif de <= 25 and sr >= 2:
            out.append("Doğu Anadolu Fayı")
        elif lon < 30.5 and 36.8 <= lat < 39.6 and "Normal" in st:
            out.append("Ege graben sistemi")
        else:
            out.append(_region_name(lat, lon))
    return pd.Series(out, index=f.index, name="system")


def _region_name(lat, lon) -> str:
    if lon < 30.5:
        return "diğer: Batı Anadolu-Ege"
    if lon < 36.0:
        return "diğer: Orta Anadolu"
    if lat >= 39.6 and lon < 41.5:
        return "diğer: Kuzeydoğu"
    if lon >= 41.0:
        return "diğer: Doğu"
    return "diğer: Güneydoğu"


# ------------------------------------------------------------------ katalog
def load_catalog(root: Path, region=TURKEY) -> tuple[pd.DataFrame, pd.DataFrame]:
    """1700'den bugüne M≥5,5 olaylar (tarihsel + aletsel) ve 1700 öncesi arka plan."""
    p = Path(root) / "data" / "processed" / "turkey_catalog_long.csv.gz"
    cat = pd.read_csv(p)
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    mcol = "mw_best" if "mw_best" in cat else "mag"
    cat["mw"] = cat[mcol].fillna(cat["mag"])
    cat = cat[(cat["time"] >= T0) & (cat["mw"] >= M0) & region.contains(cat)].copy()
    cat = cat.sort_values("time").reset_index(drop=True)
    # aynı ana düşen olaylar: 1 sn ayır
    t = np.array(days(cat["time"]), dtype=float)
    for i in range(1, len(t)):
        if t[i] <= t[i - 1]:
            t[i] = t[i - 1] + 1e-5
    cat["t"] = t
    bgp = p.with_name("historical_background.csv")
    old = pd.read_csv(bgp) if bgp.exists() else pd.DataFrame(columns=["latitude", "longitude", "mag"])
    if len(old):
        old = old[(old["mag"] >= M0 - 0.5) & region.contains(old)]
    return cat, old


# ------------------------------------------------------------------ ana veri
@dataclass
class WalkData:
    cat: pd.DataFrame
    lat_c: np.ndarray          # hücre merkezleri
    lon_c: np.ndarray
    area: np.ndarray           # km²
    unit_of_cell: np.ndarray   # hücre -> birim
    units: pd.DataFrame        # birim tablosu (segment / fay dışı), system kodu
    systems: list
    major: np.ndarray          # sistem -> büyük fay mı
    ev_cell: np.ndarray
    ev_unit: np.ndarray
    ev_sys: np.ndarray
    t_ev: np.ndarray
    m_ev: np.ndarray
    d2: np.ndarray             # (N, C, SUB²) olay → hücre alt noktaları uzaklık² (km²)
    cfs: np.ndarray            # (N, C) Coulomb işareti tanh(ΔCFF / 0,1 bar)
    fmap: np.ndarray           # (C,) fay moment haritası (toplam 1)
    smap: np.ndarray           # (Y, C) yıl başına geçmiş sismisite haritası (toplam 1)
    years: np.ndarray          # smap yılları
    tau: np.ndarray            # adım sınırları (S+1,)
    step_year: np.ndarray      # adım -> smap satırı
    step_era: np.ndarray       # adım -> dönem (0..len(ERAS); son = tam)
    L: np.ndarray              # (S, U) log doluluk
    ev_step: np.ndarray        # olay i -> olaydan hemen önceki adım
    ev_era: np.ndarray
    ev_bstep: np.ndarray       # olay i'nin başlattığı adım (t_i ile başlayan)
    extra: dict = field(default_factory=dict)


def build(root: Path, end: str = "2026-01-01", cache: Path | None = None) -> WalkData:
    if cache is not None and Path(cache).exists():
        import pickle
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    root = Path(root)
    cat, old = load_catalog(root)
    cat = cat[cat["time"] < pd.Timestamp(end, tz="UTC")].reset_index(drop=True)
    log.info("Katalog: %d olay M≥%.1f (%s → %s)", len(cat), M0, cat["time"].min().date(), cat["time"].max().date())

    # hücreler
    lat_e = np.arange(TURKEY.min_lat, TURKEY.max_lat + 1e-9, CELL)
    lon_e = np.arange(TURKEY.min_lon, TURKEY.max_lon + 1e-9, CELL)
    nla, nlo = len(lat_e) - 1, len(lon_e) - 1
    LA, LO = np.meshgrid(lat_e[:-1] + CELL / 2, lon_e[:-1] + CELL / 2, indexing="ij")
    lat_c, lon_c = LA.ravel(), LO.ravel()
    C = lat_c.size
    area = (CELL * 111.195) ** 2 * np.cos(np.radians(lat_c))
    ii = np.clip(((cat["latitude"] - TURKEY.min_lat) // CELL).astype(int), 0, nla - 1)
    jj = np.clip(((cat["longitude"] - TURKEY.min_lon) // CELL).astype(int), 0, nlo - 1)
    ev_cell = (ii * nlo + jj).to_numpy()

    # faylar
    f = gem.load_faults(gem.download(root / "data" / "ext" / "gem_active_faults_harmonized.geojson"), TURKEY)
    f = add_geometry(f).reset_index(drop=True)
    f["slip_rate_mm_yr"] = f["slip_rate_mm_yr"].fillna(0.2)
    f["m0_rate_Nm_yr"] = f["m0_rate_Nm_yr"].fillna(f["m0_rate_Nm_yr"].median())
    f["mid_lon"] = [float(np.asarray(c)[:, 0].mean()) for c in f["coords"]]
    f["mid_lat"] = [float(np.asarray(c)[:, 1].mean()) for c in f["coords"]]
    f["system"] = fault_systems(f)
    dist = np.stack([point_to_polyline_km(lat_c, lon_c, c) for c in f["coords"]], 1)  # C × F

    # fay moment haritası: iz boyunca 2 km'lik noktalar, moment hızı uzunluğa paylaştırılır
    fmap = np.zeros(C)
    for c, m0r, L in zip(f["coords"], f["m0_rate_Nm_yr"], f["length_km"]):
        c = np.asarray(c, float)
        pts = []
        for a, b in zip(c[:-1], c[1:]):
            n = max(int(np.hypot((b[0] - a[0]) * 85, (b[1] - a[1]) * 111) / 2), 1)
            s = np.linspace(0, 1, n, endpoint=False)
            pts.append(a + s[:, None] * (b - a))
        pts = np.concatenate(pts)
        pi = np.clip(((pts[:, 1] - TURKEY.min_lat) // CELL).astype(int), -1, nla)
        pj = np.clip(((pts[:, 0] - TURKEY.min_lon) // CELL).astype(int), -1, nlo)
        ok = (pi >= 0) & (pi < nla) & (pj >= 0) & (pj < nlo)
        np.add.at(fmap, pi[ok] * nlo + pj[ok], m0r / len(pts))
    fmap = fmap / fmap.sum()
    fmap = 0.9 * fmap + 0.1 / C

    # birimler: hücrenin 25 km içindeki en yüksek moment hızlı segmenti
    m0r = f["m0_rate_Nm_yr"].to_numpy()
    score = np.where(dist <= 25.0, m0r[None, :], -1.0)
    best = score.argmax(1)
    has = score.max(1) > 0
    seg_used = np.unique(best[has])
    unit_rows = []
    seg_to_unit = {}
    for k in seg_used:
        seg_to_unit[k] = len(unit_rows)
        unit_rows.append({"unit": f["fault_id"][k], "kind": "fay", "system": f["system"][k],
                          "seg": k, "slip_rate": f["slip_rate_mm_yr"][k], "m0_rate": m0r[k],
                          "mmax": f["mmax_area"][k], "area_km2": f["area_km2"][k],
                          "lat": f["mid_lat"][k], "lon": f["mid_lon"][k], "length_km": f["length_km"][k]})
    off_names = sorted({_region_name(a, b) for a, b in zip(lat_c[~has], lon_c[~has])})
    off_to_unit = {}
    for nm in off_names:
        off_to_unit[nm] = len(unit_rows)
        unit_rows.append({"unit": "fay dışı: " + nm.split(": ")[1], "kind": "fay dışı", "system": nm,
                          "seg": -1, "slip_rate": 0.0, "m0_rate": 0.0, "mmax": np.nan, "area_km2": np.nan,
                          "lat": np.nan, "lon": np.nan, "length_km": 0.0})
    units = pd.DataFrame(unit_rows)
    unit_of_cell = np.array([seg_to_unit[best[c]] if has[c] else off_to_unit[_region_name(lat_c[c], lon_c[c])]
                             for c in range(C)])
    systems = sorted(units["system"].unique(), key=lambda s: (s.startswith("diğer"), s))
    sys_code = {s: i for i, s in enumerate(systems)}
    units["sys"] = units["system"].map(sys_code)
    major = np.array([not s.startswith("diğer") for s in systems])
    ev_unit = unit_of_cell[ev_cell]
    ev_sys = units["sys"].to_numpy()[ev_unit]
    log.info("%d fay birimi, %d fay dışı birim; sistemler: %s", (units["kind"] == "fay").sum(),
             (units["kind"] != "fay").sum(), ", ".join(systems))

    # olay → hücre alt noktaları uzaklıkları
    ex, ey = xy_km(cat["latitude"], cat["longitude"])
    off = (np.arange(SUB) + 0.5) / SUB - 0.5
    sla = (lat_c[:, None, None] + off[:, None] * CELL + 0 * off[None, :]).reshape(C, -1)
    slo = (lon_c[:, None, None] + 0 * off[:, None] + off[None, :] * CELL).reshape(C, -1)
    cx, cy = xy_km(sla, slo)
    d2 = ((ex[:, None, None] - cx[None]) ** 2 + (ey[:, None, None] - cy[None]) ** 2).astype(np.float32)

    # Coulomb: her olayın, her hücrenin bağlı olduğu (veya en yakın) fay geometrisinde ΔCFF'si
    cfs = coulomb_signs(cat, f, lat_c, lon_c, dist)

    # yıl haritaları (yalnız o yıl öncesi)
    years = np.arange(1700, 2027)
    smap = np.zeros((years.size, C))
    ox, oy = xy_km(old["latitude"], old["longitude"]) if len(old) else (np.zeros(0), np.zeros(0))
    ccx, ccy = xy_km(lat_c, lon_c)
    sig = 30.0

    def kern(px, py, w):
        if len(px) == 0:
            return np.zeros(C)
        k = np.exp(-((px[:, None] - ccx) ** 2 + (py[:, None] - ccy) ** 2) / (2 * sig ** 2))
        return (w[:, None] * k).sum(0)

    base_old = kern(ox, oy, np.full(len(ox), 0.3))
    eyears = cat["time"].dt.year.to_numpy()
    kev = np.exp(-((ex[:, None] - ccx) ** 2 + (ey[:, None] - ccy) ** 2) / (2 * sig ** 2))
    cum = np.zeros(C)
    k = 0
    for yi, y in enumerate(years):
        while k < len(cat) and eyears[k] < y:
            cum += kev[k]
            k += 1
        m = base_old + cum
        m = m / m.sum() if m.sum() > 0 else np.full(C, 1 / C)
        smap[yi] = 0.95 * m + 0.05 / C

    # zaman adımları: yıl başları + olay anları
    t_ev = cat["t"].to_numpy()
    yb = np.array([year_start_days(y) for y in years])
    tau = np.unique(np.concatenate([yb, t_ev, [year_start_days(2026)]]))
    tau = tau[(tau >= 0) & (tau <= year_start_days(2026))]
    S = tau.size - 1
    step_year = np.searchsorted(yb, tau[:-1], side="right") - 1
    step_era = np.full(S, len(ERAS))
    for e, (a, b) in enumerate(ERAS):
        step_era[(tau[:-1] >= year_start_days(a)) & (tau[:-1] < year_start_days(b))] = e
    ev_bstep = np.searchsorted(tau, t_ev)          # tau[ev_bstep] == t_ev
    ev_step = ev_bstep - 1                          # olaydan hemen önceki adım
    ev_era = step_era[np.clip(ev_step, 0, S - 1)]

    L = fill_matrix(cat, f, units, tau)
    W = WalkData(cat=cat, lat_c=lat_c, lon_c=lon_c, area=area, unit_of_cell=unit_of_cell, units=units,
                 systems=systems, major=major, ev_cell=ev_cell, ev_unit=ev_unit, ev_sys=ev_sys,
                 t_ev=t_ev, m_ev=cat["mw"].to_numpy(), d2=d2, cfs=cfs, fmap=fmap, smap=smap, years=years,
                 tau=tau, step_year=step_year, step_era=step_era, L=L, ev_step=ev_step, ev_era=ev_era,
                 ev_bstep=ev_bstep)
    W.extra["faults"] = f.drop(columns=["coords"])
    if cache is not None:
        import pickle
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "wb") as fh:
            pickle.dump(W, fh)
    return W


def rupture_len_km(mw):
    return 10 ** (-2.44 + 0.59 * np.asarray(mw, float))


def fill_matrix(cat, f, units, tau) -> np.ndarray:
    """Her adım başında her fay biriminin log doluluğu.

    doluluk = kayma hızı × son kırılmadan beri geçen süre / karakteristik kayma.
    Kırılma: M≥6 bir olay, segmente (kırık yarı boyu + 10 km) içindeyse segmenti sıfırlar.
    1700'den beri kırılmamış segmentlerde saat 1700'de başlar (alt sınır; raporda belirtilir).
    """
    U = len(units)
    S = tau.size - 1
    L = np.zeros((S, U))
    big = cat[cat["mw"] >= 6.0]
    fu = units[units["kind"] == "fay"]
    for u, r in fu.iterrows():
        c = f["coords"][r["seg"]]
        d = point_to_polyline_km(big["latitude"].to_numpy(), big["longitude"].to_numpy(), c)
        hit = d <= rupture_len_km(big["mw"].to_numpy()) / 2 + 10.0
        rt = np.sort(big["t"].to_numpy()[hit])
        # her adım başı için son kırılma (adım başından önce)
        idx = np.searchsorted(rt, tau[:-1], side="left") - 1
        last = np.where(idx >= 0, rt[np.clip(idx, 0, None)], 0.0) if rt.size else np.zeros(S)
        el_yr = (tau[:-1] - last) / 365.25
        mmax = r["mmax"] if np.isfinite(r["mmax"]) else 7.0
        A = r["area_km2"] if np.isfinite(r["area_km2"]) and r["area_km2"] > 0 else 300.0
        d_char = mw_to_m0(mmax) / (3e10 * A * 1e6)               # m
        fill = r["slip_rate"] * 1e-3 * np.maximum(el_yr, 0.5) / max(d_char, 0.1)
        L[:, u] = np.clip(np.log(fill), -5.0, 2.0)
    return L


def coulomb_signs(cat, f, lat_c, lon_c, dist, mu_eff=0.4, scale_bar=0.1) -> np.ndarray:
    """(N, C): tanh(ΔCFF / 0,1 bar). Alıcı: hücrenin 60 km içindeki en yüksek moment hızlı
    segmentin geometrisi (yoksa en yakın segment), 10 km derinlikte."""
    from .faults.coulomb import build_sources, coulomb_on_receivers, stress_from_sources, to_local_m

    C = lat_c.size
    m0r = f["m0_rate_Nm_yr"].to_numpy()
    score = np.where(dist <= 60.0, m0r[None, :], -1.0)
    rec = np.where(score.max(1) > 0, score.argmax(1), dist.argmin(1))
    strike = np.zeros(C)
    for c in range(C):
        co = np.asarray(f["coords"][rec[c]], float)
        s = np.argmin(np.hypot(co[:-1, 0] - lon_c[c], co[:-1, 1] - lat_c[c])) if len(co) > 1 else 0
        a, b = co[s], co[min(s + 1, len(co) - 1)]
        dx = (b[0] - a[0]) * np.cos(np.radians(lat_c[c]))
        strike[c] = (np.degrees(np.arctan2(dx, b[1] - a[1])) + 360) % 360
    dip = f["dip"].to_numpy()[rec]
    rake = np.nan_to_num(f["rake"].to_numpy()[rec], nan=180.0)
    obs = to_local_m(lat_c, lon_c, np.full(C, 10.0))
    ev = cat.assign(mw_best=cat["mw"], depth_km=cat.get("depth_km", pd.Series(10.0, index=cat.index)).fillna(10.0))
    src = build_sources(ev, f, mw_col="mw_best", min_mw=M0, max_dist_km=30.0)
    out = np.zeros((len(cat), C), np.float32)
    # build_sources olay sırasını korur (hepsi M≥5,5)
    for j in range(len(src)):
        st = stress_from_sources(obs, src.iloc[[j]])
        dcff, _, _ = coulomb_on_receivers(st, strike, dip, rake, mu_eff)
        out[j] = np.tanh(dcff / 1e5 / scale_bar)
    log.info("Coulomb işaretleri: %d kaynak × %d hücre", len(src), C)
    return out
