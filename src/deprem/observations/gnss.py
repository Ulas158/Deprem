"""GNSS (GPS) zaman serileri: hız alanı, gerinim hızı ve geçici (transient) sinyaller.

Kaynak: Nevada Geodetic Laboratory (NGL), Blewitt vd. 2018, Eos 99.
Günlük çözümler, Avrasya sabit çerçevesinde (``IGS20/tenv3/EU``).
Ekipman değişimi ve deprem adımları NGL ``steps.txt`` dosyasından okunur.

Fiziksel gerekçe
----------------
* **Gerinim hızı:** Kabukta biriken elastik gerinim, depremin enerji
  kaynağıdır. Yüksek gerinim hızlı bölgelerde uzun dönem deprem oranı
  daha yüksektir (Bird & Kreemer 2015, BSSA; Kreemer vd. 2014 GSRM).
* **Geçici sinyal:** Yavaş kayma ve büyük depremlerden önceki hızlanmış
  yer hareketi bazı olaylarda gözlenmiştir (Tohoku 2011, Iquique 2014;
  Ruiz vd. 2014, Science). Türkiye için kanıt sınırlıdır; bu özellik bu
  yüzden yalnız geçmiş veriyle hesaplanır ve ETAS'a karşı test edilir.

Sızıntı önlemi: ``t0`` yılındaki özellik yalnızca ``t0``'dan önceki günlük
konumlardan hesaplanır.
"""

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ..region import Region
from .net import download

log = logging.getLogger(__name__)

NGL = "https://geodesy.unr.edu"
HOLDINGS_URL = f"{NGL}/NGLStationPages/DataHoldings.txt"
STEPS_URL = f"{NGL}/NGLStationPages/steps.txt"
TENV3_URL = NGL + "/gps_timeseries/IGS20/tenv3/EU/{site}.EU.tenv3"


# ------------------------------------------------------------------ indirme
def station_list(raw: Path, region: Region, pad_deg: float = 1.0, min_years: float = 2.5) -> pd.DataFrame:
    p = download(HOLDINGS_URL, raw / "DataHoldings.txt")
    rows = []
    for line in p.read_text(errors="ignore").splitlines()[1:]:
        f = line.split()
        if len(f) < 11:
            continue
        try:
            lat, lon = float(f[1]), float(f[2])
        except ValueError:
            continue
        lon = (lon + 180) % 360 - 180
        rows.append({"site": f[0], "lat": lat, "lon": lon, "t_begin": f[7], "t_end": f[8],
                     "n_sol": int(f[10]) if f[10].isdigit() else 0})
    df = pd.DataFrame(rows)
    df = df[df["lat"].between(region.min_lat - pad_deg, region.max_lat + pad_deg)
            & df["lon"].between(region.min_lon - pad_deg, region.max_lon + pad_deg)]
    df["t_begin"] = pd.to_datetime(df["t_begin"], errors="coerce")
    df["t_end"] = pd.to_datetime(df["t_end"], errors="coerce")
    span = (df["t_end"] - df["t_begin"]).dt.days / 365.25
    return df[span >= min_years].reset_index(drop=True)


def fetch_all(raw: Path, sites: list[str], workers: int = 6) -> dict[str, Path]:
    download(STEPS_URL, raw / "steps.txt")
    out = {}

    def one(s):
        return s, download(TENV3_URL.format(site=s), raw / "tenv3" / f"{s}.EU.tenv3")

    with ThreadPoolExecutor(workers) as ex:
        for s, p in ex.map(one, sites):
            if p is not None:
                out[s] = p
    log.info("GNSS: %d/%d istasyon indirildi", len(out), len(sites))
    return out


def read_tenv3(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=0)
    df.columns = ["site", "date", "t", "mjd", "week", "d", "reflon", "e0", "e", "n0", "n", "u0", "u",
                  "ant", "se", "sn", "su", "cen", "ceu", "cnu", "lat", "lon", "h"][:df.shape[1]]
    out = pd.DataFrame({
        "t": df["t"].astype(float),
        "e": (df["e0"] + df["e"]) * 1000.0,  # mm
        "n": (df["n0"] + df["n"]) * 1000.0,
        "u": (df["u0"] + df["u"]) * 1000.0,
        "se": df["se"] * 1000.0, "sn": df["sn"] * 1000.0, "su": df["su"] * 1000.0,
    })
    # birinci gözleme göre göreli (sayısal kararlılık)
    for c in ("e", "n", "u"):
        out[c] -= out[c].iloc[0]
    return out


def read_steps(path: Path) -> pd.DataFrame:
    rows = []
    for line in Path(path).read_text(errors="ignore").splitlines():
        f = line.split()
        if len(f) < 3:
            continue
        try:
            t = pd.to_datetime(f[1], format="%y%b%d")
        except ValueError:
            continue
        kind = int(f[2]) if f[2].isdigit() else 0
        mag = float(f[5]) if kind == 2 and len(f) > 5 else np.nan
        rows.append({"site": f[0], "time": t, "kind": kind, "mag": mag})
    df = pd.DataFrame(rows)
    df["t"] = df["time"].dt.year + (df["time"].dt.dayofyear - 0.5) / 365.25
    return df


# ------------------------------------------------------------------ yörünge modeli
def _design(t, t_ref, steps):
    cols = [np.ones_like(t), t - t_ref,
            np.sin(2 * np.pi * t), np.cos(2 * np.pi * t),
            np.sin(4 * np.pi * t), np.cos(4 * np.pi * t)]
    for ts in steps:
        cols.append((t >= ts).astype(float))
    return np.column_stack(cols)


def _robust_fit(t, y, s, steps, t_ref, n_iter: int = 3):
    keep = np.isfinite(y)
    for _ in range(n_iter):
        A = _design(t[keep], t_ref, steps)
        w = 1 / np.clip(s[keep], 0.5, None)
        beta, *_ = np.linalg.lstsq(A * w[:, None], y[keep] * w, rcond=None)
        res = y - _design(t, t_ref, steps) @ beta
        mad = 1.4826 * np.nanmedian(np.abs(res[keep] - np.nanmedian(res[keep])))
        keep = np.isfinite(y) & (np.abs(res) < 5 * max(mad, 0.5))
    return beta, res, keep


def station_year_features(ts: pd.DataFrame, steps_t: np.ndarray, t0: float,
                          long_years: float = 5.0, recent_years: float = 1.0,
                          min_span: float = 2.5) -> dict | None:
    """``t0`` (ondalık yıl) için tek istasyon özellikleri; yalnız t<t0 verisi.

    * ``ve, vn, vu``: [t0-5, t0) aralığında yörünge modelinden hız (mm/yıl)
    * ``tr_e, tr_n``: son 1 yılın hızı − önceki dört yılın hızı (mm/yıl)
    * ``tr_z``: yatay geçici sinyalin gürültüye oranı (beyaz+renkli gürültü
      için kaba ölçek: ortalama artığın standart hatası × 3)
    """
    w = ts[(ts["t"] >= t0 - long_years) & (ts["t"] < t0)]
    if len(w) < 300 or w["t"].iloc[-1] - w["t"].iloc[0] < min_span or w["t"].iloc[-1] < t0 - 0.25:
        return None
    t = w["t"].to_numpy()
    st = [s for s in steps_t if t[0] < s < t[-1]]
    out = {}
    for c in ("e", "n", "u"):
        beta, _, _ = _robust_fit(t, w[c].to_numpy(), w["s" + c].to_numpy(), st, t0)
        out["v" + c] = beta[1]
    # geçici: önceki dönemin modelini son yıla uzat, artığın eğimine bak
    early = w["t"] < t0 - recent_years
    late = ~early
    if early.sum() < 200 or late.sum() < 100:
        out.update(tr_e=np.nan, tr_n=np.nan, tr_z=np.nan)
        return out
    te = t[early.to_numpy()]
    st_e = [s for s in steps_t if te[0] < s < te[-1]]
    trz = []
    for c in ("e", "n"):
        y = w[c].to_numpy()
        beta, res_e, keep = _robust_fit(te, y[early.to_numpy()], w["s" + c].to_numpy()[early.to_numpy()],
                                        st_e, t0)
        # son yıldaki adımları (ekipman, deprem) ayrıca çıkar
        tl = t[late.to_numpy()]
        st_l = [s for s in steps_t if tl[0] < s < tl[-1]]
        A_full = _design(t, t0, st_e)
        res = y - A_full @ beta
        rl = res[late.to_numpy()]
        Al = np.column_stack([np.ones_like(tl), tl - t0] + [(tl >= s).astype(float) for s in st_l])
        g, *_ = np.linalg.lstsq(Al, rl, rcond=None)
        sl = g[1]
        noise = 1.4826 * np.nanmedian(np.abs(res_e[keep] - np.nanmedian(res_e[keep])))
        # renkli gürültü için kabaca 3× şişirilmiş eğim belirsizliği
        se_slope = 3 * noise / max(np.sqrt(late.sum()) * np.std(tl), 1e-6)
        out["tr_" + c] = sl
        trz.append(sl / max(se_slope, 0.3))
    out["tr_z"] = float(np.hypot(*trz))
    return out


def station_year_table(raw: Path, sites: pd.DataFrame, years: list[int]) -> pd.DataFrame:
    steps = read_steps(raw / "steps.txt")
    by_site = {s: g["t"].to_numpy() for s, g in steps.groupby("site")}
    rows = []
    for _, s in sites.iterrows():
        p = raw / "tenv3" / f"{s['site']}.EU.tenv3"
        if not p.exists():
            continue
        try:
            ts = read_tenv3(p)
        except Exception as e:  # bozuk dosya
            log.warning("%s okunamadı: %s", s["site"], e)
            continue
        for y in years:
            f = station_year_features(ts, by_site.get(s["site"], np.array([])), float(y))
            if f:
                rows.append({"site": s["site"], "lat": s["lat"], "lon": s["lon"], "year": y, **f})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ hücre özellikleri
def _xy(lat, lon, lat0=39.0):
    return (np.asarray(lon) * 111.195 * np.cos(np.radians(lat0)), np.asarray(lat) * 111.195)


def strain_rate_at(points_lat, points_lon, sta: pd.DataFrame, sigma_km: float = 40.0,
                   max_km: float = 120.0, min_sta: int = 4) -> pd.DataFrame:
    """Gauss ağırlıklı yerel doğrusal hız gradyanından gerinim hızı (nanogerinim/yıl).

    Yöntem: Shen vd. 2015 (JGR) türü ağırlıklı en küçük kareler; hız =
    v0 + L·Δx, L gradyan tensörü. İkinci değişmez ve hacim değişimi verilir.
    """
    px, py = _xy(points_lat, points_lon)
    sx, sy = _xy(sta["lat"], sta["lon"])
    ve, vn = sta["ve"].to_numpy(), sta["vn"].to_numpy()
    out = np.full((len(px), 4), np.nan)
    for k in range(len(px)):
        dx, dy = sx - px[k], sy - py[k]
        d = np.hypot(dx, dy)
        sel = d <= max_km
        if sel.sum() < min_sta:
            continue
        w = np.exp(-(d[sel] / sigma_km) ** 2)
        A = np.column_stack([np.ones(sel.sum()), dx[sel], dy[sel]])
        W = np.sqrt(w)[:, None]
        ge, *_ = np.linalg.lstsq(A * W, ve[sel] * W[:, 0], rcond=None)
        gn, *_ = np.linalg.lstsq(A * W, vn[sel] * W[:, 0], rcond=None)
        exx, eyy = ge[1], gn[2]
        exy = 0.5 * (ge[2] + gn[1])
        # mm/yıl/km = 1e-6 /yıl = 1000 nanogerinim/yıl
        out[k] = [np.sqrt(exx ** 2 + eyy ** 2 + 2 * exy ** 2) * 1e3, (exx + eyy) * 1e3,
                  np.sqrt(((exx - eyy) / 2) ** 2 + exy ** 2) * 1e3, sel.sum()]
    return pd.DataFrame(out, columns=["strain_2inv", "dilatation", "max_shear", "n_sta"])


def cell_features(sty: pd.DataFrame, cells: pd.DataFrame, year: int,
                  transient_km: float = 75.0) -> pd.DataFrame:
    s = sty[sty["year"] == year]
    idx = cells["cell_id"]
    if s.empty:
        return pd.DataFrame(index=idx)
    sr = strain_rate_at(cells["lat_c"], cells["lon_c"], s)
    px, py = _xy(cells["lat_c"], cells["lon_c"])
    sx, sy = _xy(s["lat"], s["lon"])
    trz = s["tr_z"].to_numpy()
    trm = np.hypot(s["tr_e"], s["tr_n"]).to_numpy()
    zmax = np.full(len(px), np.nan)
    tmag = np.full(len(px), np.nan)
    for k in range(len(px)):
        d = np.hypot(sx - px[k], sy - py[k])
        sel = (d <= transient_km) & np.isfinite(trz)
        if sel.any():
            zmax[k] = np.nanmax(trz[sel])
            tmag[k] = np.nanmedian(trm[sel])
    return pd.DataFrame({
        "gnss_strain_rate_ns": sr["strain_2inv"].to_numpy(),
        "gnss_max_shear_ns": sr["max_shear"].to_numpy(),
        "gnss_dilatation_ns": sr["dilatation"].to_numpy(),
        "gnss_transient_z": zmax,
        "gnss_transient_mm_yr": tmag,
    }, index=idx)
