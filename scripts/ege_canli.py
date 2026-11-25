"""Ege canlı (ileriye dönük) tahmin kaydı ve değerlendirmesi.

Ön kayıt: docs/onkayit/2026-11-25_ege_ustunluk.md, bölüm 8.

    python scripts/ege_canli.py tahmin                 # ertesi UTC günü için tahmini yaz
    python scripts/ege_canli.py degerlendir --katalog data/processed/turkey_catalog.csv.gz

``tahmin``: son günlerin depremlerini USGS + EMSC (+ erişilirse AFAD) ile çeker, durum
kataloğunu (forecasts/ege/katalog_canli.csv.gz) günceller ve ertesi UTC günü için Ege
kutusundaki her 0,5° hücrenin M≥3 ve M≥4 beklenen sayısını iki model için yazar:

* ETAS: 2026 ileri yürüyen ETAS'ın donmuş parametreleri; arka plan + bilinen depremlerin
  ertesi 24 saatteki doğrudan artçıları (Omori integrali; uzayda ETAS çekirdeğiyle,
  sabit tohumlu Monte Carlo ile hücrelere dağıtılır).
* Model: aynı günlük ETAS × hücrenin 2026 için dondurulmuş model çarpanı
  (forecasts/ege/model_carpani_2026.csv; ön kayıtlı sınamada üretildi, değişmez).

Dosyalar yalnız eklenir (geçmiş tahmin değiştirilmez); git commit zamanı, tahminin
pencere başlamadan önce verildiğinin kanıtıdır.

HİÇBİR ÇIKTI DEPREM TAHMİNİ (UYARI) DEĞİLDİR: olasılıksal bir hipotezin sınamasıdır.
"""

from __future__ import annotations

import argparse
import json
import logging
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from deprem.catalog.clean import basic_clean, merge_sources
from deprem.catalog.sources import FetchSpec, fetch
from deprem.forecast.etas import LAT0, LON0, _days, _omori_integral, _xy_km
from deprem.region import TURKEY

ROOT = Path(__file__).resolve().parents[1]
LIVE = ROOT / "forecasts/ege"
STATE = LIVE / "katalog_canli.csv.gz"
DAILY = LIVE / "gunluk"
EGE_NAME = "Ege / Batı Anadolu"
CELL = 0.5
KEEP = ["time", "latitude", "longitude", "depth_km", "mag", "mag_type", "source", "event_id"]
log = logging.getLogger("ege_canli")


def grid_edges():
    lat = np.arange(TURKEY.min_lat, TURKEY.max_lat + 1e-9, CELL)
    lon = np.arange(TURKEY.min_lon, TURKEY.max_lon + 1e-9, CELL)
    return lat, lon


def cell_of(lat, lon):
    le, oe = grid_edges()
    i = np.floor((np.asarray(lat) - le[0]) / CELL).astype(int)
    j = np.floor((np.asarray(lon) - oe[0]) / CELL).astype(int)
    ok = (i >= 0) & (i < le.size - 1) & (j >= 0) & (j < oe.size - 1)
    return np.where(ok, i * (oe.size - 1) + j, -1)


def load_state():
    cat = pd.read_csv(STATE)
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    return cat


# ------------------------------------------------------------------ katalog güncelleme
def update_catalog(now: pd.Timestamp, days: int = 5) -> tuple[pd.DataFrame, list]:
    cat = load_state()
    t0 = (now - pd.Timedelta(days=days)).floor("D")
    frames, used = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for src, mm in (("afad", 0.0), ("emsc", 2.5), ("usgs", 2.5)):
            try:
                f = fetch(src, FetchSpec(TURKEY, t0, now, mm, chunk="1D"), Path(tmp), refresh=True)
                frames.append(basic_clean(f, TURKEY))
                used.append(f"{src}:{len(f)}")
            except Exception as e:  # bir kaynak çökerse diğerleriyle devam
                log.warning("%s çekilemedi: %s", src, e)
                used.append(f"{src}:hata")
    fresh = [f for f in frames if not f.empty]
    if fresh:
        new = merge_sources(fresh, priority=("afad", "emsc", "usgs"))
        new = new[new["mag"] >= 3.0]
        old = cat[cat["time"] < t0]
        cat = pd.concat([old, new[[c for c in KEEP if c in new]]], ignore_index=True)
        cat = cat.sort_values("time").reset_index(drop=True)
    return cat, used


# ------------------------------------------------------------------ günlük ETAS
def daily_etas(cat, d0, d1, issued, etas, U, seed):
    """[d0, d1) için hücre başına beklenen M≥m0 sayısı (tüm Türkiye ızgarası)."""
    P, m0, b = etas["params"], etas["m0"], etas["b"]
    le, oe = grid_edges()
    n_cells = (le.size - 1) * (oe.size - 1)
    T = (d1 - d0).total_seconds() / 86400.0
    out = P["mu"] * U * T
    par = cat[(cat["mag"] >= m0) & (cat["time"] < issued)
              & (cat["time"] >= d0 - pd.Timedelta(days=etas["dt_max_days"]))]
    if par.empty:
        return out
    tp = _days(par["time"], d0)  # negatif: pencereden önce
    lo = -tp
    hi = np.minimum(T - tp, etas["dt_max_days"])
    m = par["mag"].to_numpy() - m0
    g = P["K"] * np.exp(P["alpha"] * m) * np.maximum(
        _omori_integral(hi, P["c"], P["p"]) - _omori_integral(lo, P["c"], P["p"]), 0.0)
    x, y = _xy_km(par["latitude"], par["longitude"])
    sig = P["D"] ** 2 * np.exp(P["gamma"] * m)
    big = g > 1e-5
    # küçük katkılar: kendi hücresine
    c_small = cell_of(par["latitude"].to_numpy()[~big], par["longitude"].to_numpy()[~big])
    ok = c_small >= 0
    out = out + np.bincount(c_small[ok], weights=g[~big][ok], minlength=n_cells)
    # büyük katkılar: ETAS uzay çekirdeğinden örnekle (bölge dışına düşen kütle kaybolur, ETAS'taki gibi)
    rng = np.random.default_rng(seed)
    S = 400
    q = P["q"]
    u = rng.uniform(size=(big.sum(), S))
    r = np.sqrt(sig[big, None] * ((1 - u) ** (1 / (1 - q)) - 1))
    th = rng.uniform(0, 2 * np.pi, size=r.shape)
    xs, ys = x[big, None] + r * np.cos(th), y[big, None] + r * np.sin(th)
    lat = ys / 111.195 + LAT0
    lon = xs / (111.195 * np.cos(np.radians(LAT0))) + LON0
    c = cell_of(lat.ravel(), lon.ravel())
    w = np.repeat(g[big] / S, S)
    ok = c >= 0
    return out + np.bincount(c[ok], weights=w[ok], minlength=n_cells)


def tahmin(args):
    now = pd.Timestamp.now(tz="UTC") if args.simdi is None else pd.Timestamp(args.simdi, tz="UTC")
    d0 = now.ceil("D") if now != now.floor("D") else now
    d1 = d0 + pd.Timedelta(days=1)
    out = DAILY / f"{d0:%Y}" / f"{d0:%Y-%m-%d}.csv"
    if out.exists():
        log.info("%s zaten var; geçmiş tahmin değiştirilmez", out)
        return
    if args.katalog_guncelleme:
        cat, used = update_catalog(now)
        cat.to_csv(STATE, index=False)
    else:
        cat, used = load_state(), ["güncelleme yok"]
    cat = cat[cat["time"] < now]
    etas = json.loads((LIVE / "etas_2026.json").read_text())
    U = pd.read_csv(LIVE / "etas_arka_plan_2026.csv").sort_values("cell_id")["arka_plan_payi"].to_numpy()
    mult = pd.read_csv(LIVE / "model_carpani_2026.csv")
    ege = mult[mult["bolge"] == EGE_NAME].sort_values("cell_id")
    lam0 = daily_etas(cat, d0, d1, now, etas, U, seed=int(d0.strftime("%Y%m%d")))
    cid = ege["cell_id"].to_numpy()
    b, m0 = etas["b"], etas["m0"]
    s3, s4 = 10 ** (-b * (3.0 - m0)), 10 ** (-b * (4.0 - m0))
    e3 = lam0[cid] * s3
    k = ege["model_carpani"].to_numpy()
    df = pd.DataFrame({"cell_id": cid, "lat_c": ege["lat_c"].to_numpy(), "lon_c": ege["lon_c"].to_numpy(),
                       "etas_m3": e3, "model_m3": e3 * k,
                       "etas_m4": e3 * s4 / s3, "model_m4": e3 * k * s4 / s3})
    out.parent.mkdir(parents=True, exist_ok=True)
    head = (f"# Ege günlük tahmini. Pencere: {d0.isoformat()} - {d1.isoformat()} (UTC). "
            f"Verildiği an: {now.isoformat()}. Son bilinen deprem: {cat['time'].max().isoformat()}. "
            f"Kaynaklar: {' '.join(used)}. Beklenen sayılar (olay/gün). Deprem uyarısı değildir.\n")
    with open(out, "w") as fh:
        fh.write(head)
        df.to_csv(fh, index=False, float_format="%.6g")
    summ = LIVE / "gunluk_ozet.csv"
    row = pd.DataFrame([{"pencere_baslangic": d0.isoformat(), "verildi": now.isoformat(),
                         "etas_m3_toplam": e3.sum(), "model_m3_toplam": (e3 * k).sum(),
                         "son_deprem": cat["time"].max().isoformat(), "kaynaklar": " ".join(used)}])
    row.to_csv(summ, mode="a", header=not summ.exists(), index=False, float_format="%.5g")
    log.info("yazıldı: %s (ETAS M≥3 toplam %.3f, model %.3f)", out, e3.sum(), (e3 * k).sum())


# ------------------------------------------------------------------ değerlendirme
def degerlendir(args):
    rows = []
    for p in sorted(DAILY.glob("*/*.csv")):
        d = pd.read_csv(p, comment="#")
        d["gun"] = pd.Timestamp(p.stem, tz="UTC")
        rows.append(d)
    if not rows:
        raise SystemExit("henüz günlük tahmin yok")
    fc = pd.concat(rows, ignore_index=True)
    cat = pd.read_csv(args.katalog)
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    cat["gun"] = cat["time"].dt.floor("D")
    cat["cell_id"] = cell_of(cat["latitude"].to_numpy(), cat["longitude"].to_numpy())
    days = set(fc["gun"])
    cat = cat[cat["gun"].isin(days) & cat["cell_id"].isin(set(fc["cell_id"]))]
    first, last = min(days), max(days)
    note = "" if last >= pd.Timestamp("2027-05-25", tz="UTC") else " (ARA BAKIŞ: ön kayıtlı 6 ay dolmadı, karar sayılmaz)"
    res = {"gun_sayisi": len(days), "donem": [first.isoformat(), last.isoformat()], "not": note.strip()}
    for mg, ce, cm in ((3.0, "etas_m3", "model_m3"), (4.0, "etas_m4", "model_m4")):
        obs = cat[cat["mag"] >= mg].groupby(["cell_id", "gun"]).size().rename("n")
        d = fc.set_index(["cell_id", "gun"]).join(obs).fillna({"n": 0}).reset_index()
        d["ay"] = d["gun"].dt.strftime("%Y-%m")
        blk = d.groupby(["cell_id", "ay"]).agg(n=("n", "sum"), lm=(cm, "sum"), le=(ce, "sum"))
        # olay başına kazanç: günlük hücre olasılıkları üzerinden (Poisson), bloklara göre topla
        d["a"] = d["n"] * np.log(d[cm] / d[ce]) - (d[cm] - d[ce])
        a = d.groupby(["cell_id", "ay"])["a"].sum().to_numpy()
        n = blk["n"].to_numpy(float)
        N = n.sum()
        g = float(a.sum() / N) if N else np.nan
        rng = np.random.default_rng(2026)
        pick = rng.integers(0, len(a), (10_000, len(a)))
        den = n[pick].sum(1)
        bs = a[pick].sum(1)[den > 0] / den[den > 0]
        p = float((bs <= 0).mean()) if N else np.nan
        res[f"M≥{mg:g}"] = {"olay": int(N), "IG": g, "GA95": [float(np.quantile(bs, .025)), float(np.quantile(bs, .975))]
                            if N else None, "p_tek_yonlu": p,
                            "etas_beklenen": float(blk["le"].sum()), "model_beklenen": float(blk["lm"].sum()),
                            # ikincil: günlük tahminler hücre × ay toplanarak aynı ölçüt
                            "IG_hucre_ay_toplam": float(((n * np.log(blk["lm"] / blk["le"])).sum()
                                                         - (blk["lm"].sum() - blk["le"].sum())) / N) if N else None}
    m3 = res["M≥3"]
    if m3["olay"] < 30:
        k = "kararsız (30'dan az olay)"
    elif m3["IG"] <= 0:
        k = "kaldı"
    elif m3["p_tek_yonlu"] >= 0.05:
        k = "kararsız"
    else:
        k = "geçti (güçlü)" if m3["IG"] >= 0.10 else "geçti (zayıf)"
    res["karar"] = k + note
    print(json.dumps(res, indent=2, ensure_ascii=False))
    if args.yaz:
        (LIVE / f"degerlendirme_{last:%Y-%m-%d}.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="adim", required=True)
    t = sub.add_parser("tahmin")
    t.add_argument("--simdi", default=None, help="test için verildiği an (UTC)")
    t.add_argument("--katalog-guncelleme", action=argparse.BooleanOptionalAction, default=True)
    d = sub.add_parser("degerlendir")
    d.add_argument("--katalog", type=Path, default=ROOT / "data/processed/turkey_catalog.csv.gz")
    d.add_argument("--yaz", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    {"tahmin": tahmin, "degerlendir": degerlendir}[args.adim](args)


if __name__ == "__main__":
    main()
