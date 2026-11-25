"""Eğitim hattı aşaması: kapasite, doluluk, kilitlenme ve zincirleme kırılma özellikleri.

Her yıl Y için yalnız Y öncesi veri kullanılır:
* GPS kilitlenme oranı (Y öncesi istasyon hızları) → elemanların birikim katsayısı
* bütçe 1900'den Y başına kadar (kilitlenmeli) → kayma açığı
* kalibrasyon: Y öncesi (1950 sonrası) M≥6 depremlerin kırılma anındaki doluluğu
* Coulomb: Y öncesi M≥6 kaynakların ΔCFF'si (μ' = ``coulomb_mu``)

Çıktılar (``runs/<ad>/``): capacity_cells.csv, capacity_segments.csv,
capacity_rupture_fills.csv, capacity_calibration.json, capacity_cascade_check.json,
locking_by_year.csv
"""

from __future__ import annotations

import json
import logging

import numpy as np
import pandas as pd

from ..features.table import assign_cell, make_cells
from ..region import TURKEY
from . import capacity as cp
from .budget import element_budget
from .coulomb import coulomb_on_receivers
from .crust import Crust1, download as crust_download
from .locking import KNOWN_CREEP, locking_table

log = logging.getLogger(__name__)

# Bilinen çok parçalı kırılmalar (zincirleme modelinin sınaması için; uzunluk km)
KNOWN_MULTI = {"19990817": ("1999 İzmit Mw7.6", 145.0), "19991112": ("1999 Düzce Mw7.2", 40.0),
               "20230206_0000008": ("2023 Pazarcık Mw7.8", 350.0),
               "20230206_0000222": ("2023 Elbistan Mw7.5", 150.0)}


def _dcff_receivers(run, year, mu):
    d = run.dir
    per_src = np.load(d / "coulomb_stress_per_source.npy", mmap_mode="r")
    src = pd.read_csv(d / "coulomb_sources.csv")
    rec = pd.read_csv(d / "coulomb_receivers.csv.gz")
    t0 = pd.Timestamp(f"{year}-01-01", tz="UTC")
    sel = (pd.to_datetime(src["time"], utc=True, format="ISO8601") < t0).to_numpy()
    stress = np.asarray(per_src[sel]).sum(0) if sel.any() else np.zeros((len(rec), 6))
    dcff, _, _ = coulomb_on_receivers(stress.astype(float), rec["strike"].to_numpy(), rec["dip"].to_numpy(),
                                      rec["rake"].fillna(180.0).to_numpy(), mu)
    return rec, np.clip(dcff / 1e5, -50, 50)  # bar


def _map_receivers(el, rec, values):
    """Her elemana aynı faydaki en yakın doğrultu indeksli alıcının değeri."""
    key = rec[["fault_id", "i_strike"]].copy()
    key["v"] = values
    key = key.sort_values(["fault_id", "i_strike"])
    e = el[["fault_id", "i_strike"]].reset_index()
    e = e.sort_values("i_strike")
    key = key.sort_values("i_strike")
    m = pd.merge_asof(e, key, on="i_strike", by="fault_id", direction="nearest")
    out = np.zeros(len(el))
    out[m["index"].to_numpy()] = m["v"].fillna(0).to_numpy()
    return out


def run_capacity(run) -> dict:
    cfg = run.cfg.get("capacity", {})
    root = run.root
    el = run.elements()
    faults = run.faults
    crust = Crust1(crust_download(root / "data" / "ext" / "crust1"))
    ph = cp.physical_capacity(el, crust)
    dphys = ph["dcap_phys_m"].to_numpy()
    reg = cp.element_region(el["lat"].to_numpy(), el["lon"].to_numpy())
    cells = make_cells(TURKEY, run.cfg["data"]["cell_deg"])
    cid = assign_cell(el.rename(columns={"lat": "latitude", "lon": "longitude"}), TURKEY, cells)
    sty_p = root / "data" / "obs" / "gnss_station_year.csv"
    sty = pd.read_csv(sty_p) if sty_p.exists() else None
    mu_c = float(cfg.get("coulomb_mu", 0.4))
    alpha = float(cfg.get("cascade_alpha", 0.5))
    thr = float(cfg.get("cascade_thr", 0.8))
    start = pd.Timestamp(run.cfg.get("budget", {}).get("start", "1900-01-01"), tz="UTC")
    years = list(run.years) + [run.years[-1] + 1]
    col_keys = cp.columns(el, np.zeros(len(el)))
    nb = cp.column_graph(col_keys)
    col_index = pd.MultiIndex.from_frame(col_keys[["fault_id", "i"]])
    el_col = col_index.get_indexer(pd.MultiIndex.from_arrays([el["fault_id"], el["i_strike"]]))
    rows, lock_rows, calib_all = [], [], {}
    last = None
    for Y in years:
        t0 = pd.Timestamp(f"{Y}-01-01", tz="UTC")
        lk = locking_table(faults, sty, Y) if sty is not None else pd.DataFrame()
        coupling = np.ones(len(el))
        lock_el = np.full(len(el), np.nan)
        creep_el = np.full(len(el), np.nan)
        if len(lk):
            m = el[["fault_id"]].merge(lk, on="fault_id", how="left")
            lock_el = m["locking"].to_numpy()
            creep_el = m["creep_mm_yr"].to_numpy()
            coupling = np.where(np.isfinite(lock_el), lock_el, 1.0)
            lock_rows.append(lk.assign(year=Y))
        res = element_budget(el, run.catalog, start, t0, snapshots=[t0],
                             min_mw=float(run.cfg.get("budget", {}).get("min_mw", 4.0)),
                             coupling=coupling, record_min_mw=6.0)
        deficit = res["snapshots"].get(t0, res["final"]["deficit_Nm"].to_numpy())
        rf = cp.rupture_fills(res["records"], el, dphys)
        calib = cp.hierarchical_calibration(rf, Y, int(cfg.get("calib_from", cp.CALIB_FROM_YEAR)))
        calib_all[Y] = calib
        scale, _ = cp.calibrated_scale(reg, calib, el["fault_id"].to_numpy())
        fill = cp.slip_m(deficit, el) / dphys / scale                 # kalibre doluluk (1 = tipik kırılma)
        rec, dcff = _dcff_receivers(run, Y, mu_c)
        dcff_el = _map_receivers(el, rec, dcff)
        fill_c = np.clip(fill + dcff_el * 1e5 / (cp.DSIGMA_PA * scale), 0, None)  # + aktarılan gerilme
        cols = cp.columns(el, fill_c)
        casc = cp.cascade_potential(cols, nb, seed_thr=thr, thr=thr, alpha=alpha)
        casc_el = casc[el_col]
        w_cap = ph["m0_cap_Nm"].to_numpy()
        a = el["area_km2"].to_numpy()
        df = pd.DataFrame({"cell_id": cid, "a": a, "wc": w_cap, "fill": fill, "fill_c": fill_c,
                           "casc": casc_el, "mu": ph["rigidity_gpa"].to_numpy(), "moho": ph["moho_km"].to_numpy(),
                           "lock": lock_el, "creep": creep_el})
        df = df[df["cell_id"] >= 0]
        g = df.groupby("cell_id")
        agg = pd.DataFrame({
            "cap_log_m0": np.log10(g["wc"].sum()),
            "rigidity_gpa": (df["mu"] * df["a"]).groupby(df["cell_id"]).sum() / g["a"].sum(),
            "moho_km": g["moho"].mean(),
            "fill_mean": (df["fill"] * df["a"]).groupby(df["cell_id"]).sum() / g["a"].sum(),
            "fill_m0w": (df["fill"] * df["wc"]).groupby(df["cell_id"]).sum() / g["wc"].sum(),
            "fill_max": g["fill"].max(),
            "fill_p90": g["fill"].quantile(0.9),
            "fillc_mean": (df["fill_c"] * df["a"]).groupby(df["cell_id"]).sum() / g["a"].sum(),
            "fillc_max": g["fill_c"].max(),
            "cascade_len_km": g["casc"].max(),
            "locking_mean": g["lock"].mean(),
            "creep_mm_yr": g["creep"].mean(),
        })
        agg["cascade_mw"] = cp.length_to_mw(agg["cascade_len_km"]).clip(min=5.0)
        agg["cascade_mw"] = np.where(agg["cascade_len_km"] > 0, agg["cascade_mw"], np.nan)
        agg["year"] = Y
        rows.append(agg.reset_index())
        last = dict(Y=Y, fill=fill, fill_c=fill_c, rf=rf, calib=calib, lock=lock_el, scale=scale,
                    records=res["records"], dcff_el=dcff_el)
        log.info("kapasite %d: kilitlenme çözülen fay %d, kalibrasyon n=%d, ortanca doluluk %.2f",
                 Y, len(lk), calib.get("n", 0), np.nanmedian(fill[el["slip_rate_mm_yr"] >= 1]))
    pd.concat(rows).to_csv(run.dir / "capacity_cells.csv", index=False)
    if lock_rows:
        lk_all = pd.concat(lock_rows)
        lk_all.to_csv(run.dir / "locking_by_year.csv", index=False)
    (run.dir / "capacity_calibration.json").write_text(json.dumps(
        {str(k): v for k, v in calib_all.items()}, indent=2))
    # ---- son yıl: kırılma anı doluluk tablosu (tüm geçmiş), segment yüzdeleri, zincirleme sınaması
    rf = last["rf"].copy()
    c_ = last["calib"]
    scale_ev = np.array([np.exp(c_.get(f"fay:{f}", c_.get(f"bölge:{r}", c_["genel"]))[0])
                         for r, f in zip(rf["bölge"], rf["fay"])])
    rf["doluluk_kalibre_%"] = 100 * rf["doluluk_fiziksel"] / scale_ev
    rf.drop(columns=["elemanlar", "ağırlık"]).to_csv(run.dir / "capacity_rupture_fills.csv", index=False)
    seg = segment_fill(el, last, cfg)
    seg.to_csv(run.dir / "capacity_segments.csv", index=False)
    check = cascade_check(el, last, nb, el_col, col_keys, dphys, alpha)
    lk_cmp = creep_comparison(faults, pd.read_csv(run.dir / "locking_by_year.csv") if lock_rows else None)
    (run.dir / "capacity_cascade_check.json").write_text(json.dumps(
        {"cascade": check, "creep_compare": lk_cmp}, indent=2, default=str))
    return {"years": len(years), "rupture_events": int(len(rf))}


def segment_fill(el, last, cfg, n_mc: int = 200, seed: int = 0) -> pd.DataFrame:
    """Fay başına yüzde doluluk (alan ağırlıklı) ve %16-84 belirsizlik aralığı.

    Belirsizlik: kalibrasyon (bölge θ sonsal sd'si ve olaydan olaya saçılım değil, tipik
    değerin belirsizliği), kilitlenme oranının sd'si ve rijitlikte %10.
    """
    rng = np.random.default_rng(seed)
    calib = last["calib"]
    reg = cp.element_region(el["lat"].to_numpy(), el["lon"].to_numpy())
    _, sd = cp.calibrated_scale(reg, calib, el["fault_id"].to_numpy())
    d = pd.DataFrame({"fault_id": el["fault_id"].to_numpy(), "a": el["area_km2"].to_numpy(),
                      "f": last["fill_c"], "sd": sd, "reg": reg, "lat": el["lat"].to_numpy(),
                      "lon": el["lon"].to_numpy(), "sr": el["slip_rate_mm_yr"].to_numpy()})
    out = []
    for fid, g in d.groupby("fault_id"):
        base = np.average(g["f"], weights=g["a"])
        z = rng.normal(0, 1, n_mc)
        s_cal = float(g["sd"].mean())
        mc = base * np.exp(-(s_cal * z) + rng.normal(0, 0.1, n_mc))
        out.append({"fault_id": fid, "bölge": g["reg"].iloc[0], "lat": g["lat"].mean(), "lon": g["lon"].mean(),
                    "kayma_hızı_mm_yıl": g["sr"].iloc[0],
                    "doluluk_%": 100 * base, "doluluk_%16": 100 * np.quantile(mc, 0.16),
                    "doluluk_%84": 100 * np.quantile(mc, 0.84),
                    "en_dolu_1km_%": 100 * g["f"].max()})
    return pd.DataFrame(out).sort_values("doluluk_%", ascending=False)


def cascade_check(el, last, nb, el_col, col_keys, dphys, alpha):
    """Bilinen çok parçalı kırılmalar: depremden hemen önceki doluluktan, episantra en yakın
    sütundan başlatılan zincirleme kırılma ne kadar uzar? (eşik 1,0 ve 0,8)"""
    out = []
    reg = cp.element_region(el["lat"].to_numpy(), el["lon"].to_numpy())
    scale, _ = cp.calibrated_scale(reg, last["calib"], el["fault_id"].to_numpy())
    xy_c = np.column_stack([col_keys["lon"].to_numpy() * 85.0, col_keys["lat"].to_numpy() * 111.2])
    for r in last["records"]:
        key = next((k for k in KNOWN_MULTI if str(r["event_id"]).startswith(k)), None)
        if key is None:
            continue
        fill = cp.slip_m(r["deficit_before"], el) / dphys / scale
        cols = cp.columns(el, fill)
        seed = int(np.argmin(np.hypot(xy_c[:, 0] - r["lon"] * 85.0, xy_c[:, 1] - r["lat"] * 111.2)))
        res = {"deprem": KNOWN_MULTI[key][0], "gözlenen_uzunluk_km": KNOWN_MULTI[key][1],
               "başlangıç_doluluğu_%": round(100 * float(cols["fill"].iloc[seed]), 1)}
        for thr in (1.0, 0.8, 0.5):
            b = cp.cascade(cols["fill"].to_numpy(), nb, seed, thr=thr, alpha=alpha, forced=True)
            res[f"model_uzunluk_km_eşik{thr:g}"] = round(float(cols["len_km"].to_numpy()[b].sum()), 1)
        out.append(res)
    return out


def creep_comparison(faults, lk_all):
    if lk_all is None or lk_all.empty:
        return []
    from ..geo import point_to_polyline_km
    last = lk_all[lk_all["year"] == lk_all["year"].max()]
    out = []
    for _, k in KNOWN_CREEP.iterrows():
        d = np.array([point_to_polyline_km([k["lat"]], [k["lon"]], c)[0] for c in faults["coords"]])
        fid = faults["fault_id"].iloc[int(d.argmin())]
        m = last[last["fault_id"] == fid]
        out.append({"yer": k["yer"], "literatür_krip_mm_yıl": k["krip_mm_yil"], "kaynak": k["kaynak"],
                    "fay": fid, "gps_krip_mm_yıl": float(m["creep_mm_yr"].iloc[0]) if len(m) else None,
                    "gps_krip_sd": float(m["creep_sd"].iloc[0]) if len(m) else None,
                    "not": None if len(m) else "fayın iki yanında yeterli açık GPS istasyonu yok"})
    return out


def write_fill_layer(run, out_dir) -> int:
    """Haritaya fay başına yüzde doluluk katmanı (belirsizlik aralığıyla)."""
    from ..external.layers import write_layer
    seg = pd.read_csv(run.dir / "capacity_segments.csv")
    f = run.faults.merge(seg, on="fault_id", how="inner")
    feats = []
    for _, r in f.iterrows():
        coords = [[round(float(x), 4), round(float(y), 4)] for x, y in np.asarray(r["coords"])[:, :2]]
        feats.append({"type": "Feature", "geometry": {"type": "LineString", "coordinates": coords},
                      "properties": {"fault_id": r["fault_id"],
                                     "doluluk": round(float(r["doluluk_%"]), 0),
                                     "doluluk_alt": round(float(r["doluluk_%16"]), 0),
                                     "doluluk_ust": round(float(r["doluluk_%84"]), 0),
                                     "en_dolu_1km": round(float(r["en_dolu_1km_%"]), 0),
                                     "kayma_hizi": round(float(r["kayma_hızı_mm_yıl"]), 1)}})
    spec = {"id": "fay_doluluk", "label": "Fay doluluk oranı (%)",
            "description": ("Her fay için biriken kayma açığının, o bölgede büyük depremlerin (M≥6,8, 1950+) "
                            "tipik olarak kırıldığı birikime oranı. %100 = geçmişte tipik kırılma düzeyi. "
                            "GPS kilitlenmesi ve Coulomb aktarımı dahil. Bütçe 1900'de sıfırdan başladığı için "
                            "göreli bir ölçüdür; geçmiş depremler çok farklı doluluklarda kırılmıştır."),
            "group": "Faylar ve enerji", "order": 12,
            "note": "Deneysel araştırma çıktısıdır, deprem tahmini değildir. Doluluk eşiği yöntemi sınamada doğrulanmadı.",
            "updated": pd.Timestamp.now().strftime("%Y-%m-%d"),
            "style": {"type": "line", "width": 3, "glow": True, "opacity": 0.9},
            "color": {"property": "doluluk", "type": "interpolate",
                      "stops": [[0, "#3f7be0"], [50, "#7fa6e8"], [100, "#ffd166"], [200, "#ff7b39"], [400, "#ff3b30"]]},
            "legend": {"title": "Doluluk", "unit": "%", "labels": ["0", "100", "≥400"]},
            "tooltip": {"title": "fault_id", "fields": [
                {"property": "doluluk", "label": "Doluluk", "unit": "%", "digits": 0},
                {"property": "doluluk_alt", "label": "Alt (%16)", "unit": "%", "digits": 0},
                {"property": "doluluk_ust", "label": "Üst (%84)", "unit": "%", "digits": 0},
                {"property": "en_dolu_1km", "label": "En dolu 1 km", "unit": "%", "digits": 0},
                {"property": "kayma_hizi", "label": "Kayma hızı", "unit": "mm/yıl", "digits": 1}]}}
    write_layer(out_dir, spec, {"type": "FeatureCollection", "metadata": {"note": spec["note"]}, "features": feats})
    return len(feats)
