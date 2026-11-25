"""İstasyon çifti dv/v pilotunun sınaması (Doğu Marmara).

    python scripts/eval_dvv_pairs.py

Önce: scripts/fetch_dvv_pairs.py ccf + dvv ve scripts/fetch_catalogs.py.
Rapor: reports/dvv_pairs_marmara.md, şekil: reports/figures/dvv_pairs_marmara.png

Sınamalar (hepsi önceden belirlendi; sonuç görüldükten sonra değiştirilmedi):

1. Ölçüm kalitesi: güvenilir ölçüm oranı (cc > 0,5) ve belirsizlik; PR 9'daki tek
   istasyon yöntemiyle (%11 güvenilir, %0,2-0,7 belirsizlik) karşılaştırma.
2. Yöntem doğrulaması (öncül değil): mevsimsel döngü ve 2019 Silivri M5,7 depremi
   sırasındaki eş-sismik hız düşüşü görülebiliyor mu?
3. Öncül sınaması: ön/artçılarından arındırılmış (Gardner-Knopoff) M≥3 ve M≥4
   depremlerinden önce, çevredeki çiftlerde dv/v değişimi rastgele günlerden farklı mı?
   İki ölçek: kısa (son 60 gün − önceki 180 gün) ve yıllık modelin özelliği
   (son 6 ay − önceki 2 yıl). Yalnız deprem anından önceki ölçümler kullanılır.
   Alarm sınaması: değişim en düşük %10'luk dilimdeyse "alarm"; depremlerin yüzde kaçı
   alarm sırasında geldi? Rastgele beklenti %10'dur. Güven aralığı bootstrap ile.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from deprem.catalog.clean import decluster_gardner_knopoff  # noqa: E402
from deprem.observations import dvv_pairs as dp  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OBS = ROOT / "data" / "obs"
REP = ROOT / "reports"
log = logging.getLogger("eval_dvv_pairs")

SILIVRI = dict(time=pd.Timestamp("2019-09-26 10:59"), lat=40.87, lon=28.19, mag=5.7)
# Doğrulama olayları: bölgedeki 2009 sonrası M≥5,5 depremler (kataloğa göre)
VALIDATION = [SILIVRI, dict(time=pd.Timestamp("2025-04-23 09:49"), lat=40.854, lon=28.226, mag=6.2)]


def split_noise(pairs: pd.DataFrame, n: int = 20, seed: int = 0) -> float:
    """Çiftler rastgele iki yarıya bölünür; iki yarının farkının std'si / 2 = bölgesel seri gürültüsü."""
    rng = np.random.default_rng(seed)
    keys = pairs[["a", "b"]].drop_duplicates().apply(tuple, axis=1).tolist()
    key = pairs[["a", "b"]].apply(tuple, axis=1)
    out = []
    for _ in range(n):
        half = set(keys[i] for i in rng.permutation(len(keys))[: len(keys) // 2])
        m = key.isin(half)
        A = dp.regional_series(pairs[m]).set_index("time")["dvv"]
        B = dp.regional_series(pairs[~m]).set_index("time")["dvv"]
        j = pd.concat([A, B], axis=1, keys=["A", "B"], sort=True).dropna()
        out.append((j["A"] - j["B"]).std() / 2)
    return round(100 * float(np.mean(out)), 4)


def local_series(pairs: pd.DataFrame, lat: float, lon: float, radius_km: float) -> pd.DataFrame:
    d = np.hypot((pairs["lon"] - lon) * 85.0, (pairs["lat"] - lat) * 111.0)
    p = pairs[(d <= radius_km)].dropna(subset=["dvv", "sigma"])
    if p.empty:
        return pd.DataFrame(columns=["time", "dvv", "sigma", "n_pairs"])
    return dp.regional_series(p)


def change(series: pd.DataFrame, t: pd.Timestamp, recent_days: int, before_days: int):
    """Son ``recent_days`` ortalaması − önceki ``before_days`` ortalaması (yalnız t öncesi)."""
    s = series[series["time"] < t]
    r = s[s["time"] >= t - pd.Timedelta(days=recent_days)]["dvv"]
    b = s[(s["time"] < t - pd.Timedelta(days=recent_days)) &
          (s["time"] >= t - pd.Timedelta(days=recent_days + before_days))]["dvv"]
    if len(r) < max(3, recent_days // 10) or len(b) < max(10, before_days // 10):
        return np.nan
    return 100 * (r.mean() - b.mean())


SCALES = {"kısa (son 60 gün − önceki 180 gün)": (60, 180), "yıllık (6 ay − önceki 2 yıl)": (182, 730)}


def precursor_test(pairs, events, rng, radius_fn, q=0.10, n_boot=2000):
    rows = []
    days = pd.date_range(pairs["time"].min() + pd.Timedelta(days=1000), pairs["time"].max(), freq="7D")
    for label, (rd, bd) in SCALES.items():
        ev_x, ctrl = [], []
        cache = {}
        for e in events.itertuples():
            R = radius_fn(e.mag)
            key = (round(e.latitude, 1), round(e.longitude, 1), R)
            if key not in cache:
                ls = local_series(pairs, e.latitude, e.longitude, R)
                cache[key] = (ls, np.array([change(ls, t, rd, bd) for t in days]))
            ls, cx = cache[key]
            x = change(ls, e.time, rd, bd)
            cx = cx[np.isfinite(cx)]
            if np.isnan(x) or len(cx) < 50:
                continue
            thr = np.quantile(cx, q)
            ev_x.append((x, (x <= thr), (cx <= x).mean(), e.mag))
        if not ev_x:
            rows.append({"ölçek": label, "olay": 0})
            continue
        a = np.array(ev_x, float)
        hits = a[:, 1]
        boot = [rng.choice(hits, len(hits)).mean() for _ in range(n_boot)]
        pct = a[:, 2]  # olayın, aynı yerdeki tüm günlere göre yüzdelik sırası (düşük = hız düşüşü)
        bootp = [rng.choice(pct, len(pct)).mean() for _ in range(n_boot)]
        rows.append({"ölçek": label, "olay": len(a),
                     "alarm içinde (%)": round(100 * hits.mean(), 1),
                     "%95 GA": f"[{100 * np.quantile(boot, .025):.1f}, {100 * np.quantile(boot, .975):.1f}]",
                     "rastgele beklenti (%)": 100 * q,
                     "olasılık kazancı": round(hits.mean() / q, 2),
                     "ortalama yüzdelik sıra": round(float(pct.mean()), 3),
                     "sıra %95 GA": f"[{np.quantile(bootp, .025):.3f}, {np.quantile(bootp, .975):.3f}]"})
    return pd.DataFrame(rows)


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    pairs = pd.read_csv(OBS / "dvv_pairs_series.csv.gz", parse_dates=["time"])
    reg = pd.read_csv(OBS / "dvv_pairs_regional.csv", parse_dates=["time"])
    cache = OBS / "dvv_pairs_all.pkl"  # kalite istatistiği için süzülmemiş hali (hesabı uzun, önbellekte)
    if cache.exists():
        allp = pd.read_pickle(cache)
    else:
        allp = dp.all_pairs_dvv(OBS / "dvv_pairs", min_snr=0)
        allp.to_pickle(cache)
    rng = np.random.default_rng(0)
    out = {}

    # 1. kalite
    n_tot = len(allp)
    rel = allp["dvv"].notna().mean()
    used = pairs[["a", "b"]].drop_duplicates()
    out["kalite"] = {
        "örnek gün": int(allp["time"].nunique()),
        "çift (toplam / kullanılan, SNR≥3)": f"{allp[['a', 'b']].drop_duplicates().shape[0]} / {len(used)}",
        "ölçüm (çift × zaman)": n_tot,
        "güvenilir ölçüm oranı (cc>0,5)": f"%{100 * rel:.0f}",
        "tek ölçüm belirsizliği, ortanca (%)": round(100 * pairs["sigma"].median(), 3),
        "bölgesel seri belirsizliği, ortanca (%)": round(100 * reg["sigma"].median(), 4),
        "bölgesel seri gürültüsü, çiftler iki yarıya bölünerek (%)": split_noise(pairs),
    }

    # 2. doğrulama: mevsim + Silivri
    r = reg.dropna()
    tt = (r["time"] - pd.Timestamp("2000-01-01")).dt.days.to_numpy() / 365.25
    X = np.column_stack([np.ones_like(tt), tt, np.sin(2 * np.pi * tt), np.cos(2 * np.pi * tt)])
    # yavaş eğilim için 2 yıllık kayan ortalama çıkar, sonra mevsim
    trend = r.set_index("time")["dvv"].rolling("730D", center=True, min_periods=50).mean().to_numpy()
    yv = r["dvv"].to_numpy() - trend
    ok = np.isfinite(yv)
    beta, *_ = np.linalg.lstsq(X[ok], yv[ok], rcond=None)
    res = yv[ok] - X[ok] @ beta
    cov = np.linalg.inv(X[ok].T @ X[ok]) * res.var() * 30  # 60 günlük yığın: ~30 örnek ilişkili
    amp = np.hypot(beta[2], beta[3])
    amp_sd = np.sqrt((beta[2] ** 2 * cov[2, 2] + beta[3] ** 2 * cov[3, 3]) / max(amp ** 2, 1e-20))
    peak_frac = (np.arctan2(beta[2], beta[3]) / (2 * np.pi)) % 1.0  # en yüksek hızın yıl içindeki yeri
    out["mevsim"] = {"yıllık genlik (%)": f"{100 * amp:.3f} ± {100 * amp_sd:.3f}",
                     "en yüksek hızın ayı": int(peak_frac * 12) + 1}
    out["silivri"] = {}
    for ev in VALIDATION:
        sil = local_series(pairs, ev["lat"], ev["lon"], 60)
        t0 = ev["time"]
        pre = sil[(sil["time"] < t0) & (sil["time"] >= t0 - pd.Timedelta(days=60))]
        post = sil[(sil["time"] >= t0 + pd.Timedelta(days=60)) & (sil["time"] < t0 + pd.Timedelta(days=90))]  # yığın tamamen sonrası
        if len(pre) and len(post):
            d = post["dvv"].mean() - pre["dvv"].mean()
            sd = float(sil.set_index("time")["dvv"].diff(45).dropna().std())  # olağan 90 günlük değişim
            out["silivri"][f"{t0:%Y-%m-%d} M{ev['mag']}"] = (
                f"{int(sil['n_pairs'].median())} çift (60 km); sonrası − öncesi %{100 * d:.3f}; "
                f"olağan 90 günlük değişim std %{100 * sd:.3f}; z = {d / sd:.2f}")
    sil = local_series(pairs, SILIVRI["lat"], SILIVRI["lon"], 60)

    # 3. öncül sınaması
    cat = pd.read_csv(ROOT / "data" / "processed" / "turkey_catalog.csv.gz")
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    cat = cat.sort_values("time").reset_index(drop=True)
    cat["independent"] = decluster_gardner_knopoff(cat)
    cat["time"] = cat["time"].dt.tz_localize(None)
    box = (cat["latitude"].between(pairs["lat"].min() - 0.3, pairs["lat"].max() + 0.3) &
           cat["longitude"].between(pairs["lon"].min() - 0.4, pairs["lon"].max() + 0.4))
    start = pairs["time"].min() + pd.Timedelta(days=1000)
    tests = []
    for m in (3.0, 4.0):
        ev = cat[box & cat["independent"] & (cat["mag"] >= m) & (cat["time"] >= start)]
        tab = precursor_test(pairs, ev, rng, lambda M: max(30.0, 10 ** (0.43 * M)))
        tab.insert(0, "eşik", f"M≥{m:g}")
        tests.append(tab)
    out["öncül"] = pd.concat(tests, ignore_index=True)

    # şekil
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    ax[0].fill_between(reg["time"], 100 * (reg["dvv"] - 2 * reg["sigma"]), 100 * (reg["dvv"] + 2 * reg["sigma"]),
                       color="0.8", lw=0)
    ax[0].plot(reg["time"], 100 * reg["dvv"], lw=0.8)
    ax[0].set_ylabel("dv/v (%) bölgesel")
    ax[1].plot(sil["time"], 100 * sil["dvv"], lw=0.8, color="C1")
    ax[1].set_ylabel("dv/v (%) Silivri 60 km")
    big = cat[box & (cat["mag"] >= 4.5) & (cat["time"] >= reg["time"].min())]
    for a_ in ax:
        for t in big["time"]:
            a_.axvline(t, color="r", lw=0.6, alpha=0.6)
        a_.grid(alpha=0.3)
    ax[0].set_title("Doğu Marmara, istasyon çifti dv/v (kırmızı: M≥4,5 depremler)")
    fig.tight_layout()
    (REP / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(REP / "figures" / "dvv_pairs_marmara.png", dpi=120)

    import json
    (OBS / "dvv_pairs_eval.json").write_text(json.dumps(
        {k: (v.to_dict("records") if isinstance(v, pd.DataFrame) else v) for k, v in out.items()},
        indent=2, ensure_ascii=False, default=str))
    for k, v in out.items():
        print(f"\n== {k}")
        print(v.to_markdown(index=False) if isinstance(v, pd.DataFrame) else
              "\n".join(f"{a}: {b}" for a, b in v.items()))


if __name__ == "__main__":
    main()
