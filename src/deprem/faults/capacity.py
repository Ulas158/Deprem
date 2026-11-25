"""Fay elemanlarının enerji biriktirme kapasitesi, doluluk oranı ve zincirleme kırılma.

Her 1 km elemanı için:

1. **Fiziksel kapasite** (kayma cinsinden): elemanın kilitli genişliği W, kayaç
   sertliği μ (CRUST1.0) ve kayacın taşıyabileceği en büyük gerilme düşümü Δσ ile

       D_kap = 2 · Δσ · W / (π · μ)          (doğrultu atımlı sonsuz fay, Knopoff 1958)

   Moment kapasitesi μ·A·D_kap. Sert kayaç aynı gerilmeyi daha az kaymayla taşır.
2. **Doluluk** = biriken kayma açığı / kapasite. Açık, ``budget.element_budget``
   ile hesaplanır (birikim − deprem boşalımı); GPS'ten kilitlenme oranı L
   bilinen faylarda birikim L ile çarpılır (sürünen kısım enerjiyi yavaş yavaş boşaltır).
3. **Coulomb aktarımı**: komşu depremlerin elemana aktardığı ΔCFF, gerilme
   kapasitesine bölünüp doluluğa eklenir (boşalan parçanın yükü komşuya geçer).
4. **Kalibrasyon** ("%70'te kırılıyorsa kapasitesi odur"): geçmiş M≥6 depremlerin
   kırıldığı andaki fiziksel doluluk oranları toplanır; bölge → genel hiyerarşik
   (kısmi havuzlama) modelle her bölge için tipik kırılma doluluğu bulunur ve
   kapasite buna ölçeklenir. Böylece kalibre doluluk %100 = o bölgede depremlerin
   tipik olarak kırıldığı birikim. Yalnız hesap yılından önceki depremler kullanılır.
5. **Zincirleme kırılma**: doğrultu boyunca sütunlar (ve uçları 5 km'den yakın faylar)
   komşudur. Eşiği aşan bir sütun kırılınca yükünün α kadarını komşularına dağıtır;
   eşiği aşan komşu da kırılır; döngü durana kadar sürer. Kırılan toplam uzunluktan
   olası büyüklük çıkarılır (Wells & Coppersmith 1994).

Sınırlar: bütçe 1900'de sıfırdan başlar (gerçek başlangıç açığı bilinmez), bu
yüzden doluluk özellikle uzun yinelenme aralıklı faylarda göreli bir ölçüdür.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from ..region import SUBREGIONS

DSIGMA_PA = 10e6        # taşınabilecek en büyük gerilme düşümü (tipik büyük depremler 3-15 MPa)
DSIGMA_LOG_SD = np.log(2.0)
CALIB_FROM_YEAR = 1950  # bütçe 1900'de sıfırdan başlıyor; ilk 50 yıl ısınma
MU_REF = 3.0e10         # bütçedeki sabit μ


def element_region(lat, lon) -> np.ndarray:
    lat, lon = np.asarray(lat), np.asarray(lon)
    out = np.full(lat.shape, "diğer", dtype=object)
    for k, r in SUBREGIONS.items():
        m = (lat >= r.min_lat) & (lat < r.max_lat) & (lon >= r.min_lon) & (lon < r.max_lon) & (out == "diğer")
        out[m] = k
    return out


def physical_capacity(el: pd.DataFrame, crust) -> pd.DataFrame:
    """Eleman başına rijitlik, kabuk kalınlığı ve fiziksel kayma kapasitesi (m)."""
    mu = crust.rigidity_pa(el["lat"].to_numpy(), el["lon"].to_numpy(),
                           el["upper_depth_km"].to_numpy(), el["lower_depth_km"].to_numpy())
    W = el["fault_wid_km"].to_numpy() * 1e3
    dcap = 2 * DSIGMA_PA * W / (np.pi * mu)
    return pd.DataFrame({"rigidity_gpa": mu / 1e9,
                         "moho_km": crust.moho_depth_km(el["lat"].to_numpy(), el["lon"].to_numpy()),
                         "dcap_phys_m": dcap,
                         "m0_cap_Nm": mu * el["area_km2"].to_numpy() * 1e6 * dcap}, index=el.index)


def slip_m(deficit_nm, el):
    return np.asarray(deficit_nm) / (MU_REF * el["area_km2"].to_numpy() * 1e6)


def rupture_fills(records, el, dcap_phys) -> pd.DataFrame:
    """Her kayıtlı depremin kırıldığı elemanlardaki fiziksel doluluk ve aynı andaki dağılımda yeri."""
    reg = element_region(el["lat"].to_numpy(), el["lon"].to_numpy())
    sr = el["slip_rate_mm_yr"].to_numpy()
    fid = el["fault_id"].to_numpy()
    rows = []
    for r in records:
        fill_all = slip_m(r["deficit_before"], el) / dcap_phys
        # kırılan fay: boşalım ağırlığının en büyük payını alan fay. Doluluk yalnız onun
        # elemanlarında (komşu yavaş faylara yayılan Gauss kuyruğu ortalamayı seyreltmesin)
        main = pd.Series(r["w"]).groupby(fid[r["idx"]]).sum().idxmax()
        on = fid[r["idx"]] == main
        wm = r["w"][on] / r["w"][on].sum()
        f = float((wm * fill_all[r["idx"][on]]).sum())
        sr_w = float((wm * sr[r["idx"][on]]).sum())
        rg = element_region([r["lat"]], [r["lon"]])[0]
        same = (reg == rg) & (sr >= 1.0)
        pct = float((fill_all[same] < f).mean()) if same.any() else np.nan
        rows.append({"event_id": r["event_id"], "time": r["time"], "mw": r["mw"], "lat": r["lat"],
                     "lon": r["lon"], "bölge": rg, "fay": main, "kayma_hızı_ağırlıklı": sr_w,
                     "doluluk_fiziksel": f,
                     "aynı_bölgedeki_yüzdelik": pct,
                     "elemanlar": r["idx"], "ağırlık": r["w"]})
    return pd.DataFrame(rows)


CALIB_MIN_MW = 6.8    # fay parçasını baştan sona kıran (karakteristik) depremler
MIN_CALIB_SLIP = 2.0  # mm/yıl: yavaş faylarda (ve fay dışı depremlerde) birikim anlamlı değil


def _pool(y, groups, prior_mean, prior_var):
    """Normal-normal kısmi havuzlama (moment yöntemi): grup ortalamalarının sonsalı."""
    sigma2 = max(np.var(np.concatenate([v - v.mean() for v in groups.values() if len(v) > 1]), ddof=1)
                 if any(len(v) > 1 for v in groups.values()) else np.var(y), 0.05)
    means = np.array([v.mean() for v in groups.values()])
    ns = np.array([len(v) for v in groups.values()])
    tau2 = max(np.var(means, ddof=1) - np.mean(sigma2 / ns), 0.01) if len(means) > 1 else prior_var
    post = {}
    for k, v in groups.items():
        prec = len(v) / sigma2 + 1 / tau2
        post[k] = ((v.sum() / sigma2 + prior_mean / tau2) / prec, np.sqrt(1 / prec))
    return post, sigma2, tau2


def hierarchical_calibration(rf: pd.DataFrame, before_year: int, from_year: int = CALIB_FROM_YEAR) -> dict:
    """log(kırılma anındaki fiziksel doluluk) için fay → bölge → genel kısmi havuzlama.

    Yalnız ``before_year`` öncesi, 1950 sonrası ve kayma hızı ≥ 2 mm/yıl faylarda kırılan
    M≥6,8 depremler. Dönen: {"genel": (θ0, sd), "bölge:<ad>": (θ, sd), "fay:<id>": (θ, sd), ...}
    """
    d = rf[(rf["time"].dt.year >= from_year) & (rf["time"].dt.year < before_year)
           & (rf["kayma_hızı_ağırlıklı"] >= MIN_CALIB_SLIP) & (rf["doluluk_fiziksel"] > 1e-4)
           & (rf["mw"] >= CALIB_MIN_MW)]
    if len(d) < 3:
        return {"genel": (0.0, 1.0), "sigma": 1.0, "tau_bölge": 0.5, "n": int(len(d))}
    y = np.log(d["doluluk_fiziksel"].to_numpy())
    theta0 = float(y.mean())
    out = {"genel": (theta0, float(y.std(ddof=1) / np.sqrt(len(y)))), "n": int(len(y)),
           "olay_olay_sd": float(y.std(ddof=1))}
    rg = {k: y[d["bölge"].to_numpy() == k] for k in np.unique(d["bölge"])}
    post_r, s2, t2 = _pool(y, rg, theta0, 0.25)
    out.update({f"bölge:{k}": (float(m), float(s)) for k, (m, s) in post_r.items()})
    out["sigma"], out["tau_bölge"] = float(np.sqrt(s2)), float(np.sqrt(t2))
    # fay düzeyi: kendi bölgesinin sonsalına doğru çekilir
    for (reg_, f), g in d.groupby(["bölge", "fay"]):
        v = np.log(g["doluluk_fiziksel"].to_numpy())
        m0, _ = post_r[reg_]
        prec = len(v) / s2 + 1 / t2
        out[f"fay:{f}"] = (float((v.sum() / s2 + m0 / t2) / prec), float(np.sqrt(1 / prec)))
    return out


def calibrated_scale(el_region: np.ndarray, calib: dict, fault_ids=None) -> tuple[np.ndarray, np.ndarray]:
    """Eleman başına kapasite ölçeği exp(θ) ve log-sd (fay > bölge > genel)."""
    g = calib["genel"]
    th = np.empty(len(el_region))
    sd = np.empty(len(el_region))
    fault_ids = np.full(len(el_region), None) if fault_ids is None else fault_ids
    cache = {}
    for i, (r, f) in enumerate(zip(el_region, fault_ids)):
        key = (r, f)
        if key not in cache:
            cache[key] = calib.get(f"fay:{f}", calib.get(f"bölge:{r}", g))
        th[i], sd[i] = cache[key]
    return np.exp(th), sd


# ------------------------------------------------------------------ sütunlar ve zincirleme
def columns(el: pd.DataFrame, values: np.ndarray) -> pd.DataFrame:
    """Elemanları doğrultu sütunlarına topla (alan ağırlıklı doluluk)."""
    d = pd.DataFrame({"fault_id": el["fault_id"].to_numpy(), "i": el["i_strike"].to_numpy(),
                      "a": el["area_km2"].to_numpy(), "v": values * el["area_km2"].to_numpy(),
                      "lat": el["lat"].to_numpy(), "lon": el["lon"].to_numpy(),
                      "len": el["elem_len_km"].to_numpy()})
    g = d.groupby(["fault_id", "i"], sort=True)
    c = pd.DataFrame({"fill": g["v"].sum() / g["a"].sum(), "lat": g["lat"].mean(), "lon": g["lon"].mean(),
                      "len_km": g["len"].first()}).reset_index()
    return c


def column_graph(cols: pd.DataFrame, junction_km: float = 5.0):
    """Komşuluk listesi: aynı fayda i±1; fay uçları ``junction_km`` içindeyse birbirine bağlı."""
    n = len(cols)
    nb = [[] for _ in range(n)]
    fid = cols["fault_id"].to_numpy()
    for k in range(n - 1):
        if fid[k] == fid[k + 1]:
            nb[k].append(k + 1)
            nb[k + 1].append(k)
    first = cols.groupby("fault_id").head(1).index.to_numpy()
    last = cols.groupby("fault_id").tail(1).index.to_numpy()
    ends = np.r_[first, last]
    xy = np.column_stack([cols["lon"].to_numpy()[ends] * 111.195 * np.cos(np.radians(39)),
                          cols["lat"].to_numpy()[ends] * 111.195])
    tree = cKDTree(xy)
    for a, b in tree.query_pairs(junction_km):
        ia, ib = ends[a], ends[b]
        if fid[ia] != fid[ib]:
            nb[ia].append(ib)
            nb[ib].append(ia)
    return nb


def cascade(fill: np.ndarray, nb, seed: int, thr: float = 1.0, alpha: float = 0.5,
            forced: bool = True) -> np.ndarray:
    """Tohumdan başlayan zincirleme kırılma; kırılan sütunların maskesi.

    Kırılan sütun yükünü sıfırlar ve düşen yükün α/k'sını her komşusuna ekler
    (k: komşu sayısı). Eşiği aşan komşu kırılır. ``forced``: tohum eşiğin altında
    olsa da kırılır (bilinen bir depremin yayılımını yeniden üretmek için).
    """
    f = fill.astype(float).copy()
    broken = np.zeros(len(f), bool)
    if not forced and f[seed] < thr:
        return broken
    stack = [seed]
    while stack:
        k = stack.pop()
        if broken[k]:
            continue
        broken[k] = True
        drop = max(f[k], 0.0)
        f[k] = 0.0
        if nb[k]:
            share = alpha * drop / len(nb[k])
            for j in nb[k]:
                if not broken[j]:
                    f[j] += share
                    if f[j] >= thr:
                        stack.append(j)
    return broken


def cascade_potential(cols: pd.DataFrame, nb, seed_thr: float = 0.8, thr: float = 0.8,
                      alpha: float = 0.5) -> np.ndarray:
    """Her sütun için: içinden geçen en uzun olası zincirleme kırılmanın uzunluğu (km).

    Doluluğu ``seed_thr`` üstündeki her sütundan kırılma başlatılır.
    """
    fill = cols["fill"].to_numpy()
    L = cols["len_km"].to_numpy()
    best = np.zeros(len(cols))
    done_seed = np.zeros(len(cols), bool)
    for s in np.argsort(-fill):
        if fill[s] < seed_thr:
            break
        if done_seed[s]:
            continue
        b = cascade(fill, nb, s, thr=thr, alpha=alpha, forced=True)
        length = L[b].sum()
        best[b] = np.maximum(best[b], length)
        done_seed[b & (fill >= seed_thr)] = True  # aynı kümeden yeniden başlatma
    return best


def length_to_mw(L_km):
    L_km = np.asarray(L_km, float)
    return np.where(L_km > 0, (np.log10(np.maximum(L_km, 1e-3)) + 2.44) / 0.59, np.nan)
