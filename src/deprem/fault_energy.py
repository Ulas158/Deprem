"""Fay fay biriken enerji, küçük depremlerden gerilme göstergesi ve komşu fay zinciri.

Ulas'ın sorusu (25 Kasım 2026): "Deprem küçük olsa bile, o fayın alabildiği enerji
bulunabilir mi? Büyük deprem için boydan boya yarık mı gerekiyor, yoksa biri sınırdaysa
yeter mi?" Bu modül her fay segmenti için, her an **yalnız o andan önceki verilerle**
üç zamanla değişen gösterge üretir ve bunları sıralı fay modeline (``faultwalk``) ek
girdi olarak verir:

1. **Moment açığı oranı R** ("fayda kalan enerji / fayın kaldırabileceği"):

       açık_k(t) = μ · A_k · v_k · (t − 1700) − Σ_{i: t_i < t} M0_i · pay_ik
       kapasite_k = μ · A_k · D_kap,k,   D_kap = 2 Δσ W / (π μ)  (Knopoff 1958)
       R_k(t) = açık_k(t) / kapasite_k

   v: GEM kayma hızı (mm/yıl), A: fay alanı, W: genişlik, Δσ = 10 MPa (kayacın
   taşıyabileceği büyük gerilme düşümü; ``faults.capacity`` ile aynı). Boşalım, **her
   büyüklükteki** depremden (küçükler dahil) hesaplanır: deprem momenti, kırık yarı
   boyu + 5 km içinde kalan segmentlere iz uzunluğu × genişlik oranında dağıtılır
   (moment korunur). 1700'deki başlangıç açığı bilinmediği için 0 alınır; R göreli bir
   ölçüdür. Biriken elastik enerji kabaca E ≈ (Δσ_ort / 2μ) · açık (Kanamori 1977;
   Δσ_ort = 3 MPa).

2. **Yerel b-değeri farkı** (küçük/büyük deprem oranı; Schorlemmer vd. 2005): segmentin
   30 km yakınındaki son 5 yılın depremleriyle b, aynı pencerede tüm Türkiye kutusunun b'si
   ile karşılaştırılır: z = (b_bölge − b_yerel) / σ. + değer = yerelde orta boy depremler
   görece fazla = daha yüksek gerilme (literatürdeki hipotez). Tamamlılık büyüklüğü Mc
   her yıl tüm kutu için en büyük eğrilik + 0,2 ile bulunur; pencere için pencere
   yıllarının en büyük Mc'si kullanılır. ≥ 50 olay yoksa 0 (etkisiz). Küçük deprem kaydı
   1998 öncesi seyrek olduğu için bu gösterge pratikte 1998 sonrası doludur.

3. **Yüklü komşu zinciri**: bir segmentin ucu diğerinin izine 5 km içinde yaklaşıyorsa
   (uç uca ya da T-kavşak) ikisi komşudur. R ≥ eşik
   (0,5) olan segmentlerin bağlı kümesi, birlikte kırılabilecek yüklü zinciri verir.
   Segmentin girdisi log10(zincir uzunluğu / 100 km); segment yüklü değilse zincir
   uzunluğu 5 km sayılır. Zincir uzunluğundan olası büyüklük Wells & Coppersmith
   (1994) ile çıkarılır.

Gerilme düşümü (stress drop) ayrı bir göstergedir ama tek tek depremlerin dalga
biçimi spektrumlarından hesaplanır; katalogdan çıkarılamaz. Bu sürümde yoktur.

Hiçbir gösterge doğrulanmış kabul edilmez; katkıları ``scripts/faultwalk.py`` ile ETAS'a
ve ana modele karşı, aynı sıralı düzende ölçülür.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .faultwalk_data import T0, WalkData, days, load_fault_table, rupture_len_km, xy_km, year_start_days
from .faults.moment import mw_to_m0
from .region import TURKEY
from .stats.gr import b_value_mle, mc_maxc

log = logging.getLogger(__name__)

MU = 3.0e10
DSIGMA_CAP = 10e6          # kapasite için gerilme düşümü (Pa)
DSIGMA_AVG = 3e6           # enerji dönüşümü için ortalama gerilme düşümü (Pa)
JUNCTION_KM = 5.0
CHAIN_THR = 0.5
CHAIN_FLOOR_KM = 5.0
B_RADIUS_KM = 30.0
B_WINDOW_YR = 5.0
B_MIN_N = 50
B_FROM_YEAR = 1905
TRACE_STEP_KM = 2.0


@dataclass
class EnergyData:
    seg: pd.DataFrame          # segment tablosu (f) + uzunluk, genişlik, kapasite
    nb: list                   # segment komşuluk listesi
    tau: np.ndarray            # adım sınırları (WalkData.tau)
    deficit: np.ndarray        # (S, F) moment açığı (N·m), adım başında
    R: np.ndarray              # (S, F) açık / kapasite
    bz: np.ndarray             # (S, F) b-değeri z farkı (veri yoksa 0)
    bn: np.ndarray             # (S, F) b hesabındaki olay sayısı
    chain_km: np.ndarray       # (S, F) yüklü zincir uzunluğu (yüklü değilse 0)
    released_small: np.ndarray # (F,) M<5,5 depremlerin boşalttığı toplam moment
    released_all: np.ndarray   # (F,) tüm depremlerin boşalttığı toplam moment
    km: "KmData | None" = None # 1 km çözünürlüklü sürüm


# ------------------------------------------------------------------ yardımcılar
def _trace_points(coords, step=TRACE_STEP_KM):
    c = np.asarray(coords, float)
    pts = []
    for a, b in zip(c[:-1], c[1:]):
        n = max(int(np.hypot((b[0] - a[0]) * 85, (b[1] - a[1]) * 111) / step), 1)
        s = np.linspace(0, 1, n, endpoint=False)
        pts.append(a + s[:, None] * (b - a))
    pts.append(c[-1:])
    return np.concatenate(pts)


def segment_graph(f: pd.DataFrame, junction_km: float = JUNCTION_KM) -> list:
    """Bir segmentin ucu başka bir segmentin izine (ucuna ya da ortasına, T-kavşak)
    ``junction_km`` içinde yaklaşıyorsa ikisi komşudur."""
    pts, owner = [], []
    for k, c in enumerate(f["coords"]):
        p = _trace_points(c, step=1.0)
        pts.append(p)
        owner.append(np.full(len(p), k))
    pts = np.concatenate(pts)
    owner = np.concatenate(owner)
    px, py = xy_km(pts[:, 1], pts[:, 0])
    tree = cKDTree(np.column_stack([px, py]))
    nb = [set() for _ in range(len(f))]
    for k, c in enumerate(f["coords"]):
        c = np.asarray(c, float)
        ends = np.array([c[0], c[-1]])
        ex, ey = xy_km(ends[:, 1], ends[:, 0])
        for h in tree.query_ball_point(np.column_stack([ex, ey]), junction_km):
            for j in np.unique(owner[h]):
                if j != k:
                    nb[k].add(int(j))
                    nb[int(j)].add(k)
    return [sorted(x) for x in nb]


def components(mask: np.ndarray, nb: list, length: np.ndarray) -> np.ndarray:
    """mask'teki segmentlerin bağlı kümelerinin toplam uzunluğu (mask dışı 0)."""
    out = np.zeros(mask.size)
    seen = np.zeros(mask.size, bool)
    for s in np.nonzero(mask)[0]:
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            k = stack.pop()
            comp.append(k)
            for j in nb[k]:
                if mask[j] and not seen[j]:
                    seen[j] = True
                    stack.append(j)
        out[comp] = length[comp].sum()
    return out


def chain_component_ids(mask: np.ndarray, nb: list) -> np.ndarray:
    lab = np.full(mask.size, -1)
    c = 0
    for s in np.nonzero(mask)[0]:
        if lab[s] >= 0:
            continue
        stack = [s]
        lab[s] = c
        while stack:
            k = stack.pop()
            for j in nb[k]:
                if mask[j] and lab[j] < 0:
                    lab[j] = c
                    stack.append(j)
        c += 1
    return lab


def length_to_mw(L_km):
    L_km = np.asarray(L_km, float)
    return np.where(L_km > 0, (np.log10(np.maximum(L_km, 1e-3)) + 2.44) / 0.59, np.nan)


def load_all_events(root: Path, end: str = "2026-01-01") -> pd.DataFrame:
    """Uzun katalog, her büyüklük (1700-1899 tarihsel M≥~5, 1900 sonrası aletsel)."""
    p = Path(root) / "data" / "processed" / "turkey_catalog_long.csv.gz"
    cat = pd.read_csv(p, low_memory=False)
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    mcol = "mw_best" if "mw_best" in cat else "mag"
    cat["mw"] = cat[mcol].fillna(cat["mag"])
    cat = cat[(cat["time"] >= T0) & (cat["time"] < pd.Timestamp(end, tz="UTC")) & cat["mw"].notna()
              & TURKEY.contains(cat)].copy()
    cat = cat.sort_values("time").reset_index(drop=True)
    cat["t"] = days(cat["time"])
    return cat


# ------------------------------------------------------------------ moment bütçesi
def release_matrix(ev: pd.DataFrame, f: pd.DataFrame, width_km: np.ndarray):
    """Her olayın momentini segmentlere dağıt. Dönen: (olay indeksi, segment, M0 payı) dizileri."""
    pts, owner = [], []
    for k, c in enumerate(f["coords"]):
        p = _trace_points(c)
        pts.append(p)
        owner.append(np.full(len(p), k))
    pts = np.concatenate(pts)
    owner = np.concatenate(owner)
    px, py = xy_km(pts[:, 1], pts[:, 0])
    tree = cKDTree(np.column_stack([px, py]))
    ex, ey = xy_km(ev["latitude"].to_numpy(), ev["longitude"].to_numpy())
    mw = ev["mw"].to_numpy()
    rr = np.maximum(rupture_len_km(mw) / 2, 1.0) + 5.0
    m0 = mw_to_m0(mw)
    I, K, V = [], [], []
    hits = tree.query_ball_point(np.column_stack([ex, ey]), rr)
    for i, h in enumerate(hits):
        if not h:
            continue
        ks, cnt = np.unique(owner[h], return_counts=True)
        area = cnt * TRACE_STEP_KM * width_km[ks]          # km²
        share = area / area.sum()
        I.append(np.full(ks.size, i))
        K.append(ks)
        V.append(m0[i] * share)
    if not I:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    return np.concatenate(I), np.concatenate(K), np.concatenate(V)


def deficit_series(load: np.ndarray, t_ev: np.ndarray, I, K, V, tau: np.ndarray, clip: bool = True) -> np.ndarray:
    """(S, F) adım başında moment açığı, F birim (segment ya da 1 km eleman) için.

    açık(t) = açık(önceki deprem sonrası) + yük hızı × geçen süre − boşalım.
    ``clip``: bir deprem fayda birikmiş olandan fazlasını boşaltamaz; açık 0'ın altına
    inmez. Böylece büyük bir kırılmadan sonra saat sıfırdan başlar ("son kırılmadan beri
    biriken"), küçük depremler ise yalnız kendi momentleri kadar düşer.

    Yalnız t ≤ tau[s] olan depremler adım s'ye girer. Adım s, olay j ile başlıyorsa
    (tau[s] = t_j) j'nin boşalımı o adıma girer; j'nin kendi puanı ise bir önceki adımın
    değeriyle hesaplanır, bu yüzden ileriye bakış yoktur.
    """
    F = load.size
    S = tau.size - 1
    ts = tau[:-1]
    D = (load[None, :] * (ts[:, None] / 365.25)).astype(np.float32)
    if I.size == 0:
        return D
    te = t_ev[I]
    o = np.lexsort((te, K))
    Ks, ts_e, Vs = K[o], te[o], V[o]
    bounds = np.searchsorted(Ks, np.arange(F + 1))
    for k in range(F):
        a, b = bounds[k], bounds[k + 1]
        if a == b:
            continue
        tt, vv = ts_e[a:b], Vs[a:b]
        after = np.empty(b - a)
        d, tp = 0.0, 0.0
        lk = load[k] / 365.25
        for n in range(b - a):
            d = d + lk * (tt[n] - tp) - vv[n]
            if clip and d < 0:
                d = 0.0
            after[n] = d
            tp = tt[n]
        idx = np.searchsorted(tt, ts, side="right") - 1
        has = idx >= 0
        D[has, k] = after[idx[has]] + lk * (ts[has] - tt[idx[has]])
    return D


def deficit_matrix(ev, f, tau, I, K, V, clip: bool = True):
    """Segment başına (S, F) moment açığı."""
    sr = f["slip_rate_mm_yr"].fillna(0.2).to_numpy() * 1e-3
    load = MU * f["area_km2"].to_numpy() * 1e6 * sr                 # N·m / yıl
    return deficit_series(load, ev["t"].to_numpy(), I, K, V, tau, clip=clip)


# ------------------------------------------------------------------ 1 km elemanlar
@dataclass
class KmData:
    el: pd.DataFrame           # eleman tablosu: seg, sıra, lat, lon, genişlik, yük hızı, kapasite
    R: np.ndarray              # (S, E) float32 doluluk
    chain_km: np.ndarray       # (S, E) float32 elemanın içinde olduğu yüklü kesintisiz parçanın boyu
    n_release: np.ndarray      # (E,) elemana düşen deprem sayısı


def km_elements(f: pd.DataFrame, step: float = 1.0) -> pd.DataFrame:
    rows = []
    sr = f["slip_rate_mm_yr"].fillna(0.2).to_numpy()
    W = f["width_km"].to_numpy()
    for k, c in enumerate(f["coords"]):
        p = _trace_points(c, step=step)
        if len(p) > 1:
            p = 0.5 * (p[1:] + p[:-1])            # eleman orta noktaları
        for n, (lon, lat) in enumerate(p):
            rows.append((k, n, lat, lon))
    el = pd.DataFrame(rows, columns=["seg", "pos", "lat", "lon"])
    el["width_km"] = W[el["seg"]]
    el["slip_rate"] = sr[el["seg"]]
    el["len_km"] = step
    el["load_m0_yr"] = MU * el["len_km"] * el["width_km"] * 1e6 * el["slip_rate"] * 1e-3
    el["cap_m0"] = MU * el["len_km"] * el["width_km"] * 1e6 * (2 * DSIGMA_CAP * el["width_km"] * 1e3 / (np.pi * MU))
    return el


def km_graph(el: pd.DataFrame, junction_km: float = JUNCTION_KM):
    """Seyrek komşuluk: aynı segmentte ardışık elemanlar; segment uç elemanı başka segmentin
    en yakın elemanına (≤ junction_km) bağlanır."""
    from scipy.sparse import coo_matrix
    seg = el["seg"].to_numpy()
    E = len(el)
    a = np.nonzero(seg[1:] == seg[:-1])[0]
    rows, cols = [a], [a + 1]
    x, y = xy_km(el["lat"].to_numpy(), el["lon"].to_numpy())
    tree = cKDTree(np.column_stack([x, y]))
    first = np.r_[0, np.nonzero(seg[1:] != seg[:-1])[0] + 1]
    last = np.r_[first[1:] - 1, E - 1]
    for e in np.r_[first, last]:
        for j in tree.query_ball_point([x[e], y[e]], junction_km):
            if seg[j] != seg[e]:
                rows.append(np.array([e]))
                cols.append(np.array([j]))
    r, c = np.concatenate(rows), np.concatenate(cols)
    A = coo_matrix((np.ones(r.size), (r, c)), shape=(E, E)).tocsr()
    return (A + A.T).tocsr()


def km_release(ev: pd.DataFrame, el: pd.DataFrame):
    x, y = xy_km(el["lat"].to_numpy(), el["lon"].to_numpy())
    tree = cKDTree(np.column_stack([x, y]))
    ex, ey = xy_km(ev["latitude"].to_numpy(), ev["longitude"].to_numpy())
    mw = ev["mw"].to_numpy()
    rr = np.maximum(rupture_len_km(mw) / 2, 1.0) + 5.0
    m0 = mw_to_m0(mw)
    area = (el["len_km"] * el["width_km"]).to_numpy()
    I, K, V = [], [], []
    for i, h in enumerate(tree.query_ball_point(np.column_stack([ex, ey]), rr)):
        if not h:
            continue
        h = np.asarray(h)
        w = area[h] / area[h].sum()
        I.append(np.full(h.size, i))
        K.append(h)
        V.append(m0[i] * w)
    return np.concatenate(I), np.concatenate(K), np.concatenate(V)


def build_km(ev: pd.DataFrame, f: pd.DataFrame, tau: np.ndarray) -> KmData:
    from scipy.sparse.csgraph import connected_components
    el = km_elements(f)
    log.info("1 km elemanlar: %d", len(el))
    I, K, V = km_release(ev, el)
    D = deficit_series(el["load_m0_yr"].to_numpy(), ev["t"].to_numpy(), I, K, V, tau)
    R = (D / el["cap_m0"].to_numpy()[None, :]).astype(np.float32)
    del D
    A = km_graph(el)
    S = tau.size - 1
    chain = np.zeros(R.shape, np.float32)
    elen = el["len_km"].to_numpy()
    for s in range(S):
        m = np.nonzero(R[s] >= CHAIN_THR)[0]
        if m.size == 0:
            continue
        _, lab = connected_components(A[m][:, m], directed=False)
        tot = np.bincount(lab, elen[m])
        chain[s, m] = tot[lab]
    return KmData(el=el, R=R, chain_km=chain, n_release=np.bincount(K, minlength=len(el)))


# ------------------------------------------------------------------ b-değeri
def yearly_mc(ev: pd.DataFrame, y0: int = 1900, y1: int = 2026) -> dict:
    out = {}
    yr = ev["time"].dt.year.to_numpy()
    m = ev["mw"].to_numpy()
    for y in range(y0, y1):
        x = m[yr == y]
        out[y] = mc_maxc(x) if x.size >= 30 else 5.0
    return out


def bvalue_matrix(ev: pd.DataFrame, f: pd.DataFrame, tau: np.ndarray):
    """(S, F) yerel b-değeri z farkı ve olay sayısı; yalnız tau[s]'ten önceki 5 yıl."""
    F = len(f)
    S = tau.size - 1
    bz = np.zeros((S, F))
    bn = np.zeros((S, F), int)
    ev = ev[ev["time"].dt.year >= 1900].reset_index(drop=True)
    mcy = yearly_mc(ev)
    t_ev = ev["t"].to_numpy()
    m_ev = ev["mw"].to_numpy()
    ex, ey = xy_km(ev["latitude"].to_numpy(), ev["longitude"].to_numpy())
    tree = cKDTree(np.column_stack([ex, ey]))
    near = []
    for c in f["coords"]:
        p = _trace_points(c)
        px, py = xy_km(p[:, 1], p[:, 0])
        h = tree.query_ball_point(np.column_stack([px, py]), B_RADIUS_KM)
        idx = np.unique(np.concatenate([np.asarray(x, int) for x in h])) if len(h) else np.zeros(0, int)
        near.append(np.sort(idx))                                  # zaman sıralı (ev zamana göre sıralı)
    win = B_WINDOW_YR * 365.25
    t0b = year_start_days(B_FROM_YEAR)
    for s in range(S):
        t = tau[s]
        if t < t0b:
            continue
        y = int(1700 + t / 365.25)                                 # kaba yıl (yalnız Mc seçimi için)
        ys = [yy for yy in range(int(y - B_WINDOW_YR), y + 1) if yy in mcy]
        if not ys:
            continue
        # cari yılın Mc'si yıl bitmeden bilinmez: bir önceki yılınki kullanılır
        mc = max(mcy.get(yy if yy < y else y - 1, 5.0) for yy in ys)
        lo, hi = np.searchsorted(t_ev, [t - win, t], side="left")
        if hi - lo < B_MIN_N:
            continue
        mreg = m_ev[lo:hi]
        breg, sreg, nreg = b_value_mle(mreg, mc)
        if not np.isfinite(breg) or nreg < 200:
            continue
        for k in range(F):
            nk = near[k]
            a, b = np.searchsorted(nk, [lo, hi], side="left")
            if b - a < B_MIN_N:
                continue
            bl, sl, nl = b_value_mle(m_ev[nk[a:b]], mc)
            if nl < B_MIN_N or not np.isfinite(bl):
                continue
            bz[s, k] = np.clip((breg - bl) / np.sqrt(sl ** 2 + sreg ** 2), -4, 4)
            bn[s, k] = nl
    return bz, bn


# ------------------------------------------------------------------ toplu hesap
def build_energy(root: Path, W: WalkData, cache: Path | None = None) -> EnergyData:
    if cache is not None and Path(cache).exists():
        import pickle
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    f = load_fault_table(root)
    width = f["width_km"].to_numpy()
    f["cap_slip_m"] = 2 * DSIGMA_CAP * width * 1e3 / (np.pi * MU)
    f["cap_m0"] = MU * f["area_km2"].to_numpy() * 1e6 * f["cap_slip_m"]
    nb = segment_graph(f)
    ev = load_all_events(root)
    log.info("Enerji: %d segment, %d olay (her büyüklük)", len(f), len(ev))
    I, K, V = release_matrix(ev, f, width)
    D = deficit_matrix(ev, f, W.tau, I, K, V)
    R = D / f["cap_m0"].to_numpy()[None, :]
    small = ev["mw"].to_numpy()[I] < 5.5
    rel_small = np.bincount(K[small], V[small], minlength=len(f))
    rel_all = np.bincount(K, V, minlength=len(f))
    log.info("b-değeri pencereleri hesaplanıyor")
    bz, bn = bvalue_matrix(ev, f, W.tau)
    L = f["length_km"].to_numpy()
    S = W.tau.size - 1
    chain = np.zeros((S, len(f)))
    for s in range(S):
        chain[s] = components(R[s] >= CHAIN_THR, nb, L)
    log.info("1 km elemanlar hesaplanıyor")
    km = build_km(ev, f, W.tau)
    E = EnergyData(seg=f, nb=nb, tau=W.tau, deficit=D, R=R,
                   bz=bz, bn=bn, chain_km=chain, released_small=rel_small, released_all=rel_all, km=km)
    if cache is not None:
        import pickle
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "wb") as fh:
            pickle.dump(E, fh)
    return E


def covariates(W: WalkData, E: EnergyData) -> dict:
    """Modele girecek (S, U) girdiler. Fay dışı birimlerde 0."""
    U = len(W.units)
    S = W.tau.size - 1
    seg = W.units["seg"].to_numpy()
    fay = W.units["kind"].to_numpy() == "fay"
    out = {k: np.zeros((S, U)) for k in ("deficit", "bval", "chain")}
    sk = seg[fay]
    out["deficit"][:, fay] = np.clip(np.log(np.maximum(E.R[:, sk], 1e-3)), -3.0, 1.5)
    out["bval"][:, fay] = E.bz[:, sk]
    ch = np.maximum(E.chain_km[:, sk], CHAIN_FLOOR_KM)
    out["chain"][:, fay] = np.log10(ch / 100.0)
    if E.km is not None:                                   # 1 km sürümü: segmentin elemanları üzerinden
        el_seg = E.km.el["seg"].to_numpy()
        o = np.argsort(el_seg, kind="stable")
        b = np.searchsorted(el_seg[o], np.arange(len(E.seg) + 1))
        dk = np.zeros((S, U))
        ck = np.zeros((S, U))
        for u in np.nonzero(fay)[0]:
            e = o[b[seg[u]]:b[seg[u] + 1]]
            Rm = np.maximum(E.km.R[:, e], 0).mean(1)
            dk[:, u] = np.clip(np.log(np.maximum(Rm, 1e-3)), -3.0, 1.5)
            ck[:, u] = np.log10(np.maximum(E.km.chain_km[:, e].max(1), CHAIN_FLOOR_KM) / 100.0)
        out["deficit_km"], out["chain_km"] = dk, ck
    return out


def energy_joule(deficit_nm):
    """Moment açığından kaba elastik enerji (J): E ≈ Δσ_ort / (2μ) · M0."""
    return DSIGMA_AVG / (2 * MU) * np.asarray(deficit_nm, float)
