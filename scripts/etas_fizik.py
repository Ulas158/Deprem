"""Fiziksel model, kümelenmeye (ETAS) bir şey katıyor mu?

1) anlik: küçük deprem aşamasının en iyi 5 ayarı × 3 tohumla Silivri ve Çınarcık yerel ağlarını
   2006-2024 boyunca koşturur; her 0,05 yılda fay boyunca (0,2 km'lik sütunlar) kilitli kutuların
   alan ağırlıklı doluluğunu kaydeder. Bir anlık görüntü yalnız kendinden önceki depremleri bilir.
   Her koşu bittiğinde kalıcı klasöre yazılır (kesintide kaldığı yerden sürer).
2) etas: fay boyunca 1 boyutlu uzay-zaman ETAS'ı ayar döneminde (2006-2017) en büyük olabilirlikle
   uydurur. Üç model karşılaştırılır, hepsi sınav döneminde (2018-2023) parametreleri donmuş olarak:
     E0  ETAS (arka plan sabit, ayar dönemi depremlerinden yumuşatılmış)
     E1  ETAS + fiziksel doluluk: arka plan × exp(beta·(g(x,t) - ort))
     E1k ETAS + aynalanmış doluluk (denetim: g fay boyunca ters çevrilmiş; zamanlama korunur,
         yer bilgisi bozulur. Zamanda döngüsel kaydırma sınav dönemi bilgisini ayara sızdırdığı
         için kullanılmadı)
   Puan: sınav olaylarının olay başına log-olabilirlik farkı (E1 - E0), 30 günlük blok
   önyüklemesiyle %95 aralık (kümelenme yüzünden olayları bağımsız saymaz).

Kullanım: python scripts/etas_fizik.py anlik --ag silivri --is 3
          python scripts/etas_fizik.py etas
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import optimize, special  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
import mekanik_ikinci_tur as M  # noqa: E402

CIKTI = oz.KALICI / "etas_fizik"
DT_ANLIK = 0.05
MC = 2.0                 # bütün kaynaklar için tam kabul edilen en küçük büyüklük (GFZ Mc <= 1,9)
T0, T_AYAR, T1 = 2006.0, 2018.0, 2024.0
AGLAR = ("silivri", "cinarcik")


# ------------------------------------------------------------------ 1) anlık görüntüler
def _profil(m, bins, nb):
    d = oz.doluluk(m)
    ok = m.vw & np.isfinite(d)
    w = np.where(ok, m.alan, 0.0)
    num = np.bincount(bins, weights=np.where(ok, d, 0.0) * w, minlength=nb)
    den = np.bincount(bins, weights=w, minlength=nb)
    return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan).astype(np.float32)


def anlik_kos(ag_ad, k, seed, th):
    f = CIKTI / f"anlik_{ag_ad}_{k}_{seed}.pkl"
    if f.exists():
        return f
    ag, ol, _ = M._hazirla("kucuk", ag_ad, True)
    a = M.AYAR["kucuk"]
    x = ag.cells.x_km.to_numpy()
    x0 = x.min() - 1e-6
    h = oz.YEREL_AGLAR[ag_ad][5]
    bins = np.floor((x - x0) / h).astype(int)
    nb = bins.max() + 1
    kayit = {"t": [], "g": []}
    asil = oz.model_kur

    def kur(*aa, **kw):
        m = asil(*aa, **kw)
        ilerle0 = m.ilerle
        sira = {"t": T0}

        def ilerle(t_son, *a2, **k2):
            while sira["t"] <= t_son + 1e-9 and sira["t"] <= T1:
                if m.t < sira["t"]:
                    ilerle0(sira["t"], *a2, **k2)
                if m.t >= T0 - 1e-9:
                    kayit["t"].append(sira["t"])
                    kayit["g"].append(_profil(m, bins, nb))
                sira["t"] = round(sira["t"] + DT_ANLIK, 6)
            if m.t < t_son:
                ilerle0(t_son, *a2, **k2)

        m.ilerle = ilerle
        return m

    oz.model_kur = kur
    try:
        kos = oz.Kosu(ag, th, ol, 900 + seed, a["isinma"], a["t_bas"], a["t_son_sinav"], a["ayar_bitis"])
        kos.calistir()
    finally:
        oz.model_kur = asil
    xs = x0 + (np.arange(nb) + 0.5) * h
    out = {"t": np.array(kayit["t"]), "g": np.vstack(kayit["g"]), "x": xs, "th": th, "k": k, "seed": seed}
    tmp = f.with_suffix(".tmp")
    tmp.write_bytes(pickle.dumps(out))
    tmp.replace(f)
    return f


def _anlik_is(arg):
    t0 = time.time()
    f = anlik_kos(*arg)
    print(f"{f.name} {time.time() - t0:.0f} s", flush=True)


def anlik(args):
    CIKTI.mkdir(parents=True, exist_ok=True)
    isler = []
    for ag_ad in args.ag:
        rs = [json.loads(s) for s in (oz.KALICI / f"ara_kucuk_{ag_ad}.jsonl").read_text().splitlines()]
        rs = sorted([r for r in rs if r["J"] > -1e8], key=lambda r: -r["J"])[:5]
        for k, r in enumerate(rs):
            for seed in range(3):
                isler.append((ag_ad, k, seed, r["th"]))
    from multiprocessing import get_context
    with get_context("fork").Pool(args.n_is) as p:
        list(p.imap_unordered(_anlik_is, isler))


# ------------------------------------------------------------------ 2) ETAS
def katalog(ag_ad):
    ag, ol, _ = M._hazirla("kucuk", ag_ad, True)
    d = pd.DataFrame([{"t": o.t, "mw": o.mw, "x": o.x} for o in ol if o.tip == "kucuk"])
    d = d[(d.mw >= MC) & (d.t >= T0) & (d.t < T1)].sort_values("t").reset_index(drop=True)
    x = ag.cells.x_km.to_numpy()
    return d, (x.min(), x.max())


def _uzay(r, s, q):
    """1B kalın kuyruklu çekirdek, R üzerinde normalize: c_q s^(2q-1) / (r^2+s^2)^q."""
    cq = special.gamma(q) / (np.sqrt(np.pi) * special.gamma(q - 0.5))
    return cq * s ** (2 * q - 1) / (r * r + s * s) ** q


def _uzay_kutle(a, b, s, q, n=200):
    """Çekirdeğin [a, b] aralığındaki kütlesi (sayısal)."""
    u = np.linspace(0, 1, n)
    r = a[:, None] + (b - a)[:, None] * u[None, :]
    return np.trapezoid(_uzay(r, s[:, None], q), r, axis=1)


class Etas:
    """lambda(t,x) = mu · b(x) · exp(beta·(g(t,x)-gbar)) / Z(t)
                     + sum_i K e^{alpha(M_i-Mc)} p'(t-t_i) f(x-x_i; s_i)
    Z(t): arka planın alan üstündeki ortalaması 1 olacak biçimde (beta yalnız yeri değiştirir,
    toplam arka plan hızını değiştirmez; böylece fizik yalnız "nerede" sorusuna katkı yapar).
    """

    AD = ["mu", "K", "alpha", "c", "p", "s0", "gamma", "q", "beta"]

    def __init__(self, d, dom, xg, tg, G, b_x, fizik=True):
        self.d, self.dom, self.xg, self.tg, self.G = d, dom, xg, tg, G
        self.b_x = b_x / np.trapezoid(b_x, xg)  # alan üstünde normalize arka plan biçimi
        self.fizik = fizik
        self.gbar = np.nanmean(G[(tg >= T0) & (tg < T_AYAR)])

    def _arka(self, beta, ti, xi):
        """Olay noktalarında normalize arka plan yoğunluğu (1/km) ve zaman ızgarasında Z."""
        Gc = np.nan_to_num(self.G - self.gbar, nan=0.0)
        e = self.b_x[None, :] * np.exp(beta * Gc)
        Z = np.trapezoid(e, self.xg, axis=1)
        j = np.clip(np.searchsorted(self.tg, ti, side="right") - 1, 0, len(self.tg) - 1)
        xb = np.clip(np.searchsorted(self.xg, xi), 0, len(self.xg) - 1)
        return e[j, xb] / Z[j]

    def nll(self, th, ta, tb, sadece_puan=False):
        mu, K, alpha, c, p, s0, gamma, q, beta = th
        if not self.fizik:
            beta = 0.0
        d = self.d
        t, mw, x = d.t.to_numpy(), d.mw.to_numpy(), d.x.to_numpy()
        hedef = np.flatnonzero((t >= ta) & (t < tb))
        s = s0 * 10 ** (gamma * (mw - MC))
        prod = K * np.exp(alpha * (mw - MC))
        lam = np.empty(len(hedef))
        bk = self._arka(beta, t[hedef], x[hedef])
        for n, j in enumerate(hedef):
            i = np.arange(j)
            i = i[t[i] > t[j] - 30.0]  # 30 yıldan eski tetiklemeler ihmal
            dtt = t[j] - t[i]
            omori = (p - 1) * c ** (p - 1) * (dtt + c) ** (-p)
            lam[n] = mu * bk[n] + np.sum(prod[i] * omori * _uzay(x[j] - x[i], s[i], q))
        loglam = np.log(np.maximum(lam, 1e-300))
        # integral: arka plan mu·(tb-ta) (Z normalizasyonu sayesinde); tetikleme analitik zaman ×
        # alan içi uzay kütlesi
        i = np.flatnonzero(t < tb)
        a_ = np.maximum(t[i], ta)
        F = lambda u: 1 - (1 + np.maximum(u, 0) / c) ** (1 - p)  # noqa: E731
        zt = F(tb - t[i]) - F(a_ - t[i])
        zx = _uzay_kutle(self.dom[0] - x[i], self.dom[1] - x[i], s[i], q)
        integ = mu * (tb - ta) + np.sum(prod[i] * zt * zx)
        if sadece_puan:
            return loglam, integ, t[hedef]
        return -(loglam.sum() - integ)


ALT = np.array([1e-2, 1e-4, 0.0, 1e-5, 1.01, 0.05, 0.0, 0.6, -20.0])
UST = np.array([1e3, 1e2, 4.0, 0.1, 2.5, 5.0, 0.8, 3.0, 20.0])
BAS = np.array([5.0, 0.3, 1.5, 1e-3, 1.1, 0.5, 0.3, 1.5, 0.0])


def uydur(e: Etas, fizik: bool, tekrar=6, seed=0):
    rng = np.random.default_rng(seed)
    lo, hi = ALT.copy(), UST.copy()
    if not fizik:
        lo[-1] = hi[-1] = 0.0
    en = None
    for r in range(tekrar):
        x0 = BAS.copy() if r == 0 else lo + (hi - lo) * rng.uniform(0.1, 0.9, len(lo))
        if r > 0:
            for k in (0, 1, 3):  # log ölçekli
                x0[k] = np.exp(np.log(lo[k]) + (np.log(hi[k]) - np.log(lo[k])) * rng.uniform(0.2, 0.8))
        if not fizik:
            x0[-1] = 0.0
        try:
            res = optimize.minimize(lambda th: e.nll(th, T0, T_AYAR), x0, method="L-BFGS-B",
                                    bounds=list(zip(lo, hi)))
        except FloatingPointError:
            continue
        if en is None or res.fun < en.fun:
            en = res
    return en


def blok_onyukleme(fark, t, n=4000, blok=30 / 365.25, seed=0):
    b = np.floor((t - t.min()) / blok).astype(int)
    ub = np.unique(b)
    toplam = np.array([fark[b == u].sum() for u in ub])
    rng = np.random.default_rng(seed)
    bs = np.array([toplam[rng.integers(0, len(ub), len(ub))].sum() for _ in range(n)])
    return np.percentile(bs, [2.5, 97.5]), float(np.mean(bs <= 0))


def etas(args):
    np.seterr(all="ignore")
    sonuc = []
    for ag_ad in AGLAR:
        d, dom = katalog(ag_ad)
        fs = sorted(CIKTI.glob(f"anlik_{ag_ad}_*.pkl"))
        assert fs, f"{ag_ad}: anlık görüntü yok"
        A = [pickle.loads(f.read_bytes()) for f in fs]
        tg, xg = A[0]["t"], A[0]["x"]
        n = min(len(a["t"]) for a in A)
        tg = tg[:n]
        G = np.nanmean(np.stack([a["g"][:n] for a in A]), 0)
        # arka plan biçimi: ayar dönemi olaylarının 1 km'lik Gauss yumuşatması + düzgün taban
        tr = d[d.t < T_AYAR]
        bx = np.exp(-0.5 * ((xg[:, None] - tr.x.to_numpy()[None, :]) / 1.0) ** 2).sum(1) + 0.05 * len(tr) / len(xg)
        # denetim: doluluk fay boyunca aynalanmış (zamanlama korunur, yer bozulur)
        Gk = G[:, ::-1].copy()
        # doluluk profili yaklaşık simetrik (ağ uçlarında yüksek: sınır etkisi), bu yüzden aynalama
        # zayıf bir denetim. Ek iki ayrıştırma:
        #   E1s: yalnız durağan profil (ayar dönemi zaman ortalaması) -> uç etkisi + sabit desen
        #   E1d: yalnız zamanla değişen kısım (g - durağan profil) -> depremlerle gelen gerilme değişimi
        ay = (tg >= T0) & (tg < T_AYAR)
        Gs = np.broadcast_to(np.nanmean(G[ay], 0), G.shape).copy()
        Gd = G - Gs
        modeller = {"E0": (G, False), "E1": (G, True), "E1k": (Gk, True), "E1s": (Gs, True), "E1d": (Gd, True)}
        sat = {"ag": ag_ad, "n_ayar": int((d.t < T_AYAR).sum()), "n_sinav": int((d.t >= T_AYAR).sum()),
               "kosu": len(A)}
        puan = {}
        for ad, (g, fz) in modeller.items():
            e = Etas(d, dom, xg, tg, g, bx, fizik=fz)
            r = uydur(e, fz)
            th = r.x
            ll, integ, tt = e.nll(th, T_AYAR, T1, sadece_puan=True)
            puan[ad] = (ll, integ, tt)
            sat[f"{ad}_param"] = dict(zip(Etas.AD, np.round(th, 4).tolist()))
            sat[f"{ad}_ayar_ll"] = float(-r.fun)
            sat[f"{ad}_sinav_ll"] = float(ll.sum() - integ)
        for ad in ("E1", "E1k", "E1s", "E1d"):
            # aynı olaylar, olay başına fark; integral farkı olay sayısına bölünüp dağıtılır
            f = (puan[ad][0] - puan["E0"][0]) - (puan[ad][1] - puan["E0"][1]) / len(puan["E0"][0])
            ga, p = blok_onyukleme(f, puan["E0"][2])
            sat[f"{ad}_eksi_E0_olay_basi"] = float(f.mean())
            sat[f"{ad}_eksi_E0_toplam"] = float(f.sum())
            sat[f"{ad}_toplam_GA95"] = ga.tolist()
            sat[f"{ad}_p_iyilesme_yok"] = p
        sonuc.append(sat)
        print(json.dumps(sat, ensure_ascii=False, indent=1), flush=True)
    (ROOT / "reports" / "tables" / "etas_fizik.json").write_text(json.dumps(sonuc, ensure_ascii=False, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("asama", choices=["anlik", "etas"])
    ap.add_argument("--ag", nargs="*", default=list(AGLAR))
    ap.add_argument("--is", dest="n_is", type=int, default=3)
    args = ap.parse_args()
    {"anlik": anlik, "etas": etas}[args.asama](args)


if __name__ == "__main__":
    main()
