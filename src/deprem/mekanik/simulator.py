"""Fay kırılma simülatörü: yükleme, yavaş kayma, sünek akma ve deprem.

Fizik (her hücre için, birimler Pa, m, yıl):

* **Yükleme.** Levha hareketi "geri kayma" yöntemiyle verilir (Savage 1983):
  dS_i/dt = Σ_j K_ij (v_j − v_pl,j). Model dışındaki fay (doğu ve batı uçlar,
  40 km altı) levha hızıyla sürekli kayıyor kabul edilir. Tamamen kilitli
  hücrede v = 0 olduğu için gerilme her yıl biraz artar; sürekli kayan
  hücrede v = v_pl olduğu için hiç artmaz.
* **Coulomb gerilmesi.** S = τ + μ' σn (σn çekme pozitif). Bir hücrenin
  kayması komşularda hem τ hem σn değiştirir; bir deprem diğerini bu yolla
  yaklaştırır ya da uzaklaştırır.
* **Yavaş kayma (krip).** Sıcaklığın sürtünmeyi hıza göre güçlendiren
  yaptığı hücreler (sığ < ~100 °C ve derin 350-450 °C) deprem üretmez,
  kararlı sürtünmeyle kayar: v = v_pl exp((S − S_ref)/A), A = (a−b) σ_eff
  (Marone 1998). Bir deprem bu hücrelere yük bindirirse hızlanır (deprem
  sonrası art kayma), gerilme gölgesinde kalırsa yavaşlar.
* **Sünek akma.** 20-40 km arası ve 450 °C üstü hücreler doğrusal ağdalı
  (viskoz) akar: v = v_pl + (S − S_ref)/(k_ii t_M), t_M sıcaklıkla değişen
  Maxwell zamanı. Büyük depremden sonra alttaki katmanın yıllarca süren
  "ittirmesi" buradan gelir.
* **Kırılma.** Hıza göre zayıflayan hücreler kilitlidir. S, statik
  dayanıma ulaşınca hücre kayar ve gerilmesi dinamik düzeye iner. Kaymanın
  aktardığı gerilme komşuları da kırabilir (çığ). Kırılma sürerken dinamik
  zayıflama ``eta_dyn`` komşuların eşiğini düşürür (sarsıntının ilerleyen
  kırığın ucunda geçici gerilme yığması; Rundle vd. 2006, Virtual Quake).
  Kırık kendiliğinden, gerilmesi yetmeyen yerde durur.

Bu bir yarı-statik hücresel simülatördür (RSQSim, Virtual Quake ailesi),
tam dinamik dalga çözümü değildir. Sünek akma, tam 3B ağdalı-elastik yarı
uzay yerine fay altındaki sünek makaslama kuşağıyla temsil edilir.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

from . import termal

MU = 3.0e10
YIL = 365.25 * 86400.0


@dataclass
class Parametreler:
    mu_eff: float = 0.4            # efektif sürtünme (Coulomb)
    dtau0: float = 4.0e6           # Pa, ortalama gerilme düşümü (derin kısım)
    sigma_h: float = 0.35          # dayanım heterojenliği (lognormal)
    eta_dyn: float = 0.5           # dinamik zayıflama (0: yok, 1: çok güçlü)
    asiri: float = 0.1             # kayma aşımı (overshoot)
    kayma_ref: float = 1.0         # m; komşu kayması bu düzeye ulaşınca dinamik zayıflama tam etkili
    sigma_grad: float = 2.5e6      # Pa/km, efektif normal gerilme artışı (yüksek gözenek basıncı)
    sigma_cap: float = 60e6        # Pa
    tM_ref: float = 20.0           # yıl, 500 °C'de Maxwell zamanı
    krip_carpan: float = 1.0       # (a-b)σ ölçeği
    isi_akisi: np.ndarray | float = 0.075  # W/m², hücre başına veya sabit
    v_pl: np.ndarray | float = 0.020       # m/yıl, hücre başına veya sabit
    sig_kirilma_mw: float = 6.8    # "büyük deprem" eşiği (kayıt için)
    mu_kayma: float = 3.0e10       # Pa, kayma modülü (moment için; ikinci tur 33 GPa)
    seed: int = 0
    extra: dict = field(default_factory=dict)


class Model:
    def __init__(self, cells: pd.DataFrame, Kt: np.ndarray, Kn: np.ndarray, p: Parametreler,
                 ek_krip: np.ndarray | None = None, K: np.ndarray | None = None):
        self.c = cells.reset_index(drop=True)
        self.p = p
        n = len(self.c)
        self.n = n
        if K is None:
            K = np.asfortranarray((Kt + np.float32(p.mu_eff) * Kn).astype(np.float32))
        self.K = K
        self.kself = -np.diag(self.K).astype(float)
        z = self.c["z_km"].to_numpy()
        q = np.broadcast_to(np.asarray(p.isi_akisi, float), (n,))
        self.T = termal.sicaklik(z, q)
        ab = termal.a_eksi_b(self.T)
        sun = (self.c["katman"].to_numpy() == "sunek") | (self.T > 450.0)
        self.viskoz = sun
        ek = np.zeros(n, bool) if ek_krip is None else np.asarray(ek_krip, bool)
        ab = np.where(ek & ~sun, np.maximum(ab, 0.004), ab)
        self.vs = (ab > 0) & ~sun                 # yavaş kayan sürtünmeli
        self.vw = (ab <= 0) & ~sun                # kilitli, deprem üreten
        self.creep = self.vs | self.viskoz
        self.ic = np.flatnonzero(self.creep)
        self.iw = np.flatnonzero(self.vw)
        sig = np.minimum(p.sigma_grad * np.maximum(z, 0.5), p.sigma_cap)
        self.A = np.where(self.vs, np.maximum(np.maximum(ab, 1e-3) * sig, 1e5) * p.krip_carpan, np.inf)
        self.tM = np.where(self.viskoz, termal.maxwell_zamani(self.T, p.tM_ref), np.inf)
        self.v_pl = np.broadcast_to(np.asarray(p.v_pl, float), (n,)).copy()
        self.alan = self.c["alan_m2"].to_numpy()
        self.x = self.c["x_km"].to_numpy()
        kab = (self.c["katman"].to_numpy() == "kabuk") if "katman" in self.c else np.ones(n, bool)
        self.h_kirilgan = float(np.median(self.c["W_km"].to_numpy()[kab])) if "W_km" in self.c else 1.0
        self.z = z
        rng = np.random.default_rng(p.seed)
        # dayanım ve gerilme düşümü (yalnız kilitli hücreler anlamlı)
        zfac = np.clip(z / 6.0, 0.25, 1.0)
        self.dtau = p.dtau0 * zfac * rng.lognormal(-0.5 * p.sigma_h ** 2, p.sigma_h, n)
        self.S_fail = self.dtau.copy()            # göreli: kırılma eşiği dtau, durma 0 civarı
        self.S = rng.uniform(0, 1, n) * self.S_fail
        self.S[self.creep] = 0.0
        self.S_ref = np.zeros(n)
        self.slip = np.zeros(n)
        self.t = 0.0
        self.rng = rng
        # sabit yükleme hızı: -K v_pl
        self.r0 = -(self.K @ self.v_pl.astype(np.float32)).astype(np.float64)
        self.Kc = np.ascontiguousarray(self.K[:, self.ic])
        self.v = np.zeros(n)
        self.v[self.ic] = self.v_pl[self.ic]
        self._kararli_baslangic()
        self.olaylar: list[dict] = []
        self.son_kendi = np.full(n, -np.inf)  # her hücrenin en son kendiliğinden kırıldığı an
        self.kayit_slip = None  # isteğe bağlı: olay kayma haritaları
        self.dt_prev = 1e-3
        self.krip_dt = 0.5
        self.tM_adim_alt = 0.0  # ikinci tur: 0,05 yıl (birinci tur davranışı: 0)

    def _kararli_baslangic(self):
        """Kilitli hücreler dururken krip hücrelerinin gerilmesini sabit tutan hızlar.

        K_CC (v_C - v_pl,C) = -K_CW (0 - v_pl,W)  ->  v_C; S_ref buna göre seçilir.
        Böylece başlangıçta krip kuşağı yapay biçimde hızlanıp kilitli kısmı aşırı yüklemez.
        """
        ic, iw = self.ic, self.iw
        Kcc = np.asarray(self.K[np.ix_(ic, ic)], float)
        rhs = np.asarray(self.K[np.ix_(ic, iw)], float) @ self.v_pl[iw]
        dv = np.linalg.solve(Kcc, rhs)
        vC = np.clip(self.v_pl[ic] + dv, 0.02 * self.v_pl[ic], 20 * self.v_pl[ic])
        vs = self.vs[ic]
        S_ref = np.zeros(len(ic))
        S_ref[vs] = -self.A[ic][vs] * np.log(vC[vs] / self.v_pl[ic][vs])
        S_ref[~vs] = -(vC[~vs] - self.v_pl[ic][~vs]) * self.kself[ic][~vs] * self.tM[ic][~vs]
        self.S_ref[ic] = S_ref
        self.v[ic] = vC

    # ------------------------------------------------------------- durum
    def durum(self):
        return {k: getattr(self, k).copy() for k in ("S", "S_ref", "slip", "v", "S_fail", "dtau", "son_kendi")
                if hasattr(self, k)} | {"t": self.t, "dt_prev": self.dt_prev}

    def yukle(self, d):
        for k, val in d.items():
            setattr(self, k, val.copy() if isinstance(val, np.ndarray) else val)

    def dayanimi_yenile(self, seed):
        """Bilinmeyen küçük ölçekli dayanım farklarını yeniden çek (topluluk üyesi çeşitliliği)."""
        rng = np.random.default_rng(seed)
        z = self.z
        zfac = np.clip(z / 6.0, 0.25, 1.0)
        frac = np.clip(self.S / np.maximum(self.S_fail, 1), -2, 1.0)
        self.dtau = self.p.dtau0 * zfac * rng.lognormal(-0.5 * self.p.sigma_h ** 2, self.p.sigma_h, self.n)
        self.S_fail = self.dtau.copy()
        self.S[self.vw] = frac[self.vw] * self.S_fail[self.vw]
        self.rng = rng

    # ------------------------------------------------------------- krip hızları
    def _krip_hizi(self, S):
        ic = self.ic
        v = np.empty(len(ic))
        vs = self.vs[ic]
        dS = S[ic] - self.S_ref[ic]
        v[vs] = self.v_pl[ic][vs] * np.exp(np.clip(dS[vs] / self.A[ic][vs], -50, 30))
        vis = ~vs
        v[vis] = self.v_pl[ic][vis] + dS[vis] / (self.kself[ic][vis] * self.tM[ic][vis])
        return v

    def _oran(self):
        vc = self._krip_hizi(self.S)
        self.v[self.ic] = vc
        rate = self.r0 + self.Kc @ vc.astype(np.float32)
        return rate, vc

    # ------------------------------------------------------------- zaman adımı
    def ilerle(self, t_son: float, max_dt: float = 1.0, deprem=True, maske_kirilma=None):
        """t_son'a kadar ilerle; kendiliğinden depremlere izin ver."""
        while self.t < t_son - 1e-9:
            rate, vc = self._oran()
            # kararlılık: sürtünmeli krip için A/(k v), ağdalı için t_M
            ic = self.ic
            vs = self.vs[ic]
            lim = np.inf
            if vs.any():
                lim = np.min(self.krip_dt * self.A[ic][vs] / (self.kself[ic][vs] * np.maximum(vc[vs], 1e-12)))
            if (~vs).any():
                # ağdalı hücrelerin kendi geri beslemesi _adim'de üstel olarak tam çözülür; çok kısa
                # Maxwell zamanlı (aşırı sıcak, fiilen anında gevşeyen) hücreler adımı sınırlamaz
                tmv = self.tM[ic][~vs]
                tmv = tmv[tmv >= self.tM_adim_alt]
                if len(tmv):
                    lim = min(lim, self.krip_dt * np.min(tmv))
            dt = min(max_dt, lim, self.dt_prev * 1.6, t_son - self.t)
            dt = max(dt, 1e-6)
            # adım içinde yükleme hızı donuk; kırılmalar analitik zamanlarında işlenir
            kalan = dt
            buyuk = False
            while kalan > 0:
                k, t_fail = (-1, np.inf)
                if deprem:
                    k, t_fail = self._ilk_kirilma(rate, maske_kirilma)
                if t_fail <= kalan:
                    self._adim(t_fail, rate, vc)
                    kalan -= t_fail
                    ev = self.kirilma(int(self.iw[k]), maske=maske_kirilma)
                    if ev is not None and ev["mw"] >= 5.0:
                        self.sig_rahatla()
                    if ev is not None and ev["mw"] >= 5.5:
                        buyuk = True  # krip hızları değişti: hızı yeniden hesapla
                        break
                else:
                    self._adim(kalan, rate, vc)
                    kalan = 0.0
            self.dt_prev = 1e-4 if buyuk else dt

    def _ilk_kirilma(self, rate, maske):
        iw = self.iw
        r = rate[iw]
        pos = r > 0
        if not pos.any():
            return -1, np.inf
        gap = self.S_fail[iw] - self.S[iw]
        tf = np.where(pos, gap / np.where(pos, r, 1), np.inf)
        if maske is not None:
            tf = np.where(maske[iw], tf, np.inf)
        k = int(np.argmin(tf))
        return k, max(tf[k], 0.0)

    def _adim(self, dt, rate, vc):
        """Kilitli hücreler: donuk hız. Krip hücreleri: kendi geri beslemesi kapalı çözümle
        (kararlılık için örtük), diğer hücrelerden gelen yük adım başındaki hızla."""
        ic = self.ic
        k = self.kself[ic]
        R = rate[ic] + k * vc               # kendi katkısı çıkarılmış dış yükleme
        S0 = self.S[ic] - self.S_ref[ic]
        vpl = self.v_pl[ic]
        dS = np.empty(len(ic))
        vs = self.vs[ic]
        if vs.any():
            A = self.A[ic][vs]
            x0 = S0[vs] / A
            a_ = R[vs] / A
            b_ = k[vs] * vpl[vs] / A
            y0 = np.exp(-np.clip(x0, -60, 60))
            ad = a_ * dt
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                y = np.where(np.abs(ad) < 1e-9, y0 + b_ * dt,
                             b_ / a_ + (y0 - b_ / a_) * np.exp(-np.clip(ad, -60, 60)))
            y = np.maximum(y, 1e-26)
            dS[vs] = A * (-np.log(y) - x0)
        vis = ~vs
        if vis.any():
            tM = self.tM[ic][vis]
            Sinf = (R[vis] - k[vis] * vpl[vis]) * tM
            Sn = Sinf + (S0[vis] - Sinf) * np.exp(-dt / tM)
            dS[vis] = Sn - S0[vis]
        dslip = (R * dt - dS) / k
        self.S += rate * dt
        self.S[ic] = self.S_ref[ic] + S0 + dS
        self.slip[ic] += dslip
        self.t += dt

    def sig_rahatla(self, v_max=1.0):
        """Kararlı (yavaş kayan) hücrelerde hız v_max'ı aşacak gerilme fazlası hemen kayarak boşalır
        (depremin sığ ve derin kenarlara sızması ve ilk saatlerdeki art kayma)."""
        ic = self.ic[self.vs[self.ic]]
        if len(ic) == 0:
            return
        lim = self.S_ref[ic] + self.A[ic] * np.log(v_max / self.v_pl[ic])
        for _ in range(50):
            ex = self.S[ic] - lim
            j = ex > 0
            if not j.any():
                break
            d = ex[j] / self.kself[ic][j]
            self.S += self.K[:, ic[j]] @ d.astype(np.float32)
            self.slip[ic[j]] += d

    # ------------------------------------------------------------- kırılma
    def kirilma(self, baslangic: int | np.ndarray, maske=None, zorla=False, kayit=True,
                tarih: str | None = None):
        """Çığ kırılması. ``maske``: kırılabilecek hücreler (zorla modunda kırık alanı).

        zorla=True: tarihsel deprem; maske içindeki tüm kilitli hücreler
        biriktirdiği gerilmeyi (eşiğe ulaşmamış olsa bile) boşaltır, kırık
        maske dışına taşmaz. Kayma ve büyüklük modelin kendi gerilme
        durumundan çıkar; gözlenen büyüklükle karşılaştırılır.
        """
        p = self.p
        S0 = self.S.copy()
        ok = self.vw.copy() if maske is None else (self.vw & maske)
        thr = self.S_fail.copy()
        d_ev = np.zeros(self.n)
        if zorla:
            fail = ok & (self.S > 0.0)
        else:
            fail = np.zeros(self.n, bool)
            fail[np.atleast_1d(baslangic)] = True
        it = 0
        from scipy.ndimage import maximum_filter
        h = self.h_kirilgan  # kilitli hücre boyu (km); komşuluk ~2 km
        kb = np.floor((self.x - self.x.min()) / h).astype(int)
        zb = np.floor(self.z / h).astype(int)
        grid = np.zeros((kb.max() + 1, zb.max() + 1))
        fs = 2 * max(int(round(2.0 / h)), 1) + 1
        target_all = -p.asiri * self.dtau
        cozuldu = 0
        while it < 600:
            if fail.any():
                idx = np.flatnonzero(fail)
                d = np.maximum((self.S[idx] - target_all[idx]) / self.kself[idx], 0.0)
                d_ev[idx] += d
                self.S += self.K[:, idx] @ d.astype(np.float32)
                np.maximum.at(grid, (kb[idx], zb[idx]), d_ev[idx])
            # kırık ucundaki dinamik gerilme yığması: 2 km komşuluktaki en büyük
            # kaymayla orantılı (büyük kırık daha güçlü iter, ince şerit itmez)
            w = np.minimum(maximum_filter(grid, size=fs)[kb, zb] / p.kayma_ref, 1.0)
            kayan = d_ev > 0
            if zorla:
                yeni = np.zeros(self.n, bool)
            else:
                yeni = ok & ~kayan & (self.S >= thr - p.eta_dyn * self.dtau * w) & (self.S > 0)
            if not yeni.any():
                # büyüme durdu: kırılan alanın tamamında gerilmeyi dinamik sürtünme düzeyine
                # indiren kaymayı tam çöz (yarı-statik çatlak; Jacobi'nin yavaş yakınsamasını atlar)
                F = np.flatnonzero(kayan & ok)
                if len(F) == 0 or cozuldu >= 3:
                    break
                for _ in range(3):
                    Kff = np.asarray(self.K[np.ix_(F, F)], dtype=np.float64)
                    x = np.linalg.solve(Kff, -(self.S[F] - target_all[F]))
                    neg = x < 0
                    if not neg.any() or neg.all():
                        break
                    if _ < 2:  # son turda F'yi kırpma: x ile F aynı boyda kalmalı
                        F = F[~neg]
                x = np.maximum(x, 0.0)
                d_ev[F] += x
                self.S += self.K[:, F] @ x.astype(np.float32)
                np.maximum.at(grid, (kb[F], zb[F]), d_ev[F])
                cozuldu += 1
                w = np.minimum(maximum_filter(grid, size=fs)[kb, zb] / p.kayma_ref, 1.0)
                kayan = d_ev > 0
                if zorla:
                    break
                yeni = ok & ~kayan & (self.S >= thr - p.eta_dyn * self.dtau * w) & (self.S > 0)
                if not yeni.any():
                    break
            # kaymaya başlamış hücre, gerilmesi dinamik düzeyin üstünde kaldıkça kaymayı sürdürür
            fail = yeni | (kayan & ok & (self.S > target_all + 0.02 * self.dtau))
            it += 1
        self.slip += d_ev
        if not zorla and hasattr(self, "son_kendi"):
            self.son_kendi[d_ev > 0] = self.t
        # kırık dışındaki krip hücreleri: S_ref sabit; art kayma doğal olarak çıkar
        m0 = p.mu_kayma * float(np.sum(d_ev * self.alan))
        if m0 <= 0:
            return None
        mw = (np.log10(m0) - 9.1) / 1.5
        hit = d_ev > 0.05 * d_ev.max()
        ev = {"t": self.t, "mw": mw, "m0": m0,
              "x0": float(self.x[hit].min()), "x1": float(self.x[hit].max()),
              "ort_kayma": float(np.sum(d_ev * self.alan) / np.sum(self.alan[hit])),
              "maks_kayma": float(d_ev.max()), "alan_km2": float(self.alan[hit].sum() / 1e6),
              "x_hipo": float(self.x[np.atleast_1d(baslangic)[0]]) if not zorla else np.nan,
              "z_hipo": float(self.z[np.atleast_1d(baslangic)[0]]) if not zorla else np.nan,
              "zorla": zorla, "tarih": tarih}
        if kayit:
            self.olaylar.append(ev)
        if self.kayit_slip is not None and mw >= self.kayit_slip[0]:
            self.kayit_slip[1].append((self.t, mw, d_ev.astype(np.float32)))
        return ev

    def katalog(self) -> pd.DataFrame:
        return pd.DataFrame(self.olaylar)

    # ------------------------------------------------------------- tanı
    def km_doluluk(self, dx=1.0):
        """Her km için kilitli hücrelerin ortalama doluluk oranı S/S_fail (0..1)."""
        frac = np.where(self.vw, self.S / np.maximum(self.S_fail, 1), np.nan)
        b = np.floor(self.x / dx).astype(int)
        df = pd.DataFrame({"km": b, "f": frac, "w": np.where(self.vw, self.alan, 0)})
        df = df[np.isfinite(df.f)]
        g = df.groupby("km")
        return pd.DataFrame({"x_km": g.size().index * dx + dx / 2,
                             "doluluk": g.apply(lambda d: np.average(d.f, weights=d.w)).to_numpy()})
