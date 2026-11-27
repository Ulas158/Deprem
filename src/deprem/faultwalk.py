"""Depremden depreme öğrenen, fay fay ayarlı model (ileriye doğru sıralı öğrenme).

Ulas'ın isteği (25 Kasım 2026): en ilk depremden başlayıp depremleri sırayla
işle; her seferinde bir sonraki depremi tahmin et, tutmazsa ayarı güncelle; her
fayı ayrı eğit; fayların (özellikle büyük fayların) birbirine etkisini, uzaktaki
fayın küçük ama sıfır olmayabilecek etkisi dahil, veriden öğren.

Bilimsel karşılığı: ileriye doğru sıralı (prequential, Dawid 1984) öğrenme ve
değerlendirme. Her olay i için:

1. Model, yalnız i'den önceki olaylarla (ve o ana kadarki boşluklarla) eğitilir.
2. Bir sonraki deprem için olasılık verir: ne zaman (bekleme süresi) ve nerede
   (0,5° hücre). Puan: log λ_{hücre_i}(t_i) − ∫_{t_{i−1}}^{t_i} Λ(t) dt.
3. Olay gerçekleşince geçmişe eklenir ve model yeniden eğitilir (sıcak başlangıç).

Puan yalnız bu, olayı görmeden önce verilmiş olasılıklardan hesaplanır; bu yüzden
1700-2025 boyunca hiçbir olay "cevap anahtarını görmüş" modelle puanlanmaz.

Model (hücre c, zaman t)
------------------------
λ_c(t) = π_dönem · [ μ · B_c(t) · exp(a_u + b_u · L_u(t))                         (arka plan, fay birimi u)
                   + Σ_j κ_j g(t − t_j) H_jc                                      (hızlı tetikleme, Omori)
                   + Σ_j κ2_j e^{−(t−t_j)/τ}/τ H2_jc ]                            (yavaş, uzun süreli etki)

* B_c = w · (geçmiş sismisite haritası) + (1 − w) · (fay moment haritası)
* L_u: fayın log doluluğu (kayma hızı × son kırılmadan beri süre / karakteristik kayma)
* a_u, b_u: **fay başına** ayar. Kısmi havuzlama: a_u = A_sistem + δ_u,
  δ_u ~ N(0, σ_δ²); A_sistem ~ N(0, σ_A²). "Ortak" sürümde δ = A = 0,
  "tamamen ayrı" sürümde σ_δ büyük ve sistem katmanı yok.
* H_jc = alan_c · f(d; D_j, q) · exp(γ_aynı·[aynı sistem] + γ_çapraz·[farklı büyük
  sistemler] + γ_cfs · tanh(ΔCFF_jc / 0,1 bar)).
  f(d) = (q−1)/(πD²) (1 + d²/D²)^−q: mesafeyle azalan, **uzaklık sınırı olmayan**
  etki; azalma hızı q veriden öğrenilir (q büyük = uzak etki hızla söner).
* κ_j = K e^{α(m_j − 5,5)}; g Omori-Utsu.
* π_dönem: 1700-1899 dönemlerinde tespit oranı (katalog eksikliği) veriden öğrenilir;
  1900 sonrası 1.

ETAS (aynı düzende): fay ayarı, büyük fay etkileşimi, Coulomb ve yavaş terim kapalı.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from .faultwalk_data import ERAS, M0, WalkData, year_start_days

log = logging.getLogger(__name__)

GLOBAL = [  # ad, öncül ortalama, öncül sd
    ("log_mu", np.log(1.0 / 365.25), 3.0),
    ("w_logit", 0.0, 1.5),
    ("logK", np.log(0.15), 1.5),
    ("alpha", 1.8, 0.7),
    ("log_c", np.log(0.1), 2.0),
    ("log_pm1", np.log(0.15), 1.0),
    ("log_D0", np.log(15.0), 1.0),
    ("gam_m", 1.5, 1.0),
    ("log_qm1", np.log(0.6), 0.7),
    ("g_same", 0.0, 1.0),
    ("g_cross", 0.0, 1.0),
    ("g_cfs", 0.0, 1.0),
    ("logK2", np.log(0.05), 2.0),
    ("alpha2", 1.8, 1.0),
    ("log_tau2", np.log(10 * 365.25), 1.0),
    ("log_D2", np.log(50.0), 1.0),
    ("log_q2m1", np.log(0.6), 0.7),
    ("g_cross2", 0.0, 1.0),
    ("g_cfs2", 0.0, 1.0),
    ("b0", 0.0, 0.5),
    ("k_cov", 0.0, 0.5),            # ek zamanla değişen girdi (ör. dv/v) katsayısı
] + [(f"era{e}", 0.0, 2.0) for e in range(len(ERAS))]


@dataclass(frozen=True)
class Variant:
    name: str
    label: str
    pooling: str = "partial"        # "pooled" | "partial" | "separate" | "none"
    interaction: bool = True        # büyük fay etkileşimi + Coulomb
    cfs: bool = True
    slow: bool = True
    rmax_km: float = 0.0            # >0: bu uzaklıktan ötesi kesilir (uzak etki kapalı)
    load: bool = True               # fay doluluğu
    sd_A: float = 1.0
    sd_delta: float = 0.5
    sd_B: float = 0.3
    sd_eps: float = 0.2
    cov: str = ""                   # W.extra["cov"][cov]: (S, U) ek girdi (eksikte 0)


VARIANTS = {
    "etas": Variant("etas", "ETAS (aynı düzende, referans)", pooling="none", interaction=False, cfs=False,
                    slow=False, load=False),
    "poisson": Variant("poisson", "Yalnız uzun dönem harita (tetiklemesiz)", pooling="none", interaction=False,
                       cfs=False, slow=False, load=False),
    "pooled": Variant("pooled", "Tüm faylar ortak ayar", pooling="pooled"),
    "partial": Variant("partial", "Fay fay ayar, ortak ayara doğru çekilen (ana model)"),
    "separate": Variant("separate", "Fay fay tamamen ayrı ayar", pooling="separate", sd_delta=3.0, sd_eps=1.0),
    "no_interact": Variant("no_interact", "Ana model − faylar arası etki (aynı/çapraz sistem, Coulomb, yavaş)",
                           interaction=False, cfs=False, slow=False),
    "no_cross": Variant("no_cross", "Ana model − büyük faylar arası çapraz etki", interaction=False),
    "no_cfs": Variant("no_cfs", "Ana model − Coulomb gerilme yönü", cfs=False),
    "no_slow": Variant("no_slow", "Ana model − yavaş (yıllar-on yıllar) etki", slow=False),
    "no_far": Variant("no_far", "Ana model, 100 km ötesi etki kesik", rmax_km=100.0),
    "no_load": Variant("no_load", "Ana model − fay doluluğu", load=False),
    "partial_dvv": Variant("partial_dvv", "Ana model + Doğu Marmara dv/v", cov="dvv"),
}


# ------------------------------------------------------------------ JAX modeli
class Model:
    def __init__(self, W: WalkData, v: Variant, cfs_override: np.ndarray | None = None,
                 unit_override: np.ndarray | None = None):
        import jax
        import jax.numpy as jnp
        jax.config.update("jax_enable_x64", True)
        self.jax, self.jnp = jax, jnp
        self.W, self.v = W, v
        N, C = W.d2.shape[:2]
        S = W.tau.size - 1
        unit_of_cell = W.unit_of_cell if unit_override is None else unit_override
        U = len(W.units)
        Sy = len(W.systems)
        self.N, self.C, self.S, self.U, self.Sy = N, C, S, U, Sy
        unit_sys = W.units["sys"].to_numpy()
        cell_sys = unit_sys[unit_of_cell]
        ev_sys = W.ev_sys if unit_override is None else unit_sys[unit_of_cell[W.ev_cell]]
        same = (ev_sys[:, None] == cell_sys[None, :]).astype(np.float32)
        maj = W.major
        cross = ((ev_sys[:, None] != cell_sys[None, :]) & maj[ev_sys][:, None] & maj[cell_sys][None, :]).astype(np.float32)
        cfs = W.cfs if cfs_override is None else cfs_override
        # birim başına harita kütleleri
        BuS = np.zeros((W.smap.shape[0], U))
        BuF = np.zeros(U)
        for u in range(U):
            m = unit_of_cell == u
            BuS[:, u] = W.smap[:, m].sum(1)
            BuF[u] = W.fmap[m].sum()
        Lm = W.L if v.load else np.zeros_like(W.L)
        Xc = W.extra.get("cov", {}).get(v.cov) if v.cov else None
        if Xc is None:
            Xc = np.zeros_like(W.L)
        dtau = np.diff(W.tau)
        era_s = W.step_era
        self.ev_unit = unit_of_cell[W.ev_cell]
        names = [g[0] for g in GLOBAL]
        self.gi = {n: i for i, n in enumerate(names)}
        G = len(names)
        self.nG = G
        self.size = G + 2 * Sy + 2 * U
        mean = np.zeros(self.size)
        sd = np.ones(self.size)
        free = np.ones(self.size, bool)
        for i, (n, m_, s_) in enumerate(GLOBAL):
            mean[i], sd[i] = m_, s_
        sA, sB = slice(G, G + Sy), slice(G + Sy, G + 2 * Sy)
        sD, sE = slice(G + 2 * Sy, G + 2 * Sy + U), slice(G + 2 * Sy + U, G + 2 * Sy + 2 * U)
        self.sl = dict(A=sA, B=sB, D=sD, E=sE)
        sd[sA], sd[sB], sd[sD], sd[sE] = v.sd_A, v.sd_B, v.sd_delta, v.sd_eps
        fix0 = []
        if v.pooling in ("pooled", "none"):
            free[sA] = free[sB] = free[sD] = free[sE] = False
        if v.pooling == "separate":
            free[sA] = free[sB] = False
        if v.pooling == "none" or not v.load:
            fix0 += ["b0"]
            free[sB] = free[sE] = False
        if not v.interaction:
            fix0 += ["g_same", "g_cross", "g_cross2"]
        if not v.cfs:
            fix0 += ["g_cfs", "g_cfs2"]
        if not v.slow:
            for n in ("logK2", "alpha2", "log_tau2", "log_D2", "log_q2m1", "g_cross2", "g_cfs2"):
                free[self.gi[n]] = False
        if not v.cov:
            fix0 += ["k_cov"]
        for n in fix0:
            free[self.gi[n]] = False
            mean[self.gi[n]] = 0.0
        self.mean, self.sd, self.free = mean, sd, free
        self.slow_on = v.slow
        self.trig_on = v.name != "poisson"
        if not self.trig_on:
            for n in ("logK", "alpha", "log_c", "log_pm1", "log_D0", "gam_m", "log_qm1"):
                free[self.gi[n]] = False

        rmax2 = (v.rmax_km ** 2) if v.rmax_km > 0 else 0.0
        # uzak çiftler: hücre merkezi; yakın çiftler (merkez < 90 km): 3×3 alt nokta ortalaması
        d2c = W.d2[:, :, W.d2.shape[2] // 2].astype(float)
        near = d2c < 90.0 ** 2
        nj, nc = np.nonzero(near)                      # j'ye göre sıralı
        self._near_cum = np.searchsorted(nj, np.arange(N + 1))
        a = dict(d2=jnp.asarray(d2c), nj=jnp.asarray(nj), nc=jnp.asarray(nc), d2n=jnp.asarray(W.d2[nj, nc].astype(float)),
                 area=jnp.asarray(W.area), same=jnp.asarray(same), cross=jnp.asarray(cross),
                 cfs=jnp.asarray(cfs), m=jnp.asarray(W.m_ev - M0), t=jnp.asarray(W.t_ev),
                 tau=jnp.asarray(W.tau), dtau=jnp.asarray(dtau), era_s=jnp.asarray(era_s),
                 ev_era=jnp.asarray(W.ev_era), ev_cell=jnp.asarray(W.ev_cell), ev_unit=jnp.asarray(self.ev_unit),
                 ev_step=jnp.asarray(np.clip(W.ev_step, 0, S - 1)), ev_year=jnp.asarray(W.step_year[np.clip(W.ev_step, 0, S - 1)]),
                 step_year=jnp.asarray(W.step_year), BuS=jnp.asarray(BuS), BuF=jnp.asarray(BuF),
                 smap=jnp.asarray(W.smap), fmap=jnp.asarray(W.fmap), L=jnp.asarray(Lm), X=jnp.asarray(Xc),
                 unit_sys=jnp.asarray(unit_sys), unit_of_cell=jnp.asarray(unit_of_cell))
        self.a = a
        gi = self.gi
        nE = len(ERAS)
        trig_on, slow_on = self.trig_on, self.slow_on

        def unpack(th):
            g = {n: th[i] for n, i in gi.items()}
            A, B, D, E = th[sA], th[sB], th[sD], th[sE]
            au = A[a["unit_sys"]] + D
            bu = g["b0"] + B[a["unit_sys"]] + E
            pi = jnp.concatenate([jax.nn.sigmoid(jnp.stack([g[f"era{e}"] for e in range(nE)])), jnp.ones(1)])
            return g, au, bu, pi

        near_cum = self._near_cum

        def kernel(a, logD, gam, logqm1, mm):
            q = 1.0 + jnp.exp(logqm1)
            D2 = jnp.exp(2 * logD + gam * mm)
            f = (q - 1) / (jnp.pi * D2[:, None]) * (1.0 + a["d2"] / D2[:, None]) ** (-q)
            if rmax2 > 0:
                f = jnp.where(a["d2"] <= rmax2, f, 0.0)
            P = int(near_cum[a["d2"].shape[0]])
            nj, nc, d2n = a["nj"][:P], a["nc"][:P], a["d2n"][:P]
            Dn = D2[nj][:, None]
            fn = ((q - 1) / (jnp.pi * Dn) * (1.0 + d2n / Dn) ** (-q)).mean(-1)
            f = f.at[nj, nc].set(fn)
            return f * a["area"][None, :]

        EV_KEYS = ("d2", "same", "cross", "cfs", "m", "t", "ev_era", "ev_cell", "ev_unit", "ev_step", "ev_year")
        ST_KEYS = ("dtau", "era_s", "step_year", "L", "X")

        def sub(nE, nS):
            b = dict(a)
            for k in EV_KEYS:
                b[k] = a[k][:nE]
            for k in ST_KEYS:
                b[k] = a[k][:nS]
            b["tau"] = a["tau"][:nS + 1]
            return b

        def parts(th, nE, nS):
            a = sub(nE, nS)  # olay ve adım dizileri ilk nE / nS ile kesilir
            g, au, bu, pi = unpack(th)
            mm = a["m"]
            mu = jnp.exp(g["log_mu"])
            w = jax.nn.sigmoid(g["w_logit"])
            # arka plan, adım başına
            ex = jnp.exp(au[None, :] + bu[None, :] * a["L"] + g["k_cov"] * a["X"])  # S × U
            Bu = w * a["BuS"][a["step_year"]] + (1 - w) * a["BuF"][None, :]       # S × U
            bg_s = mu * a["dtau"] * (Bu * ex).sum(1)
            # arka plan, olay başına
            exu = jnp.exp(au[a["ev_unit"]] + bu[a["ev_unit"]] * a["L"][a["ev_step"], a["ev_unit"]]
                          + g["k_cov"] * a["X"][a["ev_step"], a["ev_unit"]])
            Bc = w * a["smap"][a["ev_year"], a["ev_cell"]] + (1 - w) * a["fmap"][a["ev_cell"]]
            bg_e = mu * Bc * exu
            trig_e = jnp.zeros_like(bg_e)
            trig_s = jnp.zeros_like(bg_s)
            dt = a["t"][None, :] - a["t"][:, None]                                 # j × i
            if trig_on:
                K = jnp.exp(g["logK"] + jnp.clip(g["alpha"], 0, 3.5) * mm)
                c = jnp.exp(g["log_c"])
                p = 1.0 + jnp.exp(g["log_pm1"])
                mod = jnp.exp(g["g_same"] * a["same"] + g["g_cross"] * a["cross"] + g["g_cfs"] * a["cfs"])
                H = kernel(a, g["log_D0"], g["gam_m"], g["log_qm1"], mm) * mod   # N × C
                gt = jnp.where(dt > 0, (p - 1) / c * (1 + jnp.clip(dt, 0) / c) ** (-p), 0.0)
                trig_e = (K[:, None] * gt * H[:, a["ev_cell"]]).sum(0)
                Gc = lambda x: 1 - (1 + jnp.clip(x, 0) / c) ** (1 - p)
                tt = a["tau"][None, :] - a["t"][:, None]                           # N × (S+1)
                Gd = Gc(tt[:, 1:]) - Gc(tt[:, :-1])
                trig_s = ((K * H.sum(1))[:, None] * Gd).sum(0)
            if slow_on:
                K2 = jnp.exp(g["logK2"] + jnp.clip(g["alpha2"], 0, 3.5) * mm)
                tau2 = jnp.exp(g["log_tau2"])
                mod2 = jnp.exp(g["g_same"] * a["same"] + g["g_cross2"] * a["cross"] + g["g_cfs2"] * a["cfs"])
                H2 = kernel(a, g["log_D2"], g["gam_m"], g["log_q2m1"], mm) * mod2
                g2 = jnp.where(dt > 0, jnp.exp(-jnp.clip(dt, 0) / tau2) / tau2, 0.0)
                trig_e = trig_e + (K2[:, None] * g2 * H2[:, a["ev_cell"]]).sum(0)
                tt = a["tau"][None, :] - a["t"][:, None]
                G2 = lambda x: 1 - jnp.exp(-jnp.clip(x, 0) / tau2)
                trig_s = trig_s + ((K2 * H2.sum(1))[:, None] * (G2(tt[:, 1:]) - G2(tt[:, :-1]))).sum(0)
            lam_e = pi[a["ev_era"]] * (bg_e + trig_e)
            comp_s = pi[a["era_s"]] * (bg_s + trig_s)
            return jnp.log(lam_e + 1e-300), comp_s

        mean_j, sd_j = jnp.asarray(mean), jnp.asarray(sd)

        def loss(th, ev_mask, st_mask, nE, nS):
            ll, comp = parts(th, nE, nS)
            prior = 0.5 * (((th - mean_j) / sd_j) ** 2).sum()
            return -(ll * ev_mask).sum() + (comp * st_mask).sum() + prior

        self._parts = jax.jit(parts, static_argnums=(1, 2))
        self._vg = jax.jit(jax.value_and_grad(loss), static_argnums=(3, 4))

        def cell_lambda(th, i):
            """i olayından hemen önce bütün hücrelerde λ_c (olasılık haritası için)."""
            g, au, bu, pi = unpack(th)
            mm = a["m"]
            s = a["ev_step"][i]
            w = jax.nn.sigmoid(g["w_logit"])
            uc = a["unit_of_cell"]
            bg = jnp.exp(g["log_mu"]) * (w * a["smap"][a["step_year"][s]] + (1 - w) * a["fmap"]) * \
                jnp.exp(au[uc] + bu[uc] * a["L"][s][uc] + g["k_cov"] * a["X"][s][uc])
            dt = a["t"][i] - a["t"]
            tr = jnp.zeros_like(bg)
            if trig_on:
                K = jnp.exp(g["logK"] + jnp.clip(g["alpha"], 0, 3.5) * mm)
                c = jnp.exp(g["log_c"])
                p = 1.0 + jnp.exp(g["log_pm1"])
                mod = jnp.exp(g["g_same"] * a["same"] + g["g_cross"] * a["cross"] + g["g_cfs"] * a["cfs"])
                H = kernel(a, g["log_D0"], g["gam_m"], g["log_qm1"], mm) * mod
                gt = jnp.where(dt > 0, (p - 1) / c * (1 + jnp.clip(dt, 0) / c) ** (-p), 0.0)
                tr = tr + ((K * gt)[:, None] * H).sum(0)
            if slow_on:
                K2 = jnp.exp(g["logK2"] + jnp.clip(g["alpha2"], 0, 3.5) * mm)
                tau2 = jnp.exp(g["log_tau2"])
                mod2 = jnp.exp(g["g_same"] * a["same"] + g["g_cross2"] * a["cross"] + g["g_cfs2"] * a["cfs"])
                H2 = kernel(a, g["log_D2"], g["gam_m"], g["log_q2m1"], mm) * mod2
                g2 = jnp.where(dt > 0, jnp.exp(-jnp.clip(dt, 0) / tau2) / tau2, 0.0)
                tr = tr + ((K2 * g2)[:, None] * H2).sum(0)
            return bg + tr

        self._cell_lambda = jax.jit(cell_lambda)

    # -------------------------------------------------------------- yardımcılar
    def init(self):
        return self.mean.copy()

    def buckets(self, n_ev: int, n_st: int):
        nE = min(self.N, int(np.ceil(max(n_ev, 1) / 64) * 64))
        nS = min(self.S, int(np.ceil(max(n_st, 1) / 128) * 128))
        return nE, nS

    def fit(self, th0, n_ev: int, n_st: int, maxiter: int = 60):
        """İlk n_ev olay ve ilk n_st adımla MAP (öncül = başlangıç ayarı)."""
        jnp = self.jnp
        nE, nS = self.buckets(n_ev, n_st)
        ev_mask = jnp.asarray((np.arange(nE) < n_ev).astype(float))
        st_mask = jnp.asarray((np.arange(nS) < n_st).astype(float))
        free = self.free
        base = th0.copy()

        def f(x):
            th = base.copy()
            th[free] = x
            val, gr = self._vg(jnp.asarray(th), ev_mask, st_mask, nE, nS)
            gr = np.asarray(gr)[free]
            if not np.isfinite(val) or not np.all(np.isfinite(gr)):
                return 1e12, np.zeros(free.sum())
            return float(val), gr

        r = minimize(f, base[free], jac=True, method="L-BFGS-B", options={"maxiter": maxiter})
        th = base.copy()
        th[free] = r.x
        return th

    def parts(self, th, n_ev: int | None = None, n_st: int | None = None):
        nE, nS = (self.N, self.S) if n_ev is None else self.buckets(n_ev, n_st)
        ll, comp = self._parts(self.jnp.asarray(th), nE, nS)
        return np.asarray(ll), np.asarray(comp)

    def cell_probs(self, th, i):
        lam = np.asarray(self._cell_lambda(self.jnp.asarray(th), int(i)))
        return lam / lam.sum()

    def named(self, th) -> dict:
        return {n: float(th[i]) for n, i in self.gi.items()}


# ------------------------------------------------------------------ sıralı yürüyüş
def walk(W: WalkData, v: Variant, refit_every: int = 1, start_i: int = 0, freeze_at: float | None = None,
         maxiter: int = 60, first_maxiter: int = 300, cfs_override=None, unit_override=None,
         keep_params_every: int = 50, cell_probs: bool = True, progress: str = "") -> tuple[pd.DataFrame, list]:
    """Her olay için: önceki olaylarla eğit → bu olayı puanla → geçmişe ekle.

    ``freeze_at`` (gün): bu andan sonra ayar güncellenmez (sabit ayarlı model).
    Dönüş: olay başına puan tablosu ve ara ayar anlık görüntüleri.
    """
    import time as _time
    M = Model(W, v, cfs_override=cfs_override, unit_override=unit_override)
    th = M.init()
    rows, snaps = [], []
    t_start = _time.time()
    fitted_once = False
    for i in range(start_i, M.N):
        n_st = int(W.ev_bstep[i - 1]) if i > 0 else 0       # adımlar t_{i−1}'e kadar
        frozen = freeze_at is not None and fitted_once and (i > 0 and W.t_ev[i - 1] >= freeze_at)
        if not frozen and (i - start_i) % refit_every == 0 and i > 0:
            th = M.fit(th, n_ev=i, n_st=n_st, maxiter=first_maxiter if not fitted_once else maxiter)
            fitted_once = True
        s_hi = int(W.ev_bstep[i])
        ll, comp = M.parts(th, i + 1, s_hi)                             # t_i ile başlayan adım
        s_lo = n_st
        time_comp = comp[s_lo:s_hi].sum()
        row = {"i": i, "time": W.cat["time"].iloc[i], "mw": W.m_ev[i], "lat": W.cat["latitude"].iloc[i],
               "lon": W.cat["longitude"].iloc[i], "cell": int(W.ev_cell[i]), "unit": W.units["unit"].iloc[M.ev_unit[i]],
               "system": W.units["system"].iloc[M.ev_unit[i]], "loglam": float(ll[i]), "comp": float(time_comp),
               "score": float(ll[i] - time_comp)}
        if cell_probs:
            p = M.cell_probs(th, i)
            pc = p[W.ev_cell[i]]
            row["p_cell"] = float(pc)
            row["rank_frac"] = float((p > pc).sum() / p.size)
            # zaman puanı = log Λ − ∫Λ; konum puanı = log p(hücre)
            row["loc_score"] = float(np.log(pc))
            row["time_score"] = row["score"] - row["loc_score"]
        rows.append(row)
        if keep_params_every and (i % keep_params_every == 0 or i == M.N - 1):
            snaps.append({"i": i, "time": str(W.cat["time"].iloc[i].date()), **M.named(th),
                          "A": th[M.sl["A"]].tolist(), "B": th[M.sl["B"]].tolist(),
                          "D": th[M.sl["D"]].tolist(), "E": th[M.sl["E"]].tolist()})
        if progress and i % 50 == 0:
            log.info("[%s] %d/%d olay, %.0f sn", progress, i, M.N, _time.time() - t_start)
    return pd.DataFrame(rows), snaps
