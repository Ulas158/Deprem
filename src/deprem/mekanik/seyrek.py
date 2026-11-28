"""Tek ağın sıkıştırılmış etkileşim işleci (scripts/tek_ag.py'nin ürettiği parçalardan).

K = T (tam bölge, yoğun, devrik saklı) + N (yakın çiftler, seyrek) + U (uzak çiftler, 2 km'lik kaba ağ).
Simülatör K'yı yalnız şu biçimlerde kullanır; işleç bunları sağlar:
  K @ v,  K[:, idx] @ d,  K[np.ix_(A, B)] (küçük bloklar),  K.diagonal().
Uzak kısım: kaynak ince hücrelerin potansiyeli (alan × kayma) kaba hücrede toplanıp kaba hücre alanına
bölünür (eşdeğer düzgün kayma), kaba matrisle çarpılır, alıcı ince hücre kendi kaba hücresinin değerini
alır. Tam bölgedeki alıcıya tam bölgedeki kaynaktan uzak katkı gelmez (T zaten tam).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


class SeyrekK:
    def __init__(self, klasor: Path, alan_ince: np.ndarray, bellek=True):
        k = Path(klasor)
        e = np.load(k / "eslem.npz")
        self.cm, self.tam = e["kaba"], e["tam"].astype(bool)
        self.n = len(self.cm)
        self.shape = (self.n, self.n)
        self.dtype = np.float32
        self.I = np.flatnonzero(self.tam)
        self.yer = np.full(self.n, -1)
        self.yer[self.I] = np.arange(len(self.I))
        T = np.load(k / "K_tam.npy", mmap_mode="r")
        self.TT = np.ascontiguousarray(T) if bellek else T  # [kaynak, alıcı]
        R, C, V = [], [], []
        for f in sorted(k.glob("yakin_*.npz")):
            z = np.load(f)
            R.append(z["r"]); C.append(z["c"]); V.append(z["v"])
        self.N = sp.csc_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=self.shape,
                               dtype=np.float32)
        self.Kf = np.load(k / "K_uzak_kaba.npy")
        kaba = pd.read_parquet(Path(klasor).parents[0] / "mekanik_mta_kaba" / "hucreler.parquet")
        self.Ac = kaba.alan_m2.to_numpy()
        self.a = np.asarray(alan_ince, float)
        self.nc = len(self.Ac)

    # ---------------------------------------------------------------- çarpımlar
    def _uzak(self, idx, d):
        """idx kaynaklarının (kayma d) uzak katkısı, bütün alıcılar için."""
        p = self.a[idx] * d
        t = self.tam[idx]
        cs, inv = np.unique(self.cm[idx], return_inverse=True)
        p_hep = np.bincount(inv, p, minlength=len(cs)) / self.Ac[cs]
        p_dis = np.bincount(inv, np.where(t, 0.0, p), minlength=len(cs)) / self.Ac[cs]
        Kfc = self.Kf[:, cs]
        s_hep = Kfc @ p_hep.astype(np.float32)
        out = s_hep[self.cm]
        if t.any():
            s_dis = Kfc @ p_dis.astype(np.float32)
            out[self.tam] = s_dis[self.cm[self.tam]]
        return out

    def sutun_carp(self, idx, d):
        idx = np.asarray(idx)
        d = np.asarray(d, np.float32)
        out = np.asarray(self.N[:, idx] @ d, dtype=np.float32).ravel() if len(idx) else np.zeros(self.n, np.float32)
        t = self.tam[idx]
        if t.any():
            out[self.I] += d[t] @ self.TT[self.yer[idx[t]]]
        out += self._uzak(idx, d.astype(float)).astype(np.float32)
        return out

    def __matmul__(self, v):
        v = np.asarray(v)
        nz = np.flatnonzero(v)
        if len(nz) < self.n // 4:
            return self.sutun_carp(nz, v[nz])
        v = v.astype(np.float32)
        out = np.asarray(self.N @ v, dtype=np.float32).ravel()
        out[self.I] += v[self.I] @ self.TT
        out += self._uzak(np.arange(self.n), v.astype(float)).astype(np.float32)
        return out

    def diagonal(self):
        d = np.asarray(self.N.diagonal(), np.float32).copy()
        d[self.I] += np.asarray(np.diagonal(self.TT), np.float32)
        return d

    def blok(self, A, B):
        """Yoğun K[A][:, B] (küçük kümeler için)."""
        A, B = np.asarray(A), np.asarray(B)
        out = np.asarray(self.N[:, B][A, :].todense(), dtype=np.float64) if len(A) and len(B) else np.zeros((len(A), len(B)))
        ta, tb = self.tam[A], self.tam[B]
        if ta.any() and tb.any():
            out[np.ix_(ta, tb)] += np.asarray(self.TT[np.ix_(self.yer[B[tb]], self.yer[A[ta]])], float).T
        # uzak: alıcı a, kaynak b: Kf[cm[a], cm[b]] * a_b / A_c(b); tam×tam ve yakın kaba çiftler hariç
        U = self.Kf[np.ix_(self.cm[A], self.cm[B])].astype(float) * (self.a[B] / self.Ac[self.cm[B]])[None, :]
        U[np.ix_(ta, tb)] = 0.0
        return out + U

    # ---------------------------------------------------------------- dizinleme
    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 2:
            a, b = key
            if isinstance(a, slice) and a == slice(None):
                return _Sutunlar(self, np.atleast_1d(np.asarray(b)))
            if isinstance(a, np.ndarray) and isinstance(b, np.ndarray) and a.ndim == 2:
                return self.blok(a.ravel(), b.ravel())
        raise TypeError("SeyrekK: desteklenmeyen dizinleme")


class _Sutunlar:
    def __init__(self, K: SeyrekK, idx):
        self.K, self.idx = K, idx
        self.shape = (K.n, len(idx))

    def __matmul__(self, d):
        return self.K.sutun_carp(self.idx, d)
