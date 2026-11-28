"""Tek ağ (Ulaş): bütün fay, hücre boyu geçmiş deprem yoğunluğuna göre 200 m - 1 km, sabit.

Hücre boyu kuralı (3 kademe; pencere döşemesindeki 1998-2017 eğitim yoğunluğundan, km başına deprem):
  >= 10 deprem/km -> 200 m; 4-10 -> 500 m; < 4 -> 1 km. 16 km altı 1 km, 20 km altı (sünek) 4 km.

Etkileşim (Ulaş: "İstanbul sonucuna etkisi az az az olsun", "sıkıştırmayı tam 15 GB olacak şekilde yap"):
  - TAM bölge: kuzey kol + Armutlu kolu, 27,30-29,60°D (Ganos doğu ucu - Hersek; İstanbul hattı ve iki
    yanı). Bu bölgedeki bütün hücre çiftleri tam (yoğun matris K_tam.npy, ~10,3 GB; bellek 15 GB, pencere
    eğitimi ve sistem için ~4 GB pay).
  - Diğer bütün çiftlerde, 2 km'lik kaba hücre merkezleri arası <= R_YAKIN km ise tam (seyrek, yakin_*.npz).
  - Geri kalan uzak çiftler 2 km'lik kaba ağın matrisiyle, kaynak potansiyeli gruplanarak (seyrek.py).
  K = K_tau + MU_EFF K_sn, kayma modülü MU_YENI'ye ölçekli. MU_EFF tek ağda sabit (literatür 0,4).
Kullanım: python scripts/tek_ag.py hucre | ucgen | tam | yakin | uzak  [--is 2]
Her adım kaldığı yerden sürer (biten satır blokları / dosyalar atlanır).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from deprem.mekanik import geometri as g, mta  # noqa: E402

OUT = ROOT / "data" / "mekanik_tek"
TAM_LON = (27.30, 29.60)
TAM_KOL = (0, 3)
R_YAKIN = 12.0
MU_EFF = 0.4


def h_kurali():
    d = pd.read_csv(ROOT / "reports" / "tables" / "pencere_dose.csv")
    d["yog"] = d.n_egitim / d.L
    d["h"] = np.where(d.yog >= 10, 0.2, np.where(d.yog >= 4, 0.5, 1.0))

    def f(kol, x):
        r = d[(d.kol == kol) & (d.x0 <= x + 1e-6) & (d.x1 > x)]
        return float(r.h.iloc[0]) if len(r) else 1.0
    return f, d


def hucre():
    OUT.mkdir(parents=True, exist_ok=True)
    iz, sinir = mta.fay_izleri_mta(ROOT / "data/mta/tdfh_marmara.geojson", ROOT / "data/raw/gem_active_faults.geojson")
    f, d = h_kurali()
    c = g.hucreler_kademeli(iz, f)
    c["kesim"] = ""
    for r in sinir:
        m = (c.kol == r[0]) & (c.lon >= r[2]) & (c.lon <= r[3])
        c.loc[m, "kesim"] = r[1]
    c.to_parquet(OUT / "hucreler.parquet")
    print("hücre", len(c), {h: round(float(d[d.h == h].L.sum())) for h in (0.2, 0.5, 1.0)})


def kaba_eslem(c):
    """Her ince hücrenin 2 km'lik kaba ağdaki (mekanik_mta_kaba) karşılığı: aynı kol, en yakın (x, z)."""
    from scipy.spatial import cKDTree
    k = pd.read_parquet(ROOT / "data" / "mekanik_mta_kaba" / "hucreler.parquet")
    out = np.zeros(len(c), int)
    for kol in c.kol.unique():
        ki = np.flatnonzero(k.kol.to_numpy() == kol)
        ci = np.flatnonzero(c.kol.to_numpy() == kol)
        tr = cKDTree(np.column_stack([k.x_km.to_numpy()[ki], k.z_km.to_numpy()[ki]]))
        _, j = tr.query(np.column_stack([c.x_km.to_numpy()[ci], c.z_km.to_numpy()[ci]]))
        out[ci] = ki[j]
    return out


def ucgen():
    c = pd.read_parquet(OUT / "hucreler.parquet")
    np.savez(OUT / "ucgen.npz", **g.ucgen_hazirla(c))
    tam = (c.kol.isin(TAM_KOL) & (c.lon >= TAM_LON[0]) & (c.lon <= TAM_LON[1])).to_numpy()
    np.savez(OUT / "eslem.npz", kaba=kaba_eslem(c), tam=tam)
    print("üçgenler", len(c), "tam bölge hücresi", int(tam.sum()), "GB", round(tam.sum() ** 2 * 4 / 1e9, 2))


_U = {}


def _u():
    if not _U:
        _U.update(dict(np.load(OUT / "ucgen.npz")))
        e = np.load(OUT / "eslem.npz")
        _U["kaba"], _U["tam"] = e["kaba"], e["tam"]
    return _U


def _olcek():
    from deprem.mekanik.ozumseme import MU_MATRIS, MU_YENI
    return np.float32(MU_YENI / MU_MATRIS)


def _tam_satir(r0):
    U = _u()
    I = np.flatnonzero(U["tam"])
    r1 = min(r0 + 128, len(I))
    K = np.load(OUT / "K_tam.npy", mmap_mode="r+")
    for s0 in range(0, len(I), 6000):
        K[r0:r1, s0:s0 + 6000] = g.etkilesim_blok(U, I[r0:r1], I[s0:s0 + 6000], MU_EFF) * _olcek()
    K.flush()
    return r0


def tam(n_is=2):
    from multiprocessing import get_context
    U = _u()
    n = int(U["tam"].sum())
    f = OUT / "K_tam.npy"
    if not f.exists():
        np.lib.format.open_memmap(f, mode="w+", dtype=np.float32, shape=(n, n)).flush()
    bf = OUT / "K_tam_bitti.txt"
    bitti = {int(x) for x in bf.read_text().split()} if bf.exists() else set()
    is_ = [r for r in range(0, n, 128) if r not in bitti]
    print("tam bölge satır blokları", len(is_), "/", len(range(0, n, 128)), flush=True)
    with get_context("fork").Pool(n_is) as pool:
        for i, r0 in enumerate(pool.imap_unordered(_tam_satir, is_)):
            with open(bf, "a") as fh:
                fh.write(f"{r0}\n")
            if i % 20 == 0:
                print("tam", i + 1, "/", len(is_), flush=True)


def _yakin(arg):
    """Kaba alıcı grup: <= R_YAKIN çiftler, tam bölge × tam bölge hariç. Dosyaya COO."""
    gi, crs = arg
    f = OUT / f"yakin_{gi:04d}.npz"
    if f.exists():
        return gi
    U = _u()
    cm, tm = U["kaba"], U["tam"]
    dk = np.load(OUT / "kaba_uzaklik.npy", mmap_mode="r")
    R, C, V = [], [], []
    for cr in crs:
        rec = np.flatnonzero(cm == cr)
        if len(rec) == 0:
            continue
        src = np.flatnonzero(np.isin(cm, np.flatnonzero(dk[cr] <= R_YAKIN)))
        for rr, ss in [(rec[tm[rec]], src[~tm[src]]), (rec[~tm[rec]], src)]:
            if len(rr) == 0 or len(ss) == 0:
                continue
            B = g.etkilesim_blok(U, rr, ss, MU_EFF) * _olcek()
            R.append(np.repeat(rr.astype(np.int32), len(ss)))
            C.append(np.tile(ss.astype(np.int32), len(rr)))
            V.append(B.ravel())
    z = (lambda a: np.concatenate(a) if a else np.zeros(0))
    tmp = f.with_suffix(".tmp.npz")
    np.savez(tmp, r=z(R).astype(np.int32), c=z(C).astype(np.int32), v=z(V).astype(np.float32))
    tmp.replace(f)
    return gi


def kaba_uzaklik():
    f = OUT / "kaba_uzaklik.npy"
    if not f.exists():
        k = pd.read_parquet(ROOT / "data" / "mekanik_mta_kaba" / "hucreler.parquet")
        kx = g.to_local_m(k.lat.to_numpy(), k.lon.to_numpy(), k.z_km.to_numpy()) / 1e3
        np.save(f, np.linalg.norm(kx[:, None, :] - kx[None, :, :], axis=-1).astype(np.float32))
    return np.load(f, mmap_mode="r")


def yakin(n_is=2):
    from multiprocessing import get_context
    dk = kaba_uzaklik()
    gr = [(i, list(range(a, min(a + 40, len(dk))))) for i, a in enumerate(range(0, len(dk), 40))]
    with get_context("fork").Pool(n_is) as pool:
        for i, _ in enumerate(pool.imap_unordered(_yakin, gr)):
            if i % 10 == 0:
                print("yakın", i + 1, "/", len(gr), flush=True)


def uzak():
    """Uzak alan: kaba matris; yakın (<= R) kaba çiftler sıfır."""
    d0 = ROOT / "data" / "mekanik_mta_kaba"
    dk = np.asarray(kaba_uzaklik())
    Kc = (np.load(d0 / "K_tau.npy") + np.float32(MU_EFF) * np.load(d0 / "K_sn.npy")) * _olcek()
    Kc[dk <= R_YAKIN] = 0.0
    np.save(OUT / "K_uzak_kaba.npy", Kc.astype(np.float32))
    print("uzak", Kc.shape)


if __name__ == "__main__":
    n_is = int(sys.argv[sys.argv.index("--is") + 1]) if "--is" in sys.argv else 2
    a = sys.argv[1]
    {"hucre": hucre, "ucgen": ucgen, "uzak": uzak}.get(a, lambda: None)() if a in ("hucre", "ucgen", "uzak") \
        else {"tam": tam, "yakin": yakin}[a](n_is)
