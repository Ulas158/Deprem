"""İkinci tur: özümsemeli, geri saran ayar ve bağımsız sınav.

Aşamalar
    kucuk  --ag silivri|cinarcik   100 m yerel ağ, GFZ küçük depremleri; ayar 2006-2017
    buyuk                           2 km kaba ağ (yan kollar dahil), 1719-1894 büyük depremler
    sinav                           donmuş ayarlarla sınav (küçük: 2018-2023; büyük: 1912-2025),
                                    bugünkü durum ve 30 yıllık ileri bakış

Arama (plan, "Arama durdurma kuralı"): 200 Latin hiperküp örneği, sonra en iyilerin çevresinde
50'lik turlar; ardışık 3 turda en iyi ayar puanı %1'den az iyileşirse durur. Her değerlendirme
sonucu kalıcı klasöre (DEPREM_KALICI) satır satır yazılır; kesintide kaldığı yerden sürer.
Parametreler ortak (bütün fay) olduğu için "geri sarma" bütün geçmişi değiştirir: her yeni
kombinasyon 1719'dan yeniden koşulur ve önceki tutmalar yeniden denetlenir; puan önce baştan
itibaren ardışık tutan deprem sayısına göre sıralanır.

Kullanım: OMP_NUM_THREADS=1 python scripts/mekanik_ikinci_tur.py kucuk --ag silivri --is 4
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from deprem.mekanik import ozumseme as oz  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
log = logging.getLogger("ikinci_tur")

GFZ_MC = {"becker": 1.9, "bentz": 1.6, "armutlu": 1.5}
KUCUK_PAR = ["dtau0", "sigma_h", "eta_dyn", "kalan", "kayma_ref", "tM_ref", "krip_carpan", "sigma_grad",
             "kalinlik_u", "v_pl_marmara"]
if oz.GEOMETRI == "mta":
    # MTA sürümü: kayma hızı ve kalınlık Chartier'den sabit (YEREL_CHARTIER), gecikmeli tetiklenme aranır
    KUCUK_PAR = ["dtau0", "sigma_h", "eta_dyn", "kalan", "kayma_ref", "tM_ref", "krip_carpan", "sigma_grad",
                 "a_rs"]
KURAL3_EK = "_k3"  # kaydıkça zayıflama (zayif) serbest; dosya adlarına eklenir  # mu_eff aynı düzlemde etkisiz; sığ krip bölgesi ağ dışında
BUYUK_PAR = oz.P_AD

AYAR = {"kucuk": dict(t_bas=2006.0, ayar_bitis=2018.0, t_son_sinav=2024.0, isinma=30.0),
        "buyuk": dict(t_bas=1719.0, ayar_bitis=1895.0, t_son_sinav=2026.9, isinma=400.0)}

_G = {}  # süreç başına önbellek


def _hazirla(asama, ag_ad, sinav=False):
    key = (asama, ag_ad, sinav)
    if key in _G:
        return _G[key]
    a = AYAR[asama]
    ag = oz.Ag(ROOT, ag_ad)
    t1 = a["t_son_sinav"] if sinav else a["ayar_bitis"]
    if asama == "kucuk":
        ol, say = oz.cizelge(ag, ROOT, a["t_bas"], t1, a["ayar_bitis"], 1.5, gfz_mw_min=GFZ_MC, dis=False)
    else:
        ol, say = oz.cizelge(ag, ROOT, a["t_bas"], t1, a["ayar_bitis"], 3.0)
    _G.clear()
    _G[key] = (ag, ol, say)
    return _G[key]


def degerlendir(asama, ag_ad, th, seed, sinav=False, tam_yama=False, sonsal=True, kayit_yolu=None,
                ileri=0):
    ag, ol, say = _hazirla(asama, ag_ad, sinav)
    a = AYAR[asama]
    t_son = a["t_son_sinav"] if sinav else a["ayar_bitis"]
    k = oz.Kosu(ag, th, ol, seed, a["isinma"], a["t_bas"], t_son, a["ayar_bitis"], sonsal=sonsal,
                tam_yama=tam_yama, kayit_yolu=kayit_yolu)
    t0 = time.time()
    m, kayit = k.calistir()
    sahte_ayar = oz.sahte_say(m, ol, a["t_bas"], a["ayar_bitis"]) if asama == "buyuk" else []
    if asama == "kucuk":
        # küçük ağda "sahte": tarihte karşılığı olmayan M>=4 kendiliğinden deprem (2006-2017 bu
        # bölgede gözlenen en büyük olaylardan büyük)
        mmax = max([o.mw for o in ol if o.t < a["ayar_bitis"]] + [3.0])
        sahte_ayar = [e for e in m.olaylar if not e.get("zorla") and e["mw"] > mmax + 0.5
                      and e["t"] < a["ayar_bitis"]]
    oz_ = oz.ozet_puan(kayit, len(sahte_ayar), "buyuk" if asama == "buyuk" else "kucuk")
    oz_["sure_s"] = round(time.time() - t0, 1)
    kat = pd.DataFrame([e for e in m.olaylar if not e.get("zorla")])
    oz_["n_kendi"] = int(len(kat))
    oz_["mw_maks_kendi"] = float(kat.mw.max()) if len(kat) else None
    out = {"ozet": oz_, "kayit": kayit}
    if sinav:
        out["m"] = m
        out["sahte_sinav"] = oz.sahte_say(m, ol, a["ayar_bitis"], t_son)
        out["sahte_ayar"] = sahte_ayar
        out["katalog"] = pd.DataFrame(m.olaylar)
        out["say"] = say
    return out


def _is(arg):
    asama, ag_ad, u, adlar, seed = arg
    th = oz.baslangic_theta()
    if "zayif" not in adlar:
        th["zayif"] = 0.0  # eski kural
    th.update(oz.birim_to_theta(np.asarray(u), adlar))
    import signal

    def _sure_doldu(*_):
        raise TimeoutError("deneme süresi aşıldı")
    # takılan (aşırı yavaş) ayar bütün aramayı bekletmesin: küçük aşamada 10 dk sınır
    signal.signal(signal.SIGALRM, _sure_doldu)
    signal.alarm(600 if asama == "kucuk" else 0)
    try:
        r = degerlendir(asama, ag_ad, th, seed)
        o = r["ozet"]
    except Exception as e:  # sayısal çöküş veya süre aşımı: en kötü puan, kayıt
        o = {"J": -1e9, "hata": repr(e)[:300]}
    finally:
        signal.alarm(0)
    return {"u": list(map(float, u)), "th": th, "seed": seed, **o}


def kucuk_par(ek=""):
    return KUCUK_PAR + (["zayif"] if ek == KURAL3_EK else [])


def buyuk_par(ek=""):
    return [a for a in oz.P_AD if a != "zayif" or ek == KURAL3_EK]


def kucuk_en_iyiler(ek=""):
    """Biten küçük deprem aramalarının en iyi ayarı (büyük aşama için başlangıç noktası)."""
    out = []
    for ag_ad in oz.YEREL_AGLAR:
        f = oz.KALICI / f"ara_kucuk_{ag_ad}{ek}.jsonl"
        if f.exists() and (oz.KALICI / f"bitti_kucuk_{ag_ad}{ek}").exists():
            rs = [json.loads(s) for s in f.read_text().splitlines()]
            b = max(rs, key=lambda r: r["J"])
            th = oz.baslangic_theta()
            th.update({k: b["th"][k] for k in kucuk_par(ek) if k in b["th"]})
            out.append(oz.theta_to_birim(th, buyuk_par(ek)))
    return out


def lhs(n, d, rng):
    u = (rng.permuted(np.tile(np.arange(n), (d, 1)), axis=1).T + rng.uniform(size=(n, d))) / n
    return u


def ara(asama, ag_ad, adlar, n_is, sabit=None, n_lhs=200, n_tur=50, ad_ek=""):
    """Latin hiperküp + CMA benzeri daralma; durdurma: ardışık 3 tur < %1 iyileşme."""
    oz.KALICI.mkdir(parents=True, exist_ok=True)
    f = oz.KALICI / f"ara_{asama}_{ag_ad}{ad_ek}.jsonl"
    sonuc = [json.loads(s) for s in f.read_text().splitlines()] if f.exists() else []
    d = len(adlar)
    rng = np.random.default_rng(12345)
    U = lhs(n_lhs, d, rng)
    U[0] = oz.theta_to_birim(oz.baslangic_theta(), adlar)  # literatür başlangıcı
    if sabit is not None:
        for j, u0 in enumerate(sabit):
            U[1 + j] = u0
    from multiprocessing import get_context
    pool = get_context("fork").Pool(n_is)

    def kos(UU, s0):
        isler = [(asama, ag_ad, u, adlar, s0 + i) for i, u in enumerate(UU)]
        with f.open("a") as fh:
            for r in pool.imap_unordered(_is, isler):
                sonuc.append(r)
                fh.write(json.dumps(r, default=float) + "\n")
                fh.flush()
                log.info("%s %d: J=%.3f %s", asama, len(sonuc), r["J"],
                         {k: r.get(k) for k in ("ardisik", "tuttu", "p_tuttu", "auc", "sahte", "sure_s", "hata", "zaman_hata")
                          if k in r})

    yapilan = {tuple(np.round(r["u"], 9)) for r in sonuc}
    kalan = [u for u in U if tuple(np.round(u, 9)) not in yapilan]
    if kalan:
        kos(kalan, 1000 + len(sonuc))
    en_iyi_gecmis = []
    tur = 0
    while True:
        if asama == "buyuk":
            # küçük deprem aşamasının en iyileri sonradan bitmişse, daralmadan önce başlangıç
            # noktası olarak eklenir (plan sırası: küçükler önce)
            ek = [u for u in kucuk_en_iyiler(ad_ek) if tuple(np.round(u, 9)) not in
                  {tuple(np.round(r["u"], 9)) for r in sonuc}]
            if ek:
                kos(ek, 3000 + len(sonuc))
        J = np.array([r["J"] for r in sonuc])
        best = float(J.max())
        en_iyi_gecmis.append(best)
        if len(en_iyi_gecmis) >= 4:
            son = en_iyi_gecmis[-4:]
            iy = [(son[i + 1] - son[i]) / max(abs(son[i]), 1.0) for i in range(3)]
            if all(x < 0.01 for x in iy):
                log.info("durdurma kuralı: ardışık 3 turda iyileşme < %%1 (%s)", np.round(iy, 4))
                (oz.KALICI / f"bitti_{asama}_{ag_ad}{ad_ek}").write_text(json.dumps(
                    {"n_deneme": len(sonuc), "en_iyi_J": best, "zaman": time.time()}))
                break
        # daralma: en iyi 20'nin ağırlıklı ortalaması ve kovaryansı
        order = np.argsort(-J)[:20]
        X = np.array([sonuc[i]["u"] for i in order])
        w = np.log(len(order) + 0.5) - np.log(np.arange(1, len(order) + 1))
        w /= w.sum()
        mu = (w[:, None] * X).sum(0)
        C = np.cov(X.T, aweights=w) + 1e-4 * np.eye(d)
        sc = 0.7 ** tur
        trng = np.random.default_rng(777 + len(sonuc))
        UU = np.clip(trng.multivariate_normal(mu, C * sc, size=n_tur), 0, 1)
        UU[0] = X[0]  # en iyiyi yeni tohumla yeniden dene (gürültü denetimi)
        kos(UU, 5000 + len(sonuc))
        tur += 1
    pool.close()
    return sonuc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("asama", choices=["kucuk", "buyuk", "sinav", "dene"])
    ap.add_argument("--ag", default="silivri")
    ap.add_argument("--is", dest="n_is", type=int, default=4)
    ap.add_argument("--n-lhs", type=int, default=200)
    ap.add_argument("--en-iyi", type=int, default=5)
    ap.add_argument("--tohum", type=int, default=3)
    ap.add_argument("--ileri", type=int, default=20)
    ap.add_argument("--surum", default="ana", choices=["ana", "tam_yama", "oncul"])
    ap.add_argument("--ek", default="")
    ap.add_argument("--olcut", default="eski", choices=["eski", "zaman"],
                    help="küçük aşama ayar ölçütü; 'zaman' = göreli zaman hatası en az")
    ap.add_argument("--sira", type=int, nargs="*", default=None, help="yalnız bu ayar sıraları")
    ap.add_argument("--cikti-ek", default="", help="sınav çıktı dosya adına ek (yeniden koşular için)")
    args = ap.parse_args()
    oz.KUCUK_OLCUT = args.olcut
    oz.KALICI.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(oz.KALICI / f"log_{args.asama}_{args.ag}.txt"),
                                  logging.StreamHandler()])
    if args.asama == "dene":
        asama = "kucuk" if args.ag in oz.YEREL_AGLAR else "buyuk"
        th = oz.baslangic_theta()
        r = degerlendir(asama, args.ag, th, 1)
        log.info("deneme %s %s: %s", asama, args.ag, r["ozet"])
        k = pd.DataFrame(r["kayit"])
        if len(k):
            log.info("\n%s", k.drop(columns=[c for c in ("ad",) if c in k]).describe().T.to_string())
        return
    if args.asama == "sinav":
        sinav(args)
        return
    if args.asama == "kucuk":
        ara("kucuk", args.ag, kucuk_par(args.ek), args.n_is, n_lhs=args.n_lhs, ad_ek=args.ek)
    elif args.asama == "buyuk":
        ara("buyuk", oz.BUYUK_AG, buyuk_par(args.ek), args.n_is, sabit=kucuk_en_iyiler(args.ek), n_lhs=args.n_lhs,
            ad_ek=args.ek)



# ------------------------------------------------------------------ sınav ve bugünkü durum
ISTANBUL_PARCALARI = {"Orta Marmara çukuru": (27.75, 28.05), "Kumburgaz": (28.05, 28.60),
                      "Adalar": (28.60, 29.10), "İstanbul önü toplam": (27.75, 29.10)}


def km_tablosu(m, ag, ileri_p=None):
    """Her kol için 1 km'lik dilimlerde kilitli kutuların ortalama doluluğu."""
    c = ag.cells
    dol = oz.doluluk(m)
    rows = []
    for kol in sorted(c.kol.unique()):
        s = ((c.kol == kol) & (c.katman == "kabuk")).to_numpy() & m.vw
        x = m.x[s]
        b = np.floor(x).astype(int)
        df = pd.DataFrame({"km": b, "f": dol[s], "w": m.alan[s], "lon": c.lon.to_numpy()[s],
                           "lat": c.lat.to_numpy()[s], "seg": c.seg_ad.to_numpy()[s]})
        for k, g in df.groupby("km"):
            rows.append({"kol": int(kol), "km": int(k), "segment": g.seg.iloc[0], "boylam": round(g.lon.mean(), 4),
                         "enlem": round(g.lat.mean(), 4), "doluluk": float(np.average(g.f, weights=g.w)),
                         "kilitli_kutu": int(len(g))})
    t = pd.DataFrame(rows)
    if ileri_p is not None:
        t = t.merge(ileri_p, on=["kol", "km"], how="left")
    return t


def ileri_bakis(m, ufuk=30.0, n=20, seed=0, mw_min=6.8):
    """Bugünden ileri: her kol ve km için ``ufuk`` yıl içinde M>=mw_min kırılmasına girme olasılığı.

    Bir km'nin "vurulması": o kol ve km'deki kilitli kutulardan en az birinin o depremde
    en büyük kaymanın %5'inden fazla kayması (hücre hücre; x0-x1 aralığı kollar arasında
    karıştığı için kullanılmaz).
    """
    st = m.durum()
    n0 = len(m.olaylar)
    kol = m.c["kol"].to_numpy()
    km = np.floor(m.x).astype(int)
    hits = {}
    ilk = []
    eski = m.kayit_slip
    for r in range(n):
        m.yukle(st)
        m.dayanimi_yenile(seed * 1000 + r)
        k0 = len(m.olaylar)
        m.kayit_slip = (mw_min, [])
        t0 = m.t
        while m.t < t0 + ufuk - 1e-9:
            m.ilerle(min(m.t + 1.0, t0 + ufuk))
        vurulan = set()
        for t, mw, d in m.kayit_slip[1]:
            h = d > 0.05 * d.max()
            kk = set(zip(kol[h].tolist(), km[h].tolist()))
            ilk.append({"tekrar": r, "t": t, "mw": mw, "kollar": sorted({a for a, _ in kk}),
                        "x0": float(m.x[h & (kol == 0)].min()) if (h & (kol == 0)).any() else np.nan,
                        "x1": float(m.x[h & (kol == 0)].max()) if (h & (kol == 0)).any() else np.nan})
            vurulan |= kk
        for k in vurulan:
            hits[k] = hits.get(k, 0) + 1
        del m.olaylar[k0:]
    m.kayit_slip = eski
    m.yukle(st)
    del m.olaylar[n0:]
    p = pd.DataFrame([{"kol": a, "km": b, "p_30yil": v / n} for (a, b), v in hits.items()],
                     columns=["kol", "km", "p_30yil"])
    return p, pd.DataFrame(ilk)


def sinav(args):
    """En iyi ayarları donmuş olarak sınav döneminden geçir (her biri birkaç başlangıç tohumu)."""
    out = oz.KALICI / "sinav"
    out.mkdir(parents=True, exist_ok=True)
    asama = "kucuk" if args.ag in oz.YEREL_AGLAR else "buyuk"
    f = oz.KALICI / f"ara_{asama}_{args.ag if asama == 'kucuk' else oz.BUYUK_AG}{args.ek}.jsonl"
    rs = [json.loads(s) for s in f.read_text().splitlines()]
    rs = sorted([r for r in rs if r["J"] > -1e8], key=lambda r: -r["J"])[:args.en_iyi]
    ag_ad = args.ag
    for k, r in enumerate(rs):
        if args.sira is not None and k not in args.sira:
            continue
        for seed in range(args.tohum):
            ad = f"{asama}_{ag_ad}{args.ek}_{k}_{seed}_{args.surum}{args.cikti_ek}"
            if (out / f"{ad}.pkl").exists():
                continue
            th = r["th"]
            res = degerlendir(asama, ag_ad, th, 900 + seed, sinav=True, tam_yama=args.surum == "tam_yama",
                              sonsal=args.surum != "oncul")
            m = res.pop("m")
            kay = {"th": th, "J_ayar_arama": r["J"], "seed": seed, "surum": args.surum, **res}
            if asama == "buyuk":
                p, ilk = ileri_bakis(m, 30.0, args.ileri, seed=seed)
                kay["km"] = km_tablosu(m, _G[(asama, ag_ad, True)][0], p)
                kay["ileri_p"], kay["ileri_ilk"] = p, ilk
                kay["durum"] = m.durum()
                kay["hucre"] = {"S": m.S.astype("float32"), "S_fail": m.S_fail.astype("float32"),
                                "vw": m.vw, "vs": m.vs, "viskoz": m.viskoz}
            (out / f"{ad}.pkl").write_bytes(__import__("pickle").dumps(kay))
            log.info("sınav %s: %s", ad, res["ozet"])


if __name__ == "__main__":
    main()
