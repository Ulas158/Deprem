"""Zaman hata payı sınavı (Ulaş'ın ölçütü): "deprem mutlaka olsun, zamanı hata payıyla".

Her büyük deprem E (M>=6,8, aynı kesimde bir önceki büyük depremi katalogda olanlar) için:
  1. Model, E dışındaki bütün gerçek depremleri (E'nin tarihinden öncekiler) işleyerek,
     aynı kesimdeki bir önceki büyük depremin (t_onceki) hemen sonrasına getirilir.
  2. Buradan N kez (her biri farklı küçük ölçekli dayanım farklarıyla) ileri koşulur; E'nin
     tarihinden önceki diğer gerçek depremler yolda işlenmeye devam eder, E işlenmez.
  3. Her koşuda modelin E kesiminde kendiliğinden ürettiği ilk benzer deprem
     (M >= E.mw - 0,4, gerçek kırık boyunun en az %50'siyle örtüşen) "tahmin zamanı"dır.
     t_onceki + 3·T içinde yoksa ıska sayılır.
  4. Puan: tahmin zamanı gerçek tarihin ±p·T içinde mi (p = %1, %10, %20, %30; T = 250 yıl
     ortalama tekrar aralığı varsayımı).
  Kıyaslar (aynı t_onceki'den başlayarak, modelsiz):
     şans      : Poisson (üstel, ortalama T)
     tekrarlama: BPT / ters Gauss (ortalama T, düzensizlik 0,5)
  Alarm ölçütü (Ulaş: "%100 dediğimiz yerde deprem olsun, zaman kayabilir"): modelin kendiliğinden
  kırdığı her M>=6,8 yer bir alarmdır. Alarm, ±p·T içinde aynı yerde (>= %50 örtüşme) gerçek bir
  büyük depremle eşleşirse doğru, eşleşmezse yanlış alarmdır. Şans: aynı alarmlar, zamanı rastgele.
  Ayrıca: t_onceki ile E arasında modelin ürettiği, gerçekte karşılığı olmayan M>=6,8 depremler
  (sahte) sayılır.

Kullanım: python scripts/zaman_payi.py kos --varsayim kilitli --ayar literatur --n 16 --is 4
          python scripts/zaman_payi.py ozet
"""
from __future__ import annotations

import argparse
import dataclasses
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
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
import mekanik_ikinci_tur as M  # noqa: E402
import uzun_vade as U  # noqa: E402

CIKTI = oz.KALICI / "zaman_payi"
T_ORT = 250.0
PAYLAR = (0.01, 0.10, 0.20, 0.30)
ALFA_BPT = 0.5
ISINMA = 400.0


def hedefler(ol, en_az_ortusme=0.5):
    """(hedef, önceki) çiftleri. Önceki: hedef kırığının en az yarısını kaplayan en yakın M>=6,8 büyük
    deprem (eskiden en ufak örtüşme yetiyordu: 1999 İzmit'in öncesi 1943 Hendek sayılıyordu, oysa
    1943 1999 kırığının yalnız %24'ünü kaplar)."""
    out = []
    for i, o in enumerate(ol):
        if o.tip != "buyuk" or o.mw < 6.8:
            continue
        L = max(o.x1 - o.x0, 1.0)
        for j in range(i - 1, -1, -1):
            q = ol[j]
            if q.tip == "buyuk" and q.mw >= 6.8 and (min(q.x1, o.x1) - max(q.x0, o.x0)) / L >= en_az_ortusme:
                out.append((i, j))
                break
    return out


def _sessiz(o):
    return dataclasses.replace(o, puan="")


def bir_hedef(arg):
    varsayim, ayar_ad, th, i, j, n = arg
    ag, ol, _ = M._hazirla("buyuk", oz.BUYUK_AG, True)
    E, P = ol[i], ol[j]
    f = CIKTI / f"{varsayim}_{ayar_ad}_{E.ad.replace(' ', '_')}.pkl"
    if f.exists():
        return
    t0 = time.time()
    asil = oz.model_kur
    if varsayim != "eski":
        oz.model_kur = lambda ag_, th_, seed_: U.kur(ag_, th_, seed_, varsayim)
    try:
        diger = [_sessiz(o) for k, o in enumerate(ol) if k != i and o.t < E.t]
        k_onceki = next(k for k, o in enumerate(diger) if o.ad == P.ad and abs(o.t - P.t) < 1e-6)
        t_son = min(P.t + 3 * T_ORT, 2026.9)
        bas = oz.Kosu(ag, th, diger[:k_onceki + 1], 4242, ISINMA, M.AYAR["buyuk"]["t_bas"], P.t + 1e-4, 0.0, sonsal=False,
                      karsi_olgu=False)
        m, _ = bas.calistir()
        st, n0 = m.durum(), len(m.olaylar)
        # büyüklük uyumu: işlenen gerçek büyük depremlerde modelin ürettiği büyüklük
        gmw = {o.ad: o.mw for o in ol if o.tip == "buyuk"}
        mw_uyum = [(e["tarih"], float(e["mw"]), gmw[e["tarih"]]) for e in m.olaylar
                   if e.get("zorla") and e.get("tarih") in gmw]
        tahmin, sahte, alarmlar = [], [], []
        for r in range(n):
            m.yukle(st)
            del m.olaylar[n0:]
            m.dayanimi_yenile(10_000 + r)
            ileri = oz.Kosu(ag, th, diger, 4242, ISINMA, M.AYAR["buyuk"]["t_bas"], t_son, 0.0, sonsal=False, karsi_olgu=False)
            ileri.m = m
            # bas koşusunun durumundan devam: model_kur'u atlamak için calistir'a durum verilir
            m2, _ = ileri.calistir(durum=(m, k_onceki + 1, []))
            tt = np.nan
            for e in m2.olaylar[n0:]:
                if e.get("zorla") or e["mw"] < E.mw - 0.4 or e["t"] <= P.t:
                    continue
                if oz._ortusme(e, E) >= 0.5:
                    tt = e["t"]
                    break
            tahmin.append(tt)
            # Ulaş'ın ölçütü: modelin "%100 dolu" deyip kendiliğinden kırdığı her yer bir alarmdır
            for e in m2.olaylar[n0:]:
                if not e.get("zorla") and e["mw"] >= 6.8 and e["t"] > P.t:
                    alarmlar.append({"tekrar": r, "t": e["t"], "mw": e["mw"], "x0": e["x0"], "x1": e["x1"]})
            gercek = [o for o in ol if o.tip == "buyuk" and P.t < o.t < E.t]
            for e in m2.olaylar[n0:]:
                if e.get("zorla") or e["mw"] < 6.8 or not (P.t < e["t"] < E.t):
                    continue
                if not any(oz._ortusme(e, g) >= 0.3 for g in gercek):
                    sahte.append({"tekrar": r, "t": e["t"], "mw": e["mw"], "x0": e["x0"], "x1": e["x1"]})
    finally:
        oz.model_kur = asil
    out = {"varsayim": varsayim, "ayar": ayar_ad, "hedef": E.ad, "t_gercek": E.t, "mw": E.mw, "onceki": P.ad,
           "t_onceki": P.t, "tahmin": np.array(tahmin), "sahte": pd.DataFrame(sahte), "n": n,
           "alarmlar": pd.DataFrame(alarmlar), "t_son": t_son, "mw_uyum": mw_uyum,
           "gercek": [{"ad": o.ad, "t": o.t, "x0": o.x0, "x1": o.x1, "uc_km": o.uc_km}
                      for o in ol if o.tip == "buyuk" and o.mw >= 6.5],
           "sn": round(time.time() - t0, 1)}
    tmp = f.with_suffix(".tmp")
    tmp.write_bytes(pickle.dumps(out))
    tmp.replace(f)
    print(f.name, out["sn"], "s", np.round(out["tahmin"], 1), flush=True)


def kos(args):
    CIKTI.mkdir(parents=True, exist_ok=True)
    ag, ol, _ = M._hazirla("buyuk", oz.BUYUK_AG, True)
    th = U.ayarlar()[args.ayar]
    isler = [(args.varsayim, args.ayar, th, i, j, args.n) for i, j in hedefler(ol)]
    print("hedefler:", [(ol[i].ad, ol[j].ad) for i, j in hedefler(ol)], flush=True)
    from multiprocessing import get_context
    with get_context("fork").Pool(args.n_is) as p:
        list(p.imap_unordered(bir_hedef, isler))


def kiyas_olasilik(dt_gercek, w):
    """t_onceki'den itibaren, olayın [dt-w, dt+w] içinde olma olasılığı."""
    a, b = max(dt_gercek - w, 0.0), dt_gercek + w
    poi = stats.expon(scale=T_ORT)
    bpt = stats.invgauss(ALFA_BPT ** 2, scale=T_ORT / ALFA_BPT ** 2)
    return poi.cdf(b) - poi.cdf(a), bpt.cdf(b) - bpt.cdf(a)


def _eslesir(a, g, w):
    """Alarm a, gerçek deprem g ile yer (>= %50 örtüşme) ve zaman (±w yıl) olarak eşleşiyor mu."""
    a0, a1 = g["x0"] - g["uc_km"], g["x1"] + g["uc_km"]
    L = max(g["x1"] - g["x0"], 1.0)
    ort = max(0.0, min(a["x1"], a1) - max(a["x0"], a0))
    return ort / L >= 0.5 and abs(a["t"] - g["t"]) <= w


def alarm_orani(r, w, rng, n_sans=200):
    """Modelin alarmlarından (kendiliğinden M>=6,8) kaçı ±w içinde gerçek bir depremle eşleşti; şans:
    aynı alarmlar aynı yerde, zamanı [t_onceki, t_son] içinde rastgele."""
    A = r.get("alarmlar")
    if A is None or len(A) == 0:
        return np.nan, np.nan, 0
    G = r["gercek"]
    kay = A.to_dict("records")
    dogru = np.mean([any(_eslesir(a, g, w) for g in G) for a in kay])
    sans = []
    for _ in range(n_sans):
        tt = rng.uniform(r["t_onceki"], r["t_son"], len(kay))
        sans.append(np.mean([any(_eslesir({**a, "t": t}, g, w) for g in G) for a, t in zip(kay, tt)]))
    return float(dogru), float(np.mean(sans)), len(kay)


def ozet(args):
    rows = []
    rng = np.random.default_rng(0)
    for f in sorted(CIKTI.glob("*.pkl")):
        r = pickle.loads(f.read_bytes())
        dt = r["t_gercek"] - r["t_onceki"]
        for p in PAYLAR:
            w = p * T_ORT
            t = r["tahmin"]
            model = float(np.mean(np.abs(t - r["t_gercek"]) <= w)) if len(t) else np.nan
            poi, bpt = kiyas_olasilik(dt, w)
            dogru, sans_dogru, n_alarm = alarm_orani(r, w, rng)
            rows.append({"varsayim": r["varsayim"], "ayar": r["ayar"], "hedef": r["hedef"],
                         "aralik_yil": round(dt, 1), "pay": f"±%{int(p * 100)} (±{w:g} yıl)",
                         "model": model, "sans_poisson": poi, "tekrarlama_bpt": bpt,
                         "uretti_oran": float(np.mean(np.isfinite(t))),
                         "tahmin_medyan_hata_yil": float(np.nanmedian(t - r["t_gercek"])) if np.isfinite(t).any() else np.nan,
                         "sahte_tekrar_basi": len(r["sahte"]) / r["n"],
                         "alarm_tekrar_basi": n_alarm / r["n"],
                         "alarm_dogru": dogru, "alarm_yanlis": 1 - dogru if np.isfinite(dogru) else np.nan,
                         "alarm_dogru_sans": sans_dogru})
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "reports" / "tables" / "zaman_payi.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 250)
    print(d.round(3).to_string())
    print(d.groupby(["varsayim", "ayar", "pay"])[["model", "sans_poisson", "tekrarlama_bpt",
                                                       "alarm_dogru", "alarm_dogru_sans"]].mean().round(3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("asama", choices=["kos", "ozet"])
    ap.add_argument("--varsayim", default="kilitli")
    ap.add_argument("--ayar", default="literatur")
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--is", dest="n_is", type=int, default=4)
    a = ap.parse_args()
    {"kos": kos, "ozet": ozet}[a.asama](a)


if __name__ == "__main__":
    main()
