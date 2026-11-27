"""İleri-geri eğitim (Ulaş'ın yöntemi): depremler en eskisinden başlayarak sırayla işlenir.

"sırayla yavaş yavaş ileri gideceksin deprem deprem"; "beklenen zamanda deprem olmazsa geriye gidip bir
önceki depreme sebep olan parametreyi yine depreme sebep olacak şekilde bulacaksın".

Düzen (bir pencere için):
  - Başlangıç: o pencerenin 100+ takımlık aramasının en iyisi (yoksa literatür değerleri).
  - Eğitim depremleri (1964-2017) zamana göre sıralanır, BLOK'ar deprem ileri gidilir.
  - Her adımda model en eski depremden o ana kadar yeniden koşulur (parametreler pencere başına olduğu
    için önceki depremlerin hâlâ oluşup oluşmadığı ancak böyle sınanır).
  - Yanlış alarmlar da düzeltilir: o ana kadar gözlenenden 0,5 büyük kendiliğinden deprem çıkarsa
    parametreler değiştirilir; bir değişiklik yanlış alarm sayısını artıramaz.
  - Yeni bloktaki depremlerden tutmayan varsa (mutlak zaman hatası > ESIK yıl; başlangıçta ±1 yıl,
    sonra kademeli daraltılır: ±6 ay, ±1 ay) parametreler küçük adımlarla
    değiştirilir. Bir değişiklik ancak (1) o ana kadar tutan deprem sayısını azaltmıyor ve (2) o ana
    kadarki ortalama zaman hatasını düşürüyorsa kabul edilir. Sonra bir sonraki bloğa geçilir.
  - 2018-2026 sınavı bu sürece hiç girmez; sonunda bir kez ölçülür (±1 yıl, ±1 ay, ±1 gün ve şans).
Kayıt: KALICI/ileri_geri_{ag}.jsonl (her blok bir satır; kesintide kaldığı bloktan sürer).

Kullanım: python scripts/ileri_geri.py --ag d0_105 [--blok 10 --esik 0.5 --tur 3 --is 4]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from deprem.mekanik import ozumseme as oz  # noqa: E402
import mekanik_ikinci_tur as M  # noqa: E402

oz.KUCUK_OLCUT = "zaman"
PAR = M.KUCUK_PAR
CEZA = 10.0  # yıl: kaçırılan deprem kadar (TAVAN) fazladan üretilen deprem de ceza alır
TAVAN = 10.0  # yıl: tahmini olmayan ya da 10 yıldan çok sapan deprem bu kadar hata sayılır
_G = {}


def _hazir(ag_ad):
    if _G.get("ad") != ag_ad:
        _G["ad"] = ag_ad
        _G["ag"], _G["ol"], _ = M._hazirla("kucuk", ag_ad, False)
    return _G["ag"], _G["ol"]


def kos(arg):
    """Modeli en eski depremden t_son'a kadar koşar; eğitim depremlerinin göreli zaman hataları."""
    ag_ad, u, t_son = arg
    ag, ol = _hazir(ag_ad)
    a = M.AYAR["kucuk"]
    th = dict(oz.baslangic_theta(), zayif=0.0)
    th.update(oz.birim_to_theta(np.asarray(u), PAR))
    try:
        m, kayit = oz.Kosu(ag, th, ol, 1, a["isinma"], a["t_bas"], t_son, a["ayar_bitis"]).calistir()
    except Exception as e:  # sayısal çöküş
        return {"hata": repr(e)[:200]}
    d = pd.DataFrame(kayit)
    d = d[(d.puan == "ayar") & (d.tip == "kucuk")]
    # mutlak zaman hatası (yıl): göreli hata × o yerdeki tekrar aralığı; ıska (göreli >= 3) = TAVAN
    mut = np.where(d.zaman_hata.fillna(3.0) < 3.0, d.zaman_hata.fillna(3.0) * d.dT, TAVAN)
    # yanlış alarm (Ulaş: "modelde büyük deprem çıktı ama gerçekte orada öyle deprem yok"): o ana kadar
    # gözlenen en büyük depremden 0,5 büyük kendiliğinden deprem
    mmax = max([o.mw for o in ol if o.t <= t_son and o.tip == "kucuk"] + [3.0])
    sahte = [e["t"] for e in m.olaylar if not e.get("zorla") and e["mw"] > mmax + 0.5 and e["t"] <= t_son]
    return {"h": dict(zip(d.t.round(6), np.minimum(mut, TAVAN))), "sahte": sahte}


def olc(r, tlar, esik):
    """(puan, tutan, yanlış alarm): puan = ortalama zaman hatası (yıl) + yanlış alarm başına CEZA yıl."""
    v = np.array([r["h"].get(round(t, 6), TAVAN) for t in tlar])
    ns = len(r["sahte"])
    return float(v.mean()) + CEZA * ns / max(len(tlar), 1), int((v <= esik).sum()), ns


def baslangic(ag_ad):
    f = oz.KALICI / f"ara_kucuk_{ag_ad}_zh.jsonl"
    if f.exists() and f.stat().st_size:
        rs = [json.loads(x) for x in f.read_text().splitlines()]
        b = max(rs, key=lambda r: r["J"])
        th = dict(oz.baslangic_theta(), **{k: b["th"][k] for k in PAR if k in b["th"]})
        return oz.theta_to_birim(th, PAR), "arama"
    return oz.theta_to_birim(oz.baslangic_theta(), PAR), "literatur"


def egit(args):
    from multiprocessing import get_context
    ag, ol = _hazir(args.ag)
    a = M.AYAR["kucuk"]
    tlar = sorted(o.t for o in ol if o.tip == "kucuk" and o.puan == "ayar")
    log = oz.KALICI / f"ileri_geri_{args.ag}{args.ek}.jsonl"
    rng = np.random.default_rng(11)
    if log.exists() and log.stat().st_size:
        son = json.loads(log.read_text().splitlines()[-1])
        u, k0, kay = np.array(son["u"]), son["k"], son["kaynak"]
    else:
        (u, kay), k0 = baslangic(args.ag), 0
    pool = get_context("fork").Pool(args.n_is)
    for k in range(k0 + args.blok, len(tlar) + args.blok, args.blok):
        k = min(k, len(tlar))
        t_son = min(tlar[k - 1] + 1e-3, a["ayar_bitis"])
        pre = tlar[:k]
        r = kos((args.ag, u, t_son))
        if "hata" in r:
            break
        F, n_tut, n_s = olc(r, pre, args.esik)
        yeni = [t for t in tlar[k - args.blok:k] if r["h"].get(round(t, 6), TAVAN) > args.esik]
        yeni_sahte = [t for t in r["sahte"] if t > (tlar[k - args.blok - 1] if k - args.blok > 0 else -1e9)]
        F0, kabul, sig = F, 0, args.adim
        for _ in range(args.tur if (yeni or yeni_sahte) else 0):
            ad = [np.clip(u + rng.normal(0, sig, len(u)), 0, 1) for _ in range(args.n_is)]
            rs = pool.map(kos, [(args.ag, x, t_son) for x in ad])
            iyi = None
            for x, rr in zip(ad, rs):
                if "hata" in rr:
                    continue
                f, n, ns = olc(rr, pre, args.esik)
                if n >= n_tut and ns <= n_s and f < F - 1e-4 and (iyi is None or f < iyi[1]):
                    iyi = (x, f, n, ns)
            if iyi is None:
                sig *= 0.6
                continue
            u, F, n_tut, n_s = iyi
            kabul += 1
        rec = {"k": k, "t": t_son, "u": list(map(float, u)), "th": oz.birim_to_theta(u, PAR), "F_once": F0,
               "F": F, "tutan": n_tut, "n": k, "yanlis_alarm": n_s, "yeni_tutmayan": len(yeni),
               "yeni_yanlis_alarm": len(yeni_sahte), "kabul": kabul, "kaynak": kay}
        with open(log, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        print(f"{args.ag} {k}/{len(tlar)} t={t_son:.2f} hata {F0:.3f}->{F:.3f} tutan {n_tut}/{k} yanlış {n_s} kabul {kabul}", flush=True)
    pool.close()
    (oz.KALICI / f"bitti_ileri_geri_{args.ag}{args.ek}").write_text(json.dumps({"k": len(tlar)}))
    sinav(args)


def sinav(args):
    """Son parametrelerle 2018-2026 sınavı: mutlak zaman hatası ±1 gün/ay/yıl içinde tutan oranı ve şans."""
    log = oz.KALICI / f"ileri_geri_{args.ag}{args.ek}.jsonl"
    son = json.loads(log.read_text().splitlines()[-1])
    out = {"ag": args.ag}
    for ad, u in [("ileri_geri", np.array(son["u"])), ("baslangic", np.array(json.loads(log.read_text().splitlines()[0])["u"]))]:
        th = dict(oz.baslangic_theta(), zayif=0.0, **oz.birim_to_theta(u, PAR))
        r = M.degerlendir("kucuk", args.ag, th, 1, sinav=True)
        d = pd.DataFrame(r["kayit"])
        d = d[(d.puan == "sinav") & (d.tip == "kucuk")]
        mut = np.where(d.zaman_hata < 3, d.zaman_hata * d.dT, np.inf)  # yıl
        dT = d.dT.to_numpy()
        for etik, w in [("1gun", 1 / 365.25), ("1ay", 1 / 12), ("1yil", 1.0)]:
            out[f"{ad}_tuttu_{etik}"] = float(np.mean(mut <= w)) if len(d) else None
            # şans: aynı aralıkta (önceki depremden sonra dT boyunca) rastgele zaman tahmini
            out[f"sans_{etik}"] = float(np.mean(np.minimum(2 * w / np.maximum(dT, 1e-9), 1.0))) if len(d) else None
        out[f"{ad}_medyan_goreli_hata"] = float(np.median(np.minimum(d.zaman_hata.fillna(3), 3))) if len(d) else None
        out["n_sinav"] = int(len(d))
    (oz.KALICI / f"sinav_ileri_geri_{args.ag}{args.ek}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ag", required=True)
    ap.add_argument("--blok", type=int, default=10)
    ap.add_argument("--esik", type=float, default=1.0, help="tutma payı, yıl (Ulaş: başlangıçta ±1 yıl)")
    ap.add_argument("--tur", type=int, default=3)
    ap.add_argument("--adim", type=float, default=0.08)
    ap.add_argument("--is", dest="n_is", type=int, default=4)
    ap.add_argument("--ek", default="", help="kayıt adı eki (ör. pay kademesi: _6ay)")
    ap.add_argument("--yalniz-sinav", action="store_true")
    a = ap.parse_args()
    sinav(a) if a.yalniz_sinav else egit(a)


if __name__ == "__main__":
    main()
