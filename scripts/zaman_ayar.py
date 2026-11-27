"""Zaman hatasını en aza indiren parametre ayarı ve birini-dışarıda-bırak sınavı.

Ulaş: "bana en az zaman hatası veren parametreleri ayarla".

Düzen
  - Hedefler: zaman_payi.hedefler() ile seçilen 4 büyük deprem (1754, 1894, 1912, 1999).
  - Her ayar takımı için 4 depremin her biri zaman_payi.bir_hedef ile koşulur (N tekrar).
  - Ayar puanı (küçük = iyi), bir deprem için:
        zaman hatası = tekrarların medyanı |tahmin - gerçek| (deprem üretilmezse 3·T = 750 yıl)
        + YANLIS_CEZA yıl × tekrar başına yanlış alarm (±%30 içinde gerçek depremle eşleşmeyen M>=6,8)
        + MW_CEZA yıl × işlenen gerçek büyük depremlerde ortalama büyüklük eksiği (model çok zayıf kalmasın)
  - Birini dışarıda bırak: her deprem için ayar, yalnız kalan 3 depremin ortalama puanıyla seçilir;
    seçilen ayarın dışarıdaki depremdeki hatası raporlanır. Dışarıdaki deprem seçimi hiç etkilemez.
  - Kıyas: aynı aralıkta şans (Poisson) ve basit tekrarlama (BPT) modelinin pencereye düşme olasılığı;
    ayrıca ayarsız "literatur" başlangıç takımı.
  - Krip varsayımı: --varsayim (varsayılan "kilitli", Chartier vd. 2021'in en iyi bulduğu).
    Bu varsayımda kayma hızı ve derinlik sabit (Chartier), bu yüzden v_pl_marmara, kalinlik_u,
    sig_krip_z aranmaz.

Kullanım: python scripts/zaman_ayar.py kos --n-ayar 48 --n 4 --is 4
          python scripts/zaman_ayar.py ozet
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
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
import zaman_payi as Z  # noqa: E402

Z.CIKTI = oz.KALICI / "zaman_ayar"
ARANAN = ["mu_eff", "dtau0", "sigma_h", "eta_dyn", "kalan", "kayma_ref", "tM_ref", "krip_carpan",
          "sigma_grad", "zayif", "a_rs"]
YANLIS_CEZA = 50.0
MW_CEZA = 200.0  # yıl / büyüklük birimi: işlenen gerçek depremlerde model büyüklüğünün eksik kalması
ISKA = 3 * Z.T_ORT


def ayar_takimlari(n, seed=7):
    tab = {p[0]: p for p in oz.PARAMETRELER}
    u = M.lhs(n, len(ARANAN), np.random.default_rng(seed))
    out = {"literatur": dict(oz.baslangic_theta(), zayif=0.0)}
    for k in range(n):
        th = dict(oz.baslangic_theta())
        for j, ad in enumerate(ARANAN):
            _, lo, hi, _, log = tab[ad]
            th[ad] = float(np.exp(np.log(lo) + u[k, j] * (np.log(hi) - np.log(lo))) if log
                           else lo + u[k, j] * (hi - lo))
        out[f"a{k:03d}"] = th
    # Ulaş'ın sırası: küçük depremlerle (zaman hatası en az) eğitilen her kesim penceresinin en iyi
    # ayarı da aday olarak denenir (küçük aşamada aranmayan mu_eff, zayif başlangıç değerinde)
    for ag_ad in oz.YEREL_AGLAR:
        f = oz.KALICI / f"ara_kucuk_{ag_ad}_zh.jsonl"
        if f.exists():
            rs = [json.loads(x) for x in f.read_text().splitlines()]
            b = max(rs, key=lambda r: r["J"])
            out[f"kucuk_{ag_ad}"] = dict(oz.baslangic_theta(), zayif=0.0, **{k: v for k, v in b["th"].items()
                                                                              if k in ARANAN})
    return out


def kos(args):
    Z.CIKTI.mkdir(parents=True, exist_ok=True)
    tak = ayar_takimlari(args.n_ayar)
    (Z.CIKTI / "takimlar.json").write_text(json.dumps(tak, indent=1))
    ag, ol, _ = M._hazirla("buyuk", oz.BUYUK_AG, True)
    hed = Z.hedefler(ol)
    isler = [(args.varsayim, ad, th, i, j, args.n) for ad, th in tak.items() for i, j in hed]
    from multiprocessing import get_context
    with get_context("fork").Pool(args.n_is) as p:
        list(p.imap_unordered(Z.bir_hedef, isler, chunksize=1))


def puanla(r):
    t = r["tahmin"]
    hata = np.where(np.isfinite(t), np.abs(t - r["t_gercek"]), ISKA)
    A = r.get("alarmlar")
    w = 0.30 * Z.T_ORT
    yanlis = 0
    if A is not None and len(A):
        yanlis = sum(not any(Z._eslesir(a, g, w) for g in r["gercek"]) for a in A.to_dict("records"))
    mu = r.get("mw_uyum") or []
    mw_eksik = float(np.mean([max(g - mm, 0.0) for _, mm, g in mu])) if mu else 0.0
    return {"hata_medyan": float(np.median(hata)), "yanlis_tekrar_basi": yanlis / r["n"], "mw_eksik": mw_eksik,
            "puan": float(np.median(hata)) + YANLIS_CEZA * yanlis / r["n"] + MW_CEZA * mw_eksik,
            **{f"tuttu_{int(p * 100)}": float(np.mean(np.abs(t - r["t_gercek"]) <= p * Z.T_ORT))
               for p in Z.PAYLAR}}


def ozet(args):
    rows = []
    for f in sorted(Z.CIKTI.glob("*.pkl")):
        r = pickle.loads(f.read_bytes())
        rows.append({"ayar": r["ayar"], "hedef": r["hedef"], "aralik": r["t_gercek"] - r["t_onceki"], **puanla(r)})
    d = pd.DataFrame(rows)
    hedefler = sorted(d.hedef.unique())
    tam = d.groupby("ayar").hedef.nunique()
    d = d[d.ayar.isin(tam[tam == len(hedefler)].index)]
    P = d.pivot(index="ayar", columns="hedef", values="puan")
    tak = json.loads((Z.CIKTI / "takimlar.json").read_text())

    # tüm 4 depremle en iyi takım (Ulaş'ın istediği "en az zaman hatası"; bu eğitim puanıdır)
    en_iyi = P.mean(axis=1).idxmin()
    # birini dışarıda bırak
    loo = []
    for h in hedefler:
        sec = P.drop(columns=h).mean(axis=1).idxmin()
        x = d[(d.ayar == sec) & (d.hedef == h)].iloc[0]
        lit = d[(d.ayar == "literatur") & (d.hedef == h)].iloc[0]
        rec = {"disarida": h, "secilen": sec, "egitim_puani": P.drop(columns=h).loc[sec].mean(),
               "sinav_hata_yil": x.hata_medyan, "sinav_yanlis_alarm": x.yanlis_tekrar_basi,
               "literatur_hata_yil": lit.hata_medyan}
        for p in Z.PAYLAR:
            poi, bpt = Z.kiyas_olasilik(x.aralik, p * Z.T_ORT)
            rec[f"tuttu_{int(p * 100)}"] = x[f"tuttu_{int(p * 100)}"]
            rec[f"sans_{int(p * 100)}"] = poi
            rec[f"bpt_{int(p * 100)}"] = bpt
        loo.append(rec)
    L = pd.DataFrame(loo)
    tab = ROOT / "reports" / "tables"
    d.to_csv(tab / "zaman_ayar_tum.csv", index=False, float_format="%.3f")
    L.to_csv(tab / "zaman_ayar_loo.csv", index=False, float_format="%.3f")
    egitim = d[d.ayar == en_iyi].set_index("hedef")[["hata_medyan", "yanlis_tekrar_basi"]]
    out = {"en_iyi_takim": en_iyi, "degerler": tak[en_iyi], "egitim_ort_puan": float(P.loc[en_iyi].mean()),
           "egitim_hatalari": egitim.to_dict(), "literatur_ort_puan": float(P.loc["literatur"].mean()),
           "n_takim": int(len(P)), "loo_ort_sinav_hata": float(L.sinav_hata_yil.mean()),
           "loo_ort_egitim_puan": float(L.egitim_puani.mean())}
    (tab / "zaman_ayar_ozet.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    pd.set_option("display.width", 250)
    print(P.mean(axis=1).sort_values().head(10).round(1))
    print(L.round(3).to_string())
    print(json.dumps(out, indent=1, ensure_ascii=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("asama", choices=["kos", "ozet"])
    ap.add_argument("--varsayim", default="kilitli")
    ap.add_argument("--n-ayar", type=int, default=48)
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--is", dest="n_is", type=int, default=4)
    a = ap.parse_args()
    {"kos": kos, "ozet": ozet}[a.asama](a)


if __name__ == "__main__":
    main()
