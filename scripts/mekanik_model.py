"""Marmara mekanik kırılma modeli: ayar, tarihi oynatma, sınama ve ileri bakış.

Aşamalar (her biri runs/mekanik/ altına yazar, tamamlanan atlanır):
    uzun     parametre ızgarası; uzun koşuların istatistikleri gözlemle karşılaştırılır
    gps      1999 sonrası TUBI ve ISTA hız değişimiyle krip ve sünek akma ayarı
    isinma   seçilen ayarla, üç katman kalınlığı senaryosunda uzun ısınma koşuları
    oynatma  1500'den bugüne tarihsel depremlerle sırayla oynatma + her büyük
             depremden önce ileri olasılık (sınama) + bugünden 30 yıllık ileri bakış
    rapor    tablolar ve şekiller

Kullanım: python scripts/mekanik_model.py [aşama ...] [--hizli]
"""
from __future__ import annotations

import argparse
import json
import logging
import pickle
import time
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd

from deprem.mekanik import deney, geometri, hazirlik, tarih
from deprem.mekanik.simulator import Model, Parametreler

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / deney.VERI_DIZINI
log = logging.getLogger("mekanik")

# gözlem hedefleri (uzun dönem): b-değeri bu çalışmanın Marmara kataloğundan (2000-2026,
# fay 15 km içi, M>=2.5: b = 1.00); büyük deprem tekrarlanması ve en uzun kırık tarihsel
# tablodan (1500-2026) hesaplanır.
HEDEF_B = 1.0


def tarihsel_hedefler(v):
    df = tarih.tablo(v.cells)
    big = df[(df.mw >= deney.BUYUK_MW) & (df.t >= 1500)]
    xg = np.arange(60, 400, 10.0)
    n = np.array([((big.x0 <= x) & (big.x1 >= x)).sum() for x in xg])
    tekrar = (2026.9 - 1500) / np.maximum(n, 0.5)
    return {"tekrar_medyan": float(np.median(tekrar)), "uzunluk_maks": float((big.x1 - big.x0).max()),
            "uzunluk_ort": float((big.x1 - big.x0).mean())}


def puan(oz, hedef):
    s = ((oz["b"] - HEDEF_B) / 0.15) ** 2
    s += (np.log(oz["tekrar_medyan"] / hedef["tekrar_medyan"]) / np.log(1.5)) ** 2
    s += ((oz["uzunluk_maks"] - hedef["uzunluk_maks"]) / 60.0) ** 2 if np.isfinite(oz["uzunluk_maks"]) else 25
    return float(s)


def asama_uzun(v, args):
    out = RUN / "uzun.csv"
    if out.exists():
        return pd.read_csv(out)
    hedef = tarihsel_hedefler(v)
    log.info("tarihsel hedefler: %s", hedef)
    yil, ilk = (500, 200) if args.hizli else (1000, 300)
    rows = []
    grid = [(d, e, s) for d in (3e6, 5e6) for e in (0.3, 0.6) for s in (0.3, 0.6)]
    for d, e, s in grid:
        p = Parametreler(dtau0=d, eta_dyn=e, sigma_h=s)
        t0 = time.time()
        m = deney.model_kur(v, p)
        m.ilerle(yil)
        oz = deney.ozet(m.katalog(), ilk, yil)
        oz.update(dtau0=d, eta_dyn=e, sigma_h=s, puan=puan(oz, hedef), sure_s=time.time() - t0)
        log.info("uzun %s", {k: round(val, 3) if isinstance(val, float) else val for k, val in oz.items()})
        rows.append(oz)
        m.katalog().to_csv(RUN / f"uzun_katalog_{d:.0e}_{e}_{s}.csv", index=False)
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    (RUN / "hedefler.json").write_text(json.dumps(hedef, indent=1))
    return df


def en_iyi(df):
    r = df.sort_values("puan").iloc[0]
    return Parametreler(dtau0=float(r.dtau0), eta_dyn=float(r.eta_dyn), sigma_h=float(r.sigma_h))


def _isindir(v, p, yil, uye=None):
    m = deney.model_kur(v, p, uye=uye)
    m.ilerle(yil)
    return m


def asama_gps(v, args, p):
    out = RUN / "gps.csv"
    if out.exists():
        return pd.read_csv(out)
    obs = deney.gps_gozlem(ROOT)
    obs.to_csv(RUN / "gps_gozlem.csv", index=False)
    G = {}
    for s, (la, lo) in deney.GPS_ISTASYON.items():
        G[s] = geometri.yuzey_yerdegistirme(v.cells, np.array([la]), np.array([lo]))[0]
    st = pickle.loads((RUN / "isinma_orta.pkl").read_bytes())[-1]
    ev = tarih.tablo(v.cells)
    e99 = ev[ev.tarih.str.startswith("1999")]
    rows = []
    izgara = [(tM, kc) for tM in (3.0, 10.0, 30.0) for kc in (0.5, 2.0)]
    if args.gps_ayar:
        izgara = [tuple(float(x) for x in args.gps_ayar.split(","))]
    for tM, kc in izgara:
        if (RUN / f"gps_seri_{tM:g}_{kc:g}.csv").exists():
            continue
        if True:
            p2 = replace(p, tM_ref=tM, krip_carpan=kc)
            m, _ = model_senaryo(v, p2, "orta")
            # kilitli kısmın gerilmesini ısınmış durumdan al, krip kuşağı yeni ayarla dengede
            m.S[m.vw] = st["S"][m.vw]
            m.S_fail, m.dtau = st["S_fail"].copy(), st["dtau"].copy()
            m.t = 1999.6 - 80
            m.ilerle(1999.6)
            mws = []
            for e in e99.itertuples():
                m.ilerle(e.t)
                r = deney.zorla(m, e, doldur=args.ozumse)
                mws.append(r["mw"] if r else np.nan)
            mod = deney.gps_model(m, G, 1999.9)
            j = obs.merge(mod, on=["site", "a", "b", "bil"], suffixes=("_g", "_m"))
            rms = float(np.sqrt(np.mean((j.anomali_g - j.anomali_m) ** 2)))
            log.info("gps tM=%g krip=%g rms=%.2f mm/yıl Mw=%s", tM, kc, rms, np.round(mws, 2))
            j.assign(tM_ref=tM, krip_carpan=kc, mw_izmit=mws[0], mw_duzce=mws[1]).to_csv(
                RUN / f"gps_seri_{tM:g}_{kc:g}.csv", index=False)
    if args.gps_ayar:
        return None
    for f in sorted(RUN.glob("gps_seri_*.csv")):
        j = pd.read_csv(f)
        rms = float(np.sqrt(np.mean((j.anomali_g - j.anomali_m) ** 2)))
        rows.append({"tM_ref": j.tM_ref[0], "krip_carpan": j.krip_carpan[0], "rms": rms,
                     "mw_izmit": j.get("mw_izmit", pd.Series([np.nan]))[0],
                     "mw_duzce": j.get("mw_duzce", pd.Series([np.nan]))[0]})
    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    return df


SENARYO = {"ince": 0.0, "orta": 0.5, "kalin": 1.0}  # D90: AFAD (ince) ... EMSC (kalın)


def model_senaryo(v, p, ad):
    """Katman kalınlığı senaryosu: D90 = AFAD + u (EMSC - AFAD)."""

    class _R:
        def uniform(self, *a, **k):
            return SENARYO[ad]

    q, d = hazirlik.isi_akisi_profili(v.d90, v.cells, rng=_R())
    p2 = replace(p, isi_akisi=q, v_pl=hazirlik.levha_hizi(v.cells), seed=p.seed + {"ince": 1, "orta": 2, "kalin": 3}[ad])
    m = Model(v.cells, None, None, p2, ek_krip=hazirlik.sig_krip_maskesi(v.cells), K=v.K(p.mu_eff))
    m.sig_rahatla()
    return m, d


def asama_isinma(v, args, p):
    out = RUN / "isinma.pkl"
    if out.exists():
        return pickle.loads(out.read_bytes())
    yil = 600 if args.hizli else 1500
    anlar = np.linspace(yil * 0.55, yil, 8)
    durumlar = {}
    for ad in SENARYO:
        m, d = model_senaryo(v, p, ad)
        sts = []
        for t in anlar:
            m.ilerle(t)
            sts.append(m.durum())
        durumlar[ad] = sts
        m.katalog().to_csv(RUN / f"isinma_katalog_{ad}.csv", index=False)
        log.info("ısınma %s: %d olay, D90 ort %.1f km", ad, len(m.olaylar), d.mean())
    out.write_bytes(pickle.dumps(durumlar))
    return durumlar


def asama_oynatma(v, args, p, ad, sts):
    """Bir katman kalınlığı senaryosu için topluluk üyeleri: 1500'den bugüne oynatma."""
    out = RUN / (f"oynatma_{ad}" + ("_ozumse" if args.ozumse else "") + ".pkl")
    if out.exists():
        return pickle.loads(out.read_bytes())
    xb = np.arange(0, 456, 1.0)
    sonuc = {"uyeler": []}
    m, _ = model_senaryo(v, p, ad)
    secim = sts[-args.uye:]
    for k, st in enumerate(secim):
        uye = {"ince": 0, "orta": 100, "kalin": 200}[ad] + k
        t0 = time.time()
        m.yukle(st)
        m.olaylar = []
        m.dayanimi_yenile(50_000 + uye)
        m.t = 1500.0
        ev = deney.tarihsel_olaylar(v, uye)
        kayit = {"senaryo": ad, "durum": k, "uye": uye, "zorla": [], "once": {}, "sahte": 0}
        for e in ev.itertuples():
            n0 = len(m.olaylar)
            m.ilerle(e.t)
            # depremden hemen önce: her km'nin doluluğu (gerilme / dayanım)
            kayit["once"][e.tarih] = m.km_doluluk()
            kendi = [o for o in m.olaylar[n0:] if o["mw"] >= deney.BUYUK_MW]
            kayit["sahte"] += len(kendi)
            r = deney.zorla(m, e, doldur=args.ozumse)
            kayit["zorla"].append({"tarih": e.tarih, "t": e.t, "mw_gozlem": e.mw, "mw_sd": e.mw_sd,
                                   "mw_model": r["mw"] if r else np.nan, "ort_kayma": r["ort_kayma"] if r else np.nan,
                                   "sahte_once": len(kendi), "x0": e.x0, "x1": e.x1})
        m.ilerle(2026.9)
        kayit["bugun"] = m.km_doluluk()
        kayit["bugun_hucre"] = {"S": m.S.astype("float32"), "S_fail": m.S_fail.astype("float32")}
        pr, ilk = deney.ileri_olasilik(m, 30.0, args.ileri, seed=uye * 100 + 99, x_bins=xb)
        kayit["ileri"] = {"p": pr, "ilk": ilk, "n": args.ileri}
        kayit["katalog"] = m.katalog()
        sonuc["uyeler"].append(kayit)
        log.info("üye %d (%s/%d): %.0f s, sahte büyük deprem %d, model Mw: %s", uye, ad, k,
                 time.time() - t0, kayit["sahte"],
                 [(z["tarih"][:4], round(z["mw_model"], 2)) for z in kayit["zorla"]])
        out.with_suffix(".parca.pkl").write_bytes(pickle.dumps(sonuc))
    out.write_bytes(pickle.dumps(sonuc))
    return sonuc


def asama_isinma_tek(v, args, p, ad):
    out = RUN / f"isinma_{ad}.pkl"
    if out.exists():
        return pickle.loads(out.read_bytes())
    yil = args.isinma
    anlar = np.linspace(yil * 0.6, yil, args.uye)
    m, d = model_senaryo(v, p, ad)
    sts = []
    for t in anlar:
        m.ilerle(t)
        sts.append(m.durum())
    m.katalog().to_csv(RUN / f"isinma_katalog_{ad}.csv", index=False)
    log.info("ısınma %s: %d olay, D90 ort %.1f km", ad, len(m.olaylar), d.mean())
    out.write_bytes(pickle.dumps(sts))
    return sts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("asama", choices=["gps", "isinma", "oynatma"])
    ap.add_argument("--senaryo", default="orta", choices=list(SENARYO))
    ap.add_argument("--isinma", type=float, default=1000.0)
    ap.add_argument("--uye", type=int, default=3)
    ap.add_argument("--ileri", type=int, default=6)
    ap.add_argument("--hizli", action="store_true")
    ap.add_argument("--ozumse", action="store_true", help="her gerçek depremde kırılan parçayı önce %%100 doluya çek")
    ap.add_argument("--gps-ayar", default=None, help="tek GPS ayarı 'tM,krip' (paralel koşu için)")
    args = ap.parse_args()
    RUN.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                        handlers=[logging.FileHandler(RUN / f"log_{args.asama}_{args.senaryo}.txt"),
                                  logging.StreamHandler()])
    v = deney.Veri(ROOT)
    pj = json.loads((RUN / "secilen_parametreler.json").read_text())
    p = Parametreler(**{k: val for k, val in pj.items() if k in Parametreler.__dataclass_fields__
                        and k not in ("isi_akisi", "v_pl", "extra")})
    log.info("parametreler: %s", p)
    if args.asama == "gps":
        g = asama_gps(v, args, p)
        if g is None:
            return
        r = g.sort_values("rms").iloc[0]
        pj.update(tM_ref=float(r.tM_ref), krip_carpan=float(r.krip_carpan))
        (RUN / "secilen_parametreler.json").write_text(json.dumps(pj, indent=1))
        log.info("seçilen krip/sünek ayarı: tM_ref=%g krip=%g", r.tM_ref, r.krip_carpan)
    else:
        sts = asama_isinma_tek(v, args, p, args.senaryo)
        if args.asama == "oynatma":
            asama_oynatma(v, args, p, args.senaryo, sts)


if __name__ == "__main__":
    main()
