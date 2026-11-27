"""Uzun vadeli sınama (Chartier vd. 2021 tarzı): model serbest koşar, gerçek depremler işlenmez.

Girdiler Chartier vd. 2021'e çekilir (Tablo 1: kayma hızları, 15 km sismojenik derinlik, Adalar
10 km) ve krip varsayımları yarıştırılır:
    kilitli   Tekirdağ-Avcılar tamamen kilitli, 18 mm/yıl
    derin     aynı kesimler 5 km'ye kadar kilitli, altı yavaş kayar
    yari      aynı kesimlerin depremle boşalan payı 8 mm/yıl (yarısı sessiz kayar)
    tam_krip  aynı kesimler tamamen yavaş kayar
    sig_krip  bizim ikinci tur varsayımı: Orta çukurda (27,75-28,05) 5,5 km'ye kadar sığ krip
    eski      ikinci tur girdileri olduğu gibi (karşılaştırma için)

Ölçütler (tarih değil):
    * kazı (paleosismoloji) M>=7,2 tekrar oranları, Chartier Tablo 3 (kesim ortasından geçen kırıklar)
    * 1700-2026 gözlenen M>=6,8 / M>=7,0 yıllık oranları (ikinci tur büyük deprem listesi)
    * en büyük deprem, İstanbul önü (27,75-29,10) için 30 yılda M>=7 olasılığı (zamandan bağımsız)

Kesim boylam sınırları MTA haritasından değil, yaklaşık olarak bizim izimiz üzerinde seçildi.
Her koşu bitince kalıcı klasöre yazılır; yeniden başlatılınca bitenler atlanır.

Kullanım: python scripts/uzun_vade.py kos --is 4 ; python scripts/uzun_vade.py ozet
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from deprem.mekanik import ozumseme as oz, termal  # noqa: E402
from deprem.mekanik.simulator import Model, Parametreler  # noqa: E402

CIKTI = oz.KALICI / "uzun_vade"  # MTA geometrisi: KALICI yeni klasör
ISINMA, SURE = 1000.0, 4000.0

# (ad, boylam0, boylam1, kayma hızı mm/yıl ort.) -- hızlar Chartier vd. 2021 Tablo 1 ortalamaları;
# sınırlar MTA 2026 kesimlerinin ana kol üzerindeki boylamları (data/mekanik_mta_kaba/kesimler.csv):
# Ganos = Saros + Ganos; İzmit = Darıca + Karamürsel + Gölcük + Tepetarla; Adapazarı = Arifiye.
# Chartier'nin "Orta çukur" kesimi MTA'da Tekirdağ ve Kumburgaz fayları içinde kalır. Chartier'nin
# "Çınarcık" kesimi (eğik normal, 3,6 mm/yıl) MTA'da ana kol üzerinde değil, havza kenarındaki
# Çınarcık Fayı'dır; modelde güney kol (GEM Armutlu-Çınarcık) olarak duruyor.
KESIM = [("Ganos", 25.0, 27.525, 18.0), ("Tekirdağ", 27.525, 27.935, 18.0), ("Kumburgaz", 27.935, 28.525, 18.0),
         ("Avcılar", 28.525, 28.875, 18.0), ("Adalar", 28.875, 29.345, 14.3), ("İzmit", 29.345, 30.315, 20.0),
         ("Adapazarı", 30.315, 30.665, 18.0), ("Karadere", 30.665, 30.995, 15.0), ("Düzce", 30.995, 32.0, 15.0)]
KRIP_KESIM = ("Tekirdağ", "Kumburgaz", "Avcılar")
# Chartier Tablo 3: kazı yeri kesimi, yıllık M>=7,2 oranı, belirsizlik
KAZI = [("Ganos", 0.00357, 0.00100), ("Ganos", 0.00300, 0.00100), ("İzmit", 0.00250, 0.00130),
        ("İzmit", 0.00320, 0.00140), ("İzmit", 0.00320, 0.00140), ("Adapazarı", 0.00180, 0.00100),
        ("Karadere", 0.00180, 0.00100), ("Düzce", 0.00250, 0.00130)]
VARSAYIM = ["kilitli", "derin", "yari", "tam_krip", "sig_krip", "eski"]
ISTANBUL = (27.75, 29.10)


def kesim_maske(c, ad):
    _, l0, l1, _ = next(k for k in KESIM if k[0] == ad)
    lon = c.lon.to_numpy()
    return (c.kol == 0).to_numpy() & (lon >= l0) & (lon < l1)


def kur(ag, th, seed, varsayim):
    if varsayim == "eski":
        return oz.model_kur(ag, th, seed)
    c = ag.cells
    v = oz.levha_hizi(ag, th["v_pl_marmara"])  # yan kollar GEM (Armutlu, İmralı ~0,8: Chartier <1)
    for ad, _, _, hiz in KESIM:
        v[kesim_maske(c, ad)] = hiz / 1000.0
    D = np.full(len(c), 15.0)
    D[kesim_maske(c, "Adalar")] = 10.0
    z = c.z_km.to_numpy()
    krip = np.zeros(len(c), bool)
    kk = np.zeros(len(c), bool)
    for ad in KRIP_KESIM:
        kk |= kesim_maske(c, ad)
    if varsayim == "derin":
        krip = kk & (z > 5.0)
    elif varsayim == "yari":
        v[kk] = 8.0 / 1000.0
    elif varsayim == "tam_krip":
        krip = kk.copy()
    elif varsayim == "sig_krip":
        lon = c.lon.to_numpy()
        krip = (c.kol == 0).to_numpy() & (lon >= 27.75) & (lon <= 28.05) & (z < 5.5)
    p = Parametreler(mu_eff=th["mu_eff"], dtau0=th["dtau0"], sigma_h=th["sigma_h"], eta_dyn=th["eta_dyn"],
                     asiri=-th["kalan"], kayma_ref=th["kayma_ref"], sigma_grad=th["sigma_grad"],
                     tM_ref=th["tM_ref"], krip_carpan=th["krip_carpan"],
                     isi_akisi=termal.isi_akisi_derinlikten(D), v_pl=v, seed=seed, mu_kayma=oz.MU_YENI)
    m = Model(c, None, None, p, ek_krip=krip, K=ag.K(th["mu_eff"]))
    m.tM_adim_alt = 0.05
    m.sig_rahatla()
    return m


_AG = {}


def tek(arg):
    varsayim, ayar_ad, th, seed = arg
    f = CIKTI / f"{varsayim}_{ayar_ad}_{seed}.pkl"
    if f.exists():
        return
    if "ag" not in _AG:
        _AG["ag"] = oz.Ag(ROOT, oz.BUYUK_AG)
    ag = _AG["ag"]
    t0 = time.time()
    m = kur(ag, th, seed, varsayim)
    m.t = 0.0
    m.ilerle(ISINMA)
    m.olaylar = []
    m.kayit_slip = (6.8, [])
    while m.t < ISINMA + SURE - 1e-9:
        m.ilerle(min(m.t + 50.0, ISINMA + SURE))
    kol = ag.cells.kol.to_numpy()
    lon = ag.cells.lon.to_numpy()
    ana = (kol == 0) & m.vw
    buyuk = []
    for t, mw, d in m.kayit_slip[1]:
        h = (d > 0.05 * d.max()) & ana
        buyuk.append({"t": t - ISINMA, "mw": mw, "lon0": float(lon[h].min()) if h.any() else np.nan,
                      "lon1": float(lon[h].max()) if h.any() else np.nan})
    kat = pd.DataFrame([{"t": e["t"] - ISINMA, "mw": e["mw"]} for e in m.olaylar if not e.get("zorla")])
    out = {"varsayim": varsayim, "ayar": ayar_ad, "seed": seed, "th": th, "sure": SURE, "buyuk": pd.DataFrame(buyuk),
           "katalog": kat, "n_vw": int(m.vw.sum()), "sn": round(time.time() - t0, 1)}
    tmp = f.with_suffix(".tmp")
    tmp.write_bytes(pickle.dumps(out))
    tmp.replace(f)
    print(f.name, out["sn"], "s", len(kat), "olay, M>=6,8:", len(buyuk), flush=True)


def ayarlar():
    out = {"literatur": dict(oz.baslangic_theta(), zayif=0.0)}
    rs = [json.loads(s) for s in (oz.KALICI_GEM / "ara_buyuk_kaba.jsonl").read_text().splitlines()]
    out["ikinci_tur"] = dict(max(rs, key=lambda r: r["J"])["th"], zayif=0.0)
    return out


def kos(args):
    CIKTI.mkdir(parents=True, exist_ok=True)
    isler = [(v, a, th, s) for v in VARSAYIM for a, th in ayarlar().items() for s in range(2)]
    from multiprocessing import get_context
    with get_context("fork").Pool(args.n_is) as p:
        list(p.imap_unordered(tek, isler))


def gozlenen():
    """1700-2026 kuzey kol + Marmara gözlenen M>=6,8 büyük depremler (ikinci tur listesi)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    import mekanik_ikinci_tur as M
    ag, ol, _ = M._hazirla("buyuk", oz.BUYUK_AG, True)
    b = [o for o in ol if o.tip == "buyuk" and o.mw >= 6.8 and o.t >= 1700]
    return pd.DataFrame([{"t": o.t, "mw": o.mw, "ad": o.ad} for o in b])


def ozet(args):
    rows = []
    g = gozlenen()
    yil = 2026.9 - 1700.0
    for f in sorted(CIKTI.glob("*.pkl")):
        r = pickle.loads(f.read_bytes())
        b, T = r["buyuk"], r["sure"]
        s = {"varsayim": r["varsayim"], "ayar": r["ayar"], "seed": r["seed"]}
        s["M6.8_yil"] = len(b) / T
        s["M7.0_yil"] = int((b.mw >= 7.0).sum()) / T if len(b) else 0.0
        s["M7.2_yil"] = int((b.mw >= 7.2).sum()) / T if len(b) else 0.0
        s["mw_maks"] = float(b.mw.max()) if len(b) else np.nan
        # denetim (Chartier: M7,7 üstü neredeyse hiç; 1912 ~ 120 km, 1999 İzmit ~ 145 km kırık)
        s["oran_M7.7_ustu"] = float((b.mw > 7.7).mean()) if len(b) else 0.0
        bb = b[(b.mw >= 7.2) & (b.mw <= 7.7)] if len(b) else b
        s["boy_km_M7.2_7.7_medyan"] = float(((bb.lon1 - bb.lon0) * 84.4).median()) if len(bb) else np.nan
        # kazı uyumu: her kazı için model oranı (kesim ortasından geçen M>=7,2), normalize hata kareleri
        ki = 0.0
        for ad, oran, sd in KAZI:
            _, l0, l1, _ = next(k for k in KESIM if k[0] == ad)
            mid = 0.5 * (l0 + min(l1, 31.4))
            n = int(((b.mw >= 7.2) & (b.lon0 <= mid) & (b.lon1 >= mid)).sum()) if len(b) else 0
            s[f"kazi_{ad}"] = n / T
            ki += ((n / T - oran) / sd) ** 2
        s["kazi_ki2"] = ki
        # İstanbul önü: 30 yıllık pencerelerde M>=7 kırık olasılığı (zamandan bağımsız, kayan pencere)
        if len(b):
            ist = b[(b.mw >= 7.0) & (b.lon0 <= ISTANBUL[1]) & (b.lon1 >= ISTANBUL[0])].t.to_numpy()
        else:
            ist = np.array([])
        w = np.arange(0, T - 30, 1.0)
        s["istanbul_P30_M7"] = float(np.mean([((ist >= a) & (ist < a + 30)).any() for a in w]))
        s["istanbul_M7_tekrar_yil"] = T / len(ist) if len(ist) else np.inf
        k = r["katalog"]
        s["M4_yil"] = float((k.mw >= 4.0).sum() / T) if len(k) else 0.0
        s["M5_yil"] = float((k.mw >= 5.0).sum() / T) if len(k) else 0.0
        rows.append(s)
    d = pd.DataFrame(rows)
    gz = {"M6.8_yil": len(g) / yil, "M7.0_yil": float((g.mw >= 7.0).sum() / yil), "mw_maks": float(g.mw.max()),
          "kazi": {ad: oran for ad, oran, _ in KAZI}}
    d.to_csv(ROOT / "reports" / "tables" / "uzun_vade.csv", index=False, float_format="%.5f")
    (ROOT / "reports" / "tables" / "uzun_vade_gozlenen.json").write_text(json.dumps(gz, ensure_ascii=False, indent=1))
    pd.set_option("display.width", 250)
    print("gözlenen:", gz)
    print(d.groupby(["varsayim", "ayar"]).mean(numeric_only=True).round(4).drop(columns="seed").to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("asama", choices=["kos", "ozet"])
    ap.add_argument("--is", dest="n_is", type=int, default=4)
    a = ap.parse_args()
    {"kos": kos, "ozet": ozet}[a.asama](a)


if __name__ == "__main__":
    main()
