"""Etki payı raporu: runs/etki/*.csv → reports/etki_payi.md ve reports/tables/etki_payi.csv

Tanım
-----
S(m): m sürümünün 1900-2025 olaylarındaki toplam ileriye dönük log olabilirliği (her olay
yalnız öncesiyle eğitilmiş modelle puanlanmış). Türkiye: bütün Türkiye kutusu. Marmara:
Marmara kutusundaki olaylarda log λ, her aralıkta Marmara'nın kendi beklenen sayısı çıkarılır.

* Toplam bilgi  T = S(ana model) − S(sıfır noktası: sabit hız, alanla orantılı)
* Parça katkısı Δ = S(ana model) − S(ana model − parça)      (çıkarılan parçalar)
               Δ = S(ana model + parça) − S(ana model)      (eklenen parçalar)
* Etki payı     = Δ / T  (yüzde). Çıkar-bir-tane katkıları birbiriyle örtüşen bilgiyi saymaz,
  bu yüzden toplamları %100 etmez; kalan "ortak/ayrıştırılamayan" satırındadır.

Belirsizlik: yıl blokları bootstrap'ı (aynı yıllar bütün sürümlerde birlikte örneklenir).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "etki"
B = 4000
Y0 = 1900

# (Ulas'ın listesindeki ad, karşılaştırma: (iyi, kötü) sürüm, açıklama)
ROWS = [
    ("Arka plan haritası (nerede deprem olur)", ("partial", "no_bmap"), "geçmiş sismisite + fay kayma hızı haritası birlikte"),
    ("  · geçmiş sismisite haritası", ("partial", "no_smap"), "yalnız fay haritasıyla kalınca kaybolan"),
    ("Kayma hızı (harita ve doluluktan)", ("partial", "no_slip"), "fay haritası ve doluluktaki kayma hızı birlikte"),
    ("  · fay kayma hızı haritası", ("partial", "no_fmap"), "yalnız geçmiş sismisite haritasıyla kalınca kaybolan"),
    ("Artçı/tetikleme (ETAS benzeri, hızlı)", ("partial", "no_fast"), "Omori, günler-aylar"),
    ("Yavaş etki (yıllar)", ("partial", "no_slow"), "τ ≈ birkaç yıl süren ek etki"),
    ("Tetiklemenin tümü (hızlı + yavaş)", ("partial", "no_trig"), ""),
    ("Uzaklıkla azalan faylar arası etki (100 km ötesi)", ("partial", "no_far"), "100 km ötesi kesilince kaybolan"),
    ("Büyük sistemler arası ayrı çarpan", ("partial", "no_cross"), ""),
    ("Coulomb gerilme aktarımı", ("partial", "no_cfs"), ""),
    ("Son kırılmadan geçen süre", ("load_time", "no_load"), "doluluk olmadan, yalnız süre"),
    ("Doluluk / kapasite (sürenin üstüne)", ("partial", "load_time"), "kayma hızı × süre / kapasite − yalnız süre"),
    ("Doluluk terimi tümü", ("partial", "no_load"), "süre + doluluk birlikte"),
    ("Fay fay ayar (ortak ayara karşı)", ("partial", "pooled"), ""),
    ("Sürekli yeniden eğitim (1900'de dondurmaya karşı)", ("partial", "fixed"), ""),
    ("Zincirleme kırılma (yüklü komşu fay zinciri)", ("partial_chain", "partial"), "eklenince"),
    ("b-değeri (yerel, son 5 yıl)", ("partial_bval", "partial"), "eklenince; 1998 öncesi boş"),
    ("Moment açığı (küçük depremler dahil)", ("partial_deficit", "partial"), "eklenince"),
    ("dv/v (Doğu Marmara)", ("partial_dvv", "partial"), "eklenince; 2009 sonrası, yalnız Doğu Marmara"),
]


def load():
    res = {}
    for p in sorted(RUN.glob("*.csv")):
        d = pd.read_csv(p)
        d["time"] = pd.to_datetime(d["time"], utc=True, format="ISO8601")
        res[p.stem] = d.sort_values("i").reset_index(drop=True)
    return res


def per_year(res, region: bool, y0=Y0, y1=2025):
    """Sürüm × yıl matrisi: yıl içindeki olay puanları toplamı (bölge veya Türkiye)."""
    base = res["partial"]
    yr = base["time"].dt.year.to_numpy()
    keep = (yr >= y0) & (yr <= y1)
    years = np.unique(yr[keep])
    out, n = {}, None
    for k, d in res.items():
        assert (d["i"].to_numpy() == base["i"].to_numpy()).all(), k
        if region:
            s = d["loglam"].to_numpy() * d["in_region"].to_numpy() - d["comp_region"].to_numpy()
        else:
            s = d["score"].to_numpy()
        out[k] = np.array([s[keep & (yr == y)].sum() for y in years])
    ev = base["in_region"].to_numpy() if region else np.ones(len(base), bool)
    n = np.array([(ev & keep & (yr == y)).sum() for y in years])
    return out, n, years


def table(res, region: bool, y0=Y0, y1=2025):
    S, n, years = per_year(res, region, y0, y1)
    rng = np.random.default_rng(1)
    idx = rng.integers(0, len(years), size=(B, len(years)))
    tot = lambda k: S[k].sum()
    boot = lambda k: S[k][idx].sum(1)
    T, Tb = tot("partial") - tot("uniform"), boot("partial") - boot("uniform")
    N = n.sum()
    rows = []
    for name, (a, b), note in ROWS:
        if a not in S or b not in S:
            rows.append(dict(değişken=name, not_=note, durum="koşulmadı"))
            continue
        d, db = tot(a) - tot(b), boot(a) - boot(b)
        dlo, dhi = np.percentile(db, [2.5, 97.5])
        # yüzde aralığı: Δ'nın aralığı ÷ toplam bilginin nokta değeri (T'nin kendi belirsizliği ayrıca verilir)
        pct, lo, hi = 100 * d / T, 100 * dlo / T, 100 * dhi / T
        p = float((db <= 0).mean())
        if dlo > 0:
            durum = "gerçek katkı"
        elif dhi < 0:
            durum = "zararlı (modeli kötüleştiriyor)"
        else:
            durum = "sıfırdan ayırt edilemez"
        rows.append(dict(değişken=name, pay=pct, pay_lo=lo, pay_hi=hi, olay_başına=d / N, d_lo=dlo / N,
                         d_hi=dhi / N, p=p, durum=durum, not_=note))
    extra = dict(T=T / N, T_lo=np.percentile(Tb, 2.5) / N, T_hi=np.percentile(Tb, 97.5) / N, N=int(N),
                 etas=(tot("partial") - tot("etas")) / N if "etas" in S else np.nan,
                 etas_share=100 * (tot("etas") - tot("uniform")) / T if "etas" in S else np.nan)
    return pd.DataFrame(rows), extra


def fmt(df):
    out = []
    for _, r in df.iterrows():
        if r.get("durum") == "koşulmadı" or pd.isna(r.get("pay", np.nan)):
            out.append(dict(değişken=r["değişken"], **{"etki payı": "koşulmadı"}))
            continue
        out.append({"değişken": r["değişken"], "etki payı": f"%{r['pay']:.1f}",
                    "%95 aralık": f"[%{r['pay_lo']:.1f}, %{r['pay_hi']:.1f}]",
                    "olay başına (log)": f"{r['olay_başına']:+.3f}", "p": f"{r['p']:.3f}", "karar": r["durum"]})
    return pd.DataFrame(out)


COMPARE = [("pruned", "Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)"),
           ("etas", "ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model)")]


def compare(res, region, y0, y1):
    """Elenmiş modeller − tam (ana) model, olay başına, yıl bootstrap GA."""
    S, n, years = per_year(res, region, y0, y1)
    rng = np.random.default_rng(2)
    idx = rng.integers(0, len(years), size=(B, len(years)))
    N = n.sum()
    rows = []
    for k, lab in COMPARE:
        if k not in S:
            continue
        d = S[k] - S["partial"]
        db = d[idx].sum(1)
        lo, hi = np.percentile(db, [2.5, 97.5])
        rows.append({"model": lab, "olay": int(N), "tam modele göre (olay başına)": f"{d.sum() / N:+.3f}",
                     "%95 GA": f"[{lo / N:+.3f}, {hi / N:+.3f}]", "p (elenmiş daha iyi değil)": f"{(db <= 0).mean():.3f}"})
    return pd.DataFrame(rows)


def main():
    res = load()
    lines = ["# Model değişkenlerinin etki payı (yüzde)", "",
             "> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Yüzdeler, modelin geçmiş depremleri "
             "ileriye dönük (cevabı görmeden) tahmin ederken topladığı bilginin hangi parçadan geldiğini ölçer.", "",
             "Yöntem: `scripts/etki_payi.py` (fay fay sıralı model, PR 14), her sürüm baştan sıralı yürütüldü. "
             "Etki payı = o parça çıkarılınca (veya eklenince) kaybedilen (kazanılan) ileriye dönük log olabilirlik "
             "÷ ana modelin, hiçbir bilgi taşımayan sıfır noktasına (sabit hız, alanla orantılı) göre toplam "
             "bilgisi. %95 aralıklar yıl blokları bootstrap'ı. Bir parçanın tek başına payı, başka parçalarla "
             "örtüşen bilgiyi saymaz; bu yüzden toplam %100 etmez.", ""]
    all_rows = []
    for region, title in ((False, "Türkiye (1900-2025, M≥5,5)"), (True, "Marmara kutusu (1900-2025, M≥5,5)")):
        df, ex = table(res, region)
        df.insert(0, "bölge", "Marmara" if region else "Türkiye")
        all_rows.append(df)
        lines += [f"## {title}", "",
                  f"Olay sayısı: {ex['N']}. Ana modelin toplam bilgisi (sıfır noktasına göre): olay başına "
                  f"{ex['T']:.2f} [{ex['T_lo']:.2f}, {ex['T_hi']:.2f}] log birimi."
                  + (f" Ana model − ETAS: {ex['etas']:+.3f} olay başına; ETAS tek başına toplam bilginin "
                     f"%{ex['etas_share']:.0f}'ini taşıyor." if np.isfinite(ex["etas"]) else ""), "",
                  fmt(df).to_markdown(index=False), ""]
    tab = pd.concat(all_rows, ignore_index=True)
    (ROOT / "reports" / "tables").mkdir(parents=True, exist_ok=True)
    tab.to_csv(ROOT / "reports" / "tables" / "etki_payi.csv", index=False)
    lines += ["## Seçim dönemi (1900-1999) ve dokunulmamış dönem (2000-2025), Türkiye", "",
              "Eleme kararı yalnız 1900-1999 sonuçlarıyla verildi; 2000-2025 bu kararda kullanılmadı.", ""]
    for y0, y1 in ((1900, 1999), (2000, 2025)):
        df, ex = table(res, False, y0, y1)
        lines += [f"### {y0}-{y1} ({ex['N']} deprem)", "", fmt(df).to_markdown(index=False), ""]
    lines += ["## Elenmiş model tam modele karşı", "",
              "Elenen parçalar (1900-1999'da katkısı sıfırdan ayırt edilemeyen ya da zararlı, ve başka bir parçayla "
              "yer değiştirmeyen): Coulomb gerilme yönü, büyük sistemler arası ayrı çarpan ve aynı sistem çarpanı, "
              "doluluk/son kırılmadan geçen süre, fay fay ayar (ortak ayar kullanıldı). Pozitif = elenmiş model daha iyi. "
              "Seçim 1900-1999 verisiyle yapıldığı için o dönemdeki üstünlük iyimserdir; asıl ölçü 2000-2025'tir.", ""]
    for region, nm in ((False, "Türkiye"), (True, "Marmara")):
        for y0, y1 in ((1900, 1999), (2000, 2025), (1900, 2025)):
            c = compare(res, region, y0, y1)
            if len(c):
                lines += [f"**{nm} {y0}-{y1}**", "", c.to_markdown(index=False), ""]
    lines += ["## Bu çalışmada ölçülmeyenler", "",
              "- **GPS kilitlenme / krip:** sıralı modelde yok. Önceki çalışmada (reports/parametre_aramasi.md) "
              "açık GPS istasyonları seyrek olduğu için krip çözülemedi; etki payı ölçülemez.",
              "- **Gaz / sıvı:** varsayılan olarak kapalı. Ortak modelde eklenince ileriye dönük puan ETAS'a göre "
              "kötüleşti (birleşik eğitim özeti); payı sıfır veya eksi.", "",
              "## Sınırlar", "",
              "- Marmara'da 1900-2025'te az sayıda M≥5,5 deprem var; aralıklar bu yüzden çok geniş.",
              "- Çıkar-bir-tane payı, parçanın 'yerine başka bir parça geçemeyen' katkısıdır. Örneğin fay haritası "
              "ile geçmiş sismisite haritası aynı bilgiyi kısmen taşıdığından her birinin tek payı küçük, ikisinin "
              "birlikte payı büyük görünür.",
              "- b-değeri yalnız 1998 sonrası, dv/v yalnız 2009 sonrası ve yalnız Doğu Marmara'da dolu; "
              "bu dönemlerde çok az deprem var, testin gücü düşük."]
    (ROOT / "reports" / "etki_payi.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
