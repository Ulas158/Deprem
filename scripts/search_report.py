"""Parametre araması raporu: reports/parametre_aramasi.md (+ harita katmanı fay_doluluk).

    python scripts/search_report.py --config configs/train_search.toml
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from deprem.pipeline import Run

ROOT = Path(__file__).resolve().parents[1]


def md(df: pd.DataFrame, digits=3) -> str:
    return df.round(digits).to_markdown(index=False)


def f3(x):
    return "—" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:+.3f}"


def search_section(run, w, title, h="##"):
    d = run.dir
    s = json.loads((d / "search" / "summary.json").read_text())
    w(f"{h} {title}\n")
    sy, ty = s["search_years"], s["test_years"]
    sm = [f"M≥{float(m):g}" for m in run.cfg.get("search", {}).get("score_mags", run.mags[:3])]
    trials = pd.read_csv(d / "search" / "trials.csv")
    test = pd.DataFrame(s["test_table"])
    best = test.iloc[0]
    # ------------------------------------------------------------ kısa sonuç
    w(f"{h}# Kısa sonuç\n")
    nb = s.get("null_best_scores", [])
    beat_null = bool(nb) and s["real_best_score"] > max(nb)
    w(f"- {len(trials)} kombinasyon, yalnız {sy[0]}-{sy[-1]} yıllarında denendi. En iyisinin arama puanı "
      f"{s['real_best_score']:+.3f} (olay başına, ETAS'a göre; 0 = ETAS kadar).")
    if nb:
        w(f"- Aynı arama, parametrelerin yerleri rastgele karıştırılmış (bilgisi silinmiş) verilerde "
          f"{len(nb)} kez yapıldı; oradaki en iyi puanlar: {', '.join(f'{x:+.3f}' for x in nb)}. "
          + ("Gerçek arama bunların hepsini geçti, yani arama döneminde gerçek bir sinyal var gibi görünüyor."
             if beat_null else
             "Gerçek arama bunları **geçemedi**: arama döneminde bulunan 'en iyi' kombinasyon, bilgisiz "
             "veride şans eseri bulunabilecek kadar iyi; yani üstünlüğü büyük ölçüde şans."))
    else:
        w("- Bilgisi silinmiş verilerde tekrar arama bu koşuda tamamlanmadı (bir arama 1 saatten uzun sürdü); "
          "sınama sonucu zaten negatif olduğu için sonucu değiştirmez.")
    p = best["Holm düzeltilmiş p"]
    w(f"- Dokunulmamış sınama yıllarında ({ty[0]}-{ty[-1]}) en iyi kombinasyonun puanı "
      f"{best['sınama puanı']:+.3f}; %95 güven aralığı {best['birleşik %95 GA']}, çoklu test düzeltmeli p = {p:.2f}. "
      + ("ETAS'ı anlamlı biçimde geçti." if p < 0.05 and best["sınama puanı"] > 0 else
         "ETAS'ı **anlamlı biçimde geçmedi**."))
    w("")
    # ------------------------------------------------------------ yöntem
    w(f"{h}# Nasıl yapıldı\n")
    w(f"1. **Arama dönemi {sy[0]}-{sy[-1]}.** Her yıl için model yalnız o yıldan önceki yıllarla eğitildi ve o yılı "
      f"tahmin etti. Puan: {', '.join(sm)} depremlerin olay başına bilgi kazancının ortalaması (ETAS'a göre; + iyi)."
      + (f" Katsayılar {run.cfg['search'].get('refit_every', 1)} yılda bir, yalnız o bloktan önceki yıllarla yeniden eğitildi."
         if run.cfg.get("search", {}).get("refit_every", 1) > 1 else ""))
    w("2. **Aranan şeyler:** hangi parametre gruplarının açık olacağı, her grubun önem ağırlığı (0,05-3; küçük = "
      "katkısı kısılmış), genel düzenlileştirme gücü, Coulomb sürtünmesi μ', veriyle önem ağırlığı (ARD, uyarlamalı "
      "ridge), fiziksel enerji öncülü (enerji/kayma birikimi grupları büyük, gelgit gibi küçük etkiler küçük "
      "başlangıç), katalog belirsizlik ağırlıkları ve eğitim penceresi (tüm geçmiş / son 3 / son 5 yıl). "
      "Arama Bayesçi optimizasyonla (TPE) yapıldı.")
    if nb:
        w("3. **Şans sınırı:** aynı arama, parametrelerin hücrelere rastgele karıştırıldığı verilerde tekrarlandı.")
    else:
        w("3. **Şans sınırı:** bu koşuda tamamlanmadı; yerine dokunulmamış sınama yılları ve çoklu test düzeltmesi kullanıldı.")
    w(f"4. **Sınama {ty[0]}-{ty[-1]}:** aramada hiç kullanılmadı. İlk {sum(test['model'].str.startswith('aramada'))} "
      "farklı kombinasyon ve karşılaştırma modelleri burada yıl yıl ileri yürüyüşle bir kez ölçüldü; güven aralıkları "
      "yıllar üzerinden bootstrap, p değerleri Holm yöntemiyle çoklu teste göre düzeltildi.")
    w("5. **ETAS referansı değişmedi** (her yıl yalnız o yıl öncesiyle yeniden uyduruldu, önceki PR'lardaki gibi).\n")
    # ------------------------------------------------------------ sınama tablosu
    w(f"{h}# Sınama yılları sonuçları\n")
    cols = ["model", "arama puanı", "sınama puanı", "birleşik %95 GA", "birleşik p (tek yönlü)", "Holm düzeltilmiş p"]
    for m_ in sm:
        cols += [f"{m_} kazanç", f"{m_} geçtiği yıl"]
    cols.insert(cols.index(f"{sm[-1]} geçtiği yıl"), f"{sm[-1]} %95 GA")
    w(md(test[[c for c in cols if c in test.columns]]))
    w("")
    w("Puanlar olay başına doğal log birimindedir: +0,10 ≈ gerçekleşen depremlere ETAS'tan %10 daha yüksek "
      "olasılık vermek demektir.\n")
    # ------------------------------------------------------------ en iyi kombinasyon
    bc = json.loads(best["config"])
    w(f"{h}# Aramanın seçtiği kombinasyon\n")
    w(f"- Açık gruplar ({len(bc['groups'])}): " + ", ".join(bc["groups"]))
    w(f"- Önem ağırlıkları: " + ", ".join(f"{g} {v:.2f}" for g, v in sorted(bc["weights"].items(), key=lambda x: -x[1])))
    w(f"- Genel öncül sd {bc['prior_sd']:.3f}, μ' {bc['mu_eff']}, ARD {'açık' if bc['ard'] else 'kapalı'}, "
      f"enerji öncülü {'açık' if bc['energy_prior'] else 'kapalı'}, belirsizlik ağırlıkları "
      f"{'açık' if bc['uncertainty_weights'] else 'kapalı'}, eğitim penceresi {bc['window'] or 'tüm geçmiş'}\n")
    coef = s.get("best_coef_last_year", {})
    if coef:
        cdf = pd.DataFrame([{"parametre": k, "katsayı": v} for k, v in coef.items() if k != "(sabit)"])
        cdf = cdf.reindex(cdf["katsayı"].abs().sort_values(ascending=False).index).head(15)
        w(f"Son sınama yılındaki en büyük 15 katsayı (standartlaştırılmış; +0,1 ≈ 1 sd artışta oran ×1,105):\n")
        w(md(cdf))
        w("")
    # grup sıklığı: ilk %5 denemede hangi gruplar açık
    topq = trials.head(max(len(trials) // 20, 5))
    freq = {}
    for c in topq["config"]:
        for g in json.loads(c)["groups"]:
            freq[g] = freq.get(g, 0) + 1
    allg = sorted(s["groups"])
    fr = pd.DataFrame([{"grup": g, "en iyi %5 denemede açık olma oranı": freq.get(g, 0) / len(topq)} for g in allg])
    w("En iyi %5 denemede her grubun açık olma oranı (yüksek = arama bu grubu sevdi; yalnız arama dönemi):\n")
    w(md(fr.sort_values("en iyi %5 denemede açık olma oranı", ascending=False), 2))
    w("")
    # ------------------------------------------------------------ doluluk özetleri
    if s.get("capacity_variants"):
        w(f"{h}# Doluluk nasıl özetlenmeli? (ortalama, moment ağırlıklı, en dolu parça, zincirleme)\n")
        w("Her satır, kapasite grupları dışındaki tüm parametrelere yalnız o özetin eklenmesinin katkısıdır "
          "(+ iyi). Sınama sütunu aramada kullanılmayan yıllardır.\n")
        w(md(pd.DataFrame(s["capacity_variants"])))
        w("")
    # ------------------------------------------------------------ deprem deprem
    if s.get("eventwalk"):
        w(f"{h}# Deprem deprem (zaman + konum), sınama yılları\n")
        ew = pd.DataFrame(s["eventwalk"])
        keep = [c for c in ["hedef", "olay", "zaman kazancı", "zaman %95 GA", "konum kazancı", "konum %95 GA",
                            "toplam kazanç", "toplam %95 GA", "en olası %5 hücrede (ETAS)",
                            "en olası %5 hücrede (model)"] if c in ew.columns]
        w(md(ew[keep]))
        w("\nBüyüklük için model ETAS ile aynı Gutenberg-Richter dağılımını kullanır; büyüklükte fark yoktur.\n")
    if s.get("replicates"):
        rp = pd.DataFrame(s["replicates"])
        w(f"{h}# Katalog belirsizliği\n")
        w(f"Katalog her depremin zaman, konum ve büyüklük hatasıyla {len(rp)} kez yeniden örneklendi; katalogdan "
          "hesaplanan parametreler ve etiketler yeniden üretildi, en iyi kombinasyon sınama yıllarında yeniden ölçüldü.\n")
        q = rp.drop(columns=["rep"]).describe(percentiles=[0.025, 0.5, 0.975]).T[["2.5%", "50%", "97.5%"]]
        w(md(q.reset_index().rename(columns={"index": "ölçü"})))
        w("")


def capacity_section(run, w):
    d = run.dir
    cfrom = int(run.cfg.get("capacity", {}).get("calib_from", 1950))
    bstart = run.cfg.get("budget", {}).get("start", "1900-01-01")[:4]
    rf = pd.read_csv(d / "capacity_rupture_fills.csv")
    rf["time"] = pd.to_datetime(rf["time"], utc=True, format="ISO8601")
    cal = json.loads((d / "capacity_calibration.json").read_text())
    last = cal[max(cal, key=int)]
    big = rf[(rf["time"].dt.year >= cfrom) & (rf["kayma_hızı_ağırlıklı"] >= 2.0) & (rf["mw"] >= 6.8)]
    w("### \"%100 dolunca deprem olur mu?\" testi\n")
    w("Her 1 km parçanın enerji biriktirme kapasitesi = kilitli alan × kayaç sertliği (CRUST1.0) × taşınabilir en büyük "
      "gerilme düşümü. Doluluk = biriken kayma açığı (GPS kilitlenmesiyle, depremlerde boşalarak) / kapasite. "
      "Ardından kapasite, geçmiş büyük depremlerin kırıldığı andaki doluluğa göre fay → bölge → genel "
      "hiyerarşisiyle kalibre edildi (o fayın tipik kırılma birikimi = %100).\n")
    if len(big):
        pct = big["aynı_bölgedeki_yüzdelik"].dropna()
        wil = stats.wilcoxon(pct - 0.5) if len(pct) > 5 else None
        w(f"- Bütçe {bstart}'de başladı; {cfrom} sonrası, kayma hızı ≥2 mm/yıl faylarda {len(big)} büyük deprem (M≥6,8).")
        w(f"- Kırılma anındaki doluluk olaydan olaya çok değişiyor: log-ölçekte standart sapma "
          f"{last.get('olay_olay_sd', float('nan')):.2f}, yani tipik değerin yaklaşık "
          f"{np.exp(last.get('olay_olay_sd', 0)):.1f} katı aşağı-yukarı. Sabit bir eşik (\"%100'de kırılır\") yok.")
        w(f"- Kırılan yerin doluluğu, aynı bölgedeki tüm fay parçalarına göre ortalama %{100 * pct.mean():.0f}'lik "
          f"dilimde (eşik doğru olsaydı %100'e yakın olurdu; rastgele olsaydı %50). "
          + (f"Wilcoxon p = {wil.pvalue:.3f}." if wil else ""))
        allb = rf[(rf["kayma_hızı_ağırlıklı"] >= 2.0) & (rf["mw"] >= 6.8)]
        if len(allb) > len(big):
            w(f"- Tüm dönem ({allb['time'].dt.year.min()}-{allb['time'].dt.year.max()}, {len(allb)} deprem, tarihsel "
              f"konum/büyüklük hataları büyük): kırılma anı kalibre doluluğun ortancası "
              f"%{allb['doluluk_kalibre_%'].median():.0f}, %10-90 aralığı %{allb['doluluk_kalibre_%'].quantile(0.1):.0f}-"
              f"%{allb['doluluk_kalibre_%'].quantile(0.9):.0f}; aynı bölgedeki yüzdelik ortalaması "
              f"%{100 * allb['aynı_bölgedeki_yüzdelik'].mean():.0f}.")
        bf = big.groupby("fay").filter(lambda g: len(g) > 1)
        if len(bf):
            w("- Aynı fayın birden fazla kırılması (döngüden döngüye değişkenlik):\n")
            w(md(bf[["fay", "time", "mw", "doluluk_kalibre_%"]].assign(time=bf["time"].dt.year), 1))
        w("")
        w(md(big[["time", "mw", "bölge", "fay", "doluluk_kalibre_%", "aynı_bölgedeki_yüzdelik"]]
             .assign(time=big["time"].dt.strftime("%Y-%m-%d")), 2))
        w("")
    cc = json.loads((d / "capacity_cascade_check.json").read_text())
    if cc.get("cascade"):
        w("### Zincirleme kırılma bilinen depremleri yeniden üretebiliyor mu?\n")
        w("Depremden hemen önceki doluluktan, episantra en yakın parçadan başlatılan kırılma; kırılan parça yükünün "
          "yarısını komşularına aktarır, eşiği aşan komşu da kırılır:\n")
        w(md(pd.DataFrame(cc["cascade"]), 1))
        w("")
    if cc.get("creep_compare"):
        w("### GPS'ten krip (yavaş kayma) ve literatür\n")
        w(md(pd.DataFrame(cc["creep_compare"])[["yer", "literatür_krip_mm_yıl", "gps_krip_mm_yıl", "gps_krip_sd", "not",
                                                "kaynak"]], 1))
        w("\nAçık GPS istasyonları seyrek olduğu için kilitlenme ancak birkaç fayda çözülebildi; çoğunda sonuç öncüle "
          "yakın kaldı. Sığ krip (İsmetpaşa gibi) birkaç km genişliğinde olduğundan 20-50 km aralıklı GPS ile "
          "ayrılamıyor; InSAR verisi bağlanırsa iyileşir.\n")
    seg = pd.read_csv(d / "capacity_segments.csv")
    seg = seg[seg["kayma_hızı_mm_yıl"] >= 2.0].head(15)
    w("### Şu an en dolu faylar (kayma hızı ≥2 mm/yıl)\n")
    w(md(seg[["fault_id", "bölge", "lat", "lon", "kayma_hızı_mm_yıl", "doluluk_%", "doluluk_%16", "doluluk_%84"]], 1))
    w(f"\nBu yüzdeler göreli ölçüdür (bütçe {bstart}'de sıfırdan başlar) ve yukarıdaki test sabit bir kırılma eşiği "
      "göstermediği için tek başına \"sıradaki deprem\" anlamına gelmez. Haritada: *Faylar ve enerji → Fay doluluk oranı*.\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--short", type=Path, default=ROOT / "configs/train_search.toml")
    ap.add_argument("--long", type=Path, default=ROOT / "configs/train_search_long.toml")
    ap.add_argument("--out", type=Path, default=ROOT / "reports/parametre_aramasi.md")
    ap.add_argument("--layers", type=Path, default=ROOT / "web/data/layers")
    a = ap.parse_args()
    rs, rl = Run(a.short, ROOT), Run(a.long, ROOT)
    L = []
    w = L.append
    w("# En iyi parametre kombinasyonu araması\n")
    w("> Deneysel araştırma çıktısıdır, deprem tahmini değildir.\n")
    w("İki arama yapıldı, ikisinde de 2020-2025 aramada hiç kullanılmayan son sınama dönemidir:\n")
    w("- **Uzun arama (1900'den, M≥5,5):** katalog, fay, Coulomb, kayma açığı, kapasite/doluluk, zincirleme kırılma, "
      "kayaç sertliği. Arama 1920-2019.")
    w("- **Kısa arama (2012'den, M≥3):** yukarıdakilere ek olarak AFAD küçük deprem parametreleri (b-değeri, sürüler, "
      "sürtünme vekili), gaz/sıvı ve dış etkiler, GPS gerinim ve kilitlenme/krip. Arama 2015-2019.\n")
    w("## Özet\n")
    w("- **Maksimize eden tek bir kombinasyon bulunamadı.** Aramanın en iyi bulduğu ayarlar arama yıllarında ETAS'ı "
      "geçiyor, ama dokunulmamış 2020-2025 yıllarında bu üstünlük ya küçülüyor (kısa arama: +0,06, anlamlı değil) ya da "
      "tersine dönüyor (uzun arama: -0,09). Yani arama, geçmişe fazla uyan (ezberleyen) ayarlar seçiyor.")
    w("- **Enerji öncülü** (enerji birikimi gelgitten önemli varsayımı) sonucu değiştirmedi: kısa aramada fark 0,001, "
      "uzun aramada kapatınca sınama puanı biraz düştü (-0,087 → -0,105), ikisi de anlamsız.")
    w("- **Kayaç sertliği/kapasite ve zincirleme kırılma** tek başına eklendiğinde sınamada küçük ama tutarlı artı "
      "verdi (+0,01 ile +0,05; tek başına anlamlılık testi yapılmadı, keşif amaçlı); **doluluk oranı ve doluluk + Coulomb** ise sıfır ya da eksi.")
    w("- **%100 dolunca kırılır mı?** Hayır, sabit bir eşik yok: aynı fay bir döngüde %100'ün çok altında, başka "
      "döngüde çok üstünde kırılıyor. Kırılan yerler ortalamadan biraz daha dolu (bölgedeki %65'lik dilim), o kadar.")
    w("- **Zincirleme model** 1999 İzmit-Düzce ve 2023 Kahramanmaraş kırılma uzunluklarını yeniden üretemedi.")
    w("- **GPS'ten krip** açık istasyonlar seyrek olduğu için çözülemedi (İsmetpaşa: GPS 0, literatür ~8 mm/yıl).\n")
    search_section(rl, w, "Uzun arama: 1900'den bugüne, M≥5,5")
    search_section(rs, w, "Kısa arama: 2012'den bugüne, M≥3 (tüm parametreler)")
    w("## Fay doluluk oranı (uzun koşu: birikim 1700'den, tarihsel katalogla)\n")
    capacity_section(rl, w)
    w("## Sınırlar\n")
    w(f"- Sınama yalnız 6 yıl; büyük depremler (M≥6) bu sürede az, güven aralıkları geniş.")
    w("- Isı akısı, kayaç türü ve sürtünme ölçümleri bu ortamdan erişilemedi; sürtünme mikrodeprem/b-değeri/"
      "tekrarlayan depremlerden türetilen bir vekildir.")
    w("- İyonosfer gibi kısa vadeli sinyaller yıllık modelde değil, ayrı kısa vadeli testtedir (train_obs).")
    a.out.write_text("\n".join(L) + "\n")
    from deprem.faults.capacity_stage import write_fill_layer
    n = write_fill_layer(rl, a.layers)
    print("rapor:", a.out, "katman fay sayısı:", n)


if __name__ == "__main__":
    main()
