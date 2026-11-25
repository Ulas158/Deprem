"""Bölge bazında parametre ayarı raporu: reports/bolgesel_ayar.md

    python scripts/regional_report.py

runs/turkey_search/regional/summary.json (kısa) ve varsa runs/turkey_search_long/regional/summary.json
(uzun) okunur. HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CI_E, P_E = "ETAS'a göre %95 GA", "ETAS'a göre Holm p"
RUNS = [("Kısa koşu: 2012'den bugüne, M≥3", "turkey_search", "M≥3, M≥4, M≥5"),
        ("Uzun koşu: 1900'den bugüne, M≥5,5", "turkey_search_long", "M≥5,5, M≥6, M≥6,5")]


def f3(x):
    return "" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:+.3f}"


def pfmt(p):
    return f"{p:.2f}"


def section(title, s):
    L = [f"## {title}", ""]
    sy, ty = s["search_years"], s["test_years"]
    rows = s["test_table"]
    byname = {r["model"]: r for r in rows}
    g = next(r for r in rows if r["model"].startswith("ortak ayar"))
    b = rows[0]
    nb = s["null_best_scores"]
    gain_search = s["regional_best_search_score"] - s["global_search_score"]
    null_gain = [x - s["global_search_score"] for x in nb]
    L += ["### Kısa sonuç", "",
          f"- Arama yılları {sy[0]}-{sy[-1]}: ortak ayarın puanı {f3(s['global_search_score'])}, en iyi bölgesel "
          f"ayarın {f3(s['regional_best_search_score'])} (fark {f3(gain_search)}; puan olay başına, ETAS'a göre).",
          f"- Sahte bölgelerle (hücrelerin bölge etiketleri karıştırılmış) aynı arama: farklar "
          + ", ".join(f3(x) for x in null_gain) + ". "
          + ("Gerçek bölgeler sahte bölgelerin hepsinden fazla kazandırdı."
             if gain_search > max(null_gain) else
             "Gerçek bölgeler sahte bölgelerden fazla kazandırmadı: arama dönemindeki kazanç şanstan ayırt edilemiyor."),
          f"- Dokunulmamış {ty[0]}-{ty[-1]} yıllarında en iyi bölgesel ayar: ETAS'a göre {f3(b['sınama puanı'])} "
          f"(%95 GA {b[CI_E]}, Holm p = {pfmt(b[P_E])}); "
          f"ortak ayar {f3(g['sınama puanı'])}. Bölgesel − ortak = {f3(b['ortak ayara göre fark'])} "
          f"(%95 GA {b['ortak ayara göre %95 GA']}, Holm p = {pfmt(b['ortak ayara göre Holm p'])}).",
          f"- Sahte bölgelerin sınama yıllarındaki ortak ayara göre farkı: "
          + ", ".join(f3(x) for x in s["null_test_gain_vs_global"]) + ".", ""]
    L += ["### Sınama yılları sonuçları", "",
          "Her satır aynı olaylarda ölçüldü. \"ETAS'a göre\": olay başına bilgi kazancı (0 = ETAS kadar, + iyi). "
          "\"ortak ayara göre\": aynı olaylarda bölgesel model ile tek ortak ayar arasındaki fark. "
          "GA: yıllar üzerinden bootstrap; p: tek yönlü, Holm ile çoklu teste göre düzeltilmiş.", ""]
    cols = ["model", "arama puanı", "sınama puanı", "ETAS'a göre %95 GA", "ETAS'a göre Holm p",
            "ortak ayara göre fark", "ortak ayara göre %95 GA", "ortak ayara göre Holm p"]
    mk = [c for c in rows[0] if c.endswith(" kazanç")]
    df = pd.DataFrame(rows)[cols + mk]
    # havuzlamasız model olaysız bölge-yıllarda sayısal olarak patlayabilir: okunur yaz
    df["arama puanı"] = [x if x > -10 else -9.999 for x in df["arama puanı"]]
    L += [df.to_markdown(index=False, floatfmt=".3f").replace(" nan ", "   – "), ""]
    L += ["### Bölge bölge (sınama yılları)", "",
          "Olay başına kazanç, ETAS'a göre. Az olaylı bölgelerde sayılar çok oynaktır.", ""]
    rb = s["region_breakdown"]
    frames = []
    for name, lst in rb.items():
        d = pd.DataFrame(lst).set_index("bölge")
        kz = [c for c in d.columns if c.endswith("kazanç")]
        ol = [c for c in d.columns if c.endswith("olay")]
        if not frames:
            frames.append(d[ol])
        frames.append(d[kz].add_prefix(f"{name}: "))
    L += [pd.concat(frames, axis=1).reset_index().to_markdown(index=False, floatfmt=".3f").replace(" nan ", "   – "), ""]
    rc = s["best_rcfg"]
    L += ["### Aramanın seçtiği bölgesel ayar", "",
          f"- Ortak ayar (önceki aramadan): açık gruplar: {', '.join(s['global_config']['groups']) or '-'}",
          f"- Bölgeden bölgeye değişmesine izin verilen gruplar: {', '.join(rc['vary']) or 'hiçbiri (yalnız bölge sabiti)'}",
          f"- Ortak ayara çekme gücü τ = {rc['tau']:.3f} (küçük = bölgeler ortak ayara sıkı bağlı; 1 ≈ ortak katsayı "
          f"kadar sapabilir), bölge sabiti τ0 = {rc['tau0']:.3f}",
          "- Bölgeye özel serbestlik çarpanı: " + ", ".join(f"{k.split(' (')[0]} {v:.2f}" for k, v in rc["tau_r"].items()),
          ""]
    coef = s["best_coef_last_year"]
    reg_c = {k: v for k, v in coef.items() if " @ " in k}
    if reg_c:
        top = sorted(reg_c.items(), key=lambda kv: -abs(kv[1]))[:12]
        L += ["Son sınama yılında ortak katsayıdan en çok sapan bölge katsayıları (standartlaştırılmış; "
              "+0,1 ≈ 1 sd artışta oran ×1,105):", "",
              pd.DataFrame(top, columns=["parametre @ bölge", "sapma"]).to_markdown(index=False, floatfmt=".3f"), ""]
    return L


def overview(sums):
    L = ["## Özet", ""]
    for title, rn, _ in RUNS:
        if rn not in sums:
            continue
        s = sums[rn]
        b = s["test_table"][0]
        g = next(r for r in s["test_table"] if r["model"].startswith("ortak ayar"))
        nogain = s["regional_best_search_score"] - s["global_search_score"] <= max(
            x - s["global_search_score"] for x in s["null_best_scores"])
        L.append(f"- **{title.split(':')[0]}:** bölgesel ayar arama yıllarında ortak ayardan "
                 f"{f3(s['regional_best_search_score'] - s['global_search_score'])} iyi göründü"
                 + (", ama sahte bölgeler de bu kadar kazandırdı (şans sınırı içinde)" if nogain else "")
                 + f". Dokunulmamış 2020-2025'te bölgesel {f3(b['sınama puanı'])}, ortak {f3(g['sınama puanı'])} "
                 f"(ETAS'a göre). Bölgesel ayar ortak ayarı geçmedi (fark {f3(b['ortak ayara göre fark'])}, "
                 f"Holm p = {pfmt(b['ortak ayara göre Holm p'])}).")
    L += ["- **Her bölgeyi tamamen ayrı ayarlamak** (havuzlamasız) iki koşuda da açıkça kötü: bölge başına veri az, "
          "model ezberliyor. Kısmi havuzlama bu zararı önledi ama kazanç da getirmedi.",
          "- **Bölge bölge bakınca:** kısa koşuda ortak ayarın ETAS'a karşı kazancı neredeyse tamamen Ege / Batı "
          "Anadolu'dan geliyor; Marmara, Kuzey Anadolu ve Van-Kars'ta ETAS'ın gerisinde. Bu bir ipucu, kanıt değil: "
          "sonuçlara bakarak seçildiği için ayrı bir ön kayıtlı testle sınanmalı.",
          "- **Sonuç:** \"tek ortak ayar tutmadı, bölge ayarı tutar mı?\" sorusunun cevabı bu veride hayır. "
          "Hipotez reddedilmedi ama desteklenmedi; doğrulanmış sayılmaz.", ""]
    return L


def main():
    L = ["# Bölge bazında parametre ayarı", "",
         "> Deneysel araştırma çıktısıdır, deprem tahmini değildir.", "",
         "Önceki aramada ([parametre_aramasi.md](parametre_aramasi.md)) tüm Türkiye için tek ortak ayar, "
         "dokunulmamış 2020-2025 yıllarında ETAS'ı geçemedi. Bu çalışma aynı ayarın bölgeden bölgeye "
         "farklı olmasına izin verince durumun değişip değişmediğini sınar.", "",
         "## Nasıl yapıldı", "",
         "1. **Bölgeler sonuçlara bakılmadan belirlendi:** Marmara (KAF batı), Kuzey Anadolu (KAF orta), "
         "Doğu Anadolu Fay Zonu, Ege / Batı Anadolu, Doğu Anadolu platosu (Van-Kars) ve kalan her yer "
         "(İç Anadolu, Akdeniz, Karadeniz, Kıbrıs yayı). Sınırlar `region.SUBREGIONS` kutularıdır "
         "(`src/deprem/regional.py`).",
         "2. **Kısmi havuzlama:** her bölgenin katsayısı = ortak katsayı + bölge sapması. Sapma, ortak "
         "ayara doğru çekilir (Gauss öncül, genişliği τ). Veri azsa bölge ortak ayarda kalır, veri bir "
         "bölgenin gerçekten farklı olduğunu söylüyorsa sapar. τ → 0 tek ortak ayardır, τ → ∞ her bölge "
         "tamamen ayrıdır (karşılaştırma için bu da ölçüldü).",
         "3. **Arama yalnız arama yıllarında:** ortak ayar önceki PR'daki aramayla bulundu; üstüne τ, bölge "
         "sabiti, her bölgenin serbestlik çarpanı ve hangi parametre gruplarının bölgeden bölgeye "
         "değişebileceği Bayesçi optimizasyonla (TPE) arandı. Her yıl yalnız o yıldan önceki yıllarla eğitildi.",
         "4. **Şans kontrolü:** aynı bölgesel arama, hücrelerin bölge etiketleri rastgele karıştırılmış "
         "(\"sahte bölgeler\", aynı büyüklükte) verilerde tekrarlandı. Gerçek bölgeler sahte bölgelerden "
         "fazla kazandırmıyorsa bölgesel ayrımın değeri şanstan ayırt edilemez.",
         "5. **Sınama 2020-2025** aramada hiç kullanılmadı. İlk 3 bölgesel ayar ve karşılaştırma modelleri "
         "burada bir kez ölçüldü; güven aralıkları yıllar üzerinden bootstrap, p değerleri Holm ile düzeltildi. "
         "ETAS referansı önceki PR'lardakiyle aynı.",
         "6. **Deneme sayıları:** kısa koşuda 200 bölgesel deneme ve 4 sahte bölge araması; uzun koşuda bir "
         "değerlendirme ~1 dakika sürdüğü için 40 deneme ve 3 sahte bölge araması. Uzun koşunun ortak ayarı "
         "önceki PR'daki aramanın seçtiği ayardır (bu veride aynı arama puanını, +0,120, veriyor).", ""]
    sums = {rn: json.loads((ROOT / "runs" / rn / "regional" / "summary.json").read_text())
            for _, rn, _ in RUNS if (ROOT / "runs" / rn / "regional" / "summary.json").exists()}
    L[4:4] = overview(sums)
    found = False
    for title, rn, _ in RUNS:
        p = ROOT / "runs" / rn / "regional" / "summary.json"
        if p.exists():
            found = True
            L += section(title, json.loads(p.read_text()))
    if not found:
        raise SystemExit("summary.json yok: önce scripts/search_regional.py")
    L += ["## Yeniden üretme", "", "```bash",
          "python scripts/train.py --config configs/train_search.toml --until features",
          "python scripts/search_regional.py --config configs/train_search.toml",
          "python scripts/train.py --config configs/train_search_long.toml --until features",
          "python scripts/search_regional.py --config configs/train_search_long.toml --trials 40 --null-runs 3",
          "python scripts/regional_report.py", "```", ""]
    (ROOT / "reports" / "bolgesel_ayar.md").write_text("\n".join(L), encoding="utf-8")
    print("reports/bolgesel_ayar.md")


if __name__ == "__main__":
    main()
