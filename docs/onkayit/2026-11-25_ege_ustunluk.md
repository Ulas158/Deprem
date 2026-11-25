# Ön kayıt: "Model Ege'de ETAS'tan iyi" ipucunun yeni veriyle sınanması

- **Kayıt tarihi:** 2026-11-25 (UTC). Bu dosya, sınamada kullanılacak 2026 verisi indirilmeden
  ve hiçbir 2026 puanı hesaplanmadan önce yazılıp depoya gönderildi. Git commit zamanı kanıttır.
- **Bu dosyadan sonra değiştirilmeyecek kurallar:** veri, model, ölçüt, geçme/kalma kuralı.
  Sonradan bir kural değişirse bu dosyaya değil, ayrı bir "sapma" notuna yazılır ve sonuç
  ön kayıtlı sayılmaz.

> Deneysel araştırma çıktısıdır, deprem tahmini değildir.

## 1. Neden

Bölge bazında ayar çalışmasında (PR #13, `reports/bolgesel_ayar.md`, kısa koşu, M≥3,
dokunulmamış 2020-2025) tek ortak ayarlı modelin ETAS'a karşı kazancı bölge bölge
ayrıldığında **Ege / Batı Anadolu**'da olay başına +0,315 (M≥3), diğer beş bölgede eksi
ya da sıfıra yakın çıktı. Bu bölge sonuçlara bakılarak öne çıktığı için (6 bölgeden en iyisi)
şans olabilir. Burada aynı iddia, o çalışmada hiç kullanılmamış olaylarla bir kez sınanır.

## 2. Hipotez

**H1:** Aşağıda tanımlanan, önceden dondurulmuş model, Ege / Batı Anadolu kutusunda
2026-01-01 ile 2026-11-25 arasındaki M≥3 depremlere ETAS'tan daha yüksek olasılık verir
(olay başına bilgi kazancı > 0).

**H0:** Kazanç ≤ 0 (Ege'deki üstünlük, bölge seçiminden doğan şanstı).

## 3. Veri (yeni, daha önce hiç puanlanmamış)

- **Bölge:** 0,5° hücreler; hücre merkezi enlem 36,5-40,0 K ve boylam 25,5-30,5 D aralığında
  olanlar (`src/deprem/regional.py` içindeki "Ege / Batı Anadolu" kutusu, aynen; Marmara
  kutusu 40,0 K'nin kuzeyidir). Türkiye ve Yunanistan tarafındaki Ege adaları ve denizi dahildir.
- **Dönem:** 2026-01-01 00:00 UTC ≤ t < 2026-11-25 00:00 UTC (328 gün). Önceki tüm çalışmalar
  2025 sonunda kesiliyordu; bu dönem hiçbir aramada, ayarda veya raporda kullanılmadı.
- **Katalog:** `scripts/fetch_catalogs.py --end 2026-11-25` ile oluşturulan birleşik katalog
  (USGS + EMSC + AFAD; çakışan olaylarda öncelik AFAD > EMSC > USGS), derinlik ≤ 70 km,
  projenin olağan temizleme adımlarıyla. AFAD'ın ISC arşivi 2026'yı henüz kapsamadığından
  2026 için AFAD açık API'sinden 2026-01-01'den itibaren tüm günler çekilir; AFAD API'sine
  erişilemezse 2026 EMSC + USGS ile puanlanır. Bu seçim erişilebilirliğe göre yapılır,
  sonuçlara göre değil; hangi kaynağın kullanıldığı raporda yazılır.
- **Büyüklük:** birincil eşik M≥3,0 (katalogdaki tercih edilen büyüklük, önceki çalışmalardaki gibi).

## 4. Model (dondurulmuş) ve referans

- **Model:** PR #13'teki "ortak ayar (tek Türkiye)" modeli:
  λ(hücre, 2026) = ETAS(hücre, 2026) · exp(β0 + Σ β_k z_k).
  Ayar (açık parametre grupları, ağırlıklar, öncül, μ' vb.) PR #13'teki prosedürle aynen
  seçilir: `scripts/search_regional.py` içindeki ortak arama (TPE, tohum 0, 300 deneme,
  yalnız 2015-2019 arama yılları, `configs/train_search.toml`). Bu arama yalnız 2015-2019
  yıllarını puanladığı için 2026 verisini göremez. Seçilen ayar, 2026 puanı hesaplanmadan önce
  `docs/onkayit/2026-11-25_ege_dondurulmus_ayar.json` olarak ayrı bir commit'le depoya yazılır.
- **Eğitim:** katsayılar 2012-2025 yıllarıyla (2026 öncesi) bir kez eğitilir, sonra donar.
  Özellikler yalnız 2026-01-01 öncesi veriden hesaplanır (projenin olağan kuralı).
- **Referans ETAS:** projenin ileri yürüyen ETAS'ı; parametreler 2010-01-01 ile 2026-01-01
  arasındaki veriyle uydurulur, 2026-01-01 → 2026-11-25 penceresi için beklenen sayı
  benzetimle (200 benzetim, tohum 0 + 2026) hesaplanır. Model aynı ETAS'ı ofset olarak kullanır.
- **Tekrar üretim kontrolü (karar vermez):** aynı kurulumla 2020-2025 Ege sonucu yeniden
  hesaplanır; PR #13'teki +0,315'e yakın çıkmazsa bu, raporda açıkça yazılır.

## 5. Ölçüt

**Birincil ölçüt: bölgeye kısıtlı olay başına bilgi kazancı (doğal log birimi).**
Bölgedeki hücreler R, her hücrede gözlenen M≥3 olay sayısı n_c, toplam n = Σ n_c, modelin ve
ETAS'ın beklenen M≥3 sayıları λm_c ve λe_c (M≥3'e projenin b-değeri ölçeklemesiyle):

    IG = [ Σ_c n_c · ln(λm_c / λe_c)  −  (Σ_c λm_c − Σ_c λe_c) ] / n      (c ∈ R)

Bu, iki modelin yalnız Ege kutusu için verdiği Poisson tahminlerinin log-olabilirlik farkıdır.
Not: PR #13'teki +0,315, sayı farkı terimini tüm Türkiye için hesaplıyordu (bölgeye kısıtlı
değildi). Karşılaştırma için o eski ölçüt de raporlanır, ama karar birincil ölçütle verilir.

**Belirsizlik:** hücre blokları bootstrap'ı (Ege hücreleri, olay sayıları ve beklenen sayılarıyla
birlikte yerine koyarak yeniden örneklenir; 10 000 tekrar, tohum 2026). Aynı hücredeki artçı
dizileri böylece tek blok sayılır. Tek yönlü p = bootstrap örneklerinde IG ≤ 0 oranı.

## 6. Karar kuralı (değişmez)

- **Geçti (ipucu desteklendi):** IG > 0 **ve** tek yönlü p < 0,05.
- **Kaldı (ipucu desteklenmedi):** IG ≤ 0.
- **Kararsız:** IG > 0 ama p ≥ 0,05; ya da bölgede 30'dan az M≥3 olay.
- Beklenen etki büyüklüğü: eski ölçütle +0,315. "Geçti" sonucu IG ≥ +0,10 ise "güçlü",
  0 < IG < +0,10 ise "zayıf" diye nitelenir.
- "Geçti" bile olsa bu yalnız bir sınamadır; **doğrulanmış** sayılması için bölüm 8'deki
  canlı (ileriye dönük) testin de aynı yönde çıkması gerekir.

## 7. İkincil sonuçlar (karar vermez, yalnız raporlanır)

1. Aynı ölçüt M≥4 eşiğinde (az olay, geniş belirsizlik).
2. Eski ölçüt (sayı farkı terimi tüm Türkiye'den), PR #13 ile kıyas için.
3. Küçük depremler: 2,5 ≤ M < 3,0 bandındaki olaylar (hiçbir çalışmada hedef olarak
   puanlanmadı), 2020-2025 ve 2026 için ayrı ayrı; tahminler M≥3 beklenen sayısından
   b-değeriyle ölçeklenir. Ege denizinde M2,5 tamlığı zayıf olabileceğinden yalnız bilgi amaçlı.
4. Türkiye'nin geri kalanında (Ege dışı) aynı ölçüt, 2026: ipucu Ege'ye özgüyse burada
   kazanç beklenmez.

## 8. Canlı (ileriye dönük) test

2026-11-26'dan itibaren her gün, ertesi 24 saat için Ege kutusundaki her hücrenin M≥3 ve M≥4
beklenen sayıları hem ETAS hem model için (model = günlük ETAS × 2026 için dondurulmuş
exp(β0 + Σ β_k z_k)) zaman damgasıyla depoya yazılır. Değerlendirme en erken 2027-05-26
(6 ay), asıl değerlendirme 2027-11-26 (12 ay) tarihinde, aynı birincil ölçütle (günlük
tahminler toplanarak) ve aynı karar kuralıyla yapılır; bootstrap blokları hücre × takvim ayıdır.
Model katsayıları ve hücre çarpanları bu sürede değiştirilmez (2027'ye geçilince de
2026 çarpanı kullanılır; yalnız ETAS kısmı her gün yeni depremlerle güncellenir). Günlük
ETAS: parametreler bölüm 4'teki 2026 ETAS'ı (donmuş); her gün, o güne kadarki katalogla
koşullu yoğunluğun ertesi 24 saat üzerindeki integrali. Yöntem ayrıntısı `scripts/ege_canli.py` ve
`forecasts/ege/BENIOKU.md` içindedir.
