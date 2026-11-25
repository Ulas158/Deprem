# Bölge bazında parametre ayarı

> Deneysel araştırma çıktısıdır, deprem tahmini değildir.

## Özet

- **Kısa koşu:** bölgesel ayar arama yıllarında ortak ayardan +0.024 iyi göründü, ama sahte bölgeler de bu kadar kazandırdı (şans sınırı içinde). Dokunulmamış 2020-2025'te bölgesel +0.036, ortak +0.049 (ETAS'a göre). Bölgesel ayar ortak ayarı geçmedi (fark -0.013, Holm p = 1.00).
- **Uzun koşu:** bölgesel ayar arama yıllarında ortak ayardan +0.019 iyi göründü, ama sahte bölgeler de bu kadar kazandırdı (şans sınırı içinde). Dokunulmamış 2020-2025'te bölgesel -0.136, ortak -0.087 (ETAS'a göre). Bölgesel ayar ortak ayarı geçmedi (fark -0.049, Holm p = 1.00).
- **Her bölgeyi tamamen ayrı ayarlamak** (havuzlamasız) iki koşuda da açıkça kötü: bölge başına veri az, model ezberliyor. Kısmi havuzlama bu zararı önledi ama kazanç da getirmedi.
- **Bölge bölge bakınca:** kısa koşuda ortak ayarın ETAS'a karşı kazancı neredeyse tamamen Ege / Batı Anadolu'dan geliyor; Marmara, Kuzey Anadolu ve Van-Kars'ta ETAS'ın gerisinde. Bu bir ipucu, kanıt değil: sonuçlara bakarak seçildiği için ayrı bir ön kayıtlı testle sınanmalı.
- **Sonuç:** "tek ortak ayar tutmadı, bölge ayarı tutar mı?" sorusunun cevabı bu veride hayır. Hipotez reddedilmedi ama desteklenmedi; doğrulanmış sayılmaz.

Önceki aramada ([parametre_aramasi.md](parametre_aramasi.md)) tüm Türkiye için tek ortak ayar, dokunulmamış 2020-2025 yıllarında ETAS'ı geçemedi. Bu çalışma aynı ayarın bölgeden bölgeye farklı olmasına izin verince durumun değişip değişmediğini sınar.

## Nasıl yapıldı

1. **Bölgeler sonuçlara bakılmadan belirlendi:** Marmara (KAF batı), Kuzey Anadolu (KAF orta), Doğu Anadolu Fay Zonu, Ege / Batı Anadolu, Doğu Anadolu platosu (Van-Kars) ve kalan her yer (İç Anadolu, Akdeniz, Karadeniz, Kıbrıs yayı). Sınırlar `region.SUBREGIONS` kutularıdır (`src/deprem/regional.py`).
2. **Kısmi havuzlama:** her bölgenin katsayısı = ortak katsayı + bölge sapması. Sapma, ortak ayara doğru çekilir (Gauss öncül, genişliği τ). Veri azsa bölge ortak ayarda kalır, veri bir bölgenin gerçekten farklı olduğunu söylüyorsa sapar. τ → 0 tek ortak ayardır, τ → ∞ her bölge tamamen ayrıdır (karşılaştırma için bu da ölçüldü).
3. **Arama yalnız arama yıllarında:** ortak ayar önceki PR'daki aramayla bulundu; üstüne τ, bölge sabiti, her bölgenin serbestlik çarpanı ve hangi parametre gruplarının bölgeden bölgeye değişebileceği Bayesçi optimizasyonla (TPE) arandı. Her yıl yalnız o yıldan önceki yıllarla eğitildi.
4. **Şans kontrolü:** aynı bölgesel arama, hücrelerin bölge etiketleri rastgele karıştırılmış ("sahte bölgeler", aynı büyüklükte) verilerde tekrarlandı. Gerçek bölgeler sahte bölgelerden fazla kazandırmıyorsa bölgesel ayrımın değeri şanstan ayırt edilemez.
5. **Sınama 2020-2025** aramada hiç kullanılmadı. İlk 3 bölgesel ayar ve karşılaştırma modelleri burada bir kez ölçüldü; güven aralıkları yıllar üzerinden bootstrap, p değerleri Holm ile düzeltildi. ETAS referansı önceki PR'lardakiyle aynı.
6. **Deneme sayıları:** kısa koşuda 200 bölgesel deneme ve 4 sahte bölge araması; uzun koşuda bir değerlendirme ~1 dakika sürdüğü için 40 deneme ve 3 sahte bölge araması. Uzun koşunun ortak ayarı önceki PR'daki aramanın seçtiği ayardır (bu veride aynı arama puanını, +0,120, veriyor).

## Kısa koşu: 2012'den bugüne, M≥3

### Kısa sonuç

- Arama yılları 2015-2019: ortak ayarın puanı +0.207, en iyi bölgesel ayarın +0.231 (fark +0.024; puan olay başına, ETAS'a göre).
- Sahte bölgelerle (hücrelerin bölge etiketleri karıştırılmış) aynı arama: farklar +0.046, -0.001, +0.005, +0.001. Gerçek bölgeler sahte bölgelerden fazla kazandırmadı: arama dönemindeki kazanç şanstan ayırt edilemiyor.
- Dokunulmamış 2020-2025 yıllarında en iyi bölgesel ayar: ETAS'a göre +0.036 (%95 GA [-0.054, 0.212], Holm p = 1.00); ortak ayar +0.049. Bölgesel − ortak = -0.013 (%95 GA [-0.034, 0.027], Holm p = 1.00).
- Sahte bölgelerin sınama yıllarındaki ortak ayara göre farkı: -0.062, -0.003, +0.024, +0.001.

### Sınama yılları sonuçları

Her satır aynı olaylarda ölçüldü. "ETAS'a göre": olay başına bilgi kazancı (0 = ETAS kadar, + iyi). "ortak ayara göre": aynı olaylarda bölgesel model ile tek ortak ayar arasındaki fark. GA: yıllar üzerinden bootstrap; p: tek yönlü, Holm ile çoklu teste göre düzeltilmiş.

| model                                                 |   arama puanı |   sınama puanı | ETAS'a göre %95 GA   |   ETAS'a göre Holm p |   ortak ayara göre fark | ortak ayara göre %95 GA   |   ortak ayara göre Holm p |   M≥3 kazanç |   M≥4 kazanç |   M≥5 kazanç |   M≥6 kazanç |
|:------------------------------------------------------|--------------:|---------------:|:---------------------|---------------------:|------------------------:|:--------------------------|--------------------------:|-------------:|-------------:|-------------:|-------------:|
| bölgesel: aramada 1. sıradaki                         |         0.231 |          0.036 | [-0.054, 0.212]      |                1.000 |                  -0.013 | [-0.034, 0.027]           |                     1.000 |        0.077 |        0.010 |        0.021 |        0.051 |
| bölgesel: aramada 2. sıradaki                         |         0.227 |          0.029 | [-0.070, 0.190]      |                1.000 |                  -0.020 | [-0.049, 0.031]           |                     1.000 |        0.064 |       -0.007 |        0.029 |        0.051 |
| bölgesel: aramada 3. sıradaki                         |         0.226 |          0.024 | [-0.086, 0.204]      |                1.000 |                  -0.025 | [-0.063, 0.029]           |                     1.000 |        0.062 |       -0.001 |        0.010 |        0.011 |
| yalnız ETAS'ın kendi düzeltmesi (özelliksiz)          |         0.022 |         -0.014 | [-0.026, 0.027]      |                1.000 |                  -0.063 | [-0.196, 0.004]           |                     1.000 |       -0.010 |       -0.023 |       -0.011 |       -0.020 |
| ortak ayar (tek Türkiye, önceki PR'daki gibi aranmış) |         0.207 |          0.049 | [-0.022, 0.185]      |                0.834 |                   0.000 | [0.000, 0.000]            |                     –     |        0.078 |        0.033 |        0.036 |        0.053 |
| bölgesel: yalnız bölge sabiti (katsayılar ortak)      |         0.208 |          0.049 | [-0.023, 0.186]      |                0.834 |                  -0.000 | [-0.003, 0.002]           |                     1.000 |        0.078 |        0.033 |        0.035 |        0.054 |
| havuzlamasız: her bölge tamamen ayrı                  |         0.077 |         -0.147 | [-0.424, 0.208]      |                1.000 |                  -0.196 | [-0.403, 0.031]           |                     1.000 |       -0.106 |       -0.182 |       -0.152 |       -0.244 |
| önceki varsayılan model + bölgesel en iyi ayar        |         0.131 |          0.043 | [-0.086, 0.262]      |                1.000 |                  -0.006 | [-0.062, 0.089]           |                     1.000 |        0.088 |        0.018 |        0.023 |        0.018 |

### Bölge bölge (sınama yılları)

Olay başına kazanç, ETAS'a göre. Az olaylı bölgelerde sayılar çok oynaktır.

| bölge                                               |   M≥3 olay |   M≥4 olay |   M≥5 olay |   M≥6 olay |   ortak ayar: M≥3 kazanç |   ortak ayar: M≥4 kazanç |   ortak ayar: M≥5 kazanç |   ortak ayar: M≥6 kazanç |   bölgesel en iyi: M≥3 kazanç |   bölgesel en iyi: M≥4 kazanç |   bölgesel en iyi: M≥5 kazanç |   bölgesel en iyi: M≥6 kazanç |   havuzlamasız: M≥3 kazanç |   havuzlamasız: M≥4 kazanç |   havuzlamasız: M≥5 kazanç |   havuzlamasız: M≥6 kazanç |
|:----------------------------------------------------|-----------:|-----------:|-----------:|-----------:|-------------------------:|-------------------------:|-------------------------:|-------------------------:|------------------------------:|------------------------------:|------------------------------:|------------------------------:|---------------------------:|---------------------------:|---------------------------:|---------------------------:|
| Marmara (KAF batı)                                  |        349 |         44 |          6 |          1 |                   -0.131 |                   -0.150 |                   -0.216 |                   -0.339 |                        -0.285 |                        -0.286 |                        -0.403 |                        -0.473 |                     -0.439 |                     -0.362 |                     -0.587 |                     -0.772 |
| Kuzey Anadolu (KAF orta)                            |        251 |         41 |          2 |          1 |                   -0.128 |                   -0.141 |                    0.093 |                   -0.211 |                        -0.295 |                        -0.272 |                        -0.054 |                        -0.475 |                     -0.361 |                     -0.458 |                     -1.212 |                     -1.839 |
| Doğu Anadolu Fay Zonu                               |       4521 |        760 |         51 |          8 |                   -0.023 |                   -0.033 |                    0.032 |                    0.091 |                        -0.067 |                        -0.101 |                        -0.020 |                         0.032 |                     -0.544 |                     -0.548 |                     -0.401 |                     -0.327 |
| Ege / Batı Anadolu                                  |       4763 |        562 |         33 |          3 |                    0.315 |                    0.315 |                    0.341 |                    0.281 |                         0.355 |                         0.323 |                         0.362 |                         0.420 |                      0.355 |                      0.320 |                      0.345 |                      0.511 |
| Doğu Anadolu platosu (Van-Kars)                     |        793 |        131 |         13 |          1 |                   -0.173 |                   -0.221 |                   -0.261 |                   -0.565 |                        -0.039 |                        -0.062 |                        -0.068 |                        -0.084 |                     -0.131 |                     -0.109 |                     -0.109 |                     -0.102 |
| diğer (İç Anadolu, Akdeniz, Karadeniz, Kıbrıs yayı) |       2384 |        358 |         25 |          1 |                   -0.068 |                   -0.137 |                   -0.150 |                    0.341 |                        -0.074 |                        -0.152 |                        -0.192 |                         0.273 |                     -0.110 |                     -0.168 |                     -0.136 |                      0.131 |

### Aramanın seçtiği bölgesel ayar

- Ortak ayar (önceki aramadan): açık gruplar: b-değeri, fay geometrisi, gelgit ve periyodik yükler, hız değişimi, insan kaynaklı sahalar, kapasite ve kayaç, sayılar, son depremden geçen süre, sürtünme vekili, zincirleme kırılma
- Bölgeden bölgeye değişmesine izin verilen gruplar: fay geometrisi, insan kaynaklı sahalar, zincirleme kırılma
- Ortak ayara çekme gücü τ = 4.872 (küçük = bölgeler ortak ayara sıkı bağlı; 1 ≈ ortak katsayı kadar sapabilir), bölge sabiti τ0 = 0.020
- Bölgeye özel serbestlik çarpanı: Marmara 0.91, Kuzey Anadolu 0.25, Doğu Anadolu Fay Zonu 4.17, Ege / Batı Anadolu 4.47, Doğu Anadolu platosu 2.03, diğer 0.62

Son sınama yılında ortak katsayıdan en çok sapan bölge katsayıları (standartlaştırılmış; +0,1 ≈ 1 sd artışta oran ×1,105):

| parametre @ bölge                                                 |   sapma |
|:------------------------------------------------------------------|--------:|
| site_gas_seep_dist_km @ Doğu Anadolu Fay Zonu                     |   1.141 |
| site_injection_dist_km @ Doğu Anadolu platosu (Van-Kars)          |   0.980 |
| site_geothermal_dist_km @ Ege / Batı Anadolu                      |   0.901 |
| site_hydrocarbon_dist_km @ Doğu Anadolu platosu (Van-Kars)        |  -0.875 |
| cascade_mw__var @ Doğu Anadolu platosu (Van-Kars)                 |   0.740 |
| fault_log_m0_rate__eksik @ Ege / Batı Anadolu                     |   0.727 |
| cascade_len_km @ Doğu Anadolu platosu (Van-Kars)                  |  -0.685 |
| site_gas_seep_dist_km @ Ege / Batı Anadolu                        |  -0.625 |
| site_injection_dist_km @ Ege / Batı Anadolu                       |  -0.525 |
| site_gas_seep_dist_km @ Doğu Anadolu platosu (Van-Kars)           |  -0.483 |
| reservoir_years_since_fill__var @ Doğu Anadolu platosu (Van-Kars) |   0.454 |
| fault_dist_km @ Doğu Anadolu platosu (Van-Kars)                   |   0.443 |

## Uzun koşu: 1900'den bugüne, M≥5,5

### Kısa sonuç

- Arama yılları 1920-2019: ortak ayarın puanı +0.120, en iyi bölgesel ayarın +0.139 (fark +0.019; puan olay başına, ETAS'a göre).
- Sahte bölgelerle (hücrelerin bölge etiketleri karıştırılmış) aynı arama: farklar +0.005, +0.031, +0.021. Gerçek bölgeler sahte bölgelerden fazla kazandırmadı: arama dönemindeki kazanç şanstan ayırt edilemiyor.
- Dokunulmamış 2020-2025 yıllarında en iyi bölgesel ayar: ETAS'a göre -0.136 (%95 GA [-0.250, -0.055], Holm p = 1.00); ortak ayar -0.087. Bölgesel − ortak = -0.049 (%95 GA [-0.079, 0.058], Holm p = 1.00).
- Sahte bölgelerin sınama yıllarındaki ortak ayara göre farkı: -0.025, -0.049, -0.018.

### Sınama yılları sonuçları

Her satır aynı olaylarda ölçüldü. "ETAS'a göre": olay başına bilgi kazancı (0 = ETAS kadar, + iyi). "ortak ayara göre": aynı olaylarda bölgesel model ile tek ortak ayar arasındaki fark. GA: yıllar üzerinden bootstrap; p: tek yönlü, Holm ile çoklu teste göre düzeltilmiş.

| model                                                 |   arama puanı |   sınama puanı | ETAS'a göre %95 GA   |   ETAS'a göre Holm p |   ortak ayara göre fark | ortak ayara göre %95 GA   |   ortak ayara göre Holm p |   M≥5.5 kazanç |   M≥6 kazanç |   M≥6.5 kazanç |
|:------------------------------------------------------|--------------:|---------------:|:---------------------|---------------------:|------------------------:|:--------------------------|--------------------------:|---------------:|-------------:|---------------:|
| bölgesel: aramada 1. sıradaki                         |         0.139 |         -0.136 | [-0.250, -0.055]     |                1.000 |                  -0.049 | [-0.079, 0.058]           |                     1.000 |         -0.140 |       -0.078 |         -0.191 |
| bölgesel: aramada 2. sıradaki                         |         0.135 |         -0.106 | [-0.228, -0.043]     |                1.000 |                  -0.020 | [-0.031, 0.042]           |                     1.000 |         -0.126 |       -0.054 |         -0.140 |
| bölgesel: aramada 3. sıradaki                         |         0.135 |         -0.162 | [-0.250, -0.083]     |                1.000 |                  -0.075 | [-0.093, 0.050]           |                     1.000 |         -0.163 |       -0.103 |         -0.220 |
| yalnız ETAS'ın kendi düzeltmesi (özelliksiz)          |        -0.159 |         -0.041 | [-0.055, 0.021]      |                1.000 |                   0.046 | [-0.006, 0.233]           |                     0.327 |         -0.036 |       -0.025 |         -0.061 |
| ortak ayar (tek Türkiye, önceki PR'daki gibi aranmış) |         0.120 |         -0.087 | [-0.230, -0.037]     |                1.000 |                   0.000 | [0.000, 0.000]            |                     –     |         -0.118 |       -0.052 |         -0.090 |
| bölgesel: yalnız bölge sabiti (katsayılar ortak)      |         0.119 |         -0.100 | [-0.229, -0.056]     |                1.000 |                  -0.014 | [-0.019, -0.000]          |                     1.000 |         -0.128 |       -0.065 |         -0.107 |
| havuzlamasız: her bölge tamamen ayrı                  |        -9.999 |         -0.385 | [-0.652, -0.322]     |                1.000 |                  -0.299 | [-0.550, -0.167]          |                     1.000 |         -0.461 |       -0.468 |         -0.227 |
| önceki varsayılan model + bölgesel en iyi ayar        |        -0.147 |         -0.051 | [-0.386, 0.111]      |                1.000 |                   0.036 | [-0.250, 0.179]           |                     1.000 |         -0.030 |        0.080 |         -0.203 |

### Bölge bölge (sınama yılları)

Olay başına kazanç, ETAS'a göre. Az olaylı bölgelerde sayılar çok oynaktır.

| bölge                                               |   M≥5.5 olay |   M≥6 olay |   M≥6.5 olay |   ortak ayar: M≥5.5 kazanç |   ortak ayar: M≥6 kazanç |   ortak ayar: M≥6.5 kazanç |   bölgesel en iyi: M≥5.5 kazanç |   bölgesel en iyi: M≥6 kazanç |   bölgesel en iyi: M≥6.5 kazanç |   havuzlamasız: M≥5.5 kazanç |   havuzlamasız: M≥6 kazanç |   havuzlamasız: M≥6.5 kazanç |
|:----------------------------------------------------|-------------:|-----------:|-------------:|---------------------------:|-------------------------:|---------------------------:|--------------------------------:|------------------------------:|--------------------------------:|-----------------------------:|---------------------------:|-----------------------------:|
| Marmara (KAF batı)                                  |            1 |          1 |            0 |                      0.084 |                    0.075 |                      –     |                           0.072 |                         0.062 |                           –     |                       -1.235 |                     -1.271 |                        –     |
| Kuzey Anadolu (KAF orta)                            |            2 |          1 |            0 |                      0.214 |                   -0.069 |                      –     |                           0.355 |                         0.225 |                           –     |                       -0.554 |                     -0.268 |                        –     |
| Doğu Anadolu Fay Zonu                               |           16 |          8 |            4 |                      0.007 |                    0.012 |                     -0.062 |                          -0.062 |                        -0.063 |                          -0.141 |                       -0.474 |                     -0.547 |                       -0.053 |
| Ege / Batı Anadolu                                  |            7 |          3 |            1 |                     -0.294 |                   -0.114 |                     -0.201 |                          -0.240 |                        -0.133 |                          -0.389 |                       -0.285 |                     -0.409 |                       -0.924 |
| Doğu Anadolu platosu (Van-Kars)                     |            5 |          1 |            0 |                     -0.093 |                   -0.046 |                      –     |                          -0.172 |                        -0.120 |                           –     |                       -0.830 |                     -0.187 |                        –     |
| diğer (İç Anadolu, Akdeniz, Karadeniz, Kıbrıs yayı) |            9 |          1 |            0 |                     -0.314 |                   -0.501 |                      –     |                          -0.319 |                        -0.430 |                           –     |                       -0.264 |                      0.306 |                        –     |

### Aramanın seçtiği bölgesel ayar

- Ortak ayar (önceki aramadan): açık gruplar: Coulomb gerilme aktarımı, kapasite ve kayaç, kümelenme, moment boşalımı, son depremden geçen süre, zincirleme kırılma
- Bölgeden bölgeye değişmesine izin verilen gruplar: Coulomb gerilme aktarımı, kümelenme
- Ortak ayara çekme gücü τ = 0.071 (küçük = bölgeler ortak ayara sıkı bağlı; 1 ≈ ortak katsayı kadar sapabilir), bölge sabiti τ0 = 0.056
- Bölgeye özel serbestlik çarpanı: Marmara 0.45, Kuzey Anadolu 2.42, Doğu Anadolu Fay Zonu 3.63, Ege / Batı Anadolu 0.79, Doğu Anadolu platosu 0.65, diğer 0.22

Son sınama yılında ortak katsayıdan en çok sapan bölge katsayıları (standartlaştırılmış; +0,1 ≈ 1 sd artışta oran ×1,105):

| parametre @ bölge                              |   sapma |
|:-----------------------------------------------|--------:|
| dcff_max_mu0 @ Kuzey Anadolu (KAF orta)        |   0.145 |
| dcff_mean_mu0 @ Kuzey Anadolu (KAF orta)       |   0.141 |
| dcff_frac_pos_mu0 @ Doğu Anadolu Fay Zonu      |  -0.125 |
| dcff_max_mu0 @ Doğu Anadolu Fay Zonu           |   0.082 |
| dcff_mean_mu0 @ Doğu Anadolu Fay Zonu          |   0.077 |
| dcff_frac_pos_mu0 @ Kuzey Anadolu (KAF orta)   |   0.072 |
| dcff_mean_mu0 @ Ege / Batı Anadolu             |  -0.037 |
| dcff_max_mu0 @ Doğu Anadolu platosu (Van-Kars) |  -0.024 |
| dcff_frac_pos_mu0 @ Ege / Batı Anadolu         |  -0.023 |
| dcff_max_mu0 @ Ege / Batı Anadolu              |  -0.023 |
| (bölge sabiti) @ Ege / Batı Anadolu            |   0.016 |
| (bölge sabiti) @ Kuzey Anadolu (KAF orta)      |  -0.013 |

## Yeniden üretme

```bash
python scripts/train.py --config configs/train_search.toml --until features
python scripts/search_regional.py --config configs/train_search.toml
python scripts/train.py --config configs/train_search_long.toml --until features
python scripts/search_regional.py --config configs/train_search_long.toml --trials 40 --null-runs 3
python scripts/regional_report.py
```
