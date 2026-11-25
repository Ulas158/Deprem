# En iyi parametre kombinasyonu araması

> Deneysel araştırma çıktısıdır, deprem tahmini değildir.

İki arama yapıldı, ikisinde de 2020-2025 aramada hiç kullanılmayan son sınama dönemidir:

- **Uzun arama (1900'den, M≥5,5):** katalog, fay, Coulomb, kayma açığı, kapasite/doluluk, zincirleme kırılma, kayaç sertliği. Arama 1920-2019.
- **Kısa arama (2012'den, M≥3):** yukarıdakilere ek olarak AFAD küçük deprem parametreleri (b-değeri, sürüler, sürtünme vekili), gaz/sıvı ve dış etkiler, GPS gerinim ve kilitlenme/krip. Arama 2015-2019.

## Özet

- **Maksimize eden tek bir kombinasyon bulunamadı.** Aramanın en iyi bulduğu ayarlar arama yıllarında ETAS'ı geçiyor, ama dokunulmamış 2020-2025 yıllarında bu üstünlük ya küçülüyor (kısa arama: +0,06, anlamlı değil) ya da tersine dönüyor (uzun arama: -0,09). Yani arama, geçmişe fazla uyan (ezberleyen) ayarlar seçiyor.
- **Enerji öncülü** (enerji birikimi gelgitten önemli varsayımı) sonucu değiştirmedi: kısa aramada fark 0,001, uzun aramada kapatınca sınama puanı biraz düştü (-0,087 → -0,105), ikisi de anlamsız.
- **Kayaç sertliği/kapasite ve zincirleme kırılma** tek başına eklendiğinde sınamada küçük ama tutarlı artı verdi (+0,01 ile +0,05; tek başına anlamlılık testi yapılmadı, keşif amaçlı); **doluluk oranı ve doluluk + Coulomb** ise sıfır ya da eksi.
- **%100 dolunca kırılır mı?** Hayır, sabit bir eşik yok: aynı fay bir döngüde %100'ün çok altında, başka döngüde çok üstünde kırılıyor. Kırılan yerler ortalamadan biraz daha dolu (bölgedeki %65'lik dilim), o kadar.
- **Zincirleme model** 1999 İzmit-Düzce ve 2023 Kahramanmaraş kırılma uzunluklarını yeniden üretemedi.
- **GPS'ten krip** açık istasyonlar seyrek olduğu için çözülemedi (İsmetpaşa: GPS 0, literatür ~8 mm/yıl).

## Uzun arama: 1900'den bugüne, M≥5,5

### Kısa sonuç

- 300 kombinasyon, yalnız 1920-2019 yıllarında denendi. En iyisinin arama puanı +0.120 (olay başına, ETAS'a göre; 0 = ETAS kadar).
- Bilgisi silinmiş verilerde tekrar arama bu koşuda tamamlanmadı (bir arama 1 saatten uzun sürdü); sınama sonucu zaten negatif olduğu için sonucu değiştirmez.
- Dokunulmamış sınama yıllarında (2020-2025) en iyi kombinasyonun puanı -0.087; %95 güven aralığı [-0.233, -0.036], çoklu test düzeltmeli p = 1.00. ETAS'ı **anlamlı biçimde geçmedi**.

### Nasıl yapıldı

1. **Arama dönemi 1920-2019.** Her yıl için model yalnız o yıldan önceki yıllarla eğitildi ve o yılı tahmin etti. Puan: M≥5.5, M≥6, M≥6.5 depremlerin olay başına bilgi kazancının ortalaması (ETAS'a göre; + iyi). Katsayılar 5 yılda bir, yalnız o bloktan önceki yıllarla yeniden eğitildi.
2. **Aranan şeyler:** hangi parametre gruplarının açık olacağı, her grubun önem ağırlığı (0,05-3; küçük = katkısı kısılmış), genel düzenlileştirme gücü, Coulomb sürtünmesi μ', veriyle önem ağırlığı (ARD, uyarlamalı ridge), fiziksel enerji öncülü (enerji/kayma birikimi grupları büyük, gelgit gibi küçük etkiler küçük başlangıç), katalog belirsizlik ağırlıkları ve eğitim penceresi (tüm geçmiş / son 3 / son 5 yıl). Arama Bayesçi optimizasyonla (TPE) yapıldı.
3. **Şans sınırı:** bu koşuda tamamlanmadı; yerine dokunulmamış sınama yılları ve çoklu test düzeltmesi kullanıldı.
4. **Sınama 2020-2025:** aramada hiç kullanılmadı. İlk 5 farklı kombinasyon ve karşılaştırma modelleri burada yıl yıl ileri yürüyüşle bir kez ölçüldü; güven aralıkları yıllar üzerinden bootstrap, p değerleri Holm yöntemiyle çoklu teste göre düzeltildi.
5. **ETAS referansı değişmedi** (her yıl yalnız o yıl öncesiyle yeniden uyduruldu, önceki PR'lardaki gibi).

### Sınama yılları sonuçları

| model                                           |   arama puanı |   sınama puanı | birleşik %95 GA   |   birleşik p (tek yönlü) |   Holm düzeltilmiş p |   M≥5.5 kazanç | M≥5.5 geçtiği yıl   |   M≥6 kazanç | M≥6 geçtiği yıl   |   M≥6.5 kazanç | M≥6.5 %95 GA     | M≥6.5 geçtiği yıl   |
|:------------------------------------------------|--------------:|---------------:|:------------------|-------------------------:|---------------------:|---------------:|:--------------------|-------------:|:------------------|---------------:|:-----------------|:--------------------|
| aramada 1. sıradaki                             |         0.12  |         -0.087 | [-0.233, -0.036]  |                    0.997 |                    1 |         -0.119 | 1/6                 |       -0.052 | 1/5               |         -0.089 | [-0.105, -0.078] | 0/2                 |
| aramada 2. sıradaki                             |         0.117 |         -0.129 | [-0.251, -0.080]  |                    0.998 |                    1 |         -0.152 | 1/6                 |       -0.087 | 0/5               |         -0.147 | [-0.172, -0.109] | 0/2                 |
| aramada 3. sıradaki                             |         0.115 |         -0.067 | [-0.287, -0.023]  |                    0.994 |                    1 |         -0.132 | 1/6                 |       -0.036 | 1/5               |         -0.033 | [-0.083, -0.000] | 0/2                 |
| aramada 4. sıradaki                             |         0.115 |         -0.08  | [-0.293, -0.031]  |                    0.995 |                    1 |         -0.142 | 1/6                 |       -0.036 | 1/5               |         -0.061 | [-0.122, -0.021] | 0/2                 |
| aramada 5. sıradaki                             |         0.114 |         -0.077 | [-0.303, -0.034]  |                    0.997 |                    1 |         -0.145 | 1/6                 |       -0.042 | 1/5               |         -0.045 | [-0.075, -0.025] | 0/2                 |
| yalnız ETAS'ın kendi düzeltmesi (özelliksiz)    |        -0.159 |         -0.041 | [-0.055, 0.021]   |                    0.911 |                    1 |         -0.036 | 3/6                 |       -0.025 | 2/5               |         -0.061 | [-0.064, -0.057] | 0/2                 |
| önceki model: tüm parametreler, varsayılan ayar |        -0.151 |         -0.055 | [-0.377, 0.106]   |                    0.603 |                    1 |         -0.025 | 2/6                 |        0.079 | 2/5               |         -0.218 | [-0.441, -0.069] | 0/2                 |
| en iyi + enerji öncülü kapalı                   |         0.095 |         -0.105 | [-0.168, -0.039]  |                    0.996 |                    1 |         -0.098 | 1/6                 |       -0.074 | 0/5               |         -0.144 | [-0.220, -0.093] | 0/2                 |
| en iyi + veriyle önem ağırlığı (ARD) açık       |         0.093 |         -0.115 | [-0.184, -0.052]  |                    0.996 |                    1 |         -0.119 | 1/6                 |       -0.068 | 0/5               |         -0.157 | [-0.160, -0.153] | 0/2                 |

Puanlar olay başına doğal log birimindedir: +0,10 ≈ gerçekleşen depremlere ETAS'tan %10 daha yüksek olasılık vermek demektir.

### Aramanın seçtiği kombinasyon

- Açık gruplar (6): Coulomb gerilme aktarımı, kapasite ve kayaç, kümelenme, moment boşalımı, son depremden geçen süre, zincirleme kırılma
- Önem ağırlıkları: Coulomb gerilme aktarımı 2.48, moment boşalımı 1.59, son depremden geçen süre 0.26, kapasite ve kayaç 0.11, zincirleme kırılma 0.09, kümelenme 0.06
- Genel öncül sd 0.320, μ' 0.0, ARD kapalı, enerji öncülü açık, belirsizlik ağırlıkları kapalı, eğitim penceresi tüm geçmiş

Son sınama yılındaki en büyük 15 katsayı (standartlaştırılmış; +0,1 ≈ 1 sd artışta oran ×1,105):

| parametre                |   katsayı |
|:-------------------------|----------:|
| log_m0_30y               |    -0.199 |
| max_mag_all__eksik       |     0.185 |
| log_m0_5y__var           |    -0.164 |
| log_m0_30y__var          |     0.151 |
| cap_log_m0               |     0.128 |
| cascade_mw               |     0.087 |
| max_mag_all              |    -0.067 |
| dcff_frac_pos_mu0        |    -0.066 |
| log_m0_5y                |     0.049 |
| dcff_mean_mu0__eksik     |     0.041 |
| dcff_max_mu0__eksik      |     0.041 |
| dcff_frac_pos_mu0__eksik |     0.041 |
| rigidity_gpa             |     0.04  |
| cascade_len_km           |     0.03  |
| dcff_mean_mu0            |    -0.023 |

En iyi %5 denemede her grubun açık olma oranı (yüksek = arama bu grubu sevdi; yalnız arama dönemi):

| grup                          |   en iyi %5 denemede açık olma oranı |
|:------------------------------|-------------------------------------:|
| Coulomb gerilme aktarımı      |                                 1    |
| moment boşalımı               |                                 1    |
| kapasite ve kayaç             |                                 1    |
| zincirleme kırılma            |                                 1    |
| kümelenme                     |                                 0.87 |
| derinlik                      |                                 0.73 |
| sürtünme vekili               |                                 0.4  |
| son depremden geçen süre      |                                 0.27 |
| fay geometrisi                |                                 0.07 |
| doluluk + Coulomb             |                                 0    |
| GPS kilitlenme / krip         |                                 0    |
| hız değişimi                  |                                 0    |
| doluluk oranı                 |                                 0    |
| kayma açığı (enerji birikimi) |                                 0    |
| sayılar                       |                                 0    |

### Doluluk nasıl özetlenmeli? (ortalama, moment ağırlıklı, en dolu parça, zincirleme)

Her satır, kapasite grupları dışındaki tüm parametrelere yalnız o özetin eklenmesinin katkısıdır (+ iyi). Sınama sütunu aramada kullanılmayan yıllardır.

| özet                                        |   arama dönemi katkısı |   sınama dönemi katkısı |   sınama M≥5.5 katkı |   sınama M≥6 katkı |   sınama M≥6.5 katkı |
|:--------------------------------------------|-----------------------:|------------------------:|---------------------:|-------------------:|---------------------:|
| doluluk: basit ortalama                     |                 -0.006 |                  -0.037 |               -0.025 |             -0.02  |               -0.066 |
| doluluk: moment ağırlıklı ortalama          |                 -0.006 |                  -0.04  |               -0.026 |             -0.022 |               -0.073 |
| doluluk: en dolu parça                      |                 -0.006 |                  -0.002 |                0.001 |             -0.003 |               -0.004 |
| doluluk + Coulomb aktarımı                  |                 -0.01  |                  -0.044 |               -0.021 |             -0.026 |               -0.085 |
| zincirleme kırılma (olası uzunluk/büyüklük) |                 -0.004 |                   0.022 |               -0.003 |              0.02  |                0.048 |
| kapasite + kayaç sertliği + kabuk           |                  0     |                   0.048 |                0.024 |              0.049 |                0.072 |
| GPS kilitlenme / krip                       |                  0     |                   0.041 |                0.025 |              0.024 |                0.075 |
| sürtünme vekili                             |                  0.001 |                   0.001 |                0     |              0.001 |                0.001 |

### Deprem deprem (zaman + konum), sınama yılları

| hedef   |   olay |   zaman kazancı | zaman %95 GA     |   konum kazancı | konum %95 GA    |   toplam kazanç | toplam %95 GA    | en olası %5 hücrede (ETAS)   | en olası %5 hücrede (model)   |
|:--------|-------:|----------------:|:-----------------|----------------:|:----------------|----------------:|:-----------------|:-----------------------------|:------------------------------|
| M≥5.5   |     40 |          -0.041 | [-0.059, -0.004] |          -0.041 | [-0.168, 0.035] |          -0.082 | [-0.198, -0.009] | %38                          | %32                           |
| M≥6     |     15 |          -0.01  | [-0.051, 0.073]  |           0.028 | [-0.112, 0.093] |           0.017 | [-0.107, 0.110]  | %13                          | %13                           |
| M≥6.5   |      5 |          -0.058 | [-0.096, -0.038] |           0.036 | [-0.096, 0.197] |          -0.022 | [-0.177, 0.102]  | %20                          | %20                           |

Büyüklük için model ETAS ile aynı Gutenberg-Richter dağılımını kullanır; büyüklükte fark yoktur.

### Katalog belirsizliği

Katalog her depremin zaman, konum ve büyüklük hatasıyla 6 kez yeniden örneklendi; katalogdan hesaplanan parametreler ve etiketler yeniden üretildi, en iyi kombinasyon sınama yıllarında yeniden ölçüldü.

| ölçü              |   2.5% |    50% |   97.5% |
|:------------------|-------:|-------:|--------:|
| M≥5.5             | -0.222 | -0.169 |  -0.139 |
| M≥6               | -0.26  | -0.183 |  -0.161 |
| M≥6.5             | -0.406 | -0.297 |  -0.239 |
| score             | -0.276 | -0.21  |  -0.202 |
| olay M≥6 toplam   | -0.186 | -0.082 |   0.011 |
| olay M≥6.5 toplam | -0.296 | -0.198 |  -0.14  |

## Kısa arama: 2012'den bugüne, M≥3 (tüm parametreler)

### Kısa sonuç

- 400 kombinasyon, yalnız 2015-2019 yıllarında denendi. En iyisinin arama puanı +0.204 (olay başına, ETAS'a göre; 0 = ETAS kadar).
- Aynı arama, parametrelerin yerleri rastgele karıştırılmış (bilgisi silinmiş) verilerde 4 kez yapıldı; oradaki en iyi puanlar: +0.112, +0.137, +0.164, +0.123. Gerçek arama bunların hepsini geçti, yani arama döneminde gerçek bir sinyal var gibi görünüyor.
- Dokunulmamış sınama yıllarında (2020-2025) en iyi kombinasyonun puanı +0.061; %95 güven aralığı [-0.014, 0.215], çoklu test düzeltmeli p = 0.65. ETAS'ı **anlamlı biçimde geçmedi**.

### Nasıl yapıldı

1. **Arama dönemi 2015-2019.** Her yıl için model yalnız o yıldan önceki yıllarla eğitildi ve o yılı tahmin etti. Puan: M≥3, M≥4, M≥5 depremlerin olay başına bilgi kazancının ortalaması (ETAS'a göre; + iyi).
2. **Aranan şeyler:** hangi parametre gruplarının açık olacağı, her grubun önem ağırlığı (0,05-3; küçük = katkısı kısılmış), genel düzenlileştirme gücü, Coulomb sürtünmesi μ', veriyle önem ağırlığı (ARD, uyarlamalı ridge), fiziksel enerji öncülü (enerji/kayma birikimi grupları büyük, gelgit gibi küçük etkiler küçük başlangıç), katalog belirsizlik ağırlıkları ve eğitim penceresi (tüm geçmiş / son 3 / son 5 yıl). Arama Bayesçi optimizasyonla (TPE) yapıldı.
3. **Şans sınırı:** aynı arama, parametrelerin hücrelere rastgele karıştırıldığı verilerde tekrarlandı.
4. **Sınama 2020-2025:** aramada hiç kullanılmadı. İlk 5 farklı kombinasyon ve karşılaştırma modelleri burada yıl yıl ileri yürüyüşle bir kez ölçüldü; güven aralıkları yıllar üzerinden bootstrap, p değerleri Holm yöntemiyle çoklu teste göre düzeltildi.
5. **ETAS referansı değişmedi** (her yıl yalnız o yıl öncesiyle yeniden uyduruldu, önceki PR'lardaki gibi).

### Sınama yılları sonuçları

| model                                           |   arama puanı |   sınama puanı | birleşik %95 GA   |   birleşik p (tek yönlü) |   Holm düzeltilmiş p |   M≥3 kazanç | M≥3 geçtiği yıl   |   M≥4 kazanç | M≥4 geçtiği yıl   |   M≥5 kazanç | M≥5 %95 GA      | M≥5 geçtiği yıl   |
|:------------------------------------------------|--------------:|---------------:|:------------------|-------------------------:|---------------------:|-------------:|:------------------|-------------:|:------------------|-------------:|:----------------|:------------------|
| aramada 1. sıradaki                             |         0.204 |          0.061 | [-0.014, 0.215]   |                    0.086 |                0.654 |        0.091 | 5/6               |        0.045 | 3/6               |        0.048 | [-0.044, 0.186] | 3/6               |
| aramada 2. sıradaki                             |         0.204 |          0.074 | [-0.008, 0.241]   |                    0.066 |                0.598 |        0.103 | 4/6               |        0.06  | 3/6               |        0.059 | [-0.042, 0.219] | 3/6               |
| aramada 3. sıradaki                             |         0.204 |          0.021 | [-0.052, 0.152]   |                    0.197 |                0.751 |        0.054 | 5/6               |        0.001 | 4/6               |        0.007 | [-0.082, 0.128] | 4/6               |
| aramada 4. sıradaki                             |         0.203 |          0.048 | [-0.032, 0.203]   |                    0.15  |                0.751 |        0.079 | 4/6               |        0.03  | 3/6               |        0.035 | [-0.058, 0.182] | 3/6               |
| aramada 5. sıradaki                             |         0.202 |          0.017 | [-0.055, 0.150]   |                    0.229 |                0.751 |        0.05  | 4/6               |       -0.002 | 3/6               |        0.003 | [-0.083, 0.120] | 4/6               |
| yalnız ETAS'ın kendi düzeltmesi (özelliksiz)    |         0.022 |         -0.014 | [-0.026, 0.027]   |                    0.765 |                0.765 |       -0.01  | 3/6               |       -0.023 | 2/6               |       -0.011 | [-0.029, 0.045] | 3/6               |
| önceki model: tüm parametreler, varsayılan ayar |         0.128 |          0.048 | [-0.056, 0.232]   |                    0.152 |                0.751 |        0.085 | 5/6               |        0.026 | 4/6               |        0.032 | [-0.089, 0.235] | 3/6               |
| en iyi + enerji öncülü açık                     |         0.203 |          0.062 | [-0.014, 0.215]   |                    0.082 |                0.654 |        0.092 | 5/6               |        0.046 | 3/6               |        0.049 | [-0.040, 0.191] | 3/6               |
| en iyi + veriyle önem ağırlığı (ARD) kapalı     |         0.199 |          0.055 | [-0.023, 0.208]   |                    0.114 |                0.681 |        0.086 | 5/6               |        0.037 | 3/6               |        0.043 | [-0.045, 0.187] | 4/6               |

Puanlar olay başına doğal log birimindedir: +0,10 ≈ gerçekleşen depremlere ETAS'tan %10 daha yüksek olasılık vermek demektir.

### Aramanın seçtiği kombinasyon

- Açık gruplar (13): b-değeri, fay geometrisi, gelgit ve periyodik yükler, hız değişimi, insan kaynaklı sahalar, kayma açığı (enerji birikimi), kümelenme, moment boşalımı, sayılar, son depremden geçen süre, sürtünme vekili, sürüler ve akışkan, zincirleme kırılma
- Önem ağırlıkları: sürüler ve akışkan 2.78, b-değeri 1.27, son depremden geçen süre 0.72, hız değişimi 0.65, zincirleme kırılma 0.64, gelgit ve periyodik yükler 0.49, kümelenme 0.34, kayma açığı (enerji birikimi) 0.31, fay geometrisi 0.20, sayılar 0.15, insan kaynaklı sahalar 0.13, moment boşalımı 0.11, sürtünme vekili 0.08
- Genel öncül sd 1.222, μ' 0.8, ARD açık, enerji öncülü kapalı, belirsizlik ağırlıkları kapalı, eğitim penceresi tüm geçmiş

Son sınama yılındaki en büyük 15 katsayı (standartlaştırılmış; +0,1 ≈ 1 sd artışta oran ×1,105):

| parametre                       |   katsayı |
|:--------------------------------|----------:|
| rate_z_1y__eksik                |     0.922 |
| fault_log_m0_rate__eksik        |    -0.4   |
| log_m0_30y__eksik               |    -0.375 |
| frac_dependent_10y__eksik       |    -0.361 |
| yrs_since_m4__eksik             |     0.237 |
| fric_velocity_strengthening     |     0.199 |
| yrs_since_m5__var               |    -0.172 |
| reservoir_years_since_fill__var |     0.166 |
| max_mag_all                     |     0.164 |
| cascade_mw__var                 |     0.147 |
| yrs_since_m6__var               |    -0.145 |
| log_m0_5y                       |    -0.145 |
| n_mc_1y                         |    -0.144 |
| rate_z_1y                       |     0.126 |
| cascade_len_km__eksik           |     0.124 |

En iyi %5 denemede her grubun açık olma oranı (yüksek = arama bu grubu sevdi; yalnız arama dönemi):

| grup                          |   en iyi %5 denemede açık olma oranı |
|:------------------------------|-------------------------------------:|
| b-değeri                      |                                 1    |
| hız değişimi                  |                                 1    |
| insan kaynaklı sahalar        |                                 1    |
| fay geometrisi                |                                 1    |
| sürüler ve akışkan            |                                 1    |
| zincirleme kırılma            |                                 1    |
| sayılar                       |                                 1    |
| gelgit ve periyodik yükler    |                                 1    |
| son depremden geçen süre      |                                 1    |
| sürtünme vekili               |                                 1    |
| kümelenme                     |                                 1    |
| kayma açığı (enerji birikimi) |                                 0.65 |
| moment boşalımı               |                                 0.5  |
| Coulomb gerilme aktarımı      |                                 0.25 |
| kapasite ve kayaç             |                                 0.1  |
| GPS kilitlenme / krip         |                                 0.05 |
| GPS                           |                                 0    |
| doluluk + Coulomb             |                                 0    |
| doluluk oranı                 |                                 0    |
| derinlik                      |                                 0    |

### Doluluk nasıl özetlenmeli? (ortalama, moment ağırlıklı, en dolu parça, zincirleme)

Her satır, kapasite grupları dışındaki tüm parametrelere yalnız o özetin eklenmesinin katkısıdır (+ iyi). Sınama sütunu aramada kullanılmayan yıllardır.

| özet                                        |   arama dönemi katkısı |   sınama dönemi katkısı |   sınama M≥3 katkı |   sınama M≥4 katkı |   sınama M≥5 katkı |
|:--------------------------------------------|-----------------------:|------------------------:|-------------------:|-------------------:|-------------------:|
| doluluk: basit ortalama                     |                 -0.001 |                   0.003 |              0.002 |              0.003 |              0.002 |
| doluluk: moment ağırlıklı ortalama          |                  0     |                   0.003 |              0.002 |              0.003 |              0.003 |
| doluluk: en dolu parça                      |                 -0.013 |                  -0.003 |             -0.002 |             -0.004 |             -0.003 |
| doluluk + Coulomb aktarımı                  |                 -0.011 |                  -0.011 |             -0.009 |             -0.013 |             -0.013 |
| zincirleme kırılma (olası uzunluk/büyüklük) |                 -0.003 |                   0.011 |              0.011 |              0.012 |              0.01  |
| kapasite + kayaç sertliği + kabuk           |                 -0.003 |                   0.019 |              0.016 |              0.021 |              0.018 |
| GPS kilitlenme / krip                       |                 -0.001 |                  -0.002 |             -0     |             -0.001 |             -0.005 |
| sürtünme vekili                             |                  0.015 |                   0.027 |              0.024 |              0.039 |              0.019 |

### Deprem deprem (zaman + konum), sınama yılları

| hedef   |   olay |   zaman kazancı | zaman %95 GA    |   konum kazancı | konum %95 GA     |   toplam kazanç | toplam %95 GA   | en olası %5 hücrede (ETAS)   | en olası %5 hücrede (model)   |
|:--------|-------:|----------------:|:----------------|----------------:|:-----------------|----------------:|:----------------|:-----------------------------|:------------------------------|
| M≥3.5   |   4872 |           0.027 | [-0.001, 0.066] |          -0.027 | [-0.042, -0.010] |          -0     | [-0.040, 0.053] | %79                          | %79                           |
| M≥4     |   1896 |           0.024 | [-0.020, 0.079] |          -0.038 | [-0.057, -0.015] |          -0.014 | [-0.068, 0.057] | %78                          | %78                           |
| M≥5     |    130 |          -0.015 | [-0.087, 0.025] |          -0.023 | [-0.067, 0.015]  |          -0.038 | [-0.129, 0.020] | %77                          | %78                           |

Büyüklük için model ETAS ile aynı Gutenberg-Richter dağılımını kullanır; büyüklükte fark yoktur.

### Katalog belirsizliği

Katalog her depremin zaman, konum ve büyüklük hatasıyla 8 kez yeniden örneklendi; katalogdan hesaplanan parametreler ve etiketler yeniden üretildi, en iyi kombinasyon sınama yıllarında yeniden ölçüldü.

| ölçü            |   2.5% |    50% |   97.5% |
|:----------------|-------:|-------:|--------:|
| M≥3             |  0.033 |  0.053 |   0.088 |
| M≥4             | -0.037 | -0.004 |   0.037 |
| M≥5             | -0.001 |  0.016 |   0.055 |
| M≥6             | -0.1   | -0.019 |   0.034 |
| score           |  0.001 |  0.02  |   0.059 |
| olay M≥4 toplam | -0.05  | -0.037 |  -0.015 |
| olay M≥5 toplam | -0.048 | -0.026 |   0.014 |

## Fay doluluk oranı (uzun koşu: birikim 1700'den, tarihsel katalogla)

### "%100 dolunca deprem olur mu?" testi

Her 1 km parçanın enerji biriktirme kapasitesi = kilitli alan × kayaç sertliği (CRUST1.0) × taşınabilir en büyük gerilme düşümü. Doluluk = biriken kayma açığı (GPS kilitlenmesiyle, depremlerde boşalarak) / kapasite. Ardından kapasite, geçmiş büyük depremlerin kırıldığı andaki doluluğa göre fay → bölge → genel hiyerarşisiyle kalibre edildi (o fayın tipik kırılma birikimi = %100).

- Bütçe 1700'de başladı; 1900 sonrası, kayma hızı ≥2 mm/yıl faylarda 21 büyük deprem (M≥6,8).
- Kırılma anındaki doluluk olaydan olaya çok değişiyor: log-ölçekte standart sapma 1.67, yani tipik değerin yaklaşık 5.3 katı aşağı-yukarı. Sabit bir eşik ("%100'de kırılır") yok.
- Kırılan yerin doluluğu, aynı bölgedeki tüm fay parçalarına göre ortalama %65'lik dilimde (eşik doğru olsaydı %100'e yakın olurdu; rastgele olsaydı %50). Wilcoxon p = 0.024.
- Tüm dönem (1702-2023, 42 deprem, tarihsel konum/büyüklük hataları büyük): kırılma anı kalibre doluluğun ortancası %77, %10-90 aralığı %4-%394; aynı bölgedeki yüzdelik ortalaması %68.
- Aynı fayın birden fazla kırılması (döngüden döngüye değişkenlik):

| fay        |   time |   mw |   doluluk_kalibre_% |
|:-----------|-------:|-----:|--------------------:|
| ME_TRCS006 |   1942 |  7   |               541.9 |
| ME_TRCS006 |   1943 |  7.5 |               540.6 |
| ME_TRCS006 |   1951 |  7   |               109.4 |
| ME_TRCS010 |   1957 |  7.1 |               426.6 |
| ME_TRCS012 |   1967 |  7.3 |               200.9 |
| ME_TRCS012 |   1999 |  7.6 |               226   |
| ME_TRCS010 |   1999 |  7.2 |               402.6 |

| time       |   mw | bölge        | fay        |   doluluk_kalibre_% |   aynı_bölgedeki_yüzdelik |
|:-----------|-----:|:-------------|:-----------|--------------------:|--------------------------:|
| 1905-12-04 | 6.8  | dogu_anadolu | ME_TRCS133 |                1    |                      0.09 |
| 1912-08-09 | 7.23 | marmara      | ME_TRCS996 |              316.78 |                      0.76 |
| 1914-10-03 | 7.06 | bati_anadolu | ME_TRCS309 |               97.25 |                      0.87 |
| 1916-01-24 | 7    | kaf_orta     | ME_TRCS003 |              567.05 |                      0.89 |
| 1926-06-26 | 6.97 | bati_anadolu | ME_GRCS965 |               14.62 |                      0.3  |
| 1939-12-26 | 7.8  | diğer        | ME_TRCS147 |               42.37 |                      0.79 |
| 1942-12-20 | 7    | kaf_orta     | ME_TRCS006 |              541.94 |                      0.77 |
| 1943-11-26 | 7.5  | kaf_orta     | ME_TRCS006 |              540.63 |                      0.77 |
| 1944-02-01 | 7.6  | kaf_orta     | ME_TRCS378 |                3.16 |                      0.15 |
| 1949-07-23 | 7.16 | bati_anadolu | ME_GRCS880 |               50.8  |                      0.68 |
| 1951-08-13 | 6.96 | diğer        | ME_TRCS006 |              109.38 |                      0.94 |
| 1953-03-18 | 7.3  | marmara      | ME_TRCS031 |               68.91 |                      0.33 |
| 1957-05-26 | 7.09 | kaf_orta     | ME_TRCS010 |              426.6  |                      0.6  |
| 1964-10-06 | 6.8  | marmara      | ME_TRCS023 |               97.1  |                      0.4  |
| 1967-07-22 | 7.29 | kaf_orta     | ME_TRCS012 |              200.94 |                      0.48 |
| 1969-03-28 | 6.8  | bati_anadolu | ME_TRCS321 |               91.79 |                      0.84 |
| 1976-11-24 | 7.3  | dogu         | ME_IR43s1  |               72.45 |                      0.97 |
| 1999-08-17 | 7.6  | marmara      | ME_TRCS012 |              225.99 |                      0.63 |
| 1999-11-12 | 7.2  | kaf_orta     | ME_TRCS010 |              402.62 |                      0.59 |
| 2020-01-24 | 6.8  | dogu_anadolu | ME_TRCS051 |              259.48 |                      0.92 |
| 2023-02-06 | 7.8  | dogu_anadolu | ME_TRCS056 |              174.85 |                      0.88 |

### Zincirleme kırılma bilinen depremleri yeniden üretebiliyor mu?

Depremden hemen önceki doluluktan, episantra en yakın parçadan başlatılan kırılma; kırılan parça yükünün yarısını komşularına aktarır, eşiği aşan komşu da kırılır:

| deprem              |   gözlenen_uzunluk_km |   başlangıç_doluluğu_% |   model_uzunluk_km_eşik1 |   model_uzunluk_km_eşik0.8 |   model_uzunluk_km_eşik0.5 |
|:--------------------|----------------------:|-----------------------:|-------------------------:|---------------------------:|---------------------------:|
| 1999 İzmit Mw7.6    |                   145 |                  244.6 |                    278.1 |                      278.1 |                      278.1 |
| 1999 Düzce Mw7.2    |                    40 |                  407.8 |                    153.9 |                      167.8 |                      179.8 |
| 2023 Pazarcık Mw7.8 |                   350 |                   56.1 |                      1   |                        2   |                       40.3 |
| 2023 Elbistan Mw7.5 |                   150 |                   10.4 |                      1   |                        1   |                        1   |

### GPS'ten krip (yavaş kayma) ve literatür

| yer                          |   literatür_krip_mm_yıl |   gps_krip_mm_yıl |   gps_krip_sd | not                                              | kaynak                                                         |
|:-----------------------------|------------------------:|------------------:|--------------:|:-------------------------------------------------|:---------------------------------------------------------------|
| KAF İsmetpaşa                |                       8 |               0   |           0.4 | nan                                              | Çakır vd. 2005; Kaneko vd. 2013 (InSAR, kripmetre ~6-9 mm/yıl) |
| KAF Destek                   |                       6 |               1.7 |           2.8 | nan                                              | Kaneko vd. 2013 (InSAR, ~5-8 mm/yıl)                           |
| Marmara İzmit (1999 sonrası) |                       5 |               0.1 |           1.7 | nan                                              | Aslan vd. 2019 (deprem sonrası sığ krip, azalan)               |
| DAF Hazar-Palu               |                       4 |             nan   |         nan   | fayın iki yanında yeterli açık GPS istasyonu yok | Cavalié & Jónsson 2014 (sığ kısmi krip)                        |

Açık GPS istasyonları seyrek olduğu için kilitlenme ancak birkaç fayda çözülebildi; çoğunda sonuç öncüle yakın kaldı. Sığ krip (İsmetpaşa gibi) birkaç km genişliğinde olduğundan 20-50 km aralıklı GPS ile ayrılamıyor; InSAR verisi bağlanırsa iyileşir.

### Şu an en dolu faylar (kayma hızı ≥2 mm/yıl)

| fault_id   | bölge        |   lat |   lon |   kayma_hızı_mm_yıl |   doluluk_% |   doluluk_%16 |   doluluk_%84 |
|:-----------|:-------------|------:|------:|--------------------:|------------:|--------------:|--------------:|
| ME_TRCS013 | marmara      |  40.8 |  28.7 |                21.4 |       770   |         664.9 |         880.2 |
| ME_TRCS003 | kaf_orta     |  40.2 |  38.1 |                22.4 |       723   |         613.4 |         798.6 |
| ME_TRCS006 | kaf_orta     |  41   |  34.9 |                25   |       637.1 |         547.1 |         738.4 |
| ME_TRCS996 | diğer        |  40.4 |  26.5 |                21.9 |       592   |         528.2 |         680.3 |
| ME_TRCS396 | bati_anadolu |  37   |  27.8 |                 5.5 |       314.1 |         275.2 |         371.9 |
| ME_TRCS322 | bati_anadolu |  38.4 |  28.1 |                 5.7 |       311.9 |         264.5 |         346.5 |
| ME_TRCS001 | diğer        |  39.4 |  40.8 |                13.5 |       273.8 |         241.6 |         316.3 |
| ME_TRCS049 | dogu_anadolu |  38.6 |  39.8 |                10   |       271.7 |         232.9 |         314.7 |
| ME_TRCS269 | bati_anadolu |  37.5 |  29.2 |                 7.9 |       247.1 |         217.4 |         286.6 |
| ME_TRCS077 | dogu         |  40   |  42.7 |                 2.5 |       246.7 |         217.2 |         282.2 |
| ME_GRCS699 | diğer        |  35   |  32.4 |                 6   |       245   |         213.7 |         281.8 |
| ME_TRCS010 | kaf_orta     |  40.7 |  31.2 |                19.3 |       242.3 |         213.5 |         278.2 |
| ME_TRCS033 | diğer        |  39.4 |  26.2 |                 3.9 |       225.4 |         196.5 |         259.1 |
| ME_TRCS331 | bati_anadolu |  37.9 |  28.2 |                 4.3 |       220.2 |         194.5 |         253.4 |
| ME_TRCS008 | kaf_orta     |  40.6 |  31.3 |                12   |       219.8 |         193.9 |         254.1 |

Bu yüzdeler göreli ölçüdür (bütçe 1700'de sıfırdan başlar) ve yukarıdaki test sabit bir kırılma eşiği göstermediği için tek başına "sıradaki deprem" anlamına gelmez. Haritada: *Faylar ve enerji → Fay doluluk oranı*.

## Sınırlar

- Sınama yalnız 6 yıl; büyük depremler (M≥6) bu sürede az, güven aralıkları geniş.
- Isı akısı, kayaç türü ve sürtünme ölçümleri bu ortamdan erişilemedi; sürtünme mikrodeprem/b-değeri/tekrarlayan depremlerden türetilen bir vekildir.
- İyonosfer gibi kısa vadeli sinyaller yıllık modelde değil, ayrı kısa vadeli testtedir (train_obs).
