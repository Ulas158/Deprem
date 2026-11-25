# Fay fay biriken enerji, küçük depremler ve komşu fay zinciri

> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Hiçbir gösterge kanıtsız doğrulanmış sayılmaz.

## 1. Küçük depremler fayı boşaltıyor mu?

1700'den 2025 sonuna kadar fay segmentlerine düşen toplam deprem momentinin yalnız **%0.60**'i M<5,5 depremlerden geliyor. Kayıtlı bütün küçük depremler (432,386 deprem) toplamda, Türkiye faylarının yalnız ≈1.3 yılda biriktirdiği yükü boşaltmış. Yani küçük depremler fayı boşaltmıyor, yük büyük depremlere kalıyor. Küçük depremler enerjiyi değil, gerilmenin nasıl değiştiğine dair ipucunu (b-değeri) taşıyabilir. (Küçük depremlerin kaydı 1998'den önce eksik; eksik olanlar eklense de bu oran birkaç kattan fazla değişmez, çünkü her büyüklük basamağı bir öncekinin ~30 katı enerji taşır.)

## 2. Fay fay biriken yük (2026 başı)

Doluluk R = biriken moment açığı / fayın kaldırabileceği (10 MPa gerilme düşümüyle). Saat 1700'de sıfırdan başlıyor; 1700'den beri kırılmamış faylarda gerçek değer daha yüksek olabilir. R kalibre edilmemiş, göreli bir ölçüdür: R = 1, fayın taşıyabileceği en büyük gerilme düşümüne denk yük demektir, fayın o an kırılacağı anlamına gelmez. Kayma hızı ≥ 2 mm/yıl, uzunluk ≥ 10 km segmentler, en dolu 20. Tamamı: `reports/tables/fault_energy_2026.csv`.

| segment    | sistem                  |   uzunluk_km |   kayma_hızı_mm_yıl |   doluluk_R |   açık_Mw_eşdeğeri |   enerji_PJ |   yüklü_zincir_km |   zincirin_olası_Mw |   b_farkı_z |
|:-----------|:------------------------|-------------:|--------------------:|------------:|-------------------:|------------:|------------------:|--------------------:|------------:|
| ME_TRCS003 | Kuzey Anadolu Fayı      |          301 |                22.4 |        1.9  |               7.88 |        41.1 |               301 |                 8.3 |       -3.04 |
| ME_TRCS006 | Kuzey Anadolu Fayı      |          301 |                25   |        1.86 |               8    |        63.2 |               301 |                 8.3 |        0.03 |
| ME_TRCS013 | Kuzey Anadolu Fayı      |          104 |                21.4 |        1.79 |               7.64 |        17.9 |               379 |                 8.5 |        1.2  |
| ME_TRCS996 | Kuzey Anadolu Fayı      |          275 |                21.9 |        1.36 |               7.92 |        47.4 |               379 |                 8.5 |       -1.55 |
| ME_TRCS396 | diğer: Batı Anadolu-Ege |           42 |                 5.5 |        1.25 |               6.26 |         0.2 |               124 |                 7.7 |       -3.3  |
| ME_TRCS001 | Kuzey Anadolu Fayı      |           52 |                13.5 |        0.98 |               7.22 |         4.2 |                52 |                 7   |       -0.49 |
| ME_TRCS393 | Ege graben sistemi      |           31 |                 9.8 |        0.95 |               7.02 |         2.1 |               124 |                 7.7 |       -2.81 |
| ME_TRCS010 | Kuzey Anadolu Fayı      |           47 |                19.3 |        0.91 |               7.29 |         5.4 |                47 |                 7   |       -2.33 |
| ME_TRCS391 | Ege graben sistemi      |           51 |                 9   |        0.78 |               7.03 |         2.2 |               124 |                 7.7 |       -3.95 |
| ME_TRCS269 | Ege graben sistemi      |           27 |                 7.9 |        0.76 |               6.94 |         1.6 |                64 |                 7.2 |       -1.67 |
| ME_TRCS273 | Ege graben sistemi      |           25 |                 7.4 |        0.72 |               6.9  |         1.4 |                25 |                 6.5 |       -3.06 |
| ME_TRCS274 | Ege graben sistemi      |           17 |                 7.3 |        0.71 |               6.79 |         1   |                17 |                 6.2 |       -3.23 |
| ME_TRCS268 | Ege graben sistemi      |           37 |                 7.4 |        0.7  |               7    |         2   |                64 |                 7.2 |       -3.6  |
| ME_TRCS049 | Doğu Anadolu Fayı       |           45 |                10   |        0.68 |               7.2  |         4   |                85 |                 7.4 |       -2.74 |
| ME_GRCS699 | diğer: Orta Anadolu     |          178 |                 6   |        0.65 |               7.38 |         7.4 |               178 |                 7.9 |        0    |
| ME_TRCS271 | Ege graben sistemi      |           40 |                 7.4 |        0.62 |               7.08 |         2.6 |                40 |                 6.9 |       -1.43 |
| ME_TRCS275 | Ege graben sistemi      |           29 |                 5.5 |        0.54 |               6.86 |         1.2 |                29 |                 6.6 |       -0.88 |
| ME_TRCS051 | Doğu Anadolu Fayı       |           41 |                10   |        0.52 |               7.09 |         2.7 |                85 |                 7.4 |       -3.32 |
| ME_TRCS046 | diğer: Güneydoğu        |           46 |                 8.1 |        0.51 |               7.06 |         2.5 |                46 |                 7   |       -0.91 |
| ME_TRCS008 | Kuzey Anadolu Fayı      |           32 |                12   |        0.47 |               6.94 |         1.6 |                 0 |               nan   |       -2.21 |

enerji_PJ: biriken elastik enerji (petajoule; 1 PJ ≈ 240 kiloton TNT), E ≈ 3 MPa/(2μ) × moment açığı. açık_Mw_eşdeğeri: biriken açığın tamamı tek depremde boşalsa büyüklüğü.

## 3. Büyük depremler (M≥6,5, 1900 sonrası): kırılmadan hemen önce fay ne durumdaydı?

Yüklü = R ≥ 0.5. Zincir: depremin düştüğü segmente uç uca bağlı (≤ 5 km) yüklü segmentlerin toplam boyu. Gerçek kırık boyu büyüklükten (Wells & Coppersmith 1994) kaba tahmin.

| tarih      |   Mw | segment     | sistem                  |   doluluk_R |   segment_boyu_km |   segmentin_tek_başına_Mw |   yüklü_zincir_km |   zincirin_olası_Mw |   gerçek_kırık_km (tahmini) |   b_farkı_z |   1km_doluluk |   1km_yüklü_boy_km |
|:-----------|-----:|:------------|:------------------------|------------:|------------------:|--------------------------:|------------------:|--------------------:|----------------------------:|------------:|--------------:|-------------------:|
| 1905-12-04 | 6.8  | ME_TRCS133  | diğer: Güneydoğu        |        0.04 |               186 |                       8   |                 0 |                     |                          37 |        0    |          0.27 |                  0 |
| 1910-06-25 | 6.51 | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 1912-08-09 | 7.23 | ME_TRCS996  | Kuzey Anadolu Fayı      |        0.94 |               275 |                       8.3 |               379 |                 8.5 |                          67 |        0    |          0.75 |                359 |
| 1912-09-13 | 6.72 | ME_TRCS996  | Kuzey Anadolu Fayı      |        0.82 |               275 |                       8.3 |               379 |                 8.5 |                          33 |        0    |          0.25 |                  0 |
| 1914-10-03 | 7.06 | ME_TRCS309  | Ege graben sistemi      |        0.25 |                61 |                       7.2 |                 0 |                     |                          53 |        0    |          0.25 |                  0 |
| 1916-01-24 | 7    | ME_TRCS003  | Kuzey Anadolu Fayı      |        1.51 |               301 |                       8.3 |               301 |                 8.3 |                          49 |        0    |          1.52 |                298 |
| 1919-11-18 | 7    | ME_TRCS210  | Ege graben sistemi      |        0.25 |                30 |                       6.6 |                 0 |                     |                          49 |        0    |          0.09 |                  0 |
| 1924-09-13 | 6.86 | ME_TRCS074  | diğer: Doğu             |        0.02 |               100 |                       7.5 |                 0 |                     |                          40 |        0    |          0.01 |                  0 |
| 1926-03-18 | 6.75 | ME_TRCS282  | diğer: Batı Anadolu-Ege |        0.07 |                42 |                       6.9 |                 0 |                     |                          35 |        0    |          0.04 |                  0 |
| 1926-06-26 | 6.97 | ME_GRCS965  | diğer: Batı Anadolu-Ege |        0.18 |                45 |                       6.9 |                 0 |                     |                          47 |        0    |          0.2  |                  0 |
| 1928-04-14 | 7.1  | EUR_BGCS004 | diğer: Batı Anadolu-Ege |        0.01 |                56 |                       7.1 |                 0 |                     |                          56 |        0    |          0.01 |                  0 |
| 1930-05-06 | 7.1  | ME_IRCS067  | diğer: Doğu             |        0.05 |                57 |                       7.1 |                 0 |                     |                          56 |        0    |          0.08 |                  0 |
| 1938-04-19 | 6.62 | ME_TRCS247  | diğer: Orta Anadolu     |        0.07 |               160 |                       7.9 |                 0 |                     |                          29 |        0    |          0.07 |                  0 |
| 1939-12-26 | 7.8  | ME_TRCS003  | Kuzey Anadolu Fayı      |        1.58 |               301 |                       8.3 |               301 |                 8.3 |                         145 |        0    |          0.13 |                  0 |
| 1942-12-20 | 7    | ME_TRCS006  | Kuzey Anadolu Fayı      |        1.52 |               301 |                       8.3 |               301 |                 8.3 |                          49 |        0    |          1.52 |                292 |
| 1943-11-26 | 7.5  | ME_TRCS379  | Kuzey Anadolu Fayı      |        0.12 |                55 |                       7.1 |                 0 |                     |                          97 |        0    |          0.09 |                  0 |
| 1944-02-01 | 7.6  | ME_TRCS377  | diğer: Orta Anadolu     |        0.15 |                27 |                       6.6 |                 0 |                     |                         111 |        0    |          0.15 |                  0 |
| 1944-10-06 | 6.74 | ME_TRCS033  | Ege graben sistemi      |        0.17 |               132 |                       7.7 |                 0 |                     |                          34 |        0    |          0.23 |                  0 |
| 1948-02-09 | 7.3  | ME_GRCS780  | diğer: Batı Anadolu-Ege |        0.05 |                64 |                       7.2 |                 0 |                     |                          74 |        0    |          0.17 |                  0 |
| 1949-07-23 | 7.16 | ME_GRCS880  | Ege graben sistemi      |        0.22 |                56 |                       7.1 |                 0 |                     |                          61 |        0    |          0.11 |                  0 |
| 1949-08-17 | 6.78 | ME_TRCS115  | diğer: Kuzeydoğu        |        0    |                78 |                       7.3 |                 0 |                     |                          36 |        0    |          0.99 |                 44 |
| 1951-08-13 | 6.96 | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 1953-03-18 | 7.3  | ME_TRCS023  | diğer: Batı Anadolu-Ege |        0.33 |                97 |                       7.5 |                 0 |                     |                          74 |        0    |          0.2  |                  0 |
| 1955-07-16 | 6.57 | ME_TRCS337  | Ege graben sistemi      |        0.12 |                56 |                       7.1 |                 0 |                     |                          27 |        0    |          0.07 |                  0 |
| 1956-07-09 | 7.7  | EUR_GRCS958 | Ege graben sistemi      |        0.04 |                45 |                       6.9 |                 0 |                     |                         127 |        0    |          0.05 |                  0 |
| 1957-04-24 | 7.1  | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 1957-04-25 | 7.3  | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 1957-05-26 | 7.09 | ME_TRCS010  | Kuzey Anadolu Fayı      |        1.16 |                47 |                       7   |               216 |                 8.1 |                          55 |        0    |          1.16 |                 53 |
| 1964-10-06 | 6.8  | ME_TRCS023  | diğer: Batı Anadolu-Ege |        0.32 |                97 |                       7.5 |                 0 |                     |                          37 |        0    |          0.01 |                  0 |
| 1966-08-19 | 6.8  | ME_TRCS081  | diğer: Doğu             |        0.12 |                78 |                       7.3 |                 0 |                     |                          37 |        0    |          0.1  |                  0 |
| 1967-07-22 | 7.29 | ME_TRCS012  | Kuzey Anadolu Fayı      |        0.5  |               137 |                       7.8 |                 0 |                     |                          73 |        0    |          0.1  |                  0 |
| 1969-03-28 | 6.8  | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 1970-03-28 | 7.23 | ME_TRCS311  | Ege graben sistemi      |        0.09 |               173 |                       7.9 |                 0 |                     |                          67 |        0    |          0.01 |                  0 |
| 1971-05-22 | 6.58 | ME_TRCS046  | diğer: Güneydoğu        |        0.44 |                46 |                       7   |                 0 |                     |                          28 |        0    |          0.57 |                 18 |
| 1975-03-27 | 6.7  | ME_TRCS996  | Kuzey Anadolu Fayı      |        1.11 |               275 |                       8.3 |               379 |                 8.5 |                          33 |        0    |          0    |                  0 |
| 1975-09-06 | 6.7  | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 1976-11-24 | 7.3  | ME_IR149s1  | diğer: Doğu             |        0.08 |               101 |                       7.5 |                 0 |                     |                          74 |        0    |          0.14 |                  0 |
| 1983-10-30 | 6.8  | ME_TRCS113  | diğer: Doğu             |        0.09 |                78 |                       7.3 |                 0 |                     |                          37 |        0    |          0.03 |                  0 |
| 1988-12-07 | 6.8  | ME_AM2s2    | diğer: Doğu             |        0.06 |                77 |                       7.3 |                 0 |                     |                          37 |        0    |          0.06 |                  0 |
| 1991-04-29 | 7    | ME_G9s1s1   | diğer: Doğu             |        0.07 |               154 |                       7.8 |                 0 |                     |                          49 |        0    |          0.03 |                  0 |
| 1992-03-13 | 6.7  | ME_TRCS003  | Kuzey Anadolu Fayı      |        1.68 |               301 |                       8.3 |               301 |                 8.3 |                          33 |        0    |          1.11 |                298 |
| 1999-08-17 | 7.6  | ME_TRCS012  | Kuzey Anadolu Fayı      |        0.45 |               137 |                       7.8 |                 0 |                     |                         111 |        0    |          0.63 |                 80 |
| 1999-11-12 | 7.2  | ME_TRCS010  | Kuzey Anadolu Fayı      |        1.07 |                47 |                       7   |                79 |                 7.4 |                          64 |        0    |          1.15 |                 42 |
| 2002-02-03 | 6.5  | ME_TRCS204  | diğer: Orta Anadolu     |        0    |                33 |                       6.7 |                 0 |                     |                          25 |        0    |          0.01 |                  0 |
| 2011-10-23 | 7.1  | ME_TRCS423  | diğer: Doğu             |        0.02 |                27 |                       6.6 |                 0 |                     |                          56 |       -0.75 |          0.02 |                  0 |
| 2017-07-20 | 6.6  | ME_GRCS965  | diğer: Batı Anadolu-Ege |        0.07 |                45 |                       6.9 |                 0 |                     |                          28 |       -3.85 |          0.16 |                  0 |
| 2020-01-24 | 6.8  | ME_TRCS133  | diğer: Güneydoğu        |        0.07 |               186 |                       8   |                 0 |                     |                          37 |        1.01 |          0.54 |                 69 |
| 2020-10-30 | 7    | ME_GRCS915  | Ege graben sistemi      |        0.12 |                33 |                       6.7 |                 0 |                     |                          49 |       -2.29 |          0.04 |                  0 |
| 2023-02-06 | 7.8  | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |
| 2023-02-06 | 6.7  | ME_TRCS056  | Doğu Anadolu Fayı       |        0    |               110 |                       7.6 |                 0 |                     |                          33 |        1.5  |          0    |                  0 |
| 2023-02-06 | 7.5  | fay dışı    | nan                     |      nan    |               nan |                     nan   |               nan |               nan   |                         nan |      nan    |        nan    |                nan |

43 fay üstü büyük depremin 10'inde depremin başladığı segment yüklüydü. 1 km parçalarla bakıldığında (merkez üssüne en yakın parça, 25 km içinde) 43 depremin 10'inde parça yüklüydü. Büyük depremlerin çoğu, bu hesapla 'boş' görünen yavaş faylarda oldu: bu faylarda deprem tekrarı binlerce yıl, saatimiz ise 1700'de başlıyor; fayın 1700'deki yükünü bilmiyoruz.

## 4. "Biri sınırdaysa yeter mi, yoksa boydan boya yüklü olmalı mı?" sınaması

Fay üstündeki her M≥5,5 deprem için (1900-2025, ilk 60'ı yalnız eğitim, 297 deprem puanlandı, 34'i M≥6,5): yalnız önceki depremlerle eğitilmiş bir model "bu deprem M≥6,5 olacak mı" olasılığı veriyor. Karşılaştırma: büyüklüğü hiç bilmeyen sabit oran (Gutenberg-Richter). Pozitif = o bilgi büyük depremleri ayırt etmeye yardım ediyor.

| bilgi                                 |   olay başına kazanç | %95 GA           |   p |   M≥6,5 olaylarda |   Holm p |
|:--------------------------------------|---------------------:|:-----------------|----:|------------------:|---------:|
| segment boyu (yükten bağımsız)        |               -0     | [-0.004, +0.003] |   1 |            -0.035 |        1 |
| yalnız kendi segmenti (yüklüyse boyu) |               -0.002 | [-0.005, +0.002] |   1 |            -0.046 |        1 |
| yüklü komşu zinciri                   |               -0.002 | [-0.006, +0.002] |   1 |            -0.047 |        1 |
| zincir + kendi doluluğu               |               -0.003 | [-0.010, +0.007] |   1 |            -0.046 |        1 |
| 1 km parçalarla yüklü kesintisiz boy  |               -0.001 | [-0.005, +0.004] |   1 |            -0.038 |        1 |

## 5. Göstergeler deprem modeline eklenince (ne zaman ve nerede)

## Yöntem ve sınırlar

- Yöntem: `src/deprem/fault_energy.py`; modele ekleme: `src/deprem/faultwalk.py` (`partial_*` sürümleri).
- Kayma hızları GEM aktif fay veritabanından (jeolojik ve jeodezik derlemesi); GPS kilitlenme oranı bu sürümde kullanılmadı, krip yapan faylarda yük olduğundan fazla görünür.
- Deprem momenti, merkez üssünden kırık yarı boyu + 5 km içindeki fay izlerine dağıtıldı; tek yönlü kırılan büyük depremlerde (ör. 1939, 2023) boşalımın bir kısmı yanlış segmente düşer.
- 1700-1899 kataloğu eksik; o dönemde boşalım az sayılır, doluluk olduğundan yüksek görünür.
- Kapasite tek bir gerilme düşümüyle (10 MPa) hesaplandı; gerçek değer fay fay 1-30 MPa arasında değişebilir. Gerilme düşümünü fay fay ölçmek için dalga biçimi spektrumları gerekir; bu sürümde yok.
- Yüklülük eşiği (R ≥ 0,5) ve b penceresi (5 yıl, 30 km, ≥ 50 olay) sonuçlara bakılmadan önceden seçildi.
