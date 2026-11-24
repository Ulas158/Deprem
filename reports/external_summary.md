# Gaz, sıvı ve dış etkiler: ilk sınama sonuçları

> **Uyarı.** Bu rapordaki hiçbir sayı bir deprem tahmini değildir. Hepsi geçmiş veride sınanan hipotezlerdir. Küçük p değerleri tek başına nedensellik göstermez; çok sayıda bölge ve faz denendiği için yanlış pozitif beklenir (Bonferroni: 0,05 / test sayısı).

Katalog: 2008–2026-11-23, Mc = 2.2, 182,806 olay, artçılar ayıklandıktan sonra (Gardner-Knopoff) 43,481 bağımsız olay.

Bilimsel arka plan, mekanizmalar ve kanıt düzeyleri: [docs/akiskan_ve_dis_etkiler.md](../docs/akiskan_ve_dis_etkiler.md).

## 1. Gelgit ve mevsim (Schuster testi)

Her olay için şu fazlar hesaplandı: **saat** (yerel güneş saati; fiziksel değil, katalog kalitesi kontrolü), **ay_12sa** (yalnızca Ay'ın yarı günlük gelgiti; asıl gelgit testi), **gelgit_birlesik** (Ay + Güneş gelgit eğrisinde faz, Tanaka 2002), **gelgit_14g** (dolunay/yeni ay 0°) ve **mevsim** (1 Ocak 0°). Fazlar rastgeleyse p büyüktür.

### Önemli bulgu: sahte gelgit sinyali

Küçük depremlerde (M≥2.2) günün saatine çok güçlü bir yığılma var (log10 p = -309); M≥3'te bu kayboluyor (log10 p = -0.4). Bu, gece/gündüz algılama farkından ve/veya katalogda kalan taş ocağı patlamalarından kaynaklanır, depremlerin fiziğinden değil. Güneş gelgiti tam 12 saatlik olduğundan aynı etki Ay+Güneş birleşik gelgit fazında da görünüyor (M≥2.2: log10 p = -26). Güneş'ten bağımsız olan **Ay gelgiti testinde ise sinyal yok** (M≥2.2: -0.06, M≥3: -0.38, M≥4: -1.20). Sonuç: Türkiye kataloğunda gelgit etkisine dair kanıt bulunamadı; birleşik gelgit fazı küçük depremlerde kullanılmamalı. Modelde gelgit özelliği Ay fazından, saat fazı ise ayrı bir kalite kontrol özelliği olarak hesaplanıyor.

Fiziksel testler (ay_12sa, gelgit_14g, mevsim) toplam 54 kez yapıldığı için anlamlılık eşiği log10 p < -3.0 alınmalıdır (Bonferroni).

| bölge                           |   M≥ |     N |   saat_log10p |   saat_faz |   ay_12sa_log10p |   ay_12sa_faz |   gelgit_birlesik_log10p |   gelgit_birlesik_faz |   gelgit_14g_log10p |   gelgit_14g_faz |   mevsim_log10p |   mevsim_faz |
|:--------------------------------|-----:|------:|--------------:|-----------:|-----------------:|--------------:|-------------------------:|----------------------:|--------------------:|-----------------:|----------------:|-------------:|
| Türkiye                         |  2.2 | 43481 |       -308.58 |       -140 |            -0.06 |            39 |                   -25.56 |                    53 |               -2.76 |             -131 |           -4.47 |          148 |
| Türkiye                         |  3   | 10088 |         -0.38 |        -42 |            -0.38 |            37 |                    -3.24 |                    33 |               -0.09 |             -128 |           -7.28 |           75 |
| Türkiye                         |  4   |  1034 |         -0.12 |       -121 |            -1.2  |           135 |                    -0.63 |                    63 |               -0.06 |              -53 |           -0.29 |           90 |
| Marmara (KAF batı)              |  2.2 |  2687 |        -37.42 |       -130 |            -0.54 |            64 |                    -2.25 |                    60 |               -0.75 |               45 |           -1.68 |           65 |
| Marmara (KAF batı)              |  3   |   467 |         -0.32 |       -136 |            -0.03 |            55 |                    -0.01 |                    81 |               -0.07 |               94 |           -0.69 |           39 |
| Marmara (KAF batı)              |  4   |    69 |         -0.24 |       -132 |            -0.04 |          -169 |                    -0.34 |                  -123 |               -0.28 |              -76 |           -0.55 |          153 |
| KAF orta (Bolu-Erzincan)        |  2.2 |  3901 |        -44.84 |       -139 |            -0.07 |             2 |                    -4.93 |                    52 |               -2.32 |             -143 |           -0.67 |          168 |
| KAF orta (Bolu-Erzincan)        |  3   |   748 |         -0.2  |       -119 |            -0.47 |           127 |                    -0.86 |                    94 |               -1.22 |             -132 |           -0.74 |          119 |
| KAF orta (Bolu-Erzincan)        |  4   |    61 |         -0.08 |         36 |            -0.64 |           113 |                    -0.7  |                    76 |               -0.49 |               -6 |           -0.39 |           73 |
| Doğu Anadolu Fay Zonu           |  2.2 |  5775 |       -122.12 |       -143 |            -2.06 |            14 |                    -9.23 |                    50 |               -0.36 |             -163 |           -1.5  |          117 |
| Doğu Anadolu Fay Zonu           |  3   |  1333 |         -2.3  |       -127 |            -0.34 |             6 |                    -1.63 |                    55 |               -0.29 |              149 |           -4.96 |           58 |
| Doğu Anadolu Fay Zonu           |  4   |   126 |         -1.15 |        150 |            -0.01 |          -162 |                    -0.33 |                    34 |               -0.17 |             -130 |           -0.33 |           49 |
| Batı Anadolu / Ege graben       |  2.2 | 10063 |        -29.82 |       -162 |            -0.03 |          -165 |                    -3.47 |                    26 |               -1.15 |              -88 |           -1.61 |          172 |
| Batı Anadolu / Ege graben       |  3   |  2234 |         -2.97 |         54 |            -0.03 |             0 |                    -0.73 |                   -10 |               -0.57 |              -34 |           -0.93 |           48 |
| Batı Anadolu / Ege graben       |  4   |   210 |         -0.14 |         47 |            -2    |           134 |                    -1.21 |                   124 |               -0.31 |               33 |           -0.31 |           34 |
| Doğu Anadolu platosu (Van-Kars) |  2.2 |  4727 |        -10.44 |       -117 |            -0.86 |           159 |                    -1.58 |                    95 |               -0.13 |              -71 |           -0.1  |           87 |
| Doğu Anadolu platosu (Van-Kars) |  3   |  1146 |         -1.11 |        -40 |            -0.35 |           119 |                    -0.53 |                   100 |               -0.07 |             -165 |           -0.91 |          124 |
| Doğu Anadolu platosu (Van-Kars) |  4   |   131 |         -1.46 |        -95 |            -0.72 |            94 |                    -0.49 |                    78 |               -0.64 |             -104 |           -0.79 |          170 |

Fiziksel testlerde eşiği geçen satır sayısı: **3** (Türkiye M≥2.2, Türkiye M≥3, Doğu Anadolu Fay Zonu M≥3). Bunlar doğrulanmış değildir; sonraki dönem verisinde tekrar etmesi ve katalog tamlığının mevsimsel değişimiyle açıklanamaması gerekir.

### Mevsim sinyali yıl yıl

Türkiye M≥3 bağımsız olaylarda mevsim testi log10 p = -7.3; 2023 çıkarılınca -6.1. Ancak sinyal yıllara eşit dağılmıyor: 19 yılın yalnızca 3 tanesinde log10 p < −2. Bu, gerçek bir mevsimsel yükleme etkisinden çok birkaç yıldaki büyük deprem dizilerinin artçı ayıklamadan kaçan kısmından ve kataloğun yıllar içinde değişen büyüklük ölçeğinden (M≥3 olay sayısı 2008'de ~1270, sonraki yıllarda ~300) kaynaklanıyor olabilir. Şimdilik **kanıt sayılmaz**; daha homojen bir katalog ve ETAS tabanlı ayıklamayla yeniden sınanmalı.

|   yıl |    N |   mevsim_log10p |   faz |
|------:|-----:|----------------:|------:|
|  2008 | 1267 |            -1.9 |    78 |
|  2009 | 1260 |            -0.1 |   105 |
|  2010 | 1197 |            -4.6 |    67 |
|  2011 |  865 |            -4.1 |    90 |
|  2012 |  614 |            -2.2 |    36 |
|  2013 |  519 |            -0.3 |   -42 |
|  2014 |  442 |            -1.1 |  -168 |
|  2015 |  398 |            -0.2 |   115 |
|  2016 |  372 |            -0.5 |    12 |
|  2017 |  349 |            -1.7 |    74 |
|  2018 |  270 |            -0.2 |    79 |
|  2019 |  284 |            -0.7 |  -107 |
|  2020 |  295 |            -0.3 |   154 |
|  2021 |  270 |            -0.5 |   -70 |
|  2022 |  281 |            -0.2 |    16 |
|  2023 |  421 |            -1.8 |    75 |
|  2024 |  331 |            -0.2 |    45 |
|  2025 |  348 |            -0.2 |    25 |
|  2026 |  305 |            -1.7 |   153 |

## 2. Deprem sürüleri ve akışkan göçü

Uzay-zaman bağlantısıyla (10 km, 5 gün) en az 5 olaylı **2,960** küme bulundu; bunların **1,726** tanesi sürü ölçütünü sağlıyor (en büyük olay başta değil ve ikinci en büyükten en fazla 0,5 büyük). Sürülerde olayların zamanla yayılmasından hidrolik yayınım D hesaplandı (r² = 4πDt). Akışkan kaynaklı sürülerde D tipik olarak 0,01–10 m²/s'dir.

* D hesaplanabilen sürü: 496; ortanca D = 36 m²/s (%25–%75: 16–92). Konum hataları (birkaç km) küçük sürülerde D'yi büyütür; değerler üst sınır gibi okunmalı.

En büyük 15 sürü:

| t_start    |    n |   duration_days |    lat |    lon |   depth_km |   mag_max |   D_m2_s |   radius_km |
|:-----------|-----:|----------------:|-------:|-------:|-----------:|----------:|---------:|------------:|
| 2025-01-27 | 3287 |         146.819 | 36.625 | 25.682 |      9.76  |       5.3 |   61.452 |      28.072 |
| 2021-04-02 | 1555 |         199.077 | 36.402 | 27.094 |      8.015 |       5.7 |    9.299 |      22.624 |
| 2017-01-12 | 1547 |         138.384 | 39.534 | 26.115 |      7.44  |       5.4 |   12.422 |      18.759 |
| 2012-04-10 |  963 |         105.31  | 39.1   | 29.14  |      6.319 |       5.3 |   10.137 |      12.027 |
| 2008-03-23 |  809 |         173.815 | 36.987 | 29.202 |      7.407 |       4.3 |    3.601 |      13.457 |
| 2012-09-07 |  633 |         204.489 | 38.748 | 43.339 |      8.028 |       4.5 |   34.176 |      38.73  |
| 2010-02-27 |  573 |         192.821 | 37.228 | 28.197 |      7.328 |       3.2 |    5.214 |      14.679 |
| 2017-03-28 |  434 |         106.359 | 38.692 | 27.703 |      7.099 |       5   |   11.755 |      28.125 |
| 2012-11-08 |  419 |          49.454 | 36.604 | 27.953 |      8.589 |       4.7 |   15.932 |      11.774 |
| 2009-02-24 |  374 |         205.745 | 39.103 | 27.565 |      7.881 |       3.6 |    1.793 |      10.579 |
| 2009-04-16 |  322 |         208.279 | 37.22  | 28.2   |      8.279 |       3.4 |   13.116 |      19.183 |
| 2025-04-19 |  304 |          83.614 | 39.23  | 28.982 |      9.209 |       4.6 |    8.662 |       7.793 |
| 2020-02-16 |  300 |          41.281 | 38.422 | 44.481 |      6.584 |       6   |   31.66  |      17.527 |
| 2009-04-25 |  261 |         145.946 | 39.674 | 29.446 |      7.052 |       3.1 |    1.083 |       5.89  |
| 2008-01-01 |  243 |          56.594 | 39.375 | 33.148 |      6.009 |       4   |   37.843 |      14.922 |

## 3. Uzak büyük depremlerin dinamik tetiklemesi

864 küresel M≥6,5 depremden, Türkiye'ye 1000 km'den uzak ve genlik vekili yeterli olanlar (159 olay) için, her birinden sonraki 5 gündeki yerel bağımsız (artçısı ayıklanmış) olay sayısı önceki 30 günle karşılaştırıldı (β istatistiği; β > 2 kabaca anlamlı artış, β < −2 azalma).

* Türkiye geneli: β = **1.26** (5 gün içinde 5486 olay, toplam 37800)
* β > 2 olan hücre sayısı (≥20 olaylı hücreler): 10 / 417 (çoklu test nedeniyle birkaç tanesi tesadüfen beklenir)

## 4. İnsan kaynaklı sahalar

78 saha (baraj, maden, jeotermal, petrol/gaz, gaz deposu, doğal gaz çıkışı) başlangıç/bitiş yılı ve durumuyla `external_data/insan_kaynakli_sahalar.csv` dosyasında. Konumların çoğu yaklaşık ve doğrulanmamış; özellik tablosuna yalnızca o yıl faal olan sahalar girer.

| tür               |   sayı |
|:------------------|-------:|
| baraj             |     30 |
| jeotermal         |     13 |
| maden_komur       |     12 |
| maden_metal       |      8 |
| petrol            |      4 |
| maden_endustriyel |      3 |
| dogalgaz          |      3 |
| gaz_depolama      |      2 |
| dogal_gaz_cikisi  |      2 |
| maden_cozelti     |      1 |

## 5. Veri kaynakları ve erişim

| veri                                                 | kaynak                                                 | erisim   | durum                   | not                                                              |
|:-----------------------------------------------------|:-------------------------------------------------------|:---------|:------------------------|:-----------------------------------------------------------------|
| Birleşik deprem kataloğu                             | USGS, EMSC, AFAD (ISC)                                 | acik     | bağlı                   | Sürü, gelgit, mevsim, dinamik tetikleme özellikleri              |
| USGS küresel M≥6,5 kataloğu                          | earthquake.usgs.gov FDSN                               | acik     | bağlı                   | Uzak tetikleme için                                              |
| WRI Küresel Santral Veritabanı                       | github.com/wri/global-power-plant-database (CC BY 4.0) | acik     | bağlı                   | Baraj ve jeotermal santral konumları                             |
| Maden, petrol/gaz, jeotermal saha listesi            | Bu depoda elle derlendi (external_data/)               | elle     | bağlı (yaklaşık)        | Konumlar ve yıllar doğrulanmalı; MAPEG/MTA ile karşılaştırılmalı |
| Küresel Isı Akısı Veritabanı (IHFC 2024)             | ihfc-iugg.org / GFZ Data Services                      | acik     | bu ortamdan erişilemedi | Ağ izni verilirse otomatik indirilebilir                         |
| MTA jeotermal envanteri ve sıcak su kaynakları       | mta.gov.tr yayınları                                   | elle     | bağlı değil             | Tablo olarak indirme sunulmuyor                                  |
| AFAD radon / öncül gözlem istasyonları               | AFAD                                                   | izin     | bağlı değil             | Zaman serileri açık değil; resmi başvuru gerekir                 |
| KOERI / üniversite radon ve gaz ölçümleri            | Yayınlar ve araştırmacılar                             | izin     | bağlı değil             | Yayınlarda grafik olarak var, ham veri yazarlardan istenmeli     |
| DSİ yeraltı suyu gözlem kuyuları ve baraj seviyeleri | DSİ                                                    | izin     | bağlı değil             | Kurumsal başvuru gerekir                                         |
| GRACE/GRACE-FO hidrolojik yük                        | NASA PO.DAAC                                           | kayit    | bağlı değil             | Ücretsiz NASA Earthdata hesabı                                   |
| ERA5 / ERA5-Land (yağış, kar, basınç)                | Copernicus CDS                                         | kayit    | bağlı değil             | Ücretsiz Copernicus hesabı                                       |
| Uydu altimetrisi göl/rezervuar seviyeleri            | DAHITI / Hydroweb                                      | kayit    | bağlı değil             | Ücretsiz hesap                                                   |
| GNSS zaman serileri (yavaş kayma)                    | Nevada Geodetic Laboratory                             | acik     | bu ortamdan erişilemedi | Ağ izni verilirse indirilebilir                                  |
| Akışkan enjeksiyon/çekim hacimleri                   | Şirketler, EPDK, MAPEG                                 | izin     | bağlı değil             | Genellikle gizli; resmi başvuru                                  |
| OpenStreetMap madenler, taş ocakları, barajlar       | Overpass API / Geofabrik                               | acik     | bu ortamdan erişilemedi | Ağ izni verilirse ~binlerce saha otomatik eklenir                |

## 6. Fay elemanlarına aktarım ve eksik değer doldurma

546,992 fay elemanının (1 km × 1 km) her biri için gaz/sıvı ve dış etki parametreleri hesaplandı. Eksik değerler `deprem.faults.fill.fill_missing` ile önce aynı fay boyunca sağ/sol komşulardan, sonra komşu faylardan, sonra bölgesel ortancadan, hiç ölçüm yoksa literatür varsayılanından dolduruldu; her değerin kaynağı `<ad>__fill`, belirsizliği `<ad>__sd` sütunundadır. Doldurulan değerler eğitimde öncül olarak kullanılıp düzeltilecektir. Hiç ölçümü olmayan ve makul bir literatür değeri bulunmayan büyüklükler (ör. radon) uydurulmaz, 'veri yok' kalır. Oranlar:

| parametre                    |   gerçek veri |   fay boyunca komşulardan |   komşu faylardan |   bölgesel ortanca |   literatür varsayılanı |   veri yok |
|:-----------------------------|--------------:|--------------------------:|------------------:|-------------------:|------------------------:|-----------:|
| radon_bq_m3                  |         0     |                     0     |             0     |              0     |                       0 |          1 |
| radon_anomaly_z              |         0     |                     0     |             0     |              0     |                       0 |          1 |
| co2_flux_g_m2_d              |         0     |                     0     |             0     |              0     |                       0 |          1 |
| he3_he4_rc_ra                |         0     |                     0     |             0     |              0     |                       0 |          1 |
| ch4_h2_ppm                   |         0     |                     0     |             0     |              0     |                       0 |          1 |
| groundwater_level_m          |         0     |                     0     |             0     |              0     |                       0 |          1 |
| groundwater_anomaly_z        |         0     |                     0     |             0     |              0     |                       0 |          1 |
| spring_temp_c                |         0     |                     0     |             0     |              0     |                       0 |          1 |
| spring_chem_anomaly_z        |         0     |                     0     |             0     |              0     |                       0 |          1 |
| thermal_spring_dist_km       |         0     |                     0     |             0     |              0     |                       0 |          1 |
| pore_pressure_lambda         |         0     |                     0     |             0     |              0     |                       1 |          0 |
| swarm_n_10y                  |         1     |                     0     |             0     |              0     |                       0 |          0 |
| swarm_frac_10y               |         0.355 |                     0.315 |             0.237 |              0.094 |                       0 |          0 |
| hydraulic_diffusivity_log10  |         0.103 |                     0.119 |             0.434 |              0.345 |                       0 |          0 |
| tidal_semidiurnal_log10p     |         0.539 |                     0.238 |             0.162 |              0.061 |                       0 |          0 |
| hour_of_day_log10p           |         0.539 |                     0.238 |             0.162 |              0.061 |                       0 |          0 |
| tidal_fortnight_log10p       |         0.539 |                     0.238 |             0.162 |              0.061 |                       0 |          0 |
| seasonal_log10p              |         0.539 |                     0.238 |             0.162 |              0.061 |                       0 |          0 |
| dyn_trigger_beta             |         0.803 |                     0.136 |             0.048 |              0.014 |                       0 |          0 |
| hydro_load_amp_kpa           |         0     |                     0     |             0     |              0     |                       0 |          1 |
| snow_load_amp_kpa            |         0     |                     0     |             0     |              0     |                       0 |          1 |
| precip_anomaly_z             |         0     |                     0     |             0     |              0     |                       0 |          1 |
| atm_pressure_amp_hpa         |         0     |                     0     |             0     |              0     |                       0 |          1 |
| lake_level_change_m          |         0     |                     0     |             0     |              0     |                       0 |          1 |
| slow_slip_rate_mm_yr         |         0     |                     0     |             0     |              0     |                       0 |          1 |
| tec_anomaly_z                |         0     |                     0     |             0     |              0     |                       0 |          1 |
| geomagnetic_anomaly_z        |         0     |                     0     |             0     |              0     |                       0 |          1 |
| thermal_ir_anomaly_z         |         0     |                     0     |             0     |              0     |                       0 |          1 |
| site_reservoir_dist_km       |         1     |                     0     |             0     |              0     |                       0 |          0 |
| site_mining_dist_km          |         1     |                     0     |             0     |              0     |                       0 |          0 |
| site_geothermal_dist_km      |         1     |                     0     |             0     |              0     |                       0 |          0 |
| site_hydrocarbon_dist_km     |         1     |                     0     |             0     |              0     |                       0 |          0 |
| site_injection_dist_km       |         1     |                     0     |             0     |              0     |                       0 |          0 |
| site_gas_seep_dist_km        |         1     |                     0     |             0     |              0     |                       0 |          0 |
| injection_volume_m3_yr       |         0     |                     0     |             0     |              0     |                       0 |          1 |
| extraction_volume_m3_yr      |         0     |                     0     |             0     |              0     |                       0 |          1 |
| reservoir_level_rate_m_month |         0     |                     0     |             0     |              0     |                       0 |          1 |
| heat_flow_mw_m2              |         0     |                     0     |             0     |              0     |                       1 |          0 |
| vp_vs_ratio                  |         0     |                     0     |             0     |              0     |                       1 |          0 |
| curie_depth_km               |         0     |                     0     |             0     |              0     |                       1 |          0 |
| crustal_thickness_km         |         0     |                     0     |             0     |              0     |                       1 |          0 |

