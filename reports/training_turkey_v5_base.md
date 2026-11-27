# Uzun ileri eğitim: turkey_v5_base

> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Model ve varsayımlar doğrulanmamıştır.

## Ne yapıldı

Model, elimizdeki en eski kayıtlardan başlayarak bugüne kadar **yıl yıl ve deprem deprem** ilerledi. Her yılın başında yalnız o güne kadarki verilerle ETAS ve model yeniden eğitildi; o yıl içindeki her deprem gerçekleşmeden önce zamanı, yeri ve büyüklüğü için olasılık verildi, deprem gelince puanlandı ve geçmişe eklendi.

* 1700-1899 arası tarihsel depremler (EPICA, SHARE, GEM, NOAA) zaman çizelgesine, daha eskileri (MÖ 2000'den beri) yalnız "nerede deprem olur" haritasına düşük ağırlıkla girdi.
* Katalog her dönemde farklı büyüklükten itibaren eksiksiz olduğu için üç koşu var: M≥4 (1995-2025, eğitim 1990'ten).
* Eski yıllar eğitimde daha az ağırlık aldı (o yılların konum/büyüklük hataları büyük).
* **Güven aralıkları eski kayıtların belirsizliğini içerir:** katalog 10 kez yeniden örneklendi; her depremin tarihi, yeri ve büyüklüğü kendi hatası kadar rastgele kaydırıldı ve yürüyüş her seferinde baştan yapıldı. Aralıklar bu örnekler ve ay blokları üzerinden bootstrap ile hesaplandı, bu yüzden eski dönemlerde daha geniştir.
* GPS, iyonosfer gibi gözlemlerin olmadığı yıllarda bu özellikler "eksik" olarak işaretlendi.

## Sonuç: ETAS'a karşı olay başına bilgi kazancı

Pozitif = model ETAS'tan iyi. Aralık sıfırı içeriyorsa fark anlamlı değil. "Toplam kazanç" en iyi tahmin katalogundan, "katalog örnekleri ortalaması" hatalarla kaydırılmış kataloglardan gelir; kaydırma gerçek hatanın üstüne ek gürültü kattığından örnekler biraz kötümserdir, aralık bu duyarlılığı gösterir. "Yalnız örnekleme GA" katalog belirsizliği olmadan çıkan aralıktır; aradaki fark eski kayıtların belirsizliğinden gelir.

| eşik   | dönem            |   olay |   toplam kazanç |   katalog örnekleri ortalaması | toplam %95 GA    | yalnız örnekleme GA   | katalog belirsizliğinin payı   |   zaman kazanç | zaman %95 GA    |   konum kazanç | konum %95 GA     | en olası %5 hücrede (ETAS)   | en olası %5 hücrede (model)   | yıllık sayı aralığı tuttu (ETAS)   | yıllık sayı aralığı tuttu (model)   |
|:-------|:-----------------|-------:|----------------:|-------------------------------:|:-----------------|:----------------------|:-------------------------------|---------------:|:----------------|---------------:|:-----------------|:-----------------------------|:------------------------------|:-----------------------------------|:------------------------------------|
| M≥4    | 1995-2025 (tümü) |   4518 |          -0.02  |                         -0.008 | [-0.041, 0.016]  | [-0.062, 0.012]       | %0                             |         -0.008 | [-0.022, 0.024] |         -0.012 | [-0.026, 0.005]  | %64 [61-64]                  | %63 [61-63]                   | %80                                | %74                                 |
| M≥4    | 1995-1999        |    479 |          -0.248 |                         -0.155 | [-0.301, -0.036] | [-0.457, -0.068]      | %0                             |         -0.109 | [-0.186, 0.024] |         -0.139 | [-0.171, -0.023] | %57 [46-56]                  | %56 [45-55]                   | %80                                | %80                                 |
| M≥4    | 2000-2025        |   4039 |           0.007 |                          0.008 | [-0.013, 0.030]  | [-0.012, 0.031]       | %0                             |          0.004 | [-0.007, 0.031] |          0.003 | [-0.015, 0.011]  | %65 [62-64]                  | %64 [62-64]                   | %80                                | %72                                 |

## Gelecek yıl için model çıktısı (%95 aralık)

Bu bir deprem tahmini değildir; modelin Türkiye geneli için yıllık beklenen sayıları ve olasılıklarıdır. Aralıklar benzetim ve katalog belirsizliğini birlikte içerir.

| eşik   |   yıl | gerçekleşen   | gerçekleşen aralığı (katalog belirsizliği)   |   ETAS beklenen | ETAS %95   |   model beklenen | model %95   | en az bir M≥6 olasılığı (model)   | en az bir M≥6 olasılığı (ETAS)   | en az bir M≥7 olasılığı (model)   | en az bir M≥7 olasılığı (ETAS)   |
|:-------|------:|:--------------|:---------------------------------------------|----------------:|:-----------|-----------------:|:------------|:----------------------------------|:---------------------------------|:----------------------------------|:---------------------------------|
| M≥4    |  2026 |               |                                              |           189.9 | 114-353    |            196.5 | 116-365     | %71                               | %70                              | %10                               | %10                              |

## Yıl yıl beklenen ve gerçekleşen sayılar

### M≥4

|   yıl |   gerçekleşen | gerçekleşen aralığı (katalog belirsizliği)   |   ETAS beklenen | ETAS %95   |   model beklenen | model %95   |
|------:|--------------:|:---------------------------------------------|----------------:|:-----------|-----------------:|:------------|
|  1995 |            77 | 50-72                                        |            62.3 | 32-120     |             54.9 | 28-106      |
|  1996 |            73 | 48-69                                        |            64   | 35-125     |             57.9 | 31-113      |
|  1997 |            61 | 38-58                                        |            63.6 | 32-128     |             58   | 29-119      |
|  1998 |            74 | 74-90                                        |            63.2 | 31-129     |             60.5 | 29-122      |
|  1999 |           194 | 176-202                                      |            64.2 | 35-115     |             69.7 | 34-141      |
|  2000 |            83 | 82-107                                       |            94.6 | 48-192     |             97.8 | 50-191      |
|  2001 |            86 | 78-94                                        |            88.8 | 48-170     |             91.3 | 50-177      |
|  2002 |            65 | 65-81                                        |            83.6 | 47-148     |             88   | 50-159      |
|  2003 |           102 | 98-116                                       |            83.7 | 45-166     |             87.5 | 47-177      |
|  2004 |           106 | 104-124                                      |            86.2 | 49-153     |             94.8 | 54-174      |
|  2005 |           137 | 133-156                                      |            92.1 | 54-176     |            103.5 | 59-198      |
|  2006 |            60 | 58-66                                        |            96.2 | 53-168     |            109.3 | 60-197      |
|  2007 |           117 | 104-124                                      |            92.4 | 53-170     |            100.9 | 57-191      |
|  2008 |           100 | 86-100                                       |           100.6 | 59-181     |            113.7 | 66-207      |
|  2009 |            56 | 54-67                                        |            94.7 | 56-162     |            106.1 | 62-183      |
|  2010 |            75 | 63-76                                        |            90.2 | 53-164     |             99.6 | 58-178      |
|  2011 |           306 | 287-305                                      |            90.8 | 54-161     |             99.8 | 60-175      |
|  2012 |           165 | 149-165                                      |           136.9 | 80-258     |            164.8 | 94-333      |
|  2013 |            88 | 77-90                                        |           124.1 | 72-236     |            149.1 | 85-289      |
|  2014 |           107 | 86-106                                       |           117   | 65-226     |            136.2 | 76-267      |
|  2015 |            86 | 74-86                                        |           111.1 | 63-213     |            128.8 | 72-250      |
|  2016 |            62 | 54-63                                        |           109.8 | 63-219     |            125.3 | 72-249      |
|  2017 |           155 | 139-155                                      |           106.7 | 58-216     |            119.1 | 64-243      |
|  2018 |            73 | 67-78                                        |           115.4 | 62-245     |            129.9 | 69-284      |
|  2019 |           114 | 98-112                                       |           108.2 | 60-224     |            119.3 | 65-247      |
|  2020 |           246 | 220-242                                      |           110.9 | 61-231     |            123.1 | 68-263      |
|  2021 |           164 | 142-162                                      |           134.8 | 76-258     |            147.7 | 83-297      |
|  2022 |           109 | 93-108                                       |           124.9 | 73-257     |            132.5 | 76-271      |
|  2023 |           783 | 730-774                                      |           126.1 | 71-259     |            132.7 | 75-273      |
|  2024 |           132 | 119-133                                      |           196.1 | 117-354    |            215   | 128-404     |
|  2025 |           462 | 415-454                                      |           164.3 | 97-310     |            177.2 | 102-335     |
