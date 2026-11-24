# Eğitim raporu: turkey_base

> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Model ve varsayımlar doğrulanmamıştır. Sonuçlar, modelin eğitimde hiç görmediği 2020–2025 yıllarında, kabul görmüş referans model ETAS'a karşı ölçülmüştür.

## Sade dille özet

* **M≥3** (13060 olay): model ETAS'tan **istatistiksel olarak daha iyi** (olay başına bilgi kazancı 0.062, p = 8.51e-36)
* **M≥4** (1896 olay): model ETAS'tan **daha kötü** (olay başına bilgi kazancı -0.027, p = 0.0318)
* **M≥5** (130 olay): model ile ETAS arasında **anlamlı fark yok** (olay başına bilgi kazancı 0.021, p = 0.608)
* **M≥6** (15 olay): model ile ETAS arasında **anlamlı fark yok** (olay başına bilgi kazancı 0.081, p = 0.531)

Bilgi kazancı pozitifse model depremlerin olduğu hücrelere ETAS'tan daha yüksek olasılık vermiştir. p < 0,05 farkın şans eseri olma ihtimalinin düşük olduğunu gösterir; tek bir sınama dönemi olduğu için bu bile bir kanıt değil, ileriye dönük test edilmesi gereken bir işarettir.

## Eşik eşik sonuçlar

| eşik   |   gözlenen_olay |   ETAS_beklenen |   model_beklenen |   ETAS_N_testi_q |   model_N_testi_q |   logL_ETAS |   logL_model |   bilgi_kazancı_olay_başına |    p_değeri |   Brier_ETAS |   Brier_model |
|:-------|----------------:|----------------:|-----------------:|-----------------:|------------------:|------------:|-------------:|----------------------------:|------------:|-------------:|--------------:|
| M≥3    |           13060 |         10463.2 |           9288   |           0      |            0      |    -34514.3 |     -33708.1 |                       0.062 | 8.51363e-36 |      0.24822 |       0.19217 |
| M≥4    |            1896 |           926.5 |            822.5 |           0      |            0      |     -6532.5 |      -6583.8 |                      -0.027 | 0.0317979   |      0.10248 |       0.10328 |
| M≥5    |             130 |            82   |             72.8 |           0      |            0      |      -568   |       -565.2 |                       0.021 | 0.607693    |      0.02022 |       0.0203  |
| M≥6    |              15 |             7.3 |              6.4 |           0.0079 |            0.0028 |       -89.6 |        -88.4 |                       0.081 | 0.530688    |      0.00424 |       0.00423 |

## Yıl yıl (M≥3)

|   yıl |   gözlenen |   ETAS |   model |   logL_ETAS |   logL_model |
|------:|-----------:|-------:|--------:|------------:|-------------:|
|  2020 |       1921 | 1473.4 |  1186.4 |     -3353.7 |      -3122.4 |
|  2021 |       1209 | 1748.4 |  1442.6 |     -2649.6 |      -2743.8 |
|  2022 |        766 | 1573   |  1416.3 |     -1142.3 |      -1039.7 |
|  2023 |       4512 | 1573   |  1625.4 |    -12635   |     -13463.3 |
|  2024 |       1033 | 2224.5 |  1844.2 |     -1410.7 |      -1380.5 |
|  2025 |       3619 | 1870.9 |  1773.1 |    -13322.9 |     -11958.4 |

## ETAS parametreleri (eğitim dönemi)

|   mu |      K |   alpha |      c |      p |      D |      q |   gamma |
|-----:|-------:|--------:|-------:|-------:|-------:|-------:|--------:|
| 0.41 | 0.0194 |  1.7217 | 0.0141 | 1.0585 | 2.8209 | 2.2232 |   0.923 |

Dallanma oranı (10 yıl kesilmiş): 0.76, b = 1.05, hedef olay sayısı 3932.

## Öğrenilen model

Seçilen efektif sürtünme μ' = 0.0 (doğrulama olabilirliği: 0.0: -2650.5, 0.4: -2668.6, 0.8: -2659.6).

En büyük 15 katsayı (standartlaştırılmış özellik başına log-oran etkisi; pozitif = o özellik yüksekken ETAS'ın öngördüğünden fazla deprem):

| özellik                       |   katsayı |
|:------------------------------|----------:|
| frac_dependent_10y            |     0.482 |
| n_mc_1y                       |    -0.473 |
| yrs_since_m4__eksik           |     0.467 |
| yrs_since_m5__eksik           |     0.447 |
| yrs_since_m6__eksik           |     0.416 |
| n_mc_30y                      |     0.329 |
| rate_z_1y                     |     0.328 |
| max_mag_all                   |     0.326 |
| log_m0_5y                     |    -0.254 |
| slip_deficit_mean__eksik      |     0.229 |
| slip_deficit_max__eksik       |     0.229 |
| log_m0_5y__eksik              |    -0.22  |
| aftershocks_per_m4_10y__eksik |    -0.196 |
| slip_deficit_max              |     0.195 |
| depth_mean_10y                |     0.195 |

## Kalibrasyon (sınama dönemi)

Her satır, modelin benzer olasılık verdiği hücre-yıllarını gruplar: iyi kalibre bir modelde ortalama olasılık ile gözlenen sıklık birbirine yakındır.

**M3 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0.009 |      0.056 | 410 |
|    0.087 |      0.13  | 409 |
|    0.24  |      0.222 | 410 |
|    0.425 |      0.269 | 409 |
|    0.621 |      0.408 | 409 |
|    0.828 |      0.512 | 410 |
|    0.972 |      0.731 | 409 |
|    1     |      0.91  | 410 |

**M4 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0.001 |      0.01  | 410 |
|    0.008 |      0.024 | 409 |
|    0.024 |      0.041 | 410 |
|    0.048 |      0.071 | 409 |
|    0.083 |      0.081 | 409 |
|    0.148 |      0.132 | 410 |
|    0.295 |      0.281 | 409 |
|    0.641 |      0.427 | 410 |

**M5 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0     |      0     | 410 |
|    0.001 |      0     | 409 |
|    0.002 |      0.01  | 410 |
|    0.004 |      0.012 | 409 |
|    0.008 |      0.007 | 409 |
|    0.014 |      0.024 | 410 |
|    0.031 |      0.037 | 409 |
|    0.107 |      0.078 | 410 |

**M6 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0     |      0     | 410 |
|    0     |      0     | 409 |
|    0     |      0     | 410 |
|    0     |      0.002 | 409 |
|    0.001 |      0.002 | 409 |
|    0.001 |      0.005 | 410 |
|    0.003 |      0.005 | 409 |
|    0.01  |      0.02  | 410 |
