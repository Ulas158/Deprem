# Eğitim raporu: turkey_external

> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Model ve varsayımlar doğrulanmamıştır. Sonuçlar, modelin eğitimde hiç görmediği 2020–2025 yıllarında, kabul görmüş referans model ETAS'a karşı ölçülmüştür.

## Sade dille özet

* **M≥3** (13060 olay): model ile ETAS arasında **anlamlı fark yok** (olay başına bilgi kazancı 0.006, p = 0.244)
* **M≥4** (1896 olay): model ETAS'tan **daha kötü** (olay başına bilgi kazancı -0.068, p = 6.45e-08)
* **M≥5** (130 olay): model ile ETAS arasında **anlamlı fark yok** (olay başına bilgi kazancı -0.015, p = 0.733)
* **M≥6** (15 olay): model ile ETAS arasında **anlamlı fark yok** (olay başına bilgi kazancı 0.086, p = 0.43)

Bilgi kazancı pozitifse model depremlerin olduğu hücrelere ETAS'tan daha yüksek olasılık vermiştir. p < 0,05 farkın şans eseri olma ihtimalinin düşük olduğunu gösterir; tek bir sınama dönemi olduğu için bu bile bir kanıt değil, ileriye dönük test edilmesi gereken bir işarettir.

## Eşik eşik sonuçlar

| eşik   |   gözlenen_olay |   ETAS_beklenen |   model_beklenen |   ETAS_N_testi_q |   model_N_testi_q |   logL_ETAS |   logL_model |   bilgi_kazancı_olay_başına |    p_değeri |   Brier_ETAS |   Brier_model |
|:-------|----------------:|----------------:|-----------------:|-----------------:|------------------:|------------:|-------------:|----------------------------:|------------:|-------------:|--------------:|
| M≥3    |           13060 |         10463.2 |           9902.4 |           0      |            0      |    -34514.3 |     -34440.1 |                       0.006 | 0.243799    |      0.24822 |       0.19228 |
| M≥4    |            1896 |           926.5 |            876.9 |           0      |            0      |     -6532.5 |      -6662   |                      -0.068 | 6.44709e-08 |      0.10248 |       0.10606 |
| M≥5    |             130 |            82   |             77.6 |           0      |            0      |      -568   |       -570   |                      -0.015 | 0.732648    |      0.02022 |       0.02043 |
| M≥6    |              15 |             7.3 |              6.9 |           0.0079 |            0.0049 |       -89.6 |        -88.3 |                       0.086 | 0.42978     |      0.00424 |       0.00423 |

## Yıl yıl (M≥3)

|   yıl |   gözlenen |   ETAS |   model |   logL_ETAS |   logL_model |
|------:|-----------:|-------:|--------:|------------:|-------------:|
|  2020 |       1921 | 1473.4 |  1132   |     -3353.7 |      -3124.4 |
|  2021 |       1209 | 1748.4 |  1553.7 |     -2649.6 |      -2821.6 |
|  2022 |        766 | 1573   |  1453.8 |     -1142.3 |      -1118.5 |
|  2023 |       4512 | 1573   |  1631.6 |    -12635   |     -13531.5 |
|  2024 |       1033 | 2224.5 |  2146.5 |     -1410.7 |      -1484.5 |
|  2025 |       3619 | 1870.9 |  1984.8 |    -13322.9 |     -12359.6 |

## ETAS parametreleri (eğitim dönemi)

|   mu |      K |   alpha |      c |      p |      D |      q |   gamma |
|-----:|-------:|--------:|-------:|-------:|-------:|-------:|--------:|
| 0.41 | 0.0194 |  1.7217 | 0.0141 | 1.0585 | 2.8209 | 2.2232 |   0.923 |

Dallanma oranı (10 yıl kesilmiş): 0.76, b = 1.05, hedef olay sayısı 3932.

## Öğrenilen model

Seçilen efektif sürtünme μ' = 0.8 (doğrulama olabilirliği: 0.0: -3014.7, 0.4: -2976.5, 0.8: -2924.0).

En büyük 15 katsayı (standartlaştırılmış özellik başına log-oran etkisi; pozitif = o özellik yüksekken ETAS'ın öngördüğünden fazla deprem):

| özellik                            |   katsayı |
|:-----------------------------------|----------:|
| yrs_since_m5__eksik                |     0.561 |
| yrs_since_m6__eksik                |     0.516 |
| n_mc_1y                            |    -0.454 |
| yrs_since_m4__eksik                |     0.452 |
| max_mag_all                        |     0.389 |
| dyn_trigger_beta__eksik            |     0.377 |
| rate_z_1y                          |     0.353 |
| frac_dependent_10y                 |     0.299 |
| log_m0_5y__eksik                   |    -0.252 |
| n_mc_30y                           |     0.25  |
| hydraulic_diffusivity_log10__eksik |    -0.206 |
| slip_deficit_max                   |     0.186 |
| log_m0_5y                          |    -0.184 |
| reservoir_years_since_fill         |     0.175 |
| yrs_since_m6                       |     0.17  |

## Kalibrasyon (sınama dönemi)

Her satır, modelin benzer olasılık verdiği hücre-yıllarını gruplar: iyi kalibre bir modelde ortalama olasılık ile gözlenen sıklık birbirine yakındır.

**M3 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0.009 |      0.056 | 410 |
|    0.091 |      0.125 | 409 |
|    0.26  |      0.224 | 410 |
|    0.43  |      0.274 | 409 |
|    0.627 |      0.386 | 409 |
|    0.823 |      0.549 | 410 |
|    0.966 |      0.741 | 409 |
|    0.999 |      0.883 | 410 |

**M4 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0.001 |      0.01  | 410 |
|    0.009 |      0.027 | 409 |
|    0.026 |      0.044 | 410 |
|    0.049 |      0.066 | 409 |
|    0.085 |      0.09  | 409 |
|    0.145 |      0.146 | 410 |
|    0.281 |      0.276 | 409 |
|    0.642 |      0.407 | 410 |

**M5 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0     |      0     | 410 |
|    0.001 |      0     | 409 |
|    0.002 |      0.01  | 410 |
|    0.004 |      0.01  | 409 |
|    0.008 |      0.012 | 409 |
|    0.014 |      0.02  | 410 |
|    0.029 |      0.037 | 409 |
|    0.115 |      0.08  | 410 |

**M6 (model)**

|   p_mean |   obs_freq |   n |
|---------:|-----------:|----:|
|    0     |      0     | 410 |
|    0     |      0     | 409 |
|    0     |      0     | 410 |
|    0     |      0     | 409 |
|    0.001 |      0.007 | 409 |
|    0.001 |      0.002 | 410 |
|    0.003 |      0.002 | 409 |
|    0.011 |      0.022 | 410 |
