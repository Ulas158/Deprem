# Gaz, sıvı ve dış etkiler: eğitim karşılaştırması

> **Uyarı.** Bu sonuçlar deprem tahmini değildir; zamana göre ayrılmış sınama döneminde ölçülen hipotez testleridir.

Kurulum: `configs/train_base.toml` ve `configs/train_external.toml` (aynı ETAS referansı, aynı Coulomb ve moment bütçesi). Eğitim 2012–2019, sınama 2020–2025, 0,5° hücre × 1 yıl, eğitim büyüklüğü M≥3. Her özellik grubu temel modele tek tek ve birlikte eklendi. Sınama dönemi log-olabilirliği (büyük = iyi) ve ETAS'a göre olay başına bilgi kazancı:

| model                                           |   logL M≥3 |   logL M≥4 |   logL M≥5 |   logL M≥6 |   bilgi kazancı M≥3 |   bilgi kazancı M≥5 |
|:------------------------------------------------|-----------:|-----------:|-----------:|-----------:|--------------------:|--------------------:|
| ETAS (referans)                                 |   -34514.3 |    -6532.5 |     -568   |      -89.6 |               0     |               0     |
| Temel (sismisite + fay + Coulomb, ETAS üzerine) |   -33708.1 |    -6583.8 |     -565.2 |      -88.4 |               0.062 |               0.021 |
| Temel + gelgit/mevsim/saat                      |   -34416   |    -6720.4 |     -570.6 |      -88.7 |               0.008 |              -0.02  |
| Temel + sürü/hidrolik yayınım                   |   -34304.7 |    -6683.1 |     -568.8 |      -88.8 |               0.016 |              -0.006 |
| Temel + uzak deprem tetiklemesi                 |   -33962.8 |    -6585   |     -566.8 |      -88.6 |               0.042 |               0.009 |
| Temel + insan kaynaklı sahalar                  |   -33923.8 |    -6619.2 |     -569.1 |      -88.5 |               0.045 |              -0.008 |
| Temel + hepsi                                   |   -34440.1 |    -6662   |     -570   |      -88.3 |               0.006 |              -0.015 |

## Sonuç

* Hiçbir gaz/sıvı veya dış etki grubu sınama döneminde temel modeli iyileştirmedi; hepsi M≥3 log-olabilirliğini düşürdü (en az zarar: insan kaynaklı sahalar ve uzak tetikleme, en çok: gelgit/mevsim/saat). Hepsi birlikte eklenince temel modelin ETAS'a karşı kazancı (M≥3'te olay başına 0,062) neredeyse sıfıra iniyor (0,006).
* Bu, katalogdan türetilebilen bu göstergelerin 1 yıllık, 0,5°'lik pencerelerde ek bilgi taşımadığını gösteriyor. Literatürle uyumlu: gelgit ve mevsim etkileri küçük ve kısa vadelidir; sürüler ve uzak tetikleme çoğunlukla küçük depremleri etkiler.
* Sınama dönemine 2023 Kahramanmaraş ve 2025 dizileri hâkim; tek bir dönemle verilen karar kesin değildir. Hipotezler kayıtta kalıyor; ölçüm verisi (radon, yeraltı suyu, ısı akısı, GNSS) bağlandığında ve daha kısa zaman pencereleriyle (gün–hafta) yeniden sınanmalı.
* Özellik sayısı arttıkça aşırı uyum riski artıyor; birleşik eğitimde bu gruplar varsayılan olarak kapalı tutulup yalnızca ileri test geçerse açılmalı.

## Hepsi birlikte modelinde katsayılar (standartlaştırılmış, ETAS'a göre log-oran)

| özellik                     |   katsayı |
|:----------------------------|----------:|
| hour_of_day_log10p          |    -0.151 |
| site_mining_dist_km         |    -0.134 |
| dyn_trigger_beta            |    -0.101 |
| site_geothermal_dist_km     |    -0.099 |
| site_injection_dist_km      |    -0.099 |
| tidal_fortnight_log10p      |    -0.046 |
| seasonal_log10p             |    -0.044 |
| hydraulic_diffusivity_log10 |    -0.037 |
| swarm_frac_10y              |    -0.02  |
| heat_flow_mw_m2             |     0     |
| radon_anomaly_z             |     0     |
| thermal_spring_dist_km      |     0     |
| groundwater_anomaly_z       |     0     |
| swarm_n_10y                 |     0.042 |
| site_hydrocarbon_dist_km    |     0.054 |
| site_gas_seep_dist_km       |     0.055 |
| site_reservoir_dist_km      |     0.07  |
| tidal_semidiurnal_log10p    |     0.099 |
| reservoir_years_since_fill  |     0.175 |

Katsayılar küçük (|β| < 0,3) ve sınamada tutmadığı için fiziksel yorum yapılmamalıdır.

