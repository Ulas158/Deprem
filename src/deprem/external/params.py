"""Gaz, sıvı ve dış etkiler için eleman parametreleri ve gizli parametreler.

Fay elemanı kaydına (``deprem.faults.parameters.REGISTRY``) üç yeni grup ekler:

* ``akışkan ve gaz``  radon, CO2, helyum, yeraltı suyu, kaynak suyu kimyası,
  sürü/difüzyon göstergeleri
* ``dış yükler``      gelgit, mevsimsel hidrolojik yük, atmosfer basıncı,
  uzak büyük depremlerin dinamik tetiklemesi, yavaş kayma
* ``insan kaynaklı``  baraj, maden, jeotermal, petrol/gaz, gaz deposu

Durum (status) değerleri: ``seismicity`` (katalogdan hesaplanır),
``external`` (bu depodaki saha listesinden hesaplanır), ``placeholder``
(veri kaynağı bağlanmadı; NaN). Veri kaynaklarının erişim durumu
``DATA_SOURCES`` tablosundadır.
"""

from __future__ import annotations

import pandas as pd

from ..faults.parameters import register
from ..forecast.latent import LATENT, LatentParam

FLUID = "akışkan ve gaz"
LOADS = "dış yükler"
HUMAN = "insan kaynaklı"

# (ad, birim, grup, durum, açıklama, kaynak)
PARAMS = [
    # ---- akışkan ve gaz: ölçüm gerektirenler
    ("radon_bq_m3", "Bq/m³", FLUID, "placeholder", "Toprak gazı / kaynak suyu radon (²²²Rn) derişimi", "AFAD/KOERI/üniversite istasyonları (izinle)"),
    ("radon_anomaly_z", "-", FLUID, "placeholder", "Radonun mevsim ve meteoroloji etkisi çıkarılmış anomali z-değeri", "Türetilecek"),
    ("co2_flux_g_m2_d", "g/m²/gün", FLUID, "placeholder", "Toprak CO2 akısı (fay boyunca derin kaynaklı gaz)", "Saha ölçümleri / literatür"),
    ("he3_he4_rc_ra", "R/Ra", FLUID, "placeholder", "³He/⁴He oranı: manto kaynaklı akışkanın fay boyunca yükselmesi", "Literatür (ör. Güleç vd. 2002; Mutlu vd. 2008)"),
    ("ch4_h2_ppm", "ppm", FLUID, "placeholder", "Toprak gazında CH4 / H2", "Saha ölçümleri"),
    ("groundwater_level_m", "m", FLUID, "placeholder", "Yeraltı suyu seviyesi (kuyu)", "DSİ gözlem kuyuları (izinle)"),
    ("groundwater_anomaly_z", "-", FLUID, "placeholder", "Yeraltı suyu seviyesinde yağış/pompaj dışı anomali", "Türetilecek"),
    ("spring_temp_c", "°C", FLUID, "placeholder", "Sıcak su kaynağı sıcaklığı", "MTA jeotermal envanteri (elle)"),
    ("spring_chem_anomaly_z", "-", FLUID, "placeholder", "Kaynak suyu kimyası (Cl, Na, B, iletkenlik) anomalisi", "Saha ölçümleri / literatür"),
    ("thermal_spring_dist_km", "km", FLUID, "placeholder", "En yakın sıcak su kaynağına uzaklık", "MTA envanteri (elle)"),
    ("pore_pressure_lambda", "-", FLUID, "placeholder", "Gözenek basıncı / litostatik basınç (λ) yerel tahmini; ölçüm yoksa gizli parametre 'pore_pressure_ratio' kullanılır", "Kuyu ölçümleri / Vp/Vs / ters çözüm"),
    ("station_dist_km", "km", FLUID, "external", "Gaz/su verisi aktarılan en yakın istasyona uzaklık (aktarımın güvenilirliği)", "Ölçüm istasyonları dosyası"),
    # ---- akışkan ve gaz: katalogdan
    ("swarm_n_10y", "-", FLUID, "seismicity", "Son 10 yılda yakında başlayan deprem sürüsü sayısı", "Katalog (swarms.py)"),
    ("swarm_frac_10y", "-", FLUID, "seismicity", "Yakındaki ≥5 olaylı kümelerin sürü oranı (akışkan sürülerini ana şok–artçıdan ayırır)", "Katalog (swarms.py)"),
    ("hydraulic_diffusivity_log10", "log10 m²/s", FLUID, "seismicity", "Sürü göçünden hidrolik yayınım D (r² = 4πDt, Shapiro 1997)", "Katalog (swarms.py)"),
    # ---- dış yükler: katalogdan
    ("tidal_semidiurnal_log10p", "log10 p", LOADS, "seismicity", "Ay'ın yarı günlük gelgit fazı Schuster testi (bağımsız olaylar)", "Katalog + gök mekaniği (periodic.py)"),
    ("hour_of_day_log10p", "log10 p", LOADS, "seismicity", "Günün saatine yığılma (kontrol: gece/gündüz algılama farkı, taş ocağı patlamaları)", "Katalog"),
    ("tidal_fortnight_log10p", "log10 p", LOADS, "seismicity", "14 günlük (kuvvetli/zayıf gelgit) faz Schuster testi", "Katalog + gök mekaniği"),
    ("seasonal_log10p", "log10 p", LOADS, "seismicity", "Yıllık faz Schuster testi (mevsimsel yük vekili)", "Katalog"),
    ("dyn_trigger_beta", "-", LOADS, "seismicity", "Uzak büyük depremler sonrası 5 gün β istatistiği (dinamik tetikleme)", "Katalog + USGS küresel M≥6,5"),
    # ---- dış yükler: veri gerekenler
    ("hydro_load_amp_kpa", "kPa", LOADS, "placeholder", "Yıllık hidrolojik yük genliği (GRACE eşdeğer su kalınlığı)", "NASA GRACE/GRACE-FO (Earthdata hesabı)"),
    ("snow_load_amp_kpa", "kPa", LOADS, "placeholder", "Kar yükü yıllık genliği", "ERA5-Land (Copernicus hesabı)"),
    ("precip_anomaly_z", "-", LOADS, "placeholder", "Yağış anomalisi (yağışla tetiklenen sismisite, Hainzl vd. 2006)", "ERA5 / MGM"),
    ("atm_pressure_amp_hpa", "hPa", LOADS, "placeholder", "Atmosfer basıncı değişim genliği", "ERA5"),
    ("lake_level_change_m", "m", LOADS, "placeholder", "Doğal göl/rezervuar seviye değişimi", "DSİ / uydu altimetrisi (DAHITI, Hydroweb)"),
    ("slow_slip_rate_mm_yr", "mm/yıl", LOADS, "placeholder", "Yavaş kayma / geçici aseismik kayma hızı", "GNSS zaman serileri (Nevada Geodetic Lab)"),
    # ---- diğer öncü sinyal adayları (kanıt zayıf; karşılaştırma için tanımlı)
    ("tec_anomaly_z", "-", "diğer öncü adayları", "placeholder", "İyonosfer toplam elektron içeriği (TEC) anomalisi", "deprem.observations.ionosphere (hücre düzeyinde hesaplanıyor)"),
    ("geomagnetic_anomaly_z", "-", "diğer öncü adayları", "placeholder", "Yerel manyetik alan / ULF anomalisi", "deprem.observations.geomag (hücre düzeyinde hesaplanıyor)"),
    ("thermal_ir_anomaly_z", "-", "diğer öncü adayları", "placeholder", "Uydu termal kızılötesi (yüzey sıcaklığı) anomalisi", "deprem.observations.thermal (hücre düzeyinde hesaplanıyor)"),
    # ---- insan kaynaklı: saha listesinden
    ("site_reservoir_dist_km", "km", HUMAN, "external", "En yakın faal baraj rezervuarına uzaklık", "external_data/insan_kaynakli_sahalar.csv"),
    ("reservoir_years_since_fill", "yıl", HUMAN, "external", "50 km içindeki en yakın barajın dolumundan bu yana geçen süre", "Aynı"),
    ("site_mining_dist_km", "km", HUMAN, "external", "En yakın faal madene uzaklık", "Aynı"),
    ("site_geothermal_dist_km", "km", HUMAN, "external", "En yakın faal jeotermal sahaya uzaklık", "Aynı + WRI GPPD"),
    ("site_hydrocarbon_dist_km", "km", HUMAN, "external", "En yakın petrol/gaz sahası veya gaz deposuna uzaklık", "Aynı"),
    ("site_injection_dist_km", "km", HUMAN, "external", "En yakın akışkan enjeksiyonu yapılan sahaya uzaklık", "Aynı"),
    ("site_gas_seep_dist_km", "km", HUMAN, "external", "En yakın bilinen doğal gaz çıkışına uzaklık", "Aynı (Géli vd. 2018)"),
    # ---- insan kaynaklı: veri gerekenler
    ("injection_volume_m3_yr", "m³/yıl", HUMAN, "placeholder", "Yakındaki toplam akışkan enjeksiyon hacmi", "Şirket/EPDK/MAPEG verisi (izinle)"),
    ("extraction_volume_m3_yr", "m³/yıl", HUMAN, "placeholder", "Yakındaki akışkan/kütle çekimi", "Şirket/MAPEG verisi (izinle)"),
    ("reservoir_level_rate_m_month", "m/ay", HUMAN, "placeholder", "Rezervuar seviye değişim hızı", "DSİ (izinle) / uydu altimetrisi"),
]

for _n, _u, _g, _s, _d, _src in PARAMS:
    register(_n, _u, _g, _s, _d, _src)

# Daha önce "henüz veri yok" olarak tanımlanan gelgit alanı artık katalogdan hesaplanıyor.
register("tidal_stress_sensitivity", "log10 p", LOADS, "seismicity",
         "Gelgit duyarlılığı (= tidal_semidiurnal_log10p; geriye uyumluluk)", "Katalog + gök mekaniği")


# ---------------------------------------------------------------- gizli parametreler
NEW_LATENT = [
    LatentParam("hydraulic_diffusivity", "Hidrolik yayınım D [m²/s] (basınç cephesi r = √(4πDt))",
                0.3, 0.005, 20.0, -0.5, 0.7, log=True),
    LatentParam("fluid_gain", "Akışkan göstergelerinin (sürü, Vp/Vs, ısı akısı) arka plan oranına etkisi (0 = etkisiz)",
                0.0, -2.0, 2.0, 0.0, 0.5),
    LatentParam("tidal_gain", "Gelgit fazına göre oran modülasyonu genliği (1 + g·cos φ); Dieterich'e göre Δτ/Aσ",
                0.0, 0.0, 0.5, 0.0, 0.05),
    LatentParam("seasonal_gain", "Mevsimsel (yıllık) oran modülasyonu genliği", 0.0, 0.0, 0.5, 0.0, 0.05),
    LatentParam("seasonal_phase_deg", "En yüksek sismisite oranının yıllık fazı (1 Ocak = 0°)", 0.0, -180.0, 180.0, 0.0, 90.0),
    LatentParam("dynamic_trigger_gain", "Uzak büyük deprem sonrası birkaç günlük oran artışı katsayısı",
                0.0, 0.0, 5.0, 0.0, 0.5),
    LatentParam("reservoir_decay_yr", "Rezervuar tetiklemesinin dolumdan sonra sönüm süresi [yıl]",
                5.0, 0.3, 50.0, 0.7, 0.5, log=True),
    LatentParam("injection_gain", "Akışkan enjeksiyonu/çekimi yakınlığının oran etkisi (0 = etkisiz)",
                0.0, -1.0, 3.0, 0.0, 0.5),
]
for _p in NEW_LATENT:
    LATENT.setdefault(_p.name, _p)


# ---------------------------------------------------------------- veri kaynakları
# erişim: acik (doğrudan indirilebilir), kayit (ücretsiz hesap gerekir),
#         izin (kurumdan izin/yazışma gerekir), elle (yayınlardan derlenmeli)
DATA_SOURCES = pd.DataFrame([
    ("Birleşik deprem kataloğu", "USGS, EMSC, AFAD (ISC)", "acik", "bağlı", "Sürü, gelgit, mevsim, dinamik tetikleme özellikleri"),
    ("USGS küresel M≥6,5 kataloğu", "earthquake.usgs.gov FDSN", "acik", "bağlı", "Uzak tetikleme için"),
    ("WRI Küresel Santral Veritabanı", "github.com/wri/global-power-plant-database (CC BY 4.0)", "acik", "bağlı", "Baraj ve jeotermal santral konumları"),
    ("Maden, petrol/gaz, jeotermal saha listesi", "Bu depoda elle derlendi (external_data/)", "elle", "bağlı (yaklaşık)", "Konumlar ve yıllar doğrulanmalı; MAPEG/MTA ile karşılaştırılmalı"),
    ("Küresel Isı Akısı Veritabanı (IHFC 2024)", "ihfc-iugg.org / GFZ Data Services", "acik", "bu ortamdan erişilemedi", "Ağ izni verilirse otomatik indirilebilir"),
    ("MTA jeotermal envanteri ve sıcak su kaynakları", "mta.gov.tr yayınları", "elle", "bağlı değil", "Tablo olarak indirme sunulmuyor"),
    ("AFAD radon / öncül gözlem istasyonları", "AFAD", "izin", "bağlı değil", "Zaman serileri açık değil; resmi başvuru gerekir"),
    ("KOERI / üniversite radon ve gaz ölçümleri", "Yayınlar ve araştırmacılar", "izin", "bağlı değil", "Yayınlarda grafik olarak var, ham veri yazarlardan istenmeli"),
    ("DSİ yeraltı suyu gözlem kuyuları ve baraj seviyeleri", "DSİ", "izin", "bağlı değil", "Kurumsal başvuru gerekir"),
    ("GRACE/GRACE-FO hidrolojik yük", "NASA PO.DAAC", "kayit", "bağlı değil", "Ücretsiz NASA Earthdata hesabı"),
    ("ERA5 / ERA5-Land (yağış, kar, basınç)", "Copernicus CDS", "kayit", "bağlı değil", "Ücretsiz Copernicus hesabı"),
    ("Uydu altimetrisi göl/rezervuar seviyeleri", "DAHITI / Hydroweb", "kayit", "bağlı değil", "Ücretsiz hesap"),
    ("GNSS zaman serileri (yavaş kayma)", "Nevada Geodetic Laboratory", "acik", "bu ortamdan erişilemedi", "Ağ izni verilirse indirilebilir"),
    ("Akışkan enjeksiyon/çekim hacimleri", "Şirketler, EPDK, MAPEG", "izin", "bağlı değil", "Genellikle gizli; resmi başvuru"),
    ("OpenStreetMap madenler, taş ocakları, barajlar", "Overpass API / Geofabrik", "acik", "bu ortamdan erişilemedi", "Ağ izni verilirse ~binlerce saha otomatik eklenir"),
], columns=["veri", "kaynak", "erisim", "durum", "not"])
