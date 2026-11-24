# Gaz, sıvı ve dış etkiler: mekanizmalar, kanıt düzeyi ve modeldeki karşılıkları

> Bu not bir literatür özeti ve modelleme planıdır. Buradaki hiçbir etki doğrulanmış bir deprem
> öncüsü değildir. Her biri modelde **sınanacak bir değişken** olarak tanımlandı; ancak eğitimde
> görülmemiş dönemde referans modellerden (düzleştirilmiş sismisite, ETAS) daha iyi tahmin
> sağlarsa anlamlı sayılacaktır.

## 1. Temel fikir: sıvılar depremi nasıl "kolaylaştırır"?

Bir fay, üzerindeki makaslama gerilmesi sürtünme direncini aştığında kayar (Coulomb ölçütü):

    τ ≥ μ · (σn − p) + C

* τ: fayı kaydırmaya çalışan makaslama gerilmesi (levha hareketiyle yavaşça birikir)
* σn: fayı sıkıştıran normal gerilme, μ: sürtünme katsayısı (~0,6), C: kohezyon
* **p: kayaç gözeneklerindeki sıvı (su, CO2, gaz) basıncı**

Gözenek basıncı p arttıkça fayı "kenetleyen" etkin gerilme (σn − p) azalır. Yani sıvı, fayı
bir hava yastığı gibi hafifletir. **Enerjiyi sıvı vermez**; enerji levha hareketiyle fayda
birikmiştir. Sıvı yalnızca kırılmanın **ne zaman** olacağını öne çekebilir (tetikler). Gazlar
(radon, helyum, CO2) ise çoğunlukla bu süreçlerin **izleyicisidir**: kayaç çatladığında veya
sıvı göç ettiğinde yüzeye çıkışları değişebilir.

Önemli sınırlama: bir faydaki gerilme kırılmaya zaten çok yakın değilse küçük bir basınç
değişimi deprem üretmez. Bu yüzden aynı etki (ör. barajın dolması) bir yerde deprem
tetiklerken başka bir yerde hiçbir şey yapmaz.

## 2. Etkiler ve kanıt düzeyi

Kanıt düzeyi: **Güçlü** (tekrarlanmış, kontrollü veya çok sayıda örnek), **Orta** (birden çok
bağımsız çalışma, etki küçük veya bölgeye bağlı), **Zayıf** (tekil/geriye dönük örnekler,
tutarsız), **Desteklenmiyor** (sistematik testlerde tutmadı).

| Etki | Mekanizma | Kanıt | Önemli çalışmalar | Modeldeki karşılığı |
|---|---|---|---|---|
| Akışkan enjeksiyonu (atık su, jeotermal, EGS) | Gözenek basıncı artışı | **Güçlü** | Rangely deneyi (Raleigh vd. 1976: basınç açılıp kapatılarak depremler açılıp kapandı); Oklahoma 2009–2016 (Ellsworth 2013; Keranen vd. 2014), Pohang 2017 M5,5 (Grigoli vd. 2018) | `site_injection_dist_km`, `injection_volume_m3_yr` (veri gerekir), gizli `injection_gain` |
| Akışkan/gaz çekimi | Rezervuar sıkışması, gerilme değişimi | **Güçlü** | Groningen gaz sahası (Hollanda) | `site_hydrocarbon_dist_km`, `site_geothermal_dist_km` |
| Baraj (rezervuar) | Su yükü (anında) + basıncın difüzyonla derine inmesi (aylar–yıllar) | **Güçlü** (varlığı), büyük depremlerde tartışmalı | Koyna 1967 M6,3; Kremasta 1966; Zipingpu–Wenchuan 2008 tartışması (Gupta 2002) | `site_reservoir_dist_km`, `reservoir_years_since_fill`, gizli `reservoir_decay_yr` |
| Madencilik | Kütle çıkarımı, gerilme yeniden dağılımı | **Güçlü** ama genellikle M<5 | Kaya patlamaları; HiQuake veritabanı (Foulger vd. 2018) | `site_mining_dist_km` |
| Doğal sıvı göçü ve deprem sürüleri | Basınç cephesi r ≈ √(4πDt) ile ilerler | **Güçlü–Orta** | Matsushiro 1965–67; Shapiro vd. 1997; Miller vd. 2004 (İtalya 1997 artçıları CO2 basıncıyla); Ross vd. 2020 | `swarm_n_10y`, `swarm_frac_10y`, `hydraulic_diffusivity_log10`, gizli `hydraulic_diffusivity`, `pore_pressure_ratio` |
| Uzak büyük depremlerin dalgaları (dinamik tetikleme) | Kısa süreli gerilme salınımı, akışkan yeniden dağılımı | **Güçlü** küçük depremler için, büyükler için nadir | Landers 1992 (Hill vd. 1993); 2012 Hint Okyanusu M8,6 sonrası küresel artış (Pollitz vd. 2012) | `dyn_trigger_beta`, gizli `dynamic_trigger_gain` |
| Gelgit (Ay–Güneş) | Kabukta ~0,001–0,01 MPa periyodik gerilme | **Orta** (etki %1–10, bölgeye bağlı) | Cochran vd. 2004; Tanaka 2012 (Tohoku öncesi duyarlılık artışı); Ide vd. 2016 | `tidal_semidiurnal_log10p`, `tidal_fortnight_log10p`, gizli `tidal_gain` |
| Mevsimsel hidrolojik yük (yağış, kar, yeraltı suyu) | Yüzey yükü ve gözenek basıncı değişimi | **Orta** | Himalaya (Bettinelli vd. 2008); Kaliforniya (Johnson vd. 2017); Japonya kar yükü (Heki 2003); yağışla tetiklenme (Hainzl vd. 2006) | `seasonal_log10p`, `hydro_load_amp_kpa`, `snow_load_amp_kpa`, `precip_anomaly_z`, gizli `seasonal_gain` |
| Yavaş kayma / öncü kayma | Fayın bir kısmı sessizce kayıp yükü kilitli kısma aktarır | **Orta** (bazı büyük depremlerde gözlendi) | Tohoku 2011 (Kato vd. 2012); İzmit 1999 çekirdeklenme evresi (Bouchon vd. 2011) | `slow_slip_rate_mm_yr` (GNSS gerekir) |
| Yeraltı suyu kimyası | Çatlaklanma ile farklı akiferlerin karışması | **Zayıf–Orta** | İzlanda Húsavík (Skelton vd. 2014) | `spring_chem_anomaly_z`, `groundwater_anomaly_z` |
| Radon | Çatlaklanma ve gaz göçü | **Zayıf** (tutarsız; meteoroloji ve mevsim etkisi büyük) | Kobe 1995 (Igarashi vd. 1995, geriye dönük); ICEF 2011 raporu: güvenilir öncü değil | `radon_bq_m3`, `radon_anomaly_z` |
| Helyum ve CO2 | Derin (manto) kökenli akışkanın fay boyunca yükselmesi | Mekânsal gösterge olarak **Güçlü**, zamansal öncü olarak **Zayıf** | Anadolu'da manto helyumu (Güleç vd. 2002) | `he3_he4_rc_ra`, `co2_flux_g_m2_d`, `ch4_h2_ppm` |
| Isı akısı, jeotermal alanlar | Sıcak, sıvıca zengin kabuk; sığ kırılgan zon | Akışkan ve sürü yoğunluğu için **Orta** | — | `heat_flow_mw_m2`, `thermal_spring_dist_km`, `site_geothermal_dist_km` |
| Dilatans-difüzyon (öncü Vp/Vs düşüşü) | Kırılma öncesi çatlak açılması ve sıvı girişi | **Desteklenmiyor** (1970'lerdeki iddialar tekrarlanamadı) | Scholz vd. 1973; sonraki testler | `vp_vs_ratio` (yalnızca mekânsal akışkan göstergesi olarak) |
| Atmosfer basıncı | Küçük yüzey yükü | **Zayıf** | — | `atm_pressure_amp_hpa` |
| İyonosfer (TEC), manyetik, termal kızılötesi | Çeşitli öneriler | **Zayıf / tartışmalı** | — | `tec_anomaly_z`, `geomagnetic_anomaly_z`, `thermal_ir_anomaly_z` |
| Güneş etkinliği, gezegen dizilimi, hayvan davranışı | — | **Desteklenmiyor** | — | Eklenmedi |

Türkiye'ye özgü örnekler: Marmara Denizi tabanındaki gaz çıkışlarının fay ve sismisiteyle
ilişkisi (Géli vd. 2018), Kuzey Anadolu Fayı boyunca manto kökenli helyum (Güleç vd. 2002),
Batı Anadolu grabenlerinde (Gediz, Büyük Menderes, Simav) yoğun jeotermal etkinlik ve sürüler,
6 Şubat 2023 depremlerinden sonra Hatay termal sularında radon değişimleri (olay sonrası).

## 3. Modele nasıl eklendi?

Kod: `src/deprem/external/` ve `src/deprem/faults/fill.py`.

1. **Her 1 km × 1 km fay elemanı için ayrı parametre** (`deprem.faults.parameters.REGISTRY`, 3 yeni
   grup: akışkan ve gaz, dış yükler, insan kaynaklı; ayrıca "diğer öncü adayları").
   Ölçüm istasyonu verisi (radon, CO2, yeraltı suyu...) gelirse `elements.stations_to_elements`
   ile 30 km içindeki istasyonlardan ters uzaklık ağırlıklı elemanlara aktarılır.
2. **Eksik değerler boş bırakılmaz** (`fill.fill_missing`): önce aynı fay boyunca sağ ve soldaki
   en yakın gerçek değerin uzaklık ağırlıklı ortalaması, fayda veri yoksa komşu faylar, o da
   yoksa bölgesel ortanca, hiç ölçüm yoksa literatür varsayılanı. Her değerin kaynağı
   (`__fill`) ve belirsizliği (`__sd`, gerçek veriden uzaklaştıkça büyür) tutulur. Eğitimde
   doldurulmuş değerler öğrenilebilir parametredir; doldurma değeri öncül ortalama, belirsizlik
   öncül genişliğidir (`trainable_vector`, `gaussian_log_prior`). Hiç ölçümü olmayan ve makul
   literatür değeri bulunmayan bir büyüklük (ör. radon) uydurulmaz.
3. **Katalogdan hemen hesaplananlar**: gelgit (yarı günlük, 14 günlük) ve mevsim fazı Schuster
   testleri (`periodic.py`, Ay–Güneş konumları `astro.py`), deprem sürüleri ve hidrolik yayınım
   (`swarms.py`), uzak büyük depremler sonrası β istatistiği (`dynamic.py`, USGS küresel M≥6,5).
4. **İnsan kaynaklı sahalar** (`external_data/insan_kaynakli_sahalar.csv`): 78 baraj, maden,
   jeotermal, petrol/gaz sahası, gaz deposu ve doğal gaz çıkışı; başlangıç/bitiş yılı ve durum
   (aktif/pasif/kapalı) ile. Özellikler yalnızca o tarihte faal olan sahalardan hesaplanır.
   Barajların konumu WRI Küresel Santral Veritabanı'ndan; madenler, petrol/gaz ve bazı
   jeotermal sahaların konum ve yılları elle derlendi ve **doğrulanmalıdır** (`verified=hayir`).
5. **8 yeni gizli parametre** (`forecast/latent.py`): hydraulic_diffusivity, fluid_gain,
   tidal_gain, seasonal_gain, seasonal_phase_deg, dynamic_trigger_gain, reservoir_decay_yr,
   injection_gain. Hepsinin öncülü "etki yok" (0) civarında başlar; veri aksini gösterirse
   değişir.
6. **Harita katmanları** (`web/data/layers/ext_*`): her biri ayrı aç/kapa.

## 4. Veri erişimi

| Veri | Durum | Nasıl alınır |
|---|---|---|
| Deprem kataloğu (USGS, EMSC, AFAD/ISC), USGS küresel M≥6,5 | Bağlı | Otomatik |
| WRI Küresel Santral Veritabanı (baraj, jeotermal) | Bağlı | Otomatik (CC BY 4.0) |
| Küresel Isı Akısı Veritabanı (IHFC 2024) | Açık ama bu ortamdan erişilemedi | ihfc-iugg.org / GFZ için ağ izni veya elle indirme |
| OpenStreetMap madenler, taş ocakları, barajlar | Açık ama erişilemedi | Overpass API için ağ izni |
| GNSS zaman serileri (yavaş kayma) | Açık ama erişilemedi | Nevada Geodetic Laboratory için ağ izni |
| GRACE, ERA5, uydu altimetrisi | Ücretsiz hesap gerekir | NASA Earthdata / Copernicus / DAHITI |
| MTA jeotermal envanteri, sıcak su kaynakları | Elle derlenmeli | MTA yayınları |
| AFAD radon istasyonları, DSİ kuyu ve baraj seviyeleri, enjeksiyon hacimleri | İzin gerekir | Kurumlara resmi başvuru |

İstasyon verisi şu biçimde bir CSV olarak verilebilir: `lat, lon, param, value, time`
(`param` kayıttaki adlardan biri, ör. `radon_bq_m3`).

## 5. Kaynaklar (seçme)

Bettinelli vd. 2008 EPSL · Bouchon vd. 2011 Science · Cochran vd. 2004 Science · Ellsworth 2013
Science · Foulger vd. 2018 Earth-Sci. Rev. · Géli vd. 2018 Sci. Rep. · Grigoli vd. 2018 Science ·
Güleç vd. 2002 Chem. Geol. · Gupta 2002 Earth-Sci. Rev. · Hainzl vd. 2006 JGR · Heki 2003 EPSL ·
Hill vd. 1993 Science · Ide vd. 2016 Nat. Geosci. · Igarashi vd. 1995 Science · ICEF 2011 Ann.
Geophys. · Johnson vd. 2017 Science · Kato vd. 2012 Science · Keranen vd. 2014 Science · Matthews &
Reasenberg 1988 · Miller vd. 2004 Nature · Pollitz vd. 2012 Nature · Raleigh vd. 1976 Science ·
Ross vd. 2020 Science · Scholz vd. 1973 Science · Shapiro vd. 1997 GJI · Skelton vd. 2014 Nat.
Geosci. · Tanaka 2012 GRL
