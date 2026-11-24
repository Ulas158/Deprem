# Tarihsel deprem kataloğu (1900 öncesi, Türkiye ve çevresi)

Aletsel dönem öncesi büyük depremleri açık kataloglardan toplar, ortak şemaya çevirir,
kataloglar arası aynı olayları gruplar ve her kayda **belirsizliğe dayalı eğitim ağırlığı** yazar.
AFAD / USGS / EMSC aletsel katalogları ayrı modüldedir; bu modül onlarla çakışmaz (`max_year = 1899`).

> **Uyarı.** Tarihsel depremlerin büyüklükleri ölçülmemiştir; hasar dağılımından (makrosismik veri),
> başka büyüklük türlerinden dönüşümle ya da uzman yargısıyla **tahmin edilmiştir**. Tarih, saat ve konum da
> çoğu kayıtta kesin değildir. Bu yüzden her kayıt belirsizliğini ve kaynağını taşır; hiçbir tarihsel değer
> aletsel ölçüm kadar güvenilir sayılmamalıdır.

## Kapsam

- Bölge: enlem 34-44° K, boylam 24-46° D (Türkiye + Ege, Kıbrıs, Kafkasya ve Suriye/Irak sınır bölgesi)
- Dönem: MÖ 2000 - 1899
- Kaynaklar: EPICA v1.1, SHARE Orta-Doğu Türkiye kataloğu, GEM Küresel Tarihsel Katalog v1, NOAA NCEI.
  Ayrıntılar, lisanslar ve atıflar: [KAYNAKLAR.md](KAYNAKLAR.md)

| | Sayı |
|---|---|
| Bütün kayıtlar (yinelemeler dahil) | 1069 |
| Ayrı olay | 897 |
| Mw ≥ 6 olay | 327 |
| Mw ≥ 7 olay | 68 |
| 1000 yılı öncesi olay | 109 |
| Büyüklüğü olmayan (yalnız hasar/şiddet bilgisi) | 301 |

## Dosyalar

| Dosya | İçerik |
|---|---|
| `data/raw/` | Kaynaklardan indirilen ham dosyalar, değiştirilmeden (yedek) |
| `data/tarihsel_tum_kayitlar.csv` | Bütün kaynakların bütün kayıtları. Hiçbir kayıt atılmaz; yinelemeler `duplicate_group` ile bağlıdır |
| `data/tarihsel_katalog.csv` | Her olay için bir kayıt (tercih edilen kaynak) |
| `data/tamlik.json` | Stepp tamlık analizi, kullanılan parametreler, sayımlar |
| `deprem_tarihsel/` | Okuyucular, yineleme eşleştirme, ağırlıklar (`build.py`), parametreler (`config.py`) |

Yeniden üretmek ve test etmek için:

```bash
cd tarihsel
pip install pandas numpy openpyxl xlrd pytest
python -m deprem_tarihsel.build
python -m pytest -q tests
```

Kodda kullanmak için:

```python
from deprem_tarihsel.io import load
hist = load()                    # tercih edilen kayıtlar, time = UTC datetime
core = load(core_only=True)      # aletsel katalogla ortak sütunlar
```

MÖ tarihler `-0281-07-02T12:00:00Z` biçimindedir; pandas ile okurken `format="ISO8601"` verilmelidir
(`io.load` bunu yapar). Nanosaniye çözünürlüklü araçlar 1677 öncesini taşımadığı için her kayıtta ayrıca
sayısal `decimal_year` sütunu vardır.

## Şema

Aletsel katalogla ortak çekirdek (1σ belirsizlikler; kaynak vermediyse boş):
`event_id, time, latitude, longitude, depth_km, mag, mag_type, source="historical", agency,
time_err_s, horiz_err_km, depth_err_km, mag_err`

Tarihsel ekler:

| Sütun | Anlamı |
|---|---|
| `time_precision` | Tarihin bilinen en ince birimi: `year`, `month`, `day`, `hour`, `minute`, `second` |
| `time_min`, `time_max` | Olayın kesinlikle içinde olduğu aralık; `time` bunun orta noktasıdır |
| `time_err_s` | Aralığın tekdüze dağılım 1σ'sı (genişlik/√12). Yalnız yıl biliniyorsa ≈ 105 gün |
| `decimal_year` | Sürekli zaman ekseni |
| `mag_is_estimated` | Büyüklüğü olan her tarihsel kayıtta `True` |
| `mag_orig`, `mag_orig_type`, `mag_method` | Kaynağın verdiği büyüklük ve türü; Ms ise Scordilis (2006) ile Mw'ye çevrildi |
| `intensity_io`, `area`, `country`, `deaths` | Episantr şiddeti, yer adı, ülke, can kaybı (NOAA) |
| `reference` | Birincil kaynak (ör. "Ambraseys and Jackson, 1998") |
| `suspect` | Konumu kuşkulu (ör. NOAA'da bölgeyle uyuşmayan ülke adı) |
| `duplicate_group`, `is_preferred`, `n_sources`, `group_agencies` | Kataloglar arası aynı olay grubu |
| `group_mag_min`, `group_mag_max` | Aynı olay için kaynakların verdiği en küçük / en büyük Mw |

## Ağırlıklar

Ağırlıklar elle verilmez; her kaydın kendi belirsizliklerinden hesaplanır. Bütün katsayılar
`deprem_tarihsel/config.py` içindedir ve geriye dönük testlerle ayarlanmak üzere parametriktir.

```
σM   = sqrt(mag_err² + ((group_mag_max - group_mag_min)/2)²)     # kaynaklar arası tutarsızlık eklenir
σloc = horiz_err_km (yoksa kaynağa göre atanan: 30-50 km)
σt   = time_err_s (gün)

σz   = depth_err_km (yalnızca kaynak verdiyse; yoksa bu terim 0)

w_unc           = max(0.001, 1 / (1 + (σM/0.3)² + (σloc/20 km)² + (σz/30 km)² + (σt/30 gün)²))
w_completeness  = 1, kayıt tamlık penceresindeyse; değilse 0.3
w_suspect       = 0, konum kuşkuluysa; değilse 1
training_weight = w_unc · w_completeness · w_suspect   (yinelenen, tercih edilmeyen kayıtlarda 0)
```

Bileşenler (`w_mag`, `w_loc`, `w_time`, `w_unc`, `w_completeness`, `w_suspect`) ve kullanılan σ değerleri
(`mag_err_used`, `horiz_err_used`, atanıp atanmadığı `*_assigned`) ayrı sütunlarda saklanır; model başka bir
birleştirme denemek isterse yeniden hesaplayabilir. Biçim ve ölçekler aletsel katalog modülündeki
`w = 1/(1+Σ(σ/s)²)` (`src/deprem/catalog/weights.py`, `WeightConfig.scales`) ile aynıdır; iki katalog
aynı ağırlık ölçeğindedir. Yalnızca yıl düzeyinde bilinen tarihler zaman teriminden belirgin ceza alır
(w_time ≈ 0.07); ay düzeyi ≈ 0.92, gün ve daha ince ≈ 1.

Sonuç: tercih edilen tarihsel kayıtların ağırlık ortancası ≈ 0.02, en yükseği ≈ 0.5. İyi belirlenmiş bir
aletsel kaydın ağırlığı (σM≈0.1, σloc≈5 km varsayımıyla ≈0.85) bunun yaklaşık 40 katıdır. Tarihsel veri böylece etkisini korur ama aletsel veriyi bastırmaz.

### Tamlık (Stepp 1972)

Bir dönemde belli bir büyüklüğün üstündeki bütün depremler kayda geçmiş mi? 1899'dan geriye 25 yıllık
pencerelerle olay oranı hesaplanır; oran, referans oranın 2σ altına düştüğü yerde katalog eksik sayılır.
Bu veriyle bulunan başlangıç yılları (bölgenin tamamı için):

| Mw ≥ | Tam kabul edilen başlangıç |
|---|---|
| 5.5 | 1825 |
| 6.0 | 1825 |
| 6.5 | 1800 |
| 7.0 | 1800 |

Bu pencere dışındaki kayıtlar atılmaz, yalnızca ağırlıkları 0.3 ile çarpılır. Tarihsel kayıtlar özellikle
Doğu Anadolu'da ve 1500 öncesinde belirgin biçimde eksiktir; oran tahmini yapan modeller bunu hesaba katmalıdır.

## Yineleme eşleştirme

İki farklı katalogdaki kayıt şu koşullarda aynı olay sayılır: tarih aralıkları çakışıyor veya en çok 15 gün
ayrık (Jülyen/Gregoryen farkı ve kaynak tarih kaymaları için), episantrlar en çok 100 km uzak, büyüklükler
(ikisi de varsa) en çok 1.5 birim farklı. Bir grupta aynı katalogdan iki kayıt olmaz. Tercih sırası:
EPICA > SHARE-CET > GEM-GHEC > NOAA (büyüklüğü olan ve konumu kuşkusuz kayıt önce gelir).

## Bilinen sınırlar

- Tarihler kaynakta nasıl verildiyse öyledir; 1582 öncesi çoğu tarih Jülyen takvimindedir, dönüştürülmedi.
- Saatler çoğu zaman yerel saattir; kaynaklar UTC'ye çevirmiş olabilir veya olmayabilir.
- NOAA'nın antik (1000 öncesi) kayıtlarının birçoğunda büyüklük yoktur ve tarih yalnızca yıl düzeyindedir.
- SHARE-CET'teki tekrar eden 5.42 / 5.82 gibi değerler şiddetten dönüştürülmüştür (`notes` içinde `TMw=cc`).
- Ambraseys (2009), Ambraseys & Finkel (1995), Guidoboni & Comastri (2005) ve KOERI tarihsel kataloğu
  makine-okunur olarak bulunamadı; içerikleri yukarıdaki kataloglar aracılığıyla dolaylı olarak var.
