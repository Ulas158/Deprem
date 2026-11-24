# Canlı deprem haritası (`web/`)

Türkiye ve çevresindeki depremleri anlık gösteren, derleme adımı gerektirmeyen statik bir site.

![önizleme](docs/preview.jpg)

## Neler var

- **Canlı akış:** EMSC WebSocket (`standing_order`) üzerinden yeni depremler saniyeler içinde düşer;
  bağlantı koparsa kendini yeniden bağlar. Ek olarak her 60 sn'de EMSC FDSN yoklaması yapılır.
  EMSC'ye ulaşılamazsa USGS FDSN'e geçer. EMSC kataloğu AFAD ve KOERI çözümlerini de içerir.
- **Depremler:** büyüklüğe göre boyut, geçen süre veya derinliğe göre renk, son 3 saattekiler için dalga animasyonu,
  yeni deprem bildirimi, ayrıntı penceresi, zaman–büyüklük çizelgesi, 24 sa / 7 g / 30 g ve en küçük büyüklük filtresi.
- **Aktif faylar:** GEM Global Active Faults veritabanından bölgedeki 558 segment;
  mekanizmaya (sağ/sol yanal, normal, ters, verev) ya da kayma hızına (mm/yıl) göre renk, kayma hızına göre kalınlık.
- **Taban harita:** Natural Earth sınırları yerel GeoJSON olarak gömülü; harici karo sunucusu veya API anahtarı gerekmez.
- Mobil uyumlu (alttan açılan panel).

> Bu harita **tahmin yapmaz**; yalnızca kaydedilmiş depremleri ve bilinen fayları gösterir.

## Yerelde çalıştırma

```bash
cd web
python3 -m http.server 8000
# http://localhost:8000
```

(ES modülleri `file://` ile açılmaz; bir yerel sunucu gerekir.)

## Yayınlama (GitHub Pages)

`.github/workflows/pages-web.yml`, `main` dalına `web/` altında bir değişiklik geldiğinde siteyi yayınlar.
Bir kez: depo **Settings → Pages → Build and deployment → Source: GitHub Actions** seçilmeli.

## Dosyalar

| Dosya | İçerik |
|---|---|
| `index.html`, `css/style.css` | Arayüz |
| `js/feeds.js` | EMSC/USGS veri çekme ve canlı WebSocket, ortak olay biçimi |
| `js/layers.js` | Harita stili, renk ölçekleri ve **katman kayıt defteri** (`OVERLAYS`) |
| `js/app.js` | Durum, çizim, panel, zaman çizelgesi, bildirimler |
| `data/*.geojson` | Faylar, ülke sınırları, şehirler (`tools/build_data.py` ile üretilir) |

Verileri yeniden üretmek için:

```bash
curl -LO https://raw.githubusercontent.com/GEMScienceTools/gem-global-active-faults/master/geojson/gem_active_faults_harmonized.geojson
curl -LO https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_0_countries.geojson
curl -LO https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_populated_places_simple.geojson
python3 web/tools/build_data.py gem_active_faults_harmonized.geojson ne_10m_admin_0_countries.geojson ne_10m_populated_places_simple.geojson
```

## Yeni bilimsel katman ekleme

`js/layers.js` içindeki `OVERLAYS` dizisine bir giriş eklemek yeterli; panelde aç/kapat düğmesi otomatik çıkar.
`probe` alanı olan katmanlar, dosya `data/` altında bulunduğu anda kendiliğinden etkinleşir:

- **`data/moment_deficit.geojson`** — fay başına moment açığı. LineString özellikleri;
  renk `mw_equiv` ya da `mw_equiv_since_ref` (birikmiş momentin eşdeğer Mw değeri) alanından gelir.
  Üzerine gelince `fault_id`/`name`, `slip_rate_mm_yr`, `length_km`, `m0_since_ref`, `last_big_event` gösterilir;
  dosyanın `metadata.note` metni katman açıklaması olarak kullanılır. 1e15'ten büyük sayılar metne çevrilerek
  yüklenir (MapLibre karo kodlaması bu büyüklükteki değerlerde katmanı çizmiyor).
- **`data/coulomb_stress.geojson`** — Coulomb gerilme değişimi. Polygon özellikleri; renk `dcfs_bar` alanından gelir.

Python analiz hattı bu dosyaları yazdığında harita onları ek kod gerektirmeden gösterir.

## Kaynaklar ve lisanslar

- Depremler: [EMSC](https://www.seismicportal.eu) (FDSN + WebSocket), yedek [USGS](https://earthquake.usgs.gov/fdsnws/event/1/)
- Faylar: [GEM Global Active Faults](https://github.com/GEMScienceTools/gem-global-active-faults), CC BY-SA 4.0
- Sınırlar ve şehirler: [Natural Earth](https://www.naturalearthdata.com), kamu malı
- Harita motoru: [MapLibre GL JS](https://maplibre.org) 4.7.1 (BSD-3)
