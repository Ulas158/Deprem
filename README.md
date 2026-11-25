# Deprem

Depremlerin öncesinde ortaya çıkabilecek fiziksel değişimleri, veri kaynakları arasındaki gizli ilişkileri
ve tekrar eden örüntüleri araştıran; zaman, konum ve büyüklük için olasılıksal tahminler üretip bunları
geçmiş ve yeni verilerle bilimsel olarak doğrulamayı hedefleyen bir araştırma projesi.

Günümüzde depremlerin zamanını güvenilir biçimde önceden belirleyen, doğrulanmış bir yöntem yoktur.
Bu depodaki hiçbir hipotez yeterli kanıt olmadan doğrulanmış kabul edilmez.

## Araştırma kodu (Python)

Bu depo, depremlerin olasılıksal olarak öngörülebilirliğini araştırmak için
veri, istatistik ve test altyapısı içerir. **Buradaki hiçbir çıktı deprem
tahmini değildir.** Bütün sonuçlar betimsel istatistik veya ileride
ileriye dönük olarak sınanacak hipotez değişkenleridir. Bir yöntem ancak
eğitimde görmediği bir dönemde, ETAS gibi kabul görmüş referans modellerden
istatistiksel olarak daha iyi olduğunu gösterdiğinde anlamlı sayılır.

### Hızlı başlangıç

```bash
pip install -e .[dev]
python scripts/fetch_catalogs.py     # USGS + EMSC + AFAD (ISC arşivi ve son 2 hafta)
python scripts/run_analysis.py       # rapor, grafikler, fay elemanları, özellik tablosu
pytest
```

Ham ve işlenmiş veri `data/` altına yazılır ve depoya girmez (yeniden
üretilebilir). Rapor: [`reports/turkey_summary.md`](reports/turkey_summary.md).

### Model eğitimi (tek komut)

```bash
python scripts/fetch_catalogs.py                                   # katalog (önbellekli)
python scripts/train.py --config configs/train_default.toml        # tüm aşamalar
python scripts/train.py --force fit evaluate report                # yalnızca seçilen aşamalar
```

Aşamalar sırayla `etas → etas_fc → coulomb → budget → features → fit → evaluate → report`.
Her aşama çıktısını `runs/<run_name>/` altına yazar ve `status.json`'a işler; komut
yeniden çalıştırıldığında tamamlanan aşamalar atlanır (kaldığı yerden devam). Yapılandırma
`runs/<run_name>/config.toml` olarak kopyalanır. Rapor `reports/training_<run_name>.md`,
harita katmanları `web/data/layers/` altına yazılır.

Model: `λ(hücre, yıl) = Λ_ETAS · exp(β0 + Σ β_k z_k)`. Özellikler (sismisite, fay, Coulomb
ΔCFF, eleman kayma açığı) yalnızca pencere başından önceki veriden hesaplanır. Katsayılar
Gauss öncüllü Poisson olabilirliğiyle eğitim yıllarında öğrenilir, efektif sürtünme μ'
eğitimin son iki yılında doğrulamayla seçilir. Başarı, eğitimde görülmemiş sınama
yıllarında ETAS'a karşı M3/M4/M5/M6 eşiklerinde ayrı ayrı ölçülür.

Yeni bir deneme için yapılandırmayı kopyalayıp `run_name`'i değiştirin.

### Jeodezik, sismik ve hipotez gözlemleri

```bash
pip install -e ".[obs]"
python scripts/fetch_observations.py all                 # GPS, iyonosfer, ULF, termal, tekrarlayan depremler
python scripts/fetch_observations.py dvv mlcat           # ham sismogramlar (yavaş; pilot)
python scripts/train.py --config configs/train_obs.toml  # gözlemlerle eğitim + katkı ölçümü
```

GPS (gerinim hızı, geçici hareket), tekrarlayan depremler, ortam gürültüsünden sismik hız
değişimi, yapay zekâ ile zenginleştirilmiş katalog, InSAR ve GONAF kuyu gerinimi **güçlü**
grupta yıllık modele girer. İyonosfer TEC, ULF manyetik ve termal anomaliler **zayıf
kanıtlı hipotez** olarak ayrı, kısa vadeli (5 günlük) bir testte sınanır. Her özelliğin
literatüre dayalı başlangıç etkisi ve belirsizliği `src/deprem/observations/registry.py`
içindedir; eğitim bunları geçmiş depremlerle günceller. Ayrıntı: `docs/gozlemler.md`.

### En eski veriden bugüne deprem deprem eğitim

```bash
python scripts/train_long.py --config configs/train_long.toml   # ~1 saat (4 çekirdek)
```

Tarihsel katalog (1700-1899 zaman çizelgesine, daha eskisi arka plan haritasına) ve aletsel
katalogla 1900'den bugüne yıl yıl ve deprem deprem ilerler. Katalog her dönemde farklı
büyüklükten eksiksiz olduğu için üç koşu vardır (M≥5,5 1900+, M≥4,5 1975+, M≥4 1995+).
Eski kayıtlar eğitimde düşük ağırlık alır; güven aralıkları, katalogun tarih/konum/büyüklük
hatalarıyla 10 kez yeniden örneklenip yürüyüşün baştan yapılmasıyla hesaplanır, bu yüzden
eski dönemlerde daha geniştir. Rapor: `reports/training_turkey_v4_long.md`. Yöntem:
`src/deprem/longwalk.py`.

### Veri kaynakları

| Kaynak | Kapsam | Not |
|---|---|---|
| USGS ComCat | 1900+, M≥4 | Uzun geçmiş, Mw bilgisi |
| EMSC | 1998+, M≥2 | Türkiye olayları büyük ölçüde KOERI/AFAD çözümleri |
| AFAD (ISC Bülteni üzerinden) | 2010–2025, tüm büyüklükler | AFAD'ın ISC'ye gönderdiği kendi çözümleri; patlamalar ayıklanır |
| AFAD açık API | yalnızca son ~2 hafta | Düzenli çekilerek arşiv biriktirilmeli |
| AFAD web kataloğu (elle) | tam arşiv | Oturum gerektirir; dışa aktarılan dosya `--afad-export` ile okunur |
| GEM Küresel Aktif Faylar (EMME) | ~430 segment | Kayma hızı, eğim, rake, sismojenik derinlik |

AFAD web kataloğundan elle indirme: deprem.afad.gov.tr/event-catalog sayfasında
tarih aralığını ve Türkiye sınırlarını seçip sonuçları CSV/Excel olarak indirin,
dosyaları `data/raw/afad_export/` altına koyun ve
`python scripts/fetch_catalogs.py --afad-export data/raw/afad_export/*.csv` çalıştırın.

### Paket yapısı

```
src/deprem/
  region.py              Türkiye kutusu ve alt bölgeler
  catalog/
    schema.py            ortak katalog şeması (+ belirsizlik alanları)
    sources.py           USGS, EMSC, AFAD, ISC istemcileri (önbellekli)
    clean.py             temizleme, kaynaklar arası tekilleştirme, Gardner-Knopoff
    weights.py           belirsizliğe dayalı eğitim ağırlıkları
  stats/gr.py            Mc (MAXC, GFT), b-değeri (Aki-Utsu), önyükleme
  analysis/seismicity.py zaman serileri, ızgara sayımları, b haritası
  faults/
    gem.py               GEM aktif fay veritabanı
    moment.py            moment birikimi/boşalımı, olay-fay geometrisi
    mesh.py              fayları 1 km × 1 km kabuk elemanlara bölme
    parameters.py        eleman parametre kaydı (51 parametre, boş alanlar dahil)
  features/
    table.py             hücre × zaman penceresi özellik tablosu (kayıt/dekoratör)
    evaluate.py          kalibrasyon eğrisi, Brier, eşikler arası taşınma
  forecast/
    base.py              ızgara tahmin arayüzü, referans modeller
    tests.py             CSEP N/L testleri, bilgi kazancı
    latent.py            10 gizli parametre ve en iyileme
  external/              gaz, sıvı ve dış etkiler (bkz. docs/akiskan_ve_dis_etkiler.md)
    astro.py             Ay-Güneş konumu, gelgit ve mevsim fazları
    periodic.py          Schuster testi (gelgit, mevsim)
    swarms.py            deprem sürüleri, hidrolik yayınım
    dynamic.py           uzak büyük depremlerin dinamik tetiklemesi
    sites.py             baraj, maden, jeotermal, petrol/gaz sahaları (yıllarıyla)
    params.py            35+ yeni eleman parametresi, 8 gizli parametre, veri kaynakları
    elements.py          her fay elemanına aktarım (istasyonlardan IDW dahil)
    features.py          hücre × pencere tablosuna yeni özellik grupları
    layers.py            harita katmanları
  faults/fill.py         eksik eleman değerlerini fay boyunca doldurma (+ bayrak, belirsizlik)
  observations/          katalog dışı gözlemler (bkz. docs/gozlemler.md)
    registry.py          özellik listesi, kanıt grubu, başlangıç etkisi ± belirsizlik
    gnss.py              GPS zaman serileri, gerinim hızı, geçici hareket
    repeaters.py         tekrarlayan deprem aileleri, derin kayma hızı
    dvv.py               ortam gürültüsünden sismik hız değişimi
    mlcatalog.py         PhaseNet ile zengin katalog (pilot)
    insar.py, gonaf.py   InSAR hız ve kuyu gerinimi okuyucuları
    ionosphere.py        TEC anomalileri (hipotez)
    geomag.py            ULF manyetik anomaliler (hipotez)
    thermal.py           yüzey sıcaklığı anomalileri (hipotez)
    shortterm.py         kısa vadeli hipotez testi
    fit.py               öncüllü Poisson (MAP + Laplace belirsizliği)
  plotting.py
```

Gaz/sıvı ve dış etkiler için: `python scripts/run_external.py` (rapor:
[`reports/external_summary.md`](reports/external_summary.md), harita katmanları `web/data/layers/ext_*`).
Saha listesi `external_data/` altındadır.

### Değerlendirme ilkeleri

1. **Zamana göre ayırma.** Eğitim < t₀ ≤ sınama; rastgele bölme yok. Özellikler
   yalnızca pencere başlangıcından önceki veriden hesaplanır.
2. **Referansa karşı.** Her model düzleştirilmiş sismisite ve (sonraki PR'da) ETAS ile
   aynı ızgarada, Poisson log-olabilirliği ve bilgi kazancıyla karşılaştırılır.
3. **Eşik eşik.** Küçük depremlerle (M3) kalibre edilen olasılıklar M4, M5, M6
   eşiklerinde ayrıca puanlanır (`transfer_report`).
4. **Ağırlıklar belirsizlikten.** Küçük depremler küçük oldukları için değil,
   yalnızca ölçüm belirsizliği büyükse düşük ağırlık alır.

### Yol haritası

- [x] Katalog, GR istatistikleri, zaman/mekân analizi
- [x] Fay veritabanı, moment bütçesi, kabuk eleman ağı ve parametre kaydı
- [x] Özellik tablosu, kalibrasyon, gizli parametre en iyileme iskeleti
- [ ] Coulomb gerilme aktarımı (Okada) eleman ağı üzerinde
- [ ] ETAS referans modeli ve stokastik artçı ayıklama
- [ ] GNSS gerinim, krip, ısı akısı vb. boş parametrelerin bağlanması
- [ ] Tarihsel katalogla (1900 öncesi) birleştirme
- [ ] İleriye dönük (prospektif) test kaydı: tahminler önceden zaman damgalı saklanır
