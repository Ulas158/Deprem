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
  plotting.py
```

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
