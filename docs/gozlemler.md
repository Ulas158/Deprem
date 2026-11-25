# Jeodezik, sismik ve hipotez gözlemleri

Bu belge, deprem oluşumunu (zaman, konum, büyüklük) anlamak için modele eklenen
katalog dışı gözlemleri, nereden geldiklerini, her birinin depreme **başlangıç etkisini**
ve nasıl test edildiğini anlatır. Risk tarafı (zemin büyütmesi, bina stoku) bilerek
kapsam dışıdır. Kuyu suyu ve gazlar ayrı çalışmadadır (`docs/akiskan_ve_dis_etkiler.md`).

> Hiçbir çıktı deprem tahmini değildir. Her gözlem, modelin eğitimde hiç görmediği
> yıllarda kabul görmüş referans model ETAS'a karşı test edilir. Başlangıç etkileri
> literatürden alınmış varsayımlardır, kanıt değildir.

## İki grup

| Grup | Ne demek | Nasıl test edilir |
|---|---|---|
| **Güçlü bilimsel temeli olanlar** | Fiziksel mekanizması ve gözlemleri hakemli literatürde tekrar tekrar gösterilmiş | Yıllık model: 0,5° hücre × 1 yıl, ETAS'ı taban alan Poisson modeline özellik olarak girer |
| **Kanıtı zayıf olanlar (hipotez)** | Öncül sinyal olarak önerilmiş ama bağımsız testlerde tutarlı doğrulanamamış | Kısa vadeli test: 2° kutu × 5 gün; son 10 günün anomalisi sonraki 5 günde bağımsız deprem olasılığını artırıyor mu |

## Gözlemler ve veri kaynakları

| Gözlem | Grup | Kaynak (açık) | Durum |
|---|---|---|---|
| GPS hız alanı, gerinim hızı | güçlü | Nevada Geodetic Lab, Türkiye ve çevresinde ~300 istasyon, günlük | çalışıyor |
| GPS geçici hareket (yavaş kayma) | güçlü | aynı | çalışıyor |
| Tekrarlayan depremler | güçlü | birleşik katalog (AFAD/EMSC/USGS) | çalışıyor (katalog adayı; dalga biçimi doğrulaması isteğe bağlı) |
| Sismik hız değişimi (dv/v) | güçlü | KOERI KO ağı ham kayıtları (FDSN) | pilot: 12 istasyon, 10 günde bir saat |
| Yapay zekâ ile zengin katalog | güçlü | KOERI ham kayıtları + PhaseNet | pilot: Marmara, 1 gün |
| Uydu radarı (InSAR) | güçlü | COMET LiCSAR | **veri bekliyor**: açık arşivde hazır hız haritası yok |
| GONAF kuyu gerinimi | güçlü | GFZ / AFAD | **veri bekliyor**: başvuru ile alınıyor |
| İyonosfer TEC anomalisi | zayıf | CODE (Bern) küresel iyonosfer haritaları, Kp indeksi | çalışıyor |
| ULF manyetik anomali | zayıf | INTERMAGNET İznik (IZN), referans PAG | çalışıyor (dakikalık veri: 2-8 mHz) |
| Termal anomali | zayıf | NASA POWER (MERRA-2) yüzey sıcaklığı | çalışıyor (uydu LST değil, yeniden analiz) |

## Başlangıç etkisi ve eğitim

Her özelliğin `src/deprem/observations/registry.py` içinde bir başlangıç etkisi
(`prior_mean`) ve belirsizliği (`prior_sd`) vardır. Anlamı: özellik 1 standart sapma
artınca deprem oranının logaritması ne kadar değişir. Örnek: GPS gerinim hızı için
+0,30 ± 0,20, yani gerinim hızı yüksek bölgede ETAS'ın beklediğinden yaklaşık %35
fazla deprem varsayılır.

Eğitim bu değerleri Gauss öncül olarak kullanır ve 2012-2019 depremleriyle günceller.
Rapor her özellik için **başlangıç → öğrenilen** değeri ve kalan belirsizliği gösterir.
Ayrıca aynı model literatür öncülü olmadan (öncül = 0) da eğitilir; böylece sonucun
ne kadarının varsayımdan, ne kadarının veriden geldiği ayrılır.

Zayıf kanıtlı hipotezlerin başlangıç etkileri bilerek küçük (≤ 0,10) tutuldu:
literatürde iddia edilen büyük etkiler bağımsız testlerde tekrarlanamadı.

## Çalıştırma

```bash
pip install -e ".[obs]"
python scripts/fetch_catalogs.py                       # katalog
python scripts/fetch_observations.py all              # GPS, TEC, ULF, termal, tekrarlayan
python scripts/fetch_observations.py dvv              # dalga biçimi, saatler sürer
python scripts/fetch_observations.py mlcat --mlcat-days 1
python scripts/train.py --config configs/train_obs.toml
```

Rapor: `reports/training_turkey_v2_obs.md`. Harita katmanları: `web/data/layers/`
altında "Jeodezi ve sismik gözlemler" grubu.

## Bilinen sınırlamalar

* **GPS kapsaması:** Nevada arşivinde Türkiye istasyonlarının çoğu 2012-2016 arası;
  2017 sonrası istasyon sayısı 20-40'a düşüyor. Sınama yıllarında (2020-2025) GPS
  özellikleri birçok hücrede eksik; model eksik bayrağıyla bunu bilir.
* **Tekrarlayan depremler:** katalog konum hatası (1-5 km) nedeniyle aileler adaydır.
  Gerçek tekrarlayanlar dalga biçimi benzerliğiyle (korelasyon ≥ 0,9) doğrulanmalıdır.
* **dv/v:** hafif sürüm (10 günde 1 saat, 12 istasyon). Tam yöntem sürekli kayıt ve
  istasyon çiftleri arası korelasyon ister.
* **ULF:** dakikalık veri klasik ULF bandını (10-100 mHz) içermez; bu hipotezin
  zayıflatılmış bir sınamasıdır.
* **Termal:** NASA POWER bir yeniden analiz ürünüdür; gerçek uydu kızılötesi (MODIS)
  NASA Earthdata hesabı gerektirir.
* **Yapay zekâ katalogu (pilot sonucu):** Marmara'da 1 Mart 2024 günü, 10 KOERI istasyonunda
  PhaseNet 513 faz seçimi yaptı; basit ızgara birleştirmesi 16 olay buldu (resmi katalogda
  aynı kutu ve günde 8 olay; oran 2,0). Resmi olayların yalnız 3'ü eşleşti; birleştirme ve
  konumlama kaba olduğundan bu sayı gerçek bir iyileşme kanıtı değildir. Eğitime girmesi için
  aylarca kesintisiz işleme ve düzgün konumlama (ör. GaMMA + NonLinLoc) gerekir.
* **InSAR:** hazır hız haritası açık arşivde yok. Ham interferogramlardan LiCSBAS ile
  üretmek çerçeve başına onlarca GB ve saatlerce hesap ister. Üretilen ya da
  yayımlanmış (ör. Weiss vd. 2020) hız GeoTIFF'leri `data/raw/insar/` altına konunca
  eğitim hattı otomatik okur.
* **GONAF:** veri GFZ'den başvuruyla alınır; biçimi `src/deprem/observations/gonaf.py`
  içinde tarif edildi.
