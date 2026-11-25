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
| Sismik hız değişimi (dv/v) | güçlü | KOERI KO ağı ham kayıtları (FDSN) | tamamlandı: 12 istasyon, 10 günde bir saat, 2012-2026; etki sıfırdan ayırt edilemiyor |
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

## Düzeltme: büyük depremli yıllardan sonra fazla tahmin (Kasım 2026)

**Sorun.** İleriye kayan eğitimde model, eğitim yıllarında ETAS'ın *yıl başındaki* tahminini
taban alıyordu. Bu tahmin o yıl içinde olacak büyük depremlerin artçılarını içermez. 2017,
2020, 2023 gibi yıllarda gözlenen sayı tabandan çok yüksek kalınca model bunu "geçen yıl çok
deprem olan hücrede daha çok deprem olur" diye öğrendi. Ertesi yıl ETAS artçıları zaten hesaba
kattığı için aynı etki iki kez sayıldı: 2018'de 652 depreme 2106, 2021'de 1209 depreme 2286
tahmin edildi (ETAS 1635 ve 1799).

**Düzeltme.** Eğitim yıllarının tabanı artık ETAS'ın gerçekleşen geçmişle hesaplanan koşullu
yoğunluk integrali (yıl içi artçılar dahil). Böylece model yalnız ETAS'ın açıklayamadığı kısmı
öğrenir; deprem deprem puanlamayla da tutarlıdır. Eski davranış `[walkforward] train_offset =
"forecast"` ile seçilebilir.

**Sonuç** (M≥3, 2015-2025 ileriye kayan, ETAS'a göre olay başına kazanç, %95 GA):

| | önce | sonra |
|---|---|---|
| 2018 tahmini (gözlenen 652) | 2106 | 1466 |
| 2021 tahmini (gözlenen 1209) | 2286 | 1548 |
| 2024 tahmini (gözlenen 1033) | 1809 | 1935 |
| Yıl yıl M≥3 | −0,133 [−0,259; 0,081], 6/11 yıl | +0,074 [−0,037; 0,173], 10/11 yıl |
| Yıl yıl M≥4 | −0,223 [−0,335; 0,014] | −0,002 [−0,121; 0,125] |
| Yıl yıl M≥5 | −0,167 [−0,337; 0,149] | +0,051 [−0,082; 0,234] |
| Deprem deprem M≥3,5 | −0,172 [−0,212; −0,122] (daha kötü) | −0,000 [−0,030; 0,032] |
| Deprem deprem M≥4 | −0,207 [−0,259; −0,125] (daha kötü) | −0,034 [−0,084; 0,023] |
| Deprem deprem M≥5 | −0,138 [−0,213; −0,013] (daha kötü) | +0,005 [−0,050; 0,087] |

Model artık ETAS'tan anlamlı olarak kötü değil, ama anlamlı olarak iyi de değil: bütün güven
aralıkları sıfırı içeriyor. Büyük dizilerden sonraki fazlalığın kalanı ETAS'ın kendisinden
geliyor (2024'te ETAS 2449 tahmin etti); bu, ETAS'ın artçı azalmasını büyük dizilerde yavaş
bulmasıyla ilgili ayrı bir konudur.

## Bilinen sınırlamalar

* **GPS kapsaması:** Nevada arşivinde Türkiye istasyonlarının çoğu 2012-2016 arası;
  2017 sonrası istasyon sayısı 20-40'a düşüyor. Sınama yıllarında (2020-2025) GPS
  özellikleri birçok hücrede eksik; model eksik bayrağıyla bunu bilir.
* **Tekrarlayan depremler:** katalog konum hatası (1-5 km) nedeniyle aileler adaydır.
  Gerçek tekrarlayanlar dalga biçimi benzerliğiyle (korelasyon ≥ 0,9) doğrulanmalıdır.
* **dv/v (sonuç):** hafif sürüm (10 günde 1 saat, 12 KOERI istasyonu, 2012-2026) tamamlandı.
  5348 ölçüm noktasının yalnız 613'ü (%11) güvenilir çıktı (korelasyon > 0,5); ölçüm gürültüsü
  istasyon başına %0,2-0,7, yani depremlerden önce beklenen ~%0,1'lik değişimlerden büyük.
  Sadece istasyonlara yakın hücrelerde değer var (hücre-yıllarının %1'inden azı). Eğitilen etki
  +0,08 ± 0,11: sıfırdan ayırt edilemiyor; sismik hız eklenince sınama sonucu değişmedi
  (M≥3: +0,029 → +0,027). İşe yaraması için sürekli kayıt ve istasyon çiftleri arası korelasyon
  gerekir. Çok seyrek gözlemlerde model "eksik" yerine "değer var" bayrağı kullanır; bu bayrağın
  pozitif katsayısı istasyonların zaten aktif bölgelerde kurulu olmasından gelir, öncü sinyal değildir.
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
