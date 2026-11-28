# Marmara fayı için fiziksel kırılma modeli: sonuçlar

Tarih: Kasım 2026. Bu rapor, Kuzey Anadolu Fayı'nın Marmara kolu için kurulan fiziksel (mekanik) modelin iki turunu, sınav sonuçlarını ve bugünkü durumunu sade bir dille özetler.

## Kısa sonuç

- **Model sınavı geçemedi.** Ayardan sonra gelen beş büyük depremin (1912 Mürefte-Şarköy, 1999 İzmit, 1999 Düzce, 2022 Düzce-Gölyaka, 2025 Silivri açığı) hiçbirini sıkı ölçütle tutturamadı. Deprem olan yerlerde model doluluğu ortalama %50 civarındaydı; "tuttu" sayılmak için %90-98 gerekiyordu.
- Büyük depremlerin yerini sıralamada model şansa çok yakın: ortalama sıra 0,54 (şans 0,50, p ≈ 0,38, yalnız 5 olay).
- Küçük depremlerde (2018-2023, Silivri ve Çınarcık) model şansın üstünde: sıra 0,60-0,66. Ama sıkı ölçütle tutma oranı yalnız %0,5-1. "Son deprem yeni olduysa yine olur" diyen basit kümelenme kuralı da en az model kadar iyi.
- 1895-2026 arasında model, gerçekte olmamış ortalama 1,5 büyük (M6,8+) deprem üretti.
- **Bu yüzden İstanbul için verilen sayılar tahmin değildir.** Modelin bugünkü iç durumunu gösterir. Proje kuralı gereği (kanıtsız hipotez doğrulanmış sayılmaz), bu model şu hâliyle İstanbul depremi için olasılık kaynağı olarak kullanılmamalıdır.

## Model ne yapıyor

Fay yüzeyi kutulara bölündü:

- kaba ağ 2 km × 2 km, 5.415 kutu;
- Silivri ve Çınarcık için ayrıca 200 m'lik yerel ağlar.

Kuzey kolun yanında güney kol (Biga-Bandırma-Gemlik-İznik-Geyve), İmralı ve Armutlu kolları da modelde. Modele giren fiziksel etkiler:

- Levhalar her yıl birkaç milimetre itiyor (Marmara'da 12,5-17 mm arası serbest), kilitli kutularda gerilme birikiyor.
- Bir kutu kırılınca, gerilmesini komşularına ve başka kollara aktarıyor (Coulomb gerilme aktarımı).
- Bazı kesimler yavaşça kayıyor (krip); bu kesimler gerilme biriktirmiyor.
- Kırılgan katmanın kalınlığı bölgeye göre değişiyor (8-24 km).
- Altındaki sıcak, sünek kuşak sıcaklığa bağlı bir hızla akıyor ve üsttekini itiyor.
- Marmara dışındaki büyük depremlerin etkisi de eklendi (1855 Bursa'dan 2014 Gökçeada'ya 9 deprem).

## Nasıl ayarlandı (ikinci tur, Ulaş'ın yöntemi)

1. Model, gerçek depremleri zaman sırasıyla "yaşadı". Her gerçek depremde, depremin olduğu yer (konum belirsizliği içinde, olasılıkla ağırlıklı) kırılma sınırına çekildi, kırıldı ve gerilmesini aktardı. 0,1 bar üstü etki yapan her olay kullanıldı. Küçük depremler GFZ, AFAD, EMSC ve USGS kataloglarından alındı: 19.176 GFZ kaydı ve öncesi.
2. 12 parametrenin hepsi serbest bırakıldı. Başlangıç değerleri literatürden alındı.
3. Ayar sayısına sınır konmadı. Arama, üst üste 3 turda %1'den az iyileşince durdu:

| Aşama | Dönem | Denenen kombinasyon |
|---|---|---|
| Küçük depremler, Silivri | 2006-2017 | 500 |
| Küçük depremler, Çınarcık | 2006-2017 | 450 |
| Büyük depremler, tüm fay, geri sararak | 1719-1894 | 402 |

4. Sınav dönemi ayara hiç girmedi. En iyi 5 ayar donduruldu, her biri 3 farklı rastgele başlangıçla sınavdan geçirildi.

"Tuttu" sayılmak için üçü birden gerekti:

- **Doğru yer.**
- **Yeterli doluluk:** tarihsel depremde %90, erken aletsel dönemde %95, modern dönemde %98.
- **Doğru zaman:** model depremi kendiliğinden ±20, ±10 ya da ±5 yıl içinde üretmeli; küçük depremlerde aralığın ±%25'i.

## Sonuçlar

### Büyük depremler (kaba ağ, 15 koşu)

| Deprem | Bölüm | Model doluluğu | Tuttu | Model sırası | Basit "son depremden beri geçen süre" sırası |
|---|---|---|---|---|---|
| 1766 Mayıs Doğu Marmara | ayar | %47 | 0/15 | 0,48 | 0,25 |
| 1766 Ağustos Ganos | ayar | %53 | 0/15 | 0,68 | 1,00 |
| 1894 Adalar-İzmit Körfezi | ayar | %47 | 0/15 | 0,48 | 0,35 |
| 1912 Mürefte-Şarköy | sınav | %48 | 0/15 | 0,56 | 0,58 |
| 1999 İzmit | sınav | %51 | 0/15 | 0,63 | 0,86 |
| 1999 Düzce | sınav | %49 | 0/15 | 0,28 | 0,90 |
| 2022 Düzce-Gölyaka | sınav | %50 | 0/15 | 0,65 | 0,03 |
| 2025 Silivri açığı | sınav | %58 | 0/15 | 0,59 | 0,06 |

Notlar:

- Sıra, deprem olan yerin modeldeki doluluğunun, aynı derinlikteki 200 rastgele kilitli kutuya göre yüzdelik konumudur. 0,5 şans demek.
- İki ek sürüm aynı sonucu verdi:
  - konum ağırlığını gözlemle güncellemeyen sürümde sıra 0,55;
  - depremin tüm kırılma yamasını %100'e çeken sürümde sıra 0,49.
- Sahte büyük depremler (gerçekte olmamış M6,8+): ayar döneminde ortalama 1,6-1,9, sınav döneminde 1,5.

### Küçük depremler (yerel 200 m ağlar, 15 koşu)

| Ağ | Sınav olayı (2018-2023) | Model sırası | Basit süre kuralı | Sıkı tutma oranı |
|---|---|---|---|---|
| Silivri (2019 dizisi dahil) | ~735 | 0,66 | 0,28 | ~%1 |
| Çınarcık | 46 | 0,60 | 0,42 | ~%0,5 |

Notlar:

- Basit süre kuralının sırası 0,5'in çok altında. Yani bu bölgede depremler "uzun süre sessiz kalan yerde" değil, "yakında deprem olmuş yerde" oluyor (kümelenme). Bu kuralın tersi (≈0,72) Silivri'de modelden iyi.
- p-değerleri olayları bağımsız sayıyor. Kümelenme yüzünden gerçek anlamlılık, tabloda görünenden daha zayıf.

### Birinci tur (karşılaştırma için)

- Birinci tur tek bir kalın kırılgan katman senaryosuyla, 3 koşuyla yapıldı (orta senaryo, makine yeniden başlayınca kayboldu).
- 8 tarihsel depremde ortalama sıra 0,55, şansa göre p ≈ 0,12.
- Model 800 yılda ortalama 11 sahte büyük deprem üretti.
- O turda İstanbul önü doluluğu %8-20 çıktı, 30 yıllık kırılma payı 0.

İkinci turdaki doluluğun daha yüksek çıkması yöntem farkından geliyor: gerçek depremler modele zorla işlendi. Bu fark, modelin daha doğru olduğunu göstermiyor.

## Bugünkü durum (2026 sonu) ve İstanbul özeti

Aşağıdaki değerler sınavı geçemeyen modelden geliyor. Tahmin olarak okunmamalı.

| Kesim | Boylam | Model doluluğu (koşu aralığı) | Model, 30 yılda M6,8+ kırılma payı (koşu aralığı) |
|---|---|---|---|
| Orta Marmara çukuru | 27,75-28,05 | %60 (%47-74) | %2 (%0-16) |
| Kumburgaz | 28,05-28,60 | %47 (%32-68) | %6 (%0-42) |
| Adalar | 28,60-29,10 | %55 (%47-66) | %7 (%0-42) |
| İstanbul önü toplam | 27,75-29,10 | %53 (%43-65) | %8 (%0-42) |

- Tüm fay sistemi için model 30 yılda ortalama 0,5 tane M6,8+ deprem üretiyor. Koşudan koşuya bu sayı 0 ile 2,6 arasında değişiyor.
- Bu kadar geniş aralık, modelin sonucunun ayarın küçük farklarına ve rastgele başlangıca çok bağlı olduğunu gösteriyor.
- Tarihte (1719-2026) bu bölgede 30 yılda yaklaşık 1 tane M6,8+ deprem olmuş.

Dosyalar:

- [km km doluluk ve 30 yıllık pay, tüm kollar](tables/mekanik2_bugun_km.csv)
- [İstanbul özeti](tables/mekanik2_istanbul.csv)
- [büyük deprem sınavı, olay olay](tables/mekanik2_buyuk_sinav.csv), [sahte depremler](tables/mekanik2_sahte.csv)
- [küçük deprem sınavı](mekanik2_kucuk_sinav.csv)
- 3B görünüm: `reports/mekanik_3b/` (yayınlanmış sayfa: https://claude.ai/artifact/JqM2BMZJzfNBaPc9Ve51Vt)

## Plandan sapmalar (açıkça)

- **Yerel ağlar 100 m yerine 200 m.** 100 m'de tek deneme 15 dakika sürüyordu. Ulaş onayıyla 200 m'ye geçildi.
- **Viskoz kuşakta zaman adımı alt sınırı (0,05 yıl).** Simülatörü 22 kat hızlandırdı, sonuç farkı ihmal edilebilir düzeyde.
- **Marmara dışı depremler yaklaşık hesaplandı.** Dikdörtgen kaynak kullanıldı ve yakın alan etkisi ±1 MPa ile kırpıldı. Mudurnu 1967'nin ham etkisi 4 MPa'ydı.
- **Kesik deneme sayısı.** Büyük aramada 6 deneme bir hesap hatasıyla çöktü. Kırılma çözücüsünde son turda kutu listesi ile kayma listesi boy tutmuyordu. Hata düzeltildi ve 6 deneme yeniden koşuldu; hiçbiri en iyi sonucu geçmedi, sonuç değişmedi.
- **30 yıllık bakışta sayım hatası.** İlk sürüm, birden çok kola yayılan kırılmanın km aralığını kollar arasında karıştırıyordu. Düzeltildi; tablolar düzeltilmiş sürümden.
- Kaba ağ 2 km'lik. Tam 1 km ağla doğrulama yapılmadı.

## Sınırlar

- Büyük deprem sayısı çok az: ayarda 3 tarihsel, sınavda 5 olay. Böyle bir az veriyle anlamlı bir başarı göstermek zaten zor. Ama model şansı bile belirgin biçimde geçemedi.
- Tarihsel depremlerin yeri ve büyüklüğü belirsiz. Ayar bu belirsizliği konum adaylarıyla taşıdı, yine de kırılma uzunlukları kaba.
- Model kutularının dayanımı rastgele dağılıyor. Aynı ayarla farklı rastgele başlangıç, çok farklı gelecek veriyor (30 yılda 0 ile 2,6 deprem).
- Küçük deprem başarısı büyük ölçüde kümelenmeden geliyor olabilir. Kümelenmeyi hesaba katan bir anlamlılık testi (blok önyükleme) henüz yapılmadı.

## Ne yapılabilir

Ayar sayısını artırmak sonucu değiştirmeyecek gibi: arama durma kuralına ulaştı. Daha umutlu yollar şunlar:

- kırılma kuralını değiştirmek (gerçek depremlerin yayılma boyunu tutturacak dinamik zayıflama);
- 1 km ağla doğrulama;
- küçük depremler için kümelenmeyi (art sarsıntı) ayrıca modelleyip fiziksel modelin bunun üstüne bir şey katıp katmadığını sınamak.
