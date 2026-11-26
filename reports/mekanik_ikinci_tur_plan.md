# Mekanik model, ikinci tur: önceden sabitlenen plan

Yazılma: 2026-11-26 15:00 UTC. Birinci tur sınavının (oynatma) sonuçları **bitmeden ve bu plan
için kullanılmadan** sabitlendi. Aşağıdaki aralıklar, başlangıç değerleri ve sınav bölmeleri
birinci tur sonucuna göre değiştirilmeyecek.

## Yöntem: veri özümseme
Her gerçek depremde kırılan kilitli parça o an %100 doluya çekilir, ardından boşaltılır ve
yükü elastik etkileşimle bütün komşulara aktarılır (`deney.zorla(..., doldur=True)`,
`mekanik_model.py oynatma --ozumse`). Model gerçek tarih boyunca ilerler; her gerçek deprem
sırası gelince işlenir. Ölçüt: bir sonraki gerçek depremden hemen önce kırılacak yerin
doluluğu fayın geri kalanına göre yüksek mi (AUC), şans düzeyine göre (kırılmaları fay boyunca
rastgele kaydırma) p-değeri.

## Sıra
1. Küçük depremler (GFZ hassas katalogları: Becker 2023/2025, Bentz 2020, Martínez-Garzón 2022).
   - 100 m'lik yerel ağ: Silivri önü (2019 dizisi) ve Çınarcık.
   - Büyük ağda 10 km dilimlerde sıklık ve derinlik uyumu; sessiz dilimlerde model küçük deprem üretmemeli.
   - Ayar 2006-2017, sınav 2018-2023 (Silivri dizisinde dizinin ilk yarısı / ikinci yarısı).
2. Büyük depremler: ayar 1509-1894 (7 olay), sınav 1912 Mürefte, 1999 İzmit, 1999 Düzce.

## Serbest parametreler ve arama aralıkları (11)
| parametre | aralık | başlangıç (tüm fayda ortak) | dayanak | birinci tur |
|---|---|---|---|---|
| etkin sürtünme mu_eff | 0,2 - 0,6 | 0,4 | King, Stein & Lin 1994; Coulomb çalışmalarının ortak değeri | 0,4 |
| dayanım fazlası dtau0 | 5 - 50 MPa | 10 MPa | ölçülen gerilme düşüşleri 1-10 MPa (Allmann & Shearer 2009) | 30 MPa |
| dayanım rastgeleliği sigma_h | 0,05 - 0,5 | 0,3 | büyük depremlerde kayma dağılımının pürüzlülüğü (Mai & Beroza 2002) | 0,1 |
| dinamik zayıflama eta_dyn | 0 - 0,6 | 0,3 | hızlı kaymada sürtünme düşüşü deneyleri (Di Toro vd. 2011) | 0,2 |
| deprem sonu kalan seviye | -%20 ile +%50 | -%10 | dinamik aşırı kayma %10-20 (Madariaga 1976) | -%10 |
| zayıflama kayma ölçeği kayma_ref | 0,2 - 4 m | 1 m | büyük depremlerde kritik kayma 0,1-1 m (Tinti vd. 2005) | 2 m |
| sünek Maxwell zamanı tM_ref | 1 - 50 yıl | 10 yıl | 1999 sonrası viskozite ~1e19 Pa·s (Hearn vd. 2009; Ergintav vd. 2009) | 3 yıl |
| krip çarpanı | 0,3 - 5 | 1 | laboratuvar a-b değerleri olduğu gibi (Blanpied vd. 1995) | 2 |
| normal gerilme artışı sigma_grad | 1 - 17 MPa/km | 5 MPa/km | fay zonunda yüksek gözenek basıncı (Rice 1992); 17 = hidrostatik | 2,5 MPa/km |
| kırılgan kalınlık (AFAD→EMSC oranı u) | 0 - 1 | 0,5 | iki kurumun ortası | 0,5 / 1 |
| Orta Marmara sığ krip derinliği | 0 - 8 km | 3 km | deniz tabanı ölçümleri ve tekrarlayan depremler (Schmittbuhl vd. 2016; Yamamoto vd. 2019) | 0-6 km |

Başlangıç değerleri birinci tur ayarından değil, yalnız literatürden alındı. Başta her parametre
bütün fay hücreleri için aynıdır; bölgesel farklar (10 km blok) yalnız ayar verisi gerektirirse eklenir.

Sabit (ölçümle bilinen): levha hızı (GPS), fay yeri ve yönü, elastik yük aktarım matrisi.

## Karşılaştırma
- A: 11 parametre serbest. B: en etkili 4-5 parametre (duyarlılık analiziyle, yalnız ayar
  verisiyle seçilir), diğerleri birinci tur değerinde.
- İkisinin de yalnız bağımsız sınav sonuçları yan yana raporlanır.
- Arama: önce kaba ve hızlı model (2 km kutu, kısa dönemler) ile Latin hiperküp + yerel
  iyileştirme, en iyi 5 aday tam modelle doğrulanır.

## Saatin başlangıcı (Ulas'ın tarifi: yeri ve tarihi kesin olan en eski depremden ileri sar)
- Büyük depremler: tarihler hepsinde kesin, ama yer (kırılma boyu) 1509, 1659 ve 1754'te
  tartışmalı (güven "düşük"). Saat **1719 İzmit** ile başlar (yer güveni "orta"). 1509 ve 1659'un
  bıraktığı yük bilinmediği için 1719'daki başlangıç durumu rastgele ısınma durumlarından bir
  toplulukla temsil edilir. 1754 uygulanır ama düşük güvenli olduğu için puanlanmaz.
  Ayarda puanlanan: 1766 Mayıs, 1766 Ağustos, 1894. Sınav: 1912, 1999 İzmit, 1999 Düzce.
- Küçük depremler: saat GFZ kataloğunun başında (2006-01-01) başlar; 2006-2017 ayar,
  2018-2023 sınav.

## Geri saran sıralı ayar (Ulas'ın tarifi)
1. Ortak başlangıç değerleriyle saat ileri sarılır; her gerçek depremden sonra model durumu kaydedilir.
2. Bir sonraki gerçek deprem bir parçada gelince kontrol: "tuttu" sayılması için modelde o parça
   gerçek tarihte en az %90 dolu olmalı ve model o parçayı gerçek tarihten 20 yıldan fazla önce
   kendiliğinden kırmamış olmalı.
3. Tutmadıysa: o parçanın önceki depreminden sonraki kayda geri dönülür, o parçanın (10 km
   bloklar) ayarı yeni bir kombinasyonla değiştirilir, oradan yeniden ileri sarılır.
   Kombinasyonlar aralık içinde Latin hiperküp ile seçilir, en çok 30 deneme.
4. Tutunca bir sonraki depreme geçilir; daha önce tutmuş depremlerin bozulup bozulmadığı her
   değişiklikten sonra yeniden kontrol edilir.
5. Sınav depremleri bu döngüye hiç girmez; ayar bitince donmuş parametrelerle bir kez koşulur.
Hız için döngü önce 2 km'lik kaba ağda koşar, bulunan ayar tam ağda doğrulanır.

### Eşik güncellemesi (2026-11-26 15:03 UTC, hiçbir ikinci tur sonucu ve birinci tur kesin sonucu görülmeden eklendi)
"Tuttu" eşiği veri kalitesine göre:
- tarihsel (1719-1894, yer ve büyüklük yazılı kaynaklardan): parça en az %90 dolu
- erken aletsel (1912, 1943, 1963): en az %95
- modern (1999 ve sonrası; 2006 sonrası küçük depremler): en az %98
Erken kendiliğinden kırılmama şartı (20 yıl; küçük depremlerde bir önceki tekrar aralığının yarısı) korunur.
Sonuçlar ayrıca tek tip %90 eşiğiyle de raporlanır (eşik duyarlılığı).

### Zaman şartı (2026-11-26 15:04 UTC, hiçbir sonuç görülmeden eklendi)
Doluluk eşiğine ek olarak, modelin o parçayı kendiliğinden kıracağı tarih gerçek tarihe yakın
olmalı. Kontrol: gerçek depremden önceki kayıttan model, o parça zorlanmadan ileri sarılır;
parçanın kendiliğinden kırıldığı ilk tarih ölçülür (erken ya da geç, ikisi de sayılır;
gerçek tarihten sonra en çok toleransın sonuna kadar beklenir).
- tarihsel (1719-1894): ±20 yıl (tekrar aralığının yaklaşık %10'u)
- erken aletsel (1912, 1943, 1963): ±10 yıl
- modern (1999 ve sonrası): ±5 yıl
- küçük depremler: aynı yerdeki iki deprem arasındaki gerçek sürenin ±%25'i
"Tuttu" = doluluk eşiği VE zaman şartı birlikte sağlanırsa.

### Yer şartı ve konum hatası (2026-11-26 15:05 UTC, hiçbir sonuç görülmeden eklendi)
"Tuttu" = doğru yer VE doluluk eşiği VE zaman toleransı. Yer şartı: modelin kendiliğinden
kırdığı alan gerçek kırılma alanıyla örtüşmeli (büyük depremlerde kırılma boyunun en az
%50'si ortak; uç belirsizliği tablodaki uc_km kadar tolere edilir).
- Büyük modern depremler (1999 İzmit, 1999 Düzce, 2022, 2025): kırılma izi yüzey kırığı,
  InSAR ve GPS çalışmalarından alınır (uç belirsizliği 3-5 km).
- Küçük depremler: yalnız GFZ'nin yeniden konumlanmış kataloglarından ve büyük eksen hatası
  2 km'nin altında olan olaylar kullanılır; ham AFAD/EMSC/USGS kataloğu kullanılmaz. Her olay
  için konum toleransı o olayın kendi hata elipsi kadardır.

### Veri kullanımı genişletildi (2026-11-26 15:06 UTC, hiçbir sonuç görülmeden eklendi)
Ulas olabildiğince çok veri kullanılmasını istedi.
- Tek tek özümsemede GFZ'nin yeniden konumlanmış bütün depremleri kullanılır (Becker 2006-2023
  ~16.000, Bentz Silivri 2019, Martínez-Garzón Armutlu 2019), hata kesimi yok.
- Konum belirsizliği yumuşak özümsemeyle işlenir: her deprem, hata elipsi içindeki hücreleri
  konum olasılığıyla ağırlıklı olarak sınıra yaklaştırır (hassas olay dar alanı tam, hatalı
  olay geniş alanı az etkiler).
- Puanlamada her olayın ağırlığı 1/(hata²) ile orantılıdır.
- Ham AFAD/EMSC/USGS kataloğu yalnız toplu istatistikte (10 km blok başına sıklık, derinlik,
  b-değeri) kullanılır.
- Sonuç ayrıca yalnız hatası 2 km'nin altındaki olaylarla da raporlanır (yan yana).

### Konum adayları üzerinde iterasyon (2026-11-26 15:07 UTC, hiçbir sonuç görülmeden eklendi)
Ulas'ın önerisi: hatası birkaç km olan depremde hata alanındaki aday kutuların her biri için
ayrı "burada oldu" varsayımıyla ileri sar, geçmiş ve sonraki ayar depremleriyle en tutarlı
olanı seç.
- Yalnız ayar döneminde uygulanır, sıralı ve açgözlü: her olayda adaylar (hata elipsi içinde
  1 km aralıklı, en çok 5-10 aday) denenir, o ana kadarki ayar olaylarıyla en tutarlı olan
  seçilir. Hesap yükü olay başına aday sayısı kadar artar; büyük ağda yalnız M≥3 olaylara,
  100 m'lik yerel ağlarda bütün olaylara uygulanır.
- Sınav döneminde konum seçimi yapılmaz: adayların olasılık ağırlıklı ortalaması kullanılır
  (modelin uygun kutuyu kendisinin seçmesi tutmayı neredeyse garanti ederdi).
- Bu yöntem ile yumuşak (ağırlıklı) özümseme yan yana koşulur, ikisinin sınav sonucu raporlanır.
- Aday sayısı kuralı (15:07 UTC): aday kutu sayısı hata payının km cinsinden yuvarlanmış
  değeri kadardır (1 km hata → 1 kutu, 4 km → 4 kutu, 5 km → 5 kutu), 1 km aralıkla hata
  elipsinin uzun ekseni boyunca dizilir; 10 km'yi aşan hatalarda en çok 10 aday.
- Olasılıklı adaylar (15:08 UTC): ayar döneminde tek aday seçilmez. Her aday kutu için
  "deprem burada oldu" olasılığı = katalog hata elipsinden gelen önsel × modelin o andaki
  doluluğundan gelen olabilirlik (doluluk sınıra ne kadar yakınsa o kadar olası), normalize.
  Özümseme güncellemeleri bu olasılıklarla ağırlıklandırılır ve olasılıklar saklanır; raporda
  belirsiz depremler için aday dağılımları gösterilir. Sınav döneminde yalnız katalog önseli
  kullanılır (model doluluğu olasılığa girmez, döngüsellik olmasın).
- Aday sayısı kuralı güncellendi (15:09 UTC; bir önceki "hata kadar aday" kuralının yerine):
  katalog hatası 1σ (yaklaşık %68 güven) kabul edilir. Adaylar uzun eksen boyunca ±2σ
  aralığında 1 km arayla dizilir (5 km hata → ±10 km, 21 aday; %95 kapsama). Öncül olasılık
  Gauss: ortadakiler yüksek, kenardakiler düşük. Olasılığı %1'in altında kalan adaylar hesap
  yükünü azaltmak için atılır ve kalanlar yeniden normalize edilir.
