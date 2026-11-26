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

## Bilimsel düzeltmeler (2026-11-26 15:10 UTC, hiçbir ikinci tur sonucu görülmeden eklendi)
Ulas "bilimsel ol, manipüle olma" dedi. Plandaki zayıf noktalar ve önlemler:
1. **"Kırılan parça %100 doluydu" varsayımı fiziksel olarak tam doğru değil.** Dinamik kırılma,
   sınırın altındaki alanları da kırar (birinci tur modelinde kırılan hücrelerin çoğu %20-60
   doluydu). Bütün parçayı %100'e çekmek komşulara aktarılan yükü abartır. Önlem: %100'e
   yalnız başlangıç noktası (bilinen ya da olasılıklı hiposantr çevresi, ~5 km) çekilir;
   parçanın geri kalanı modeldeki gerçek gerilmesiyle kırılır. Bütün parçayı %100'e çeken
   sürüm de yan yana koşulur ve raporlanır.
2. **Serbestlik derecesi sayılır.** 11 ortak parametre + 10 km blok farkları + geri sarma +
   konum seçimi, ayarda puanlanan 3 büyük depremden çok daha fazla serbestlik verir; ayarın
   tutması bu yüzden kanıt sayılmaz. Blok farkı yalnız ayar verisi içinde birini-dışarıda-bırak
   çapraz doğrulamasıyla iyileşme gösterirse eklenir. Raporda toplam serbest sayı yazılır.
3. **Sınavın gücü düşük.** 3 büyük sınav depremi tek başına hiçbir şeyi kanıtlayamaz (±20 yıllık
   pencereye rastgele bir modelin düşme olasılığı bile ~%15-25). Asıl sınav 2018-2023 küçük
   depremleridir (yüzlerce olay). Her sonuç iki kıyasla verilir: (a) şans (rastgele kaydırma),
   (b) basit kestirme "son depremden beri geçen süre × kayma hızı". Model bu ikisini geçemezse
   sonuç "model işe yaramıyor" diye yazılır.
4. Ulas'ın veya başka birinin önerileri, bu bölümle çelişirse uygulanmadan önce gerekçesiyle
   itiraz edilir.

### Zaman sırası ve komşu depremler (2026-11-26 15:14 UTC, hiçbir sonuç görülmeden eklendi)
Bütün depremler (büyük, küçük, ana kol, yan kol) tek bir zaman çizelgesinde eskiden yeniye
işlenir; hiçbiri atlanmaz. Her biri sırası gelince yükünü elastik matrisle bütün hücrelere
aktarır. Modellenmiş herhangi bir fay koluna (ikinci sürüm ağında güney/orta kol, İmralı,
Armutlu dahil) 10 km'den yakın olaylar özümsenir. Daha uzaktaki olaylar modelde fay olmadığı
için yerleştirilemez; sayıları ve toplam momentleri raporda ayrıca listelenir.

### Model dışı büyük depremler (2026-11-26 15:15 UTC, hiçbir sonuç görülmeden eklendi)
Coulomb gerilmesi kırılma boyuna göre azalır (uzakta kabaca (L/r)^3). Bu yüzden 10 km kesimi
küçük depremler için yeterli ama büyük depremler için değil: M5 (L~3 km) 10 km'de birkaç
yüzde, M6 (L~10 km) 10-20 km'de ~%5-15, M7 uzaklarda bile 0,1 bar'ı (tetikleme eşiği olarak
anılan 0,01 MPa) aşabilir. Kural: modellenmiş faylara 10 km'den uzak, M≥5 olaylar dış kaynak
olarak eklenir. Bilinen kırılma düzlemi (odak mekanizması, yayın) ya da büyüklük-boy
bağıntısından yaklaşık düzlem kullanılır; modele binen gerilme cutde ile hesaplanır ve olay
zamanında uygulanır. Hücrelerde 0,01 MPa'nın altında etki üretenler atlanır. Adaylar en az:
1953 Yenice-Gönen, 1957 Abant, 1964 Manyas, 1967 Mudurnu (birinci tur tablosunda yok;
ana kola çok yakın, kontrol edilecek), 1855 Bursa.
- Etkiye dayalı kural (15:16 UTC, yukarıdaki M≥5 / 10 km kuralının yerine): büyüklük ve
  uzaklıktan bağımsız olarak, herhangi bir model hücresinde Coulomb gerilmesini 0,01 MPa'dan
  (0,1 bar) fazla değiştiren her olay eklenir. Ulas %1 doluluk eşiğini önerdi; %1 doluluk
  dayanım fazlasına bağlı olduğu için (başlangıçta derinde 0,1 MPa, sığda 0,025 MPa) hem
  0,1 bar'dan 2,5-10 kat daha kaba hem de ayarlanan bir parametreye bağlı (ayar değiştikçe
  hangi depremlerin girdiği değişirdi). Bu yüzden eşik sabit fiziksel birimde (0,01 MPa)
  tutulur; bu, %1 kuralının seçeceği bütün olayları da kapsar. Eklenen olay sayısı raporlanır.
- Düzeltme (15:18 UTC): yukarıdaki "M5 10 km'de birkaç yüzde" tahmini yanlıştı. Elastik
  matristen ölçülen fay boyu Coulomb değişimi (kırılma ucundan, 10 km derinlik): M5 5 km'de
  0,027 MPa, 10 km'de 0,004 MPa; M6 10 km'de 0,063, 20 km'de 0,014 MPa; M7 (50 km kırılma)
  20 km'de 0,17, 50 km'de 0,008 MPa. Etkiye dayalı kural (0,01 MPa) değişmez.

## Birinci turdan sızıntı önlemi (2026-11-26 15:19 UTC)
Birinci tur 1912 ve 1999 depremlerini de puanladı; bunlar ikinci turun sınav depremleri.
Bu yüzden birinci tur sonuçları ikinci turda aralık daraltmak, senaryo seçmek ya da başlangıç
değeri belirlemek için kullanılmaz; yalnız "ayarsız karşılaştırma" olarak raporlanır.
Açıklama: bu plan yazılmadan önce birinci turun ilk 1-2 koşusunun ara sonuçları görüldü
(13:53 ve 14:56 UTC; 8 test depreminde AUC ortalaması ~0,50 ve ~0,57). Başlangıç değerleri
bu yüzden bilerek birinci tur ayarından değil literatürden alındı; aralıklar fiziksel olarak
makul geniş tutuldu. Okuyucu bu bilgiyi göz önünde tutmalı.
- Uzlaşma (15:20 UTC, Ulas'ın isteği): birinci tur parametreleri ikinci tura alternatif
  başlangıç noktası olarak girebilir, şu şartla: karar yalnız ikinci turun ayar depremlerindeki
  (1766 Mayıs, 1766 Ağustos, 1894) puana göre verilir. Literatür başlangıcı ile birinci tur
  başlangıcı bu üç depremde karşılaştırılır, iyi olan başlangıç alınır. 1912 ve 1999 sonuçları
  hiçbir aşamada kullanılmaz. Aralıklar değişmez. Not: birinci tur ayarı genel tarih
  istatistiklerine (tekrarlanma, büyüklük) göre yapıldığı için sınav depremlerinden zayıf bir
  bilgi taşır; bu raporda belirtilir.

### Doğu boşluğu (2026-11-26 15:22 UTC, hiçbir sonuç görülmeden eklendi)
GFZ kataloğu 30,65°D'de bittiği için İzmit-Sapanca doğusu ve Düzce'de hassas veri yok. Bu
bölgede 2006 sonrası ham katalogda ~1.300 kayıt var (AFAD 495, EMSC 798; çakışanlar ayıklanır),
ama kayıtlarda hata payı bilgisi yok. Karar: ham AFAD kayıtları (EMSC ile çakışanlar tek sayılır)
bu bölgede σ = 5 km varsayılarak ±2σ aday kutularıyla özümsemeye girer, puanlamada ağırlıkları
GFZ olaylarından düşük kalır (1/σ²). 1999 kırılma alanı için yayınlanmış yeniden konumlanmış
bir artçı kataloğu bulunursa o tercih edilir. Sonuç bu bölge dahil ve hariç ayrı raporlanır.
