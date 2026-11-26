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
