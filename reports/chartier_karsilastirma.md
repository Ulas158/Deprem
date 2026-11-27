# Bizim model ile Chartier ve ekibinin (2021) Marmara modeli: ne aynı, ne farklı

Kaynak: Chartier vd. 2021, NHESS 21, 2733-2751 (doi:10.5194/nhess-21-2733-2021). Değerleri makalenin metninden ve Tablo 1-3'ten aldım. Tablo 1'in PDF'ten çıkarılan metninde sütunlar kaymıştı, bu yüzden bazı satırların eşleşmesi en iyi okumamdır. Makalenin ek dosyasını (tam fay tablosu) açmadım.

## Kısa sonuç

- **Aynı yaptığımız:** fizik tabanlı, hücrelere bölünmüş fay; levha hareketiyle yükleme; deprem gerilmesinin komşu hücrelere ve faylara aktarılması. Ganos ve İzmit kayma hızları da hemen hemen aynı.
- **En önemli farklar:**
  1. Amaç ve sınav. Onlar tarih değil, uzun vadeli oranları sınamış.
  2. Ölçüt. Paleosismoloji (kazı) verisini ölçüt olarak kullanmışlar, biz hiç kullanmadık.
  3. Krip varsayımı. Onlar birkaç varsayımı yan yana koyup veriye göre puanlamış; bizim varsayımımız, onların en iyi bulduğunun neredeyse tersi.
  4. Orta Marmara kayma hızı ve sismojenik derinlik.

## Girdiler

| | Biz | Chartier vd. 2021 | Aynı mı? |
|---|---|---|---|
| Fay haritası | GEM / EMME bölgesel veritabanı. Kuzey kol 4 parça, yan kollar 5 parça | MTA diri fay haritası (Emre vd. 2018). Kuzey kol 11 kesim (Tekirdağ, Orta çukur, Kumburgaz, Avcılar, Adalar, Çınarcık…) | **Farklı:** onlarınki daha ayrıntılı |
| Fay eğimi ve kayma yönü | Hepsi düşey, saf sağ yanal | Adalar (-168°) ve Çınarcık (-147°) eğik, normal bileşenli; diğerleri sağ yanal | **Farklı** |
| Hücre boyu | 2 km (yerel ağlar 200 m) | 1 km üçgen | Benzer |
| Ganos kayma hızı | 18,5 mm/yıl | 17-19 (ort. 18) | **Aynı** |
| Orta Marmara (Tekirdağ-Kumburgaz-Avcılar) | 12,5-17 mm/yıl, ayarla seçildi | Varsayıma göre: kilitli 17-19; yarı krip 7,5-8,5; tam krip 0 | **Farklı:** bizimki daha düşük ve ayara bırakıldı |
| Adalar | Orta Marmara ile aynı (12,5-17) | 12,2-16,5 (ort. 14,3) | Benzer |
| Çınarcık | Orta Marmara ile aynı (12,5-17) | 2,2-5 (ort. 3,6) | **Çok farklı:** onlar hareketi paralel kollar arasında paylaştırmış |
| İzmit | 22 | 18-22 (ort. 20) | Yakın |
| Sapanca-Akyazı-Karadere | 22 | Adapazarı 18, Karadere 15 | **Farklı:** bizimki yüksek |
| Düzce | 18 | 11,8-18,2 (ort. 15) | Biraz yüksek |
| Kayma hızı kaynağı | Kendi GPS profil uyumumuz, Reilinger 2006, Ergintav 2014; yan kollar GEM | Hergert ve Heidbach 2010 (GPS'ten mekanik model); Reilinger 2006'yı "jeolojik hızları abartıyor" diye reddetmişler | **Farklı** |
| Armutlu kolu | 0,8 mm/yıl | 1'den az | **Aynı** |
| Sismojenik derinlik | Depremlerin derinliğinden: Orta Marmara 10-12, Adalar 10-18, batı 18-24 km | Her yerde 15 km, Adalar 10 km | **Farklı:** bizimki batıda daha derin, ortada daha sığ |
| Krip | Orta Marmara'da (27,75-28,05°D) sığ krip (3-8 km) + derinde sıcaklığa bağlı yavaş kayma | Dört varsayım yan yana: tam krip, yarı krip, 5 km'ye kadar kilitli (altı kayıyor), tamamen kilitli. Sonuç: tam krip ve yarı krip kötü; 5 km'ye kadar kilitli ve tamamen kilitli varsayımları en iyi | **Farklı:** biz tek varsayım seçtik ve onların zayıf bulduğu yöne yakın |
| Katalog | Tarihsel büyük depremler (1719-1894), USGS M≥5, AFAD/EMSC, GFZ küçük depremler | Kadirioğlu vd. 2018 (1900-2012, Mw'ye çevrilmiş), SHEEC (1000-1899) | **Farklı** |
| Paleosismoloji | Kullanılmadı | Ganos, İzmit, Adapazarı, Karadere, Düzce kazılarından M≥7,2 tekrar oranları (her biri yaklaşık 300-550 yılda bir) | **Bizde yok** |

## Yöntem ve kurallar

| | Biz | Chartier vd. 2021 | Aynı mı? |
|---|---|---|---|
| Simülatör türü | Kendi yarı-statik hücresel simülatörümüz (RSQSim ve Virtual Quake ailesinden) | RSQSim | Aynı aile |
| Sürtünme kuralı | Dayanım eşiği + rastgele dayanım farkları + kırık ucunda geçici zayıflama; krip için sıcaklığa bağlı (a-b) | Hız ve duruma bağlı sürtünme, a = 0,01, b = 0,015 | **Farklı** |
| Normal gerilme | Derinlikle artan, 1-17 MPa/km (ayarla) | Sabit ~100 MPa (özetten, doğrulanmadı) | **Farklı** |
| Gerçek depremleri modele işleme | Evet: her gerçek depremde o yer kırıldı | Hayır: model serbest koştu | **Farklı** |
| Koşu uzunluğu | 1719-2026 (~300 yıl), ısınma 400 yıl | 10.000 yıl, ilk 1.000 yıl atılmış | **Farklı** |
| Ayar ve sınav ölçütü | Belirli depremlerin yeri, doluluğu ve zamanı ("tuttu") | Uzun vadeli büyüklük-sıklık dağılımı ve paleosismolojik tekrar oranlarıyla uyum | **Farklı:** onlar olasılık ve oran sınamış |
| Varsayımları yarıştırma | Tek model, parametre araması | Mantık ağacı: krip varsayımları, kırık kümeleri, büyüklük dağılımı biçimleri ayrı ayrı puanlanmış | **Farklı** |
| En büyük deprem | Model en çok ~M7,1 üretti | M7,7'ye kadar (Küme 1), M8,0'a kadar (Küme 2) izin verilmiş | **Farklı** |
| Küçük depremler | Sınandı (Silivri, Çınarcık) | Sınanmadı. Fizik modeli Gutenberg-Richter dengesini kendiliğinden üretemedi | Biz daha ileri gittik ama sonuç benzer |
| İstanbul olasılığı | Model %8 (koşudan koşuya 0-42) | Bu makalede yok, sonraki çalışmaya bırakılmış | — |

## Bundan ne çıkarıyorum

1. **Kayma hızlarını ve derinliği onlarınkine çekmek ucuz ve doğrudan sınanabilir.** Değiştirilecek dört şey var:
   - Orta Marmara'yı kilitli kabul edip 17-19 mm/yıl vermek;
   - Çınarcık'ı 3,6 mm/yıl yapıp hareketi paralel kollarla paylaştırmak;
   - Sapanca-Karadere ve Düzce'yi düşürmek;
   - derinliği 15 km (Adalar 10 km) yapmak.
2. **Krip varsayımını tek seçmek yerine yarıştırmak.** Onların yaptığı gibi tam kilitli, 5 km'ye kadar kilitli, yarı krip ve bizim sığ krip varsayımı yan yana koşulup aynı ölçütle puanlanabilir.
3. **Ölçütü değiştirmek.** Tarih tutturmak yerine iki şeyi ölçmek: paleosismolojik tekrar oranları (Tablo 3) ve 1700'den beri büyüklük-sıklık dağılımı. Bu, hedefin olasılık olmasıyla da uyumlu.
4. **Fay haritasını MTA'ya geçirmek daha büyük iş.** Yük aktarım hesaplarının baştan yapılması gerekir (kaba ağda ~8 dk, tam ağda saatler). Kayma hızı ve krip denemeleri bir şey gösterirse yapılmaya değer.
5. Kadirioğlu 2018 kataloğu şu an elimizde yok. Erişilebilirse eklenir, yoksa bunu sınır olarak yazarım.
