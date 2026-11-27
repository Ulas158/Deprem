# Fiziksel model, kümelenmeye (ETAS) bir şey katıyor mu?

Tarih: Kasım 2026. Program: `scripts/etas_fizik.py`. Ham sonuçlar: [tables/etas_fizik.json](tables/etas_fizik.json).

## Kısa sonuç

**Hayır, fiziğin kümelenmenin ötesinde bir şey kattığına dair kanıt yok.**

Doluluk haritası eklenince sınav puanı biraz artıyor. Ama artış iki kaynaktan geliyor, ikisi de fizik değil:

- **Ağın ucunda kalan sabit bir desen.** Bu desen büyük ölçüde modelin kenar etkisinden geliyor.
- **"Yakında deprem olmuş yerde yine olur" bilgisi.** Bu bilgi haritaya, gerçek depremler modele işlendiği için sızıyor. İşareti de fiziğin beklediğinin tersi: depremler, modelin gerilmesi azalmış dediği yerlerde oluyor.

## Ne yapıldı

- **Bölgeler:** Silivri ve Çınarcık yerel ağları, fay boyunca 14 km.
- **Veri:** M≥2,0 depremler (bütün kaynaklar için tam).
- **Ayar ve sınav dönemleri:**

| Bölge | Ayar 2006-2017 | Sınav 2018-2023 |
|---|---|---|
| Silivri | 89 deprem | 320 deprem (2019 dizisi dahil) |
| Çınarcık | 76 deprem | 34 deprem |

- **Temel model (ETAS):** Her depremin ardından, yakınında, zamanla azalan yeni deprem olasılığı. Buna ek olarak, ayar dönemindeki depremlerden çıkarılmış sabit bir arka plan.
- **Fizik eklenmiş model:** Arka planın yeri, fiziksel modelin o anki doluluk haritasına göre değişebilir. Toplam arka plan hızı aynı kalır; fizik yalnız "nerede" sorusuna katkı yapabilir.
- **Doluluk haritası:** Küçük deprem aşamasının en iyi 5 ayarı × 3 tohum, 15 koşunun ortalaması. Her 0,05 yılda bir alındı. Her görüntü yalnız kendinden önceki depremleri bilir.
- **Ayar ve puan:** Bütün modeller yalnız ayar döneminde ayarlandı, sınavda donduruldu.
  - Puan: sınav depremlerinde olabilirlik farkı (büyük = daha iyi).
  - Belirsizlik: 30 günlük bloklarla önyükleme, böylece kümelenmiş depremler bağımsız sayılmadı.

## Sonuçlar (sınav, ETAS'a göre toplam log-olabilirlik farkı)

| Model | Silivri | Çınarcık |
|---|---|---|
| ETAS + doluluk haritası | +8,3 [-0,9, +22,6], p = 0,06 | +7,3 [-1,2, +18,2], p = 0,06 |
| Denetim: harita fay boyunca ters çevrilmiş | +15,6 [-3,2, +46,4], p = 0,10 | +6,6 [-2,3, +17,5], p = 0,09 |
| Yalnız sabit desen (ayar dönemi ortalaması) | +16,0 [+2,1, +36,6], p = 0,007 | +5,5 [-3,1, +15,8], p = 0,13 |
| Yalnız depremlerle değişen kısım | +18,1 [+3,0, +38,2], p = 0,004, **işaret ters** | -1,2 [-2,5, +0,1] |

Köşeli parantez %95 aralık; p, "hiç iyileşme yok" olasılığı.

## Nasıl okunmalı

1. **Ters çevrilmiş harita da aynı kadar iyileştiriyor.** Haritanın kendisi neredeyse simetrik (ağın iki ucunda yüksek, ortada düşük), bu yüzden çevirmek pek bir şey değiştirmiyor. Yani iyileşme, depremin tam nerede olacağına dair özel bir bilgiden gelmiyor.
2. **Sabit desen Silivri'de anlamlı katkı yapıyor.** Ama bu desen büyük ölçüde ağın sınırında yüklemenin yığılmasından geliyor. Bu bir hesap ağı etkisi, fay fiziği değil.
3. **Depremlerle değişen kısım Silivri'de anlamlı, ama işareti ters.**
   - Fizik derdi ki: gerilmesi artan yerde deprem olur.
   - Veri diyor ki: modelde gerilmesi düşen yerde deprem oluyor.
   - Model, gerçek bir depremi işlediğinde o noktanın gerilmesini boşaltıyor. Bu yüzden "gerilmesi düşmüş yer" aslında "yakında deprem olmuş yer" demek. Bu, ETAS'ın zaten yakaladığı kümelenmenin bir kopyası. Çınarcık'ta bu katkı hiç yok.

## Sınırlar

- Yalnız iki küçük bölge ve 6 yıllık sınav. Silivri sınavına 2019 dizisi hâkim.
- ETAS fay boyunca tek boyutlu kuruldu. Derinlik ve faydan uzaklık kullanılmadı.
- Doluluk haritası, gerçek depremlerin modele işlendiği sürümden alındı. İşlemesiz bir sürümle sınama yapılmadı.
