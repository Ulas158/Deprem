# Model değişkenlerinin etki payı (yüzde)

> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Yüzdeler, modelin geçmiş depremleri ileriye dönük (cevabı görmeden) tahmin ederken topladığı bilginin hangi parçadan geldiğini ölçer.

Yöntem: `scripts/etki_payi.py` (fay fay sıralı model, PR 14), her sürüm baştan sıralı yürütüldü. Etki payı = o parça çıkarılınca (veya eklenince) kaybedilen (kazanılan) ileriye dönük log olabilirlik ÷ ana modelin, hiçbir bilgi taşımayan sıfır noktasına (sabit hız, alanla orantılı) göre toplam bilgisi. %95 aralıklar yıl blokları bootstrap'ı. Bir parçanın tek başına payı, başka parçalarla örtüşen bilgiyi saymaz; bu yüzden toplam %100 etmez.

## Türkiye (1900-2025, M≥5,5)

Olay sayısı: 421. Ana modelin toplam bilgisi (sıfır noktasına göre): olay başına 1.75 [1.14, 2.48] log birimi. Ana model − ETAS: -0.011 olay başına; ETAS tek başına toplam bilginin %101'ini taşıyor.

| değişken                                          | etki payı   | %95 aralık      |   olay başına (log) |     p | karar                           |
|:--------------------------------------------------|:------------|:----------------|--------------------:|------:|:--------------------------------|
| Arka plan haritası (nerede deprem olur)           | %0.0        | [%-1.7, %1.7]   |               0.001 | 0.484 | sıfırdan ayırt edilemez         |
| · geçmiş sismisite haritası                       | %0.2        | [%-1.4, %1.9]   |               0.004 | 0.405 | sıfırdan ayırt edilemez         |
| Kayma hızı (harita ve doluluktan)                 | %0.4        | [%-1.1, %1.9]   |               0.007 | 0.317 | sıfırdan ayırt edilemez         |
| · fay kayma hızı haritası                         | %0.6        | [%-0.8, %2.0]   |               0.01  | 0.214 | sıfırdan ayırt edilemez         |
| Artçı/tetikleme (ETAS benzeri, hızlı)             | %20.8       | [%11.6, %32.9]  |               0.366 | 0     | gerçek katkı                    |
| Yavaş etki (yıllar)                               | %0.8        | [%-0.1, %1.9]   |               0.014 | 0.051 | sıfırdan ayırt edilemez         |
| Tetiklemenin tümü (hızlı + yavaş)                 | %84.9       | [%53.5, %123.0] |               1.49  | 0     | gerçek katkı                    |
| Uzaklıkla azalan faylar arası etki (100 km ötesi) | %4.7        | [%1.7, %8.1]    |               0.083 | 0     | gerçek katkı                    |
| Büyük sistemler arası ayrı çarpan                 | %-0.3       | [%-0.7, %0.2]   |              -0.005 | 0.852 | sıfırdan ayırt edilemez         |
| Coulomb gerilme aktarımı                          | %-0.1       | [%-0.2, %0.0]   |              -0.001 | 0.907 | sıfırdan ayırt edilemez         |
| Son kırılmadan geçen süre                         | %-0.1       | [%-0.6, %0.5]   |              -0.001 | 0.58  | sıfırdan ayırt edilemez         |
| Doluluk / kapasite (sürenin üstüne)               | %-0.7       | [%-1.2, %-0.3]  |              -0.012 | 0.999 | zararlı (modeli kötüleştiriyor) |
| Doluluk terimi tümü                               | %-0.8       | [%-1.4, %-0.1]  |              -0.013 | 0.99  | zararlı (modeli kötüleştiriyor) |
| Fay fay ayar (ortak ayara karşı)                  | %-1.5       | [%-3.4, %0.2]   |              -0.027 | 0.958 | sıfırdan ayırt edilemez         |
| Sürekli yeniden eğitim (1900'de dondurmaya karşı) | %12.8       | [%7.1, %21.5]   |               0.225 | 0     | gerçek katkı                    |
| Zincirleme kırılma (yüklü komşu fay zinciri)      | %0.0        | [%-0.1, %0.1]   |               0     | 0.491 | sıfırdan ayırt edilemez         |
| b-değeri (yerel, son 5 yıl)                       | %-0.1       | [%-0.8, %0.5]   |              -0.001 | 0.543 | sıfırdan ayırt edilemez         |
| Moment açığı (küçük depremler dahil)              | %0.0        | [%-0.1, %0.1]   |               0     | 0.448 | sıfırdan ayırt edilemez         |
| dv/v (Doğu Marmara)                               | %0.0        | [%-0.0, %0.0]   |               0     | 0.145 | sıfırdan ayırt edilemez         |

## Marmara kutusu (1900-2025, M≥5,5)

Olay sayısı: 29. Ana modelin toplam bilgisi (sıfır noktasına göre): olay başına 1.29 [-0.01, 2.87] log birimi. Ana model − ETAS: +0.059 olay başına; ETAS tek başına toplam bilginin %95'ini taşıyor.

| değişken                                          | etki payı   | %95 aralık      |   olay başına (log) |     p | karar                           |
|:--------------------------------------------------|:------------|:----------------|--------------------:|------:|:--------------------------------|
| Arka plan haritası (nerede deprem olur)           | %-5.6       | [%-11.9, %0.9]  |              -0.072 | 0.958 | sıfırdan ayırt edilemez         |
| · geçmiş sismisite haritası                       | %-1.6       | [%-6.6, %2.9]   |              -0.021 | 0.759 | sıfırdan ayırt edilemez         |
| Kayma hızı (harita ve doluluktan)                 | %2.7        | [%-2.3, %7.9]   |               0.035 | 0.146 | sıfırdan ayırt edilemez         |
| · fay kayma hızı haritası                         | %2.2        | [%-2.9, %7.5]   |               0.028 | 0.2   | sıfırdan ayırt edilemez         |
| Artçı/tetikleme (ETAS benzeri, hızlı)             | %19.0       | [%3.6, %37.6]   |               0.245 | 0.002 | gerçek katkı                    |
| Yavaş etki (yıllar)                               | %-0.1       | [%-2.0, %1.7]   |              -0.001 | 0.518 | sıfırdan ayırt edilemez         |
| Tetiklemenin tümü (hızlı + yavaş)                 | %91.1       | [%22.9, %174.4] |               1.177 | 0     | gerçek katkı                    |
| Uzaklıkla azalan faylar arası etki (100 km ötesi) | %2.0        | [%-2.2, %5.6]   |               0.025 | 0.166 | sıfırdan ayırt edilemez         |
| Büyük sistemler arası ayrı çarpan                 | %-0.6       | [%-4.6, %1.9]   |              -0.008 | 0.568 | sıfırdan ayırt edilemez         |
| Coulomb gerilme aktarımı                          | %-0.1       | [%-0.8, %0.4]   |              -0.002 | 0.628 | sıfırdan ayırt edilemez         |
| Son kırılmadan geçen süre                         | %-0.7       | [%-4.9, %3.2]   |              -0.009 | 0.622 | sıfırdan ayırt edilemez         |
| Doluluk / kapasite (sürenin üstüne)               | %1.0        | [%-0.3, %2.6]   |               0.013 | 0.08  | sıfırdan ayırt edilemez         |
| Doluluk terimi tümü                               | %0.4        | [%-2.6, %3.1]   |               0.005 | 0.391 | sıfırdan ayırt edilemez         |
| Fay fay ayar (ortak ayara karşı)                  | %-0.1       | [%-6.1, %5.7]   |              -0.001 | 0.528 | sıfırdan ayırt edilemez         |
| Sürekli yeniden eğitim (1900'de dondurmaya karşı) | %20.1       | [%12.3, %27.3]  |               0.26  | 0     | gerçek katkı                    |
| Zincirleme kırılma (yüklü komşu fay zinciri)      | %-0.1       | [%-1.5, %1.3]   |              -0.001 | 0.521 | sıfırdan ayırt edilemez         |
| b-değeri (yerel, son 5 yıl)                       | %-2.1       | [%-5.4, %-0.1]  |              -0.027 | 0.983 | zararlı (modeli kötüleştiriyor) |
| Moment açığı (küçük depremler dahil)              | %-0.1       | [%-0.9, %0.7]   |              -0.001 | 0.571 | sıfırdan ayırt edilemez         |
| dv/v (Doğu Marmara)                               | %0.0        | [%-0.0, %0.0]   |               0     | 0.515 | sıfırdan ayırt edilemez         |

## Seçim dönemi (1900-1999) ve dokunulmamış dönem (2000-2025), Türkiye

Eleme kararı yalnız 1900-1999 sonuçlarıyla verildi; 2000-2025 bu kararda kullanılmadı.

### 1900-1999 (324 deprem)

| değişken                                          | etki payı   | %95 aralık      |   olay başına (log) |     p | karar                           |
|:--------------------------------------------------|:------------|:----------------|--------------------:|------:|:--------------------------------|
| Arka plan haritası (nerede deprem olur)           | %-1.7       | [%-3.9, %0.5]   |              -0.024 | 0.931 | sıfırdan ayırt edilemez         |
| · geçmiş sismisite haritası                       | %-1.1       | [%-3.4, %1.3]   |              -0.017 | 0.835 | sıfırdan ayırt edilemez         |
| Kayma hızı (harita ve doluluktan)                 | %0.0        | [%-2.1, %2.2]   |               0     | 0.491 | sıfırdan ayırt edilemez         |
| · fay kayma hızı haritası                         | %0.7        | [%-1.4, %2.7]   |               0.01  | 0.268 | sıfırdan ayırt edilemez         |
| Artçı/tetikleme (ETAS benzeri, hızlı)             | %21.1       | [%13.2, %29.7]  |               0.307 | 0     | gerçek katkı                    |
| Yavaş etki (yıllar)                               | %1.4        | [%-0.0, %2.8]   |               0.021 | 0.026 | sıfırdan ayırt edilemez         |
| Tetiklemenin tümü (hızlı + yavaş)                 | %87.6       | [%59.4, %118.6] |               1.278 | 0     | gerçek katkı                    |
| Uzaklıkla azalan faylar arası etki (100 km ötesi) | %6.1        | [%2.3, %10.9]   |               0.089 | 0     | gerçek katkı                    |
| Büyük sistemler arası ayrı çarpan                 | %-0.3       | [%-1.0, %0.3]   |              -0.005 | 0.824 | sıfırdan ayırt edilemez         |
| Coulomb gerilme aktarımı                          | %-0.1       | [%-0.2, %0.1]   |              -0.001 | 0.857 | sıfırdan ayırt edilemez         |
| Son kırılmadan geçen süre                         | %-0.2       | [%-0.9, %0.6]   |              -0.003 | 0.709 | sıfırdan ayırt edilemez         |
| Doluluk / kapasite (sürenin üstüne)               | %-1.3       | [%-1.9, %-0.6]  |              -0.019 | 1     | zararlı (modeli kötüleştiriyor) |
| Doluluk terimi tümü                               | %-1.5       | [%-2.3, %-0.6]  |              -0.022 | 1     | zararlı (modeli kötüleştiriyor) |
| Fay fay ayar (ortak ayara karşı)                  | %-3.0       | [%-5.4, %-0.6]  |              -0.044 | 0.994 | zararlı (modeli kötüleştiriyor) |
| Sürekli yeniden eğitim (1900'de dondurmaya karşı) | %11.9       | [%8.5, %15.4]   |               0.174 | 0     | gerçek katkı                    |
| Zincirleme kırılma (yüklü komşu fay zinciri)      | %-0.0       | [%-0.2, %0.1]   |              -0.001 | 0.666 | sıfırdan ayırt edilemez         |
| b-değeri (yerel, son 5 yıl)                       | %0.0        | [%-0.0, %0.0]   |               0     | 0.283 | sıfırdan ayırt edilemez         |
| Moment açığı (küçük depremler dahil)              | %0.0        | [%-0.1, %0.2]   |               0     | 0.328 | sıfırdan ayırt edilemez         |
| dv/v (Doğu Marmara)                               | %0.0        | [%-0.0, %0.0]   |               0     | 0.283 | sıfırdan ayırt edilemez         |

### 2000-2025 (97 deprem)

| değişken                                          | etki payı   | %95 aralık      |   olay başına (log) |     p | karar                   |
|:--------------------------------------------------|:------------|:----------------|--------------------:|------:|:------------------------|
| Arka plan haritası (nerede deprem olur)           | %3.1        | [%1.2, %5.3]    |               0.085 | 0.002 | gerçek katkı            |
| · geçmiş sismisite haritası                       | %2.6        | [%0.9, %4.3]    |               0.07  | 0     | gerçek katkı            |
| Kayma hızı (harita ve doluluktan)                 | %1.0        | [%-0.6, %2.7]   |               0.027 | 0.125 | sıfırdan ayırt edilemez |
| · fay kayma hızı haritası                         | %0.4        | [%-1.2, %2.2]   |               0.012 | 0.305 | sıfırdan ayırt edilemez |
| Artçı/tetikleme (ETAS benzeri, hızlı)             | %20.4       | [%1.5, %49.0]   |               0.56  | 0.011 | gerçek katkı            |
| Yavaş etki (yıllar)                               | %-0.2       | [%-1.1, %0.6]   |              -0.007 | 0.71  | sıfırdan ayırt edilemez |
| Tetiklemenin tümü (hızlı + yavaş)                 | %80.1       | [%17.7, %168.4] |               2.196 | 0     | gerçek katkı            |
| Uzaklıkla azalan faylar arası etki (100 km ötesi) | %2.3        | [%-1.1, %7.8]   |               0.063 | 0.229 | sıfırdan ayırt edilemez |
| Büyük sistemler arası ayrı çarpan                 | %-0.2       | [%-0.8, %0.3]   |              -0.004 | 0.696 | sıfırdan ayırt edilemez |
| Coulomb gerilme aktarımı                          | %-0.1       | [%-0.2, %0.1]   |              -0.002 | 0.778 | sıfırdan ayırt edilemez |
| Son kırılmadan geçen süre                         | %0.2        | [%-0.7, %1.2]   |               0.006 | 0.342 | sıfırdan ayırt edilemez |
| Doluluk / kapasite (sürenin üstüne)               | %0.3        | [%-0.2, %0.7]   |               0.008 | 0.108 | sıfırdan ayırt edilemez |
| Doluluk terimi tümü                               | %0.5        | [%-0.3, %1.4]   |               0.014 | 0.125 | sıfırdan ayırt edilemez |
| Fay fay ayar (ortak ayara karşı)                  | %1.1        | [%-1.0, %3.4]   |               0.029 | 0.175 | sıfırdan ayırt edilemez |
| Sürekli yeniden eğitim (1900'de dondurmaya karşı) | %14.4       | [%0.3, %35.9]   |               0.393 | 0.019 | gerçek katkı            |
| Zincirleme kırılma (yüklü komşu fay zinciri)      | %0.1        | [%-0.2, %0.3]   |               0.002 | 0.27  | sıfırdan ayırt edilemez |
| b-değeri (yerel, son 5 yıl)                       | %-0.2       | [%-2.3, %1.5]   |              -0.005 | 0.554 | sıfırdan ayırt edilemez |
| Moment açığı (küçük depremler dahil)              | %-0.0       | [%-0.3, %0.2]   |              -0.001 | 0.559 | sıfırdan ayırt edilemez |
| dv/v (Doğu Marmara)                               | %0.0        | [%-0.0, %0.0]   |               0     | 0.054 | sıfırdan ayırt edilemez |

## Elenmiş model tam modele karşı

Elenen parçalar (1900-1999'da katkısı sıfırdan ayırt edilemeyen ya da zararlı, ve başka bir parçayla yer değiştirmeyen): Coulomb gerilme yönü, büyük sistemler arası ayrı çarpan ve aynı sistem çarpanı, doluluk/son kırılmadan geçen süre, fay fay ayar (ortak ayar kullanıldı). Pozitif = elenmiş model daha iyi. Seçim 1900-1999 verisiyle yapıldığı için o dönemdeki üstünlük iyimserdir; asıl ölçü 2000-2025'tir.

**Türkiye 1900-1999**

| model                                                                         |   olay |   tam modele göre (olay başına) | %95 GA           |   p (elenmiş daha iyi değil) |
|:------------------------------------------------------------------------------|-------:|--------------------------------:|:-----------------|-----------------------------:|
| Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)              |    324 |                           0.055 | [+0.014, +0.095] |                        0.003 |
| ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model) |    324 |                           0.017 | [-0.026, +0.062] |                        0.229 |

**Türkiye 2000-2025**

| model                                                                         |   olay |   tam modele göre (olay başına) | %95 GA           |   p (elenmiş daha iyi değil) |
|:------------------------------------------------------------------------------|-------:|--------------------------------:|:-----------------|-----------------------------:|
| Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)              |     97 |                          -0.022 | [-0.093, +0.043] |                        0.742 |
| ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model) |     97 |                          -0.006 | [-0.065, +0.051] |                        0.574 |

**Türkiye 1900-2025**

| model                                                                         |   olay |   tam modele göre (olay başına) | %95 GA           |   p (elenmiş daha iyi değil) |
|:------------------------------------------------------------------------------|-------:|--------------------------------:|:-----------------|-----------------------------:|
| Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)              |    421 |                           0.037 | [+0.003, +0.074] |                        0.018 |
| ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model) |    421 |                           0.011 | [-0.023, +0.048] |                        0.261 |

**Marmara 1900-1999**

| model                                                                         |   olay |   tam modele göre (olay başına) | %95 GA           |   p (elenmiş daha iyi değil) |
|:------------------------------------------------------------------------------|-------:|--------------------------------:|:-----------------|-----------------------------:|
| Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)              |     26 |                           0.003 | [-0.085, +0.094] |                        0.475 |
| ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model) |     26 |                          -0.028 | [-0.122, +0.088] |                        0.721 |

**Marmara 2000-2025**

| model                                                                         |   olay |   tam modele göre (olay başına) | %95 GA           |   p (elenmiş daha iyi değil) |
|:------------------------------------------------------------------------------|-------:|--------------------------------:|:-----------------|-----------------------------:|
| Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)              |      3 |                          -0.156 | [-0.494, +0.126] |                        0.827 |
| ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model) |      3 |                          -0.332 | [-0.447, -0.234] |                        1     |

**Marmara 1900-2025**

| model                                                                         |   olay |   tam modele göre (olay başına) | %95 GA           |   p (elenmiş daha iyi değil) |
|:------------------------------------------------------------------------------|-------:|--------------------------------:|:-----------------|-----------------------------:|
| Elenmiş model (ortak ayar; Coulomb, sistem çarpanı, doluluk yok)              |     29 |                          -0.013 | [-0.097, +0.076] |                        0.635 |
| ETAS (aynı düzende; elenmiş model + yavaş etki de yok ile birebir aynı model) |     29 |                          -0.059 | [-0.142, +0.049] |                        0.886 |

## Bu çalışmada ölçülmeyenler

- **GPS kilitlenme / krip:** sıralı modelde yok. Önceki çalışmada (reports/parametre_aramasi.md) açık GPS istasyonları seyrek olduğu için krip çözülemedi; etki payı ölçülemez.
- **Gaz / sıvı:** varsayılan olarak kapalı. Ortak modelde eklenince ileriye dönük puan ETAS'a göre kötüleşti (birleşik eğitim özeti); payı sıfır veya eksi.

## Sınırlar

- Marmara'da 1900-2025'te az sayıda M≥5,5 deprem var; aralıklar bu yüzden çok geniş.
- Çıkar-bir-tane payı, parçanın 'yerine başka bir parça geçemeyen' katkısıdır. Örneğin fay haritası ile geçmiş sismisite haritası aynı bilgiyi kısmen taşıdığından her birinin tek payı küçük, ikisinin birlikte payı büyük görünür.
- b-değeri yalnız 1998 sonrası, dv/v yalnız 2009 sonrası ve yalnız Doğu Marmara'da dolu; bu dönemlerde çok az deprem var, testin gücü düşük.
