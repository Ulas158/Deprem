# Ege'deki üstünlük: ön kayıtlı sınama sonucu (2026 verisi)

> Deneysel araştırma çıktısıdır, deprem tahmini değildir.

Ön kayıt: [`docs/onkayit/2026-11-25_ege_ustunluk.md`](../docs/onkayit/2026-11-25_ege_ustunluk.md)
(commit `9a4b4ec`, 2026-11-25 20:57 UTC, 2026 verisi indirilmeden önce).
Dondurulmuş ayar: [`docs/onkayit/2026-11-25_ege_dondurulmus_ayar.json`](../docs/onkayit/2026-11-25_ege_dondurulmus_ayar.json)
(commit `5206493`, 22:10 UTC, 2026 puanından önce). Ham sonuç:
[`docs/onkayit/2026-11-25_ege_sonuc.json`](../docs/onkayit/2026-11-25_ege_sonuc.json).

## Sonuç: ipucu desteklenmedi (ön kayıtlı karar: "kaldı")

| | olay | olay başına kazanç (ETAS'a göre) | %95 GA | tek yönlü p | karar |
|---|---:|---:|---|---:|---|
| **Ege, M≥3, 2026-01-01 → 2026-11-25 (birincil)** | 312 | **−0,056** | [−0,336, +0,250] | 0,69 | **kaldı** |
| Ege, M≥4 (ikincil) | 24 | −0,111 | [−0,738, +0,453] | 0,72 | – |
| Ege, 2,5 ≤ M < 3 (ikincil, küçük depremler) | 559 | −0,055 | [−0,332, +0,283] | 0,67 | – |

Model, hiç görmediği 2026 Ege depremlerine ETAS'tan biraz **daha düşük** olasılık verdi.
Fark sıfırdan ayırt edilemiyor, ama ön kayıtlı kurala göre "kazanç ≤ 0" olduğu için ipucu
desteklenmedi. Küçük depremlerde de aynı yön çıktı.

## Eski sayı neden yanıltıcıydı?

PR #13'teki +0,315, puanın "toplam beklenen sayı" kısmını tüm Türkiye için hesaplıyordu.
Model Ege dışında toplam sayıyı ETAS'tan düşük tahmin ettiği için bu terim her olaya (Ege'dekilere
de) artı yazıyordu. Bölgeye kısıtlı doğru ölçütle aynı 2020-2025 verisi:

| 2020-2025 (eski veri, yeniden üretim) | olay | eski ölçüt | bölgeye kısıtlı ölçüt | %95 GA | p |
|---|---:|---:|---:|---|---:|
| Ege M≥3 | 4763 | +0,264 (PR #13: +0,315) | +0,098 | [−0,027, +0,173] | 0,07 |
| Ege M≥4 | 562 | +0,236 | +0,113 | [−0,059, +0,179] | 0,11 |

Yani üstünlüğün büyük kısmı ölçüm biçiminden geliyordu; doğru ölçütle eski veride bile anlamlı
değildi. 2026'da eski ölçüt yine +0,396 veriyor, bölgeye kısıtlı ölçüt −0,056: aynı yanılsama.

## Ege dışı (ikincil)

Ege dışındaki Türkiye'de 2026 için bölgeye kısıtlı kazanç +0,43 (482 olay, p < 0,001). Bunun
neredeyse tamamı modelin toplam sayıyı ETAS'tan düşük tahmin etmesinden geliyor (beklenen 1141'e
karşı 1481; gerçekleşen 482): 2025'teki büyük diziler (Santorini, Sındırgı) sonrası ETAS 2026'yı
fazla tahmin ediyor, model bu fazlayı kısmen düzeltiyor. Bu bilinen bir ETAS sorunudur (PR #10,
#11) ve ön kayıtlı bir hipotez değildi; ayrıca sınanmadan bulgu sayılmamalı.

## Sınırlar ve sapmalar

- **Katalog:** AFAD'ın ISC arşivi ve açık API'si 2026'nın Ocak-Ekim aylarını vermedi; 2026 M≥3
  olayları büyük ölçüde EMSC'den (5884 EMSC, 846 AFAD, 1 USGS kaydı, tüm büyüklükler). Önceki
  yıllarda AFAD ağırlıklıydı. İki model aynı olaylarla puanlandığı için karşılaştırma adil, ama
  büyüklük ölçeği farkı küçük bir etki yapabilir. Ön kayıt bu durumu öngörüyordu.
- **Ayar yeniden seçildi:** PR #13'ün ayar dosyası depoda yoktu; aynı arama (300 deneme, tohum 0,
  2015-2019) yeniden yapıldı. Olay sayıları birebir aynı çıktı (4763 / 562), ama seçilen ayar
  biraz farklı (2020-2025 tüm Türkiye puanı −0,059; PR #13'te +0,049).
- Tek dönem (11 ay) ve geniş güven aralığı: sonuç "Ege'de üstünlük yok" demekten çok
  "üstünlüğe dair kanıt yok" demektir.

## Canlı test

Ön kaydın 8. bölümündeki günlük canlı tahmin kuruldu: [`forecasts/ege/BENIOKU.md`](../forecasts/ege/BENIOKU.md).
İlk tahmin (2026-11-26 günü) bu PR'da; sonrakiler PR birleşince her gün otomatik yazılır.
Değerlendirme 2027-05-26 (ara) ve 2027-11-26 (asıl).

## Yeniden üretme

```bash
python scripts/fetch_catalogs.py --end 2026-11-25 --afad-days 328
python scripts/train.py --config configs/train_ege_2026.toml --until features
python scripts/ege_onkayit.py secim
python scripts/ege_onkayit.py sinama
```
