# Ege canlı (ileriye dönük) tahmin kaydı

Ön kayıt: [`docs/onkayit/2026-11-25_ege_ustunluk.md`](../../docs/onkayit/2026-11-25_ege_ustunluk.md), bölüm 8.

> Deneysel araştırma çıktısıdır, deprem tahmini veya uyarısı değildir.

## Ne kaydediliyor?

Aynı düzenle Marmara için de kayıt tutulur (`forecasts/marmara/`).

Her gün 21:17 UTC'de (Türkiye saatiyle 00:17) GitHub Actions
(`.github/workflows/ege_canli.yml`) ertesi UTC günü için Ege / Batı Anadolu kutusundaki
(enlem 36,5-40 K, boylam 25,5-30,5 D) 70 hücrenin her biri için beklenen deprem sayısını
iki model için yazar:

| sütun | anlamı |
|---|---|
| `etas_m3`, `etas_m4` | ETAS: bilinen en iyi basit model (artçı kümelenmesi + arka plan) |
| `model_m3`, `model_m4` | Projenin modeli: aynı günlük ETAS × hücrenin dondurulmuş çarpanı |

Dosyalar `ege-canli-kayit` dalında `forecasts/ege/gunluk/YYYY/YYYY-MM-DD.csv` olarak tutulur;
her dosyanın ilk satırı tahminin verildiği anı ve kullanılan son depremi yazar. Dosyalar
yalnız eklenir, sonradan değiştirilmez. Git commit zamanı, tahminin pencere başlamadan önce
verildiğinin kanıtıdır.

## Dondurulmuş parçalar (değişmez)

- `model_carpani_2026.csv`: hücre başına model / ETAS oranı, ön kayıtlı 2026 sınamasında
  2012-2025 verisiyle eğitilen modelden. 2027'de de aynısı kullanılır.
- `etas_2026.json`, `etas_arka_plan_2026.csv`: 2010-2025 verisiyle uydurulan ETAS parametreleri
  ve arka plan dağılımı.
- `katalog_canli.csv.gz`: tahmin verirken kullanılan deprem listesi (M≥3, Türkiye kutusu).
  Her gün son 5 günün depremleri AFAD (erişilirse), EMSC ve USGS'den yeniden çekilip eklenir.

## Günlük ETAS nasıl hesaplanıyor?

Arka plan hızı + o ana kadar bilinen her depremin ertesi 24 saatteki doğrudan artçılarının
beklenen sayısı (Omori-Utsu integrali). Artçılar uzayda ETAS'ın kendi çekirdeğiyle hücrelere
dağıtılır (sabit tohumlu Monte Carlo, 400 örnek). Pencere içinde doğacak depremlerin kendi
artçıları sayılmaz; bu iki model için aynı olduğundan karşılaştırmayı etkilemez.

## Değerlendirme

```bash
python scripts/fetch_catalogs.py --end <tarih> --afad-days <2026-11-26'dan beri gün>
python scripts/ege_canli.py degerlendir --bolge ege --yaz
```

- En erken 2027-05-26 (6 ay), asıl değerlendirme 2027-11-26 (12 ay).
- Puanlama kataloğu: değerlendirme günü projenin birleşik kataloğu (AFAD > EMSC > USGS).
- Birincil ölçüt: olay başına bilgi kazancı, gün × hücre Poisson olasılıklarıyla
  (Σ n·ln(λ_model/λ_ETAS) − Σ(λ_model − λ_ETAS)) / Σ n; M≥3.
  İkincil: günlük tahminler hücre × ay toplanarak aynı ölçüt; M≥4.
- Belirsizlik: hücre × takvim ayı blokları bootstrap'ı (10 000).
- Karar kuralı ön kayıttaki gibi: kazanç ≤ 0 → kaldı; > 0 ve p < 0,05 → geçti; arası kararsız.
- 6 aydan önceki her bakış "ara bakış" diye işaretlenir ve karar sayılmaz.
