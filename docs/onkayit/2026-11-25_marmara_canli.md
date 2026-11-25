# Ön kayıt: Marmara canlı (ileriye dönük) tahmin testi

- **Kayıt tarihi:** 2026-11-25 (UTC), ilk Marmara tahmininden önce. Git commit zamanı kanıttır.
- İstek: Ulas, "Marmara'yı da yap" (Ege canlı testinin aynısı Marmara için).

> Deneysel araştırma çıktısıdır, deprem tahmini değildir.

## Hipotez

**H1:** Projenin dondurulmuş modeli, Marmara kutusunda ileriye dönük günlük tahminlerde M≥3
depremlere ETAS'tan daha yüksek olasılık verir (olay başına bilgi kazancı > 0).
**H0:** Kazanç ≤ 0.

**Önceden beklenti (dürüst not):** PR #13'te 2020-2025 Marmara'da model ETAS'ın gerisindeydi
(eski ölçütle M≥3 −0,131). Ayrıca 2026 ön kayıtlı Ege sınamasında, Ege dışı Türkiye'nin (Marmara
dahil) toplu sonucu görüldü; Marmara'ya ait ayrı bir 2026 sayısı hesaplanmadı. Bu yüzden Marmara
için geriye dönük bir 2026 sınaması ön kayıtlı sayılmaz; karar yalnız bu canlı testle verilir.

## Bölge, modeller, veri

- **Bölge:** 0,5° hücreler; hücre merkezi enlem 40,0-41,3 K, boylam 25,5-30,5 D
  (`src/deprem/regional.py`, "Marmara (KAF batı)"; 30 hücre).
- **Modeller:** Ege ön kaydının (docs/onkayit/2026-11-25_ege_ustunluk.md) 4. ve 8. bölümlerindeki
  ile birebir aynı: günlük ETAS (2026 parametreleri, donmuş) ve model = günlük ETAS × hücrenin
  donmuş 2026 çarpanı (`forecasts/ege/model_carpani_2026.csv`, Marmara satırları). Hiçbir
  parametre bu bölge için ayrıca ayarlanmadı ve dönem boyunca değiştirilmez.
- **Kayıt:** ilk pencere 2026-11-26; sonrakiler her gün GitHub Actions ile
  (`forecasts/marmara/gunluk/`, `ege-canli-kayit` dalı).
- **Puanlama kataloğu:** değerlendirme günündeki birleşik proje kataloğu (AFAD > EMSC > USGS).

## Ölçüt ve karar kuralı (Ege ile aynı)

- Birincil: M≥3, gün × hücre Poisson olasılıklarıyla olay başına bilgi kazancı
  (Σ n·ln(λ_model/λ_ETAS) − Σ(λ_model − λ_ETAS)) / Σ n.
- Belirsizlik: hücre × takvim ayı blokları bootstrap'ı (10 000, tohum 2026); tek yönlü p = P(IG ≤ 0).
- **Geçti:** IG > 0 ve p < 0,05 (IG ≥ 0,10 ise "güçlü"). **Kaldı:** IG ≤ 0.
  **Kararsız:** IG > 0 ve p ≥ 0,05, ya da 30'dan az M≥3 olay.
- Ara değerlendirme 2027-05-26, asıl değerlendirme 2027-11-26. Öncesindeki bakışlar "ara bakış"tır.
- Marmara'da yılda M≥3 olay sayısı Ege'den az olduğu için 12 ayda "kararsız" çıkması olasıdır;
  bu durumda kayıt aynı kurallarla sürdürülür ve 24. ayda (2028-11-26) yeniden bakılır.
- İkincil (karar vermez): M≥4; günlük tahminlerin hücre × ay toplamıyla aynı ölçüt.

```bash
python scripts/ege_canli.py degerlendir --bolge marmara --yaz
```
