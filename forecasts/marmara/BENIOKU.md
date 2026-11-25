# Marmara canlı (ileriye dönük) tahmin kaydı

Ön kayıt: [`docs/onkayit/2026-11-25_marmara_canli.md`](../../docs/onkayit/2026-11-25_marmara_canli.md).

> Deneysel araştırma çıktısıdır, deprem tahmini veya uyarısı değildir.

Ege kaydıyla aynı düzen ([`forecasts/ege/BENIOKU.md`](../ege/BENIOKU.md)): her gün ertesi UTC
günü için Marmara kutusundaki (enlem 40,0-41,3 K, boylam 25,5-30,5 D) 30 hücrenin M≥3 ve M≥4
beklenen sayıları, ETAS ve model için, `gunluk/YYYY/YYYY-MM-DD.csv` dosyalarına yazılır.
Model parçaları (ETAS parametreleri, hücre çarpanları, katalog) `forecasts/ege/` altında
ortaktır ve donmuştur.

Değerlendirme: `python scripts/ege_canli.py degerlendir --bolge marmara --yaz`
(ara 2027-05-26, asıl 2027-11-26; olay azsa 2028-11-26).
