"""Deprem oluşumuna dair yeni gözlem türleri.

Bu paket, katalog dışı gözlemleri (jeodezi, sismik dalga biçimi türevleri,
atmosfer/iyonosfer) indirir, işler ve modelin kullanacağı özelliklere çevirir.

İki kanıt grubu vardır (bkz. ``registry.OBS_FEATURES``):

* ``güçlü``  fiziksel mekanizması ve deprem öncesi/sırası gözlemleri hakemli
  literatürde tekrar tekrar gösterilmiş olanlar (GNSS, InSAR, kuyu gerinimi,
  tekrarlayan depremler, sismik hız değişimi, zenginleştirilmiş katalog).
* ``zayıf``  öncül sinyal olarak önerilmiş ama bağımsız, ileriye dönük
  testlerde tutarlı biçimde doğrulanamamış olanlar (iyonosfer TEC, ULF
  manyetik, termal kızılötesi). Bunlar **test edilecek hipotez** olarak ayrı
  tutulur ve ayrı raporlanır.

Her özellik için literatüre dayalı bir **başlangıç etkisi** (öncül ortalama)
ve **belirsizlik** (öncül standart sapma) tanımlıdır; model bunları geçmiş
depremlerle eğitirken günceller. Öncül, kanıt değildir.
"""
