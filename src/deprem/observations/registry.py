"""Gözlem özellikleri kaydı: kanıt grubu, dönüşüm ve literatüre dayalı öncül etki.

Her özellik için:

* ``evidence``  "güçlü" veya "zayıf" (zayıf = test edilecek hipotez)
* ``horizon``   "yıllık" (hücre × yıl modeli, ETAS ofsetli) veya "kısa"
                (bölge × 5 gün modeli; öncül sinyal hipotezleri)
* ``transform`` ham değere uygulanan dönüşüm (log10, log1p, clip3 = işaretli log)
* ``prior_mean`` başlangıç etkisi: dönüştürülmüş ve standartlaştırılmış
                özellik 1 standart sapma arttığında deprem oranının
                logaritmasındaki beklenen değişim. Örn. +0,3 → oran ×1,35.
* ``prior_sd``  bu başlangıç tahmininin belirsizliği. Geniş = veriye bırak.

Öncüller, eğitimde Gauss öncül olarak kullanılır (ridge'in sıfır yerine
literatür değerine çekilen hali). Eğitim verisi güçlüyse öncülü ezer;
zayıfsa sonuç öncüle yakın kalır. Bu yüzden asıl ölçüt öncül değil,
**sınama yıllarında ETAS'a karşı bilgi kazancıdır**. Zayıf kanıtlı
özelliklerin öncülleri bilerek küçük tutulmuştur: literatürdeki iddia edilen
etkiler (olasılık kazancı 2-5) bağımsız testlerde tekrarlanamamıştır.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ObsFeature:
    name: str
    group: str
    evidence: str       # "güçlü" | "zayıf"
    horizon: str        # "yıllık" | "kısa"
    unit: str
    transform: str      # "none" | "log10" | "log1p" | "clip3"
    prior_mean: float
    prior_sd: float
    status: str         # "veri var" | "pilot" | "veri bekliyor"
    description: str
    rationale: str
    source: str


OBS_FEATURES: dict[str, ObsFeature] = {}


def reg(*a, **k):
    f = ObsFeature(*a, **k)
    OBS_FEATURES[f.name] = f


# ---------------------------------------------------------------- GÜÇLÜ: jeodezi
reg("gnss_strain_rate_ns", "GPS", "güçlü", "yıllık", "nanogerinim/yıl (log10)", "log10", 0.30, 0.20,
    "veri var", "GPS hızlarından ikinci değişmez gerinim hızı (son 5 yıl)",
    "Uzun dönem deprem oranı gerinim hızıyla orantılıdır (Bird & Kreemer 2015). ETAS arka planı "
    "bunun bir kısmını zaten içerdiğinden kalan etki orta düzeyde varsayıldı.",
    "Nevada Geodetic Lab (Blewitt vd. 2018)")
reg("gnss_max_shear_ns", "GPS", "güçlü", "yıllık", "nanogerinim/yıl (log10)", "log10", 0.15, 0.20,
    "veri var", "En büyük makaslama gerinim hızı",
    "Doğrultu atımlı faylarda (KAF, DAF) yüklenmeyi doğrudan ölçer; gerinim hızıyla ilişkili olduğundan "
    "öncül yarıya bölündü.", "Nevada Geodetic Lab")
reg("gnss_dilatation_ns", "GPS", "güçlü", "yıllık", "nanogerinim/yıl", "clip3", 0.0, 0.15,
    "veri var", "Hacim değişim hızı (+ genişleme, − sıkışma)",
    "İşaret bölgeye göre değişir (Ege açılması, Doğu sıkışması); öncül nötr.", "Nevada Geodetic Lab")
reg("gnss_transient_z", "GPS", "güçlü", "yıllık", "-", "log1p", 0.10, 0.15,
    "veri var", "Son 1 yılda uzun dönem hızdan sapmanın gürültüye oranı (75 km içindeki en büyük)",
    "Yavaş kayma/öncü hızlanma bazı büyük depremlerden önce gözlendi (Ruiz vd. 2014, Kato vd. 2012); "
    "genel geçerliliği tartışmalı, öncül küçük.", "Nevada Geodetic Lab")
reg("gnss_transient_mm_yr", "GPS", "güçlü", "yıllık", "mm/yıl", "log1p", 0.05, 0.15,
    "veri var", "Son 1 yıldaki hız sapmasının büyüklüğü (ortanca)", "Yukarıdakiyle aynı.",
    "Nevada Geodetic Lab")
reg("insar_los_rate_mm_yr", "InSAR", "güçlü", "yıllık", "mm/yıl", "clip3", 0.0, 0.20,
    "veri bekliyor", "Uydu radarı görüş doğrultusu hızı",
    "Krip eden fay kesimlerini ve yerel gerinimi gösterir.", "COMET LiCSAR / LiCSBAS")
reg("insar_strain_rate_ns", "InSAR", "güçlü", "yıllık", "nanogerinim/yıl (log10)", "log10", 0.20, 0.20,
    "veri bekliyor", "InSAR+GPS birleşik gerinim hızı",
    "GPS'ten daha yüksek çözünürlüklü yüklenme haritası (Weiss vd. 2020).", "COMET")
reg("gonaf_strain_anom", "Kuyu gerinimi", "güçlü", "yıllık", "-", "log1p", 0.10, 0.20,
    "veri bekliyor", "GONAF kuyu gerinimölçer anomalisi (Marmara)",
    "Nanogerinim duyarlılığında yavaş deformasyon; Marmara'da yavaş kayma olayları kaydedildi "
    "(Martínez-Garzón vd. 2019).", "GFZ / ICDP GONAF (erişim başvurusu gerekir)")
# ---------------------------------------------------------------- GÜÇLÜ: sismik
reg("rep_slip_rate_cm_yr", "Tekrarlayan depremler", "güçlü", "yıllık", "cm/yıl", "log1p", 0.15, 0.20,
    "veri var", "Tekrarlayan deprem ailelerinden derin kayma hızı (son 5 yıl, Nadeau-Johnson)",
    "Fayın krip eden kısmı çevresindeki kilitli yamaları yükler (Uchida & Bürgmann 2019).",
    "Birleşik katalog (dalga biçimi doğrulaması isteğe bağlı)")
reg("rep_accel_log", "Tekrarlayan depremler", "güçlü", "yıllık", "log oran", "clip3", 0.10, 0.15,
    "veri var", "Son 1 yıl / önceki 4 yıl tekrarlayan kayma oranı",
    "Tohoku öncesi tekrarlayan depremlerde hızlanma (Uchida vd. 2016).", "Birleşik katalog")
reg("rep_n_families", "Tekrarlayan depremler", "güçlü", "yıllık", "-", "log1p", 0.05, 0.15,
    "veri var", "30 km içindeki tekrarlayan aile sayısı", "Krip eden fay alanının göstergesi.",
    "Birleşik katalog")
reg("dvv_change_pct", "Sismik hız", "güçlü", "yıllık", "%", "clip3", -0.10, 0.15,
    "pilot", "Son 6 ay dv/v − önceki 2 yıl (ortam gürültüsü, tek istasyon)",
    "Hız düşüşü hasar/akışkan basıncı artışını gösterir (Brenguier vd. 2008); negatif öncül "
    "= hız düştükçe oran artar.", "KOERI KO ağı (FDSN)")
reg("mlcat_rate_ratio", "Yapay zekâ kataloğu", "güçlü", "yıllık", "-", "log1p", 0.10, 0.20,
    "pilot", "PhaseNet ile bulunan olay sayısı / AFAD olay sayısı",
    "Küçük depremler gerilme durumunu daha ayrıntılı gösterir (Ross vd. 2019).", "KOERI + PhaseNet")
# ---------------------------------------------------------------- ZAYIF: hipotezler (kısa vade)
reg("tec_z_max", "İyonosfer TEC", "zayıf", "kısa", "z", "clip3", 0.10, 0.15, "veri var",
    "Önceki 10 günde günlük en büyük pozitif TEC anomalisi (yerel − kuşak)",
    "Liu vd. 2006 öncesinde pozitif anomali bildirdi; Dautermann vd. 2007 doğrulayamadı.", "CODE GIM")
reg("tec_z_min", "İyonosfer TEC", "zayıf", "kısa", "z", "clip3", -0.05, 0.15, "veri var",
    "Önceki 10 günde en büyük negatif TEC anomalisi",
    "Negatif anomali raporları da var (Liu vd. 2009); negatif öncül = daha düşük TEC'de oran artar.",
    "CODE GIM")
reg("tec_frac_pos", "İyonosfer TEC", "zayıf", "kısa", "oran", "none", 0.05, 0.15, "veri var",
    "Önceki 10 günde z>2 olan saatlerin oranı", "Süreklilik ölçütü.", "CODE GIM")
reg("ulf_zh_rel_z", "Manyetik ULF", "zayıf", "kısa", "z", "clip3", 0.10, 0.15, "veri var",
    "Önceki 10 günde gece Z/H polarizasyon oranı anomalisinin en büyüğü (İznik/PAG)",
    "Hayakawa vd. 1996; Masci 2011 jeomanyetik etkinlikle açıkladı.", "INTERMAGNET IZN, PAG")
reg("ulf_z_rel_z", "Manyetik ULF", "zayıf", "kısa", "z", "clip3", 0.05, 0.15, "veri var",
    "Önceki 10 günde gece Z bileşeni gücü anomalisinin en büyüğü", "Aynı hipotez.", "INTERMAGNET")
reg("thermal_dts_z", "Termal", "zayıf", "kısa", "z", "clip3", 0.10, 0.15, "veri var",
    "Önceki 10 günde gece yüzey sıcaklığı anomalisinin en büyüğü (RETIRA benzeri)",
    "Tronin 2006, Tramutoli vd. 2005; bulut ve hava etkisiyle karışır.", "NASA POWER (MERRA-2)")
reg("thermal_dsa_z", "Termal", "zayıf", "kısa", "z", "clip3", 0.05, 0.15, "veri var",
    "Önceki 10 günde yüzey − hava sıcaklık farkı anomalisinin en büyüğü",
    "Hava sıcaklığından bağımsız yüzey ısınmasını ayırmak için.", "NASA POWER")


def table() -> pd.DataFrame:
    return pd.DataFrame([vars(f) for f in OBS_FEATURES.values()])


def apply_transform(name: str, x):
    f = OBS_FEATURES.get(name)
    x = np.asarray(x, dtype=float)
    if f is None or f.transform == "none":
        return x
    if f.transform == "log10":
        return np.log10(np.clip(np.abs(x), 1e-3, None))
    if f.transform == "log1p":
        return np.log1p(np.clip(x, 0, None))
    if f.transform == "clip3":
        # işaretli log: uç değerleri bastırır, satırdan bağımsızdır (sızıntı yok)
        return np.sign(x) * np.log1p(np.abs(x))
    return x


def prior_for(column: str) -> tuple[float, float] | None:
    """Tasarım matrisi sütunu için (ortalama, sd); bilinmeyen sütun → None."""
    if column.endswith("__eksik"):
        return None
    f = OBS_FEATURES.get(column)
    return (f.prior_mean, f.prior_sd) if f else None
