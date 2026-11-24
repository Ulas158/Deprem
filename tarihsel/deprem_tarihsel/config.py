"""Bütün ayarlanabilir parametreler tek yerde.

Buradaki varsayılanlar bilimsel olarak doğrulanmış değerler değildir; ilk
yaklaşım olarak seçilmiştir ve ileride geriye dönük testlerle ayarlanmalıdır.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Config:
    # --- Bölge ve dönem ---------------------------------------------------
    lat_min: float = 34.0
    lat_max: float = 44.0
    lon_min: float = 24.0
    lon_max: float = 46.0
    # Bu yıldan itibaren aletsel kataloglar (AFAD vb.) kullanılır.
    max_year: int = 1899

    # --- Bilinmeyen belirsizlikler için atanan 1σ değerleri ----------------
    # Kaynak belirsizlik vermediğinde kullanılır; *_assigned sütunu True olur.
    default_mag_err: dict = field(
        default_factory=lambda: {
            "EPICA": 0.5,
            "SHARE_CET": 0.5,
            "GEM_GHEC": 0.4,  # Ambraseys & Jackson (1998) sınıf genişliği ±0.4
            "NOAA_NCEI": 0.5,
        }
    )
    default_horiz_err_km: dict = field(
        default_factory=lambda: {
            "EPICA": 30.0,
            "SHARE_CET": 50.0,
            "GEM_GHEC": 30.0,
            "NOAA_NCEI": 50.0,
        }
    )
    # Ms -> Mw dönüşümünün kendi saçılması (Scordilis 2006 civarı), karesel eklenir.
    ms_to_mw_err: float = 0.17
    # Türü belirsiz büyüklük Mw sayıldığında eklenen ek belirsizlik.
    unknown_mag_type_err: float = 0.3

    # --- Yineleme (kataloglar arası aynı olay) eşleştirmesi ----------------
    dup_max_km: float = 100.0
    # Tarih aralıkları çakışmıyorsa izin verilen en büyük boşluk (gün).
    # Jülyen/Gregoryen farkı (10-12 gün) ve kaynaklar arası tarih kaymaları için 15 gün.
    dup_time_slack_days: float = 15.0
    # Büyüklükler ikisi de varsa izin verilen en büyük fark. Tarihsel büyüklük tahminleri
    # kaynaklar arasında 1 birimden fazla ayrışabiliyor (ör. 1668 Anadolu: 6.8 - 7.9).
    dup_max_dmag: float = 1.5
    # Tercih sırası: önce gelen kaynak birleşik katalogda tutulur.
    source_priority: tuple = ("EPICA", "SHARE_CET", "GEM_GHEC", "NOAA_NCEI")

    # --- Ağırlık formülü ---------------------------------------------------
    # w_unc = max(w_floor, 1 / (1 + (σM/sM)² + (σloc/sL)² + (σz/sZ)² + (σt/sT)²))
    # Aletsel katalog modülüyle (src/deprem/catalog/weights.py, WeightConfig.scales) aynı
    # biçim ve ölçekler. Derinlik terimi yalnızca kaynak derinlik belirsizliği verdiyse eklenir.
    scale_mag: float = 0.3
    scale_horiz_km: float = 20.0
    scale_depth_km: float = 30.0
    scale_time_days: float = 30.0
    w_floor: float = 1e-3
    # Büyüklüğü olmayan kayıtlarda büyüklük terimi için kullanılan σM.
    missing_mag_err: float = 1.0
    # Tamlık penceresi dışındaki kayıtların çarpanı.
    incomplete_factor: float = 0.3
    # Konumu kuşkulu kayıtların çarpanı (ör. NOAA'da bölge dışı ülke adı).
    suspect_factor: float = 0.0

    # --- Tamlık (Stepp 1972) ----------------------------------------------
    completeness_mag_bins: tuple = (5.5, 6.0, 6.5, 7.0)
    # Stepp başarısız olursa kullanılacak muhafazakâr varsayım (Mw alt sınırı: başlangıç yılı).
    default_completeness: dict = field(
        default_factory=lambda: {5.5: 1880, 6.0: 1850, 6.5: 1700, 7.0: 1500}
    )
    stepp_window_years: int = 25
    stepp_tolerance_sigma: float = 2.0
    stepp_min_events: int = 5


DEFAULT = Config()
