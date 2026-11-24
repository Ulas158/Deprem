"""Periyodik yüklerle (gelgit, mevsim) deprem zamanlaması ilişkisi.

Schuster testi: N olayın fazları φ_i için D² = (Σcos φ)² + (Σsin φ)²;
fazlar rastgeleyse p = exp(−D²/N). Küçük p, olayların belli bir faza
yığıldığını gösterir. **Artçı kümeleri p'yi yapay olarak küçültür**; bu
yüzden test artçıları ayıklanmış (bağımsız) olaylarla yapılmalıdır.

Literatür: Tanaka vd. 2002, 2012 (Sumatra ve Tohoku öncesi gelgit
duyarlılığı artışı); Cochran vd. 2004 (sığ bindirmelerde gelgit
tetiklemesi); Ide vd. 2016 (büyük gelgit genliğinde büyük deprem oranı);
Christiansen vd. 2007, Bettinelli vd. 2008, Johnson vd. 2017 (mevsimsel
hidrolojik yükleme ile mevsimsel sismisite). Etkiler küçüktür (oran
modülasyonu tipik olarak %1–%10) ve çoğu bölgede istatistiksel olarak
anlamlı değildir.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .astro import (annual_phase, fortnightly_phase, local_hour_phase, lunar_semidiurnal_phase,
                    semidiurnal_phase)


def schuster(phase_deg) -> tuple[float, float, int]:
    """(log10 p, ortalama faz [derece], N). N < 10 ise (nan, nan, N)."""
    ph = np.radians(np.asarray(phase_deg, dtype=float))
    ph = ph[np.isfinite(ph)]
    n = len(ph)
    if n < 10:
        return np.nan, np.nan, n
    c, s = np.cos(ph).sum(), np.sin(ph).sum()
    d2 = c * c + s * s
    return float(-d2 / n / np.log(10)), float(np.degrees(np.arctan2(s, c))), n


def event_phases(events: pd.DataFrame, tidal: bool = True) -> pd.DataFrame:
    """Her olay için gelgit, yıllık ve günlük (kontrol) fazlar.

    * ``phase_lunar``: yalnızca Ay'ın yarı günlük fazı. Ana gelgit testi budur.
    * ``phase_semidiurnal``: Ay + Güneş birleşik gelgit eğrisinde faz (Tanaka 2002).
      Güneş bileşeni tam 12 saatlik olduğundan gece/gündüz algılama farkı
      (küçük depremlerin gece daha iyi algılanması, gündüz taş ocağı patlamaları)
      bu fazda **sahte gelgit sinyali** üretir; yalnızca Mc'nin yeterince üstündeki
      olaylarda kullanılmalıdır.
    * ``phase_hour``: yerel güneş saati. Gelgitle ilgisi yoktur; katalog kalitesi
      kontrolüdür. Burada güçlü sinyal varsa birleşik gelgit sonucu güvenilmez.
    """
    out = pd.DataFrame(index=events.index)
    out["phase_lunar"] = lunar_semidiurnal_phase(events["time"], events["longitude"].to_numpy())
    out["phase_hour"] = local_hour_phase(events["time"], events["longitude"].to_numpy())
    out["phase_annual"] = annual_phase(events["time"])
    out["phase_fortnight"] = fortnightly_phase(events["time"])
    if tidal:
        out["phase_semidiurnal"] = semidiurnal_phase(
            events["time"], events["latitude"].to_numpy(), events["longitude"].to_numpy())
    else:
        out["phase_semidiurnal"] = np.nan
    return out


def grouped_schuster(phases: pd.Series, groups: pd.Series, prefix: str) -> pd.DataFrame:
    rows = {}
    for g, ph in phases.groupby(groups):
        lp, mean, n = schuster(ph.to_numpy())
        rows[g] = {f"{prefix}_log10p": lp, f"{prefix}_phase": mean, f"{prefix}_n": n}
    return pd.DataFrame.from_dict(rows, orient="index")


def modulation_factor(phase_deg, gain: float, preferred_deg: float = 0.0):
    """Kısa vadeli modeller için oran çarpanı: 1 + gain·cos(φ − φ0) (≥ 0)."""
    return np.clip(1 + gain * np.cos(np.radians(np.asarray(phase_deg) - preferred_deg)), 0, None)
