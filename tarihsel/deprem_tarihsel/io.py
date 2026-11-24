"""Üretilmiş tarihsel kataloğu okumak için yardımcılar."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"

CORE_COLUMNS = [
    "event_id", "time", "latitude", "longitude", "depth_km", "mag", "mag_type", "source", "agency",
    "time_err_s", "horiz_err_km", "depth_err_km", "mag_err",
]


def load(preferred_only: bool = True, core_only: bool = False, path: Path | None = None) -> pd.DataFrame:
    """Kataloğu okur; ``time`` UTC datetime olarak döner.

    preferred_only=False ise bütün kaynaklardaki bütün kayıtlar (yinelemeler dahil) gelir.
    core_only=True ise yalnızca aletsel katalogla ortak sütunlar döner (doğrudan birleştirme için).
    """
    if path is None:
        path = DATA / ("tarihsel_katalog.csv" if preferred_only else "tarihsel_tum_kayitlar.csv")
    df = pd.read_csv(path)
    for c in ("time", "time_min", "time_max"):
        if c in df:
            # MÖ yıllar ("-0281-...") nedeniyle ISO8601 biçimi açıkça verilir.
            df[c] = pd.to_datetime(df[c], utc=True, format="ISO8601")
    return df[CORE_COLUMNS] if core_only else df
