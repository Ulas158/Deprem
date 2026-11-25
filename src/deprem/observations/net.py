"""Önbellekli indirme yardımcıları (tüm gözlem kaynakları için ortak)."""

from __future__ import annotations

import logging
import time
from pathlib import Path

import requests

log = logging.getLogger(__name__)
UA = {"User-Agent": "deprem-research-prototype/0.1"}


def download(url: str, dest: Path, retries: int = 4, timeout: int = 180,
             params: dict | None = None, min_bytes: int = 1) -> Path | None:
    """``url``'yi ``dest``'e indir; dosya zaten varsa indirme. 404'te None döner."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size >= min_bytes:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    for k in range(retries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=timeout)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(r.content)
            tmp.rename(dest)
            return dest
        except requests.RequestException as e:
            wait = 2 ** k
            log.warning("%s indirilemedi (%s), %d s sonra tekrar", url, e, wait)
            time.sleep(wait)
    return None
