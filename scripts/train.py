"""Eğitim hattını çalıştır.

    python scripts/train.py --config configs/train_default.toml
    python scripts/train.py --config configs/train_default.toml --force fit evaluate report
    python scripts/train.py --config configs/train_default.toml --force all

Önce katalog hazır olmalı: python scripts/fetch_catalogs.py
Çıktılar runs/<run_name>/ altında; rapor reports/training_<run_name>.md.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from deprem.pipeline import STAGES, Run

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "configs/train_default.toml")
    ap.add_argument("--force", nargs="*", default=[], choices=STAGES + ["all"],
                    help="yeniden hesaplanacak aşamalar")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    Run(args.config, ROOT, force=args.force).run()


if __name__ == "__main__":
    main()
