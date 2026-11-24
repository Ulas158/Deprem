"""En eski veriden bugüne deprem deprem ileri eğitim.

    python scripts/train_long.py --config configs/train_long.toml
    python scripts/train_long.py --config configs/train_long.toml --summary-only

Ayrıntı ve yöntem: src/deprem/longwalk.py. Çıktılar runs/<run_name>/longwalk/,
özet runs/<run_name>/longwalk.json ve reports/training_<run_name>.md.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from deprem import longwalk
from deprem.pipeline import Run

ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=ROOT / "configs/train_long.toml")
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--replicates", type=int, default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
    run = Run(args.config, ROOT)
    cat_path = ROOT / run.cfg["data"]["catalog"]
    if not cat_path.exists():
        longwalk.build_long_catalog(ROOT, cat_path)
    if args.replicates is not None:
        run.cfg["longwalk"]["replicates"] = args.replicates
    if not args.summary_only:
        for st in ("coulomb", "budget", "observations"):
            if not run.done(st):
                logging.info("[%s] başlıyor", st)
                run.mark(st, **(getattr(run, f"stage_{st}")() or {}))
        longwalk.run_all(run, ROOT, workers=run.cfg["longwalk"].get("workers", 3))
    longwalk.summarize(run.dir, run.cfg["longwalk"])
    from deprem.report import write_long_report
    write_long_report(run, ROOT)


if __name__ == "__main__":
    main()
