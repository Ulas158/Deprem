"""İstasyon çifti gürültü korelasyonları ve dv/v (Doğu Marmara pilotu).

    python scripts/fetch_dvv_pairs.py ccf  --start 2009-01-01 --end 2026-11-20   # indirme (uzun, sürdürülebilir)
    python scripts/fetch_dvv_pairs.py dvv                                        # dv/v serileri

Çıktılar: data/obs/dvv_pairs/ (günlük korelasyonlar), data/obs/dvv_pairs_series.csv.gz (çift başına),
data/obs/dvv_pairs_regional.csv (bölgesel ortalama).
"""
import argparse
import logging
from pathlib import Path

from deprem.observations import dvv_pairs as dp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "obs" / "dvv_pairs"

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("task", choices=["ccf", "dvv"])
    ap.add_argument("--start", default="2009-01-01")
    ap.add_argument("--end", default="2026-11-20")
    ap.add_argument("--step-days", type=int, default=2)
    ap.add_argument("--hours", type=float, default=4.0)
    ap.add_argument("--hour-ut", type=int, default=20)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--stack-days", type=int, default=60)
    a = ap.parse_args()
    if a.task == "ccf":
        dp.run_ccfs(OUT, dp.MARMARA_EAST, a.start, a.end, a.step_days, a.hour_ut, a.hours, a.workers)
    else:
        pairs = dp.all_pairs_dvv(OUT, a.stack_days)
        pairs.to_csv(OUT.parent / "dvv_pairs_series.csv.gz", index=False)
        reg = dp.regional_series(pairs)
        reg.to_csv(OUT.parent / "dvv_pairs_regional.csv", index=False)
        logging.info("%d çift, %d zaman; bölgesel ortanca belirsizlik %%%.3f",
                     pairs[["a", "b"]].drop_duplicates().shape[0], len(reg), 100 * reg["sigma"].median())
