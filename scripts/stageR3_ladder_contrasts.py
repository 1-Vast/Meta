"""Stage R4 gates T1/T2: paired component bootstraps along the ablation ladder.

Reads the rows already produced by `stageR3_compare_arms.py` so the arms are the
exact ones scored on the exact same bank. Seeds are averaged inside a target
before components are resampled, so every interval is conditional on the trained
seeds.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.stageR0_retrieval_falsification import component_bootstrap

SUPPORT_SIZES = (0, 1, 2, 3, 5)
FIELDS = ("full_mse_pk", "full_mse_pk_lt40", "full_ci", "full_spearman",
          "zero_shot_mse_pk")


def per_target(rows, arm: str, field: str, k: int) -> dict[tuple[str, str], float]:
    collected: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row["arm"] == arm and row["k"] == k and row.get(field) is not None:
            collected[(row["component"], row["target"])].append(float(row[field]))
    return {key: float(np.mean(v)) for key, v in collected.items()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--pair", action="append", required=True,
                        help="new_vs_old, e.g. B2:B1")
    parser.add_argument("--draws", type=int, default=9999)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = [json.loads(line) for line in
            args.rows.read_text(encoding="utf-8").splitlines() if line.strip()]
    result: dict[str, dict] = {}
    for pair in args.pair:
        new, old = pair.split(":")
        result[pair] = {}
        for k in SUPPORT_SIZES:
            entry = {}
            for field in FIELDS:
                a = per_target(rows, new, field, k)
                b = per_target(rows, old, field, k)
                shared = sorted(set(a) & set(b))
                # old minus new: positive favours `new` for an error metric and
                # disfavours it for a ranking metric. Stated per field below.
                entry[field] = component_bootstrap(
                    [(component, target, b[(component, target)]
                      - a[(component, target)]) for component, target in shared],
                    args.draws, 20260815)
            result[pair][str(k)] = entry

    payload = {"schema": "MetaSieve.StageR3LadderContrast.v1",
               "rows": str(args.rows),
               "sign": "each value is old minus new; positive means the newer "
                       "arm has the lower number, which favours it for MSE and "
                       "disfavours it for CI/Spearman",
               "contrasts": result}
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    for pair, table in result.items():
        print(f"\n{pair}  (old minus new)")
        for k in SUPPORT_SIZES:
            entry = table[str(k)]
            print("  k=%d  MSE %+.4f [%+.4f,%+.4f]   lt40 %+.4f [%+.4f,%+.4f]   "
                  "CI %+.4f [%+.4f,%+.4f]" % (
                      k, entry["full_mse_pk"]["mean"], entry["full_mse_pk"]["lo"],
                      entry["full_mse_pk"]["hi"],
                      entry["full_mse_pk_lt40"]["mean"],
                      entry["full_mse_pk_lt40"]["lo"],
                      entry["full_mse_pk_lt40"]["hi"],
                      entry["full_ci"]["mean"], entry["full_ci"]["lo"],
                      entry["full_ci"]["hi"]))


if __name__ == "__main__":
    main()
