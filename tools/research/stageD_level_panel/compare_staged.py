"""Stage E paired arm contrast with component bootstrap.

Reads two evaluation row files (evaluate_staged.py output) and reports, for
each k and condition, the mean pairwise difference (arm A minus arm B) per
target with a component-level bootstrap interval (paired per target across
draws). Positive MSE diff = A worse than B.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.stageR0_retrieval_falsification import component_bootstrap

SUPPORT_SIZES = (0, 1, 2, 3, 5)


def load(path: Path):
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arm_a", type=Path)
    parser.add_argument("arm_b", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    a_rows, b_rows = load(args.arm_a), load(args.arm_b)
    a_name = args.arm_a.stem
    b_name = args.arm_b.stem
    out = {"arm_a": a_name, "arm_b": b_name, "by_k": {}}
    for k in SUPPORT_SIZES:
        block = {}
        conditions = sorted({r["condition"] for r in a_rows if r["k"] == k})
        for condition in conditions:
            a_map = {(r["target"], r["draw"]): r for r in a_rows
                     if r["k"] == k and r["condition"] == condition}
            b_map = {(r["target"], r["draw"]): r for r in b_rows
                     if r["k"] == k and r["condition"] == condition}
            pairs = []
            for key, ar in a_map.items():
                br = b_map.get(key)
                if br is None:
                    continue
                for field in ("mse_pk", "rmse_pk", "level_squared",
                              "centered_mse_pk", "spearman", "pearson",
                              "r_squared", "ci", "cliff_sign"):
                    pairs.append((ar["component"], ar["target"], field,
                                  ar[field] - br[field]))
            fields = sorted({p[2] for p in pairs})
            metrics = {}
            for field in fields:
                interval = component_bootstrap(
                    [(c, t, v) for c, t, f, v in pairs if f == field],
                    9999, 20260816)
                metrics[field] = interval
            block[condition] = metrics
        out["by_k"][str(k)] = block
    args.output.write_text(json.dumps(out, indent=1), encoding="utf-8")
    for k in SUPPORT_SIZES:
        block = out["by_k"][str(k)].get("correct", {})
        for field in ("mse_pk", "level_squared", "centered_mse_pk",
                      "spearman", "pearson", "ci"):
            interval = block.get(field)
            if interval:
                print(f"k={k} {field:<16} {interval['mean']:+.4f} "
                      f"[{interval['lo']:+.4f}, {interval['hi']:+.4f}]")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
