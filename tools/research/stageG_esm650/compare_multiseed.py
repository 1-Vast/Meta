"""Multi-seed paired contrast: pool per-seed per-target differences across
seeds, then component bootstrap (9999 draws)."""
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
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arm_a_rows", type=Path, nargs="+")
    parser.add_argument("--arm_b_rows", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = {"by_k": {}}
    for k in SUPPORT_SIZES:
        block = {}
        conditions = sorted({r["condition"] for r in load(args.arm_a_rows[0])
                             if r["k"] == k})
        for condition in conditions:
            diffs = {}  # (component, target) -> list over seeds of mean over draws
            for a_path, b_path in zip(args.arm_a_rows, args.arm_b_rows):
                a_rows = [r for r in load(a_path)
                          if r["k"] == k and r["condition"] == condition]
                b_rows = [r for r in load(b_path)
                          if r["k"] == k and r["condition"] == condition]
                for field in ("mse_pk", "level_squared", "centered_mse_pk",
                              "spearman", "pearson", "ci", "cliff_sign"):
                    a_map, b_map = {}, {}
                    for r in a_rows:
                        a_map.setdefault((r["component"], r["target"]), [])
                        a_map[(r["component"], r["target"])].append(r[field])
                    for r in b_rows:
                        b_map.setdefault((r["component"], r["target"]), [])
                        b_map[(r["component"], r["target"])].append(r[field])
                    for key in a_map:
                        if key not in b_map:
                            continue
                        diff = float(np.mean(a_map[key]) - np.mean(b_map[key]))
                        diffs.setdefault((key[0], key[1], field), []).append(diff)
            metrics = {}
            for (component, target, field), values in diffs.items():
                pass
            by_field = {}
            for (component, target, field), values in diffs.items():
                by_field.setdefault(field, []).append((component, target,
                                                      float(np.mean(values))))
            for field, pairs in by_field.items():
                metrics[field] = component_bootstrap(pairs, 9999, 20260816)
                metrics[field]["seeds"] = int(len(args.arm_a_rows))
            block[condition] = metrics
        out["by_k"][str(k)] = block
    args.output.write_text(json.dumps(out, indent=1), encoding="utf-8")
    for k in SUPPORT_SIZES:
        for field in ("mse_pk", "centered_mse_pk", "spearman", "ci"):
            interval = out["by_k"][str(k)].get("correct", {}).get(field)
            if interval:
                print(f"k={k} {field:<16} {interval['mean']:+.4f} "
                      f"[{interval['lo']:+.4f}, {interval['hi']:+.4f}] "
                      f"seeds={interval.get('seeds')}")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
