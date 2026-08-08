"""Apply the gauge-separated component statistic without reading query labels.

Input is JSONL.  Each row must contain ``id``, ``biological_surface`` and
``support_residual``.  Output contains the same identifier, the scalar location
and the resulting prediction vector.  This is an audit/deployment-interface
utility, not evidence that the biological surface has passed an affinity Gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from model.component_statistic import component_prediction


def transform(input_path: Path, output_path: Path, *, ridge: float, bound: float) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with input_path.open("r", encoding="utf-8") as source, output_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as destination:
        for line_number, line in enumerate(source, start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            if "query_label" in row or "query_labels" in row:
                raise ValueError(f"query labels are forbidden (line {line_number})")
            prediction, location = component_prediction(
                row["biological_surface"], row["support_residual"],
                ridge=ridge, bound=bound,
            )
            destination.write(json.dumps({
                "id": row["id"],
                "location": location,
                "prediction": prediction.tolist(),
            }, sort_keys=True) + "\n")
            count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, required=True)
    parser.add_argument("--bound", type=float, default=0.5)
    args = parser.parse_args()
    print(json.dumps({"rows": transform(
        args.input, args.output, ridge=args.ridge, bound=args.bound
    )}, sort_keys=True))


if __name__ == "__main__":
    main()
