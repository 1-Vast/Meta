"""Assemble a comparison table from retained and new nested-k result files.

Reads only `RESULT.json`-shaped artifacts and prints/writes a single table, so
that narrative numbers cannot drift from the numerical authority.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

KS = ("0", "1", "2", "3", "5")
FIELDS = ("full_mse_pk", "zero_shot_mse_pk", "sar_cut_mse_pk",
          "permuted_mse_pk", "foreign_code_state_mse_pk",
          "wrong_protein_state_mse_pk", "wrong_protein_zero_shot_gap_mse_pk",
          "zero_shot_query_spread_pk")


def read(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return payload["test"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", action="append", required=True,
                        help="name=path/to/RESULT.json (repeatable)")
    parser.add_argument("--field", default="full_mse_pk")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = {}
    for item in args.label:
        name, _, path = item.partition("=")
        rows[name] = read(Path(path))

    lines = ["| configuration | " + " | ".join(f"k={k}" for k in KS) + " |",
             "|---" * (len(KS) + 1) + "|"]
    for name, test in rows.items():
        values = []
        for k in KS:
            value = test.get(k, {}).get(args.field)
            values.append("—" if value is None else f"{value:.3f}")
        lines.append(f"| {name} | " + " | ".join(values) + " |")
    table = "\n".join(lines)
    print(f"field: {args.field}\n{table}")
    payload = {
        "field": args.field,
        "table_markdown": table,
        "values": {name: {k: {f: test.get(k, {}).get(f) for f in FIELDS}
                          for k in KS} for name, test in rows.items()},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                               encoding="utf-8")


if __name__ == "__main__":
    main()
