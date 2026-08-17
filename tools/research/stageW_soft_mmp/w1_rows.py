"""Materialize W1 D-row banks for the KIBA local-interaction training."""
from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

HERE = Path(__file__).resolve().parent


def load_observations():
    rows = []
    with gzip.open(HERE / "W1_OBSERVATIONS.jsonl.gz", "rt", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def build_rows():
    data = json.loads((HERE / "W1_DATA.json").read_text(encoding="utf-8"))
    fit = set(data["split"]["fit"])
    heldout = set(data["split"]["heldout"])
    observations = load_observations()

    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for item in observations:
        groups[(item["family_key"], item["target"])].append(item)
    effects = []
    for (family, target), items in sorted(groups.items()):
        values = [x["delta_y"] for x in items]
        effects.append({
            "family": family, "target": target,
            "component": items[0]["component"],
            "delta_y": float(np.median(values)),
            "murcko_core": items[0]["murcko_core"],
            "category_a": items[0]["category_a"],
            "category_b": items[0]["category_b"],
        })
    by_family: dict[str, list[dict]] = defaultdict(list)
    for effect in effects:
        by_family[effect["family"]].append(effect)

    fit_rows, heldout_rows, heldout_repeated = [], [], []
    fit_families = {e["family"] for e in effects if e["component"] in fit}
    for family, fam_effects in sorted(by_family.items()):
        ordered = sorted(fam_effects, key=lambda x: x["target"])
        for i, left in enumerate(ordered):
            for right in ordered[i + 1:]:
                row = {
                    "family": family,
                    "target_left": left["target"],
                    "target_right": right["target"],
                    "component_left": left["component"],
                    "component_right": right["component"],
                    "value": left["delta_y"] - right["delta_y"],
                    "murcko_core": left["murcko_core"],
                    "category_a": left["category_a"],
                    "category_b": left["category_b"],
                }
                if left["component"] in fit and right["component"] in fit:
                    if left["component"] != right["component"]:
                        fit_rows.append(row)
                elif left["component"] in heldout and right["component"] in heldout:
                    heldout_rows.append(row)
                    if family in fit_families:
                        heldout_repeated.append(row)

    rng = np.random.default_rng(data["seed"])
    order = rng.permutation(len(fit_rows))
    holdout_index = set(order[:max(1, int(len(fit_rows) * 0.10))])
    train = [r for i, r in enumerate(fit_rows) if i not in holdout_index]
    fit_unsampled = [r for i, r in enumerate(fit_rows) if i in holdout_index]

    return {
        "train": train, "fit_unsampled": fit_unsampled,
        "heldout_repeated": heldout_repeated, "heldout_all": heldout_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "W1_ROWS.jsonl.gz")
    args = parser.parse_args()
    banks = build_rows()
    with gzip.open(args.output, "wt", encoding="utf-8") as handle:
        for bank, rows in banks.items():
            for row in rows:
                row = dict(row)
                row["bank"] = bank
                handle.write(json.dumps(row, sort_keys=True) + "\n")
    print(json.dumps({bank: len(rows) for bank, rows in banks.items()}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
