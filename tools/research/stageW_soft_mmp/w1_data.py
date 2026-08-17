"""Stage W W1 — KIBA data construction and frozen split admission.

Builds soft-family target effects and D rows for KIBA, applies the frozen
component split, and checks the W1 split admission gate before any training.
Run:
    python -m tools.research.stageW_soft_mmp.w1_data
"""
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

from tools.research.stageW_soft_mmp.w0_census import (
    KIBA, build_observations, cdhit_components, read_dataset,
)

HERE = Path(__file__).resolve().parent
PREREG_SHA = "038f4d97f74841023c48a2e9b3bab5592a0bad2bb9fa54a464d5290641549082"
SEED = 20260821
HELDOUT_COMPONENTS = 24
FIT_UNSAMPLED_FRACTION = 0.10
ADMISSION = {"heldout_repeated_rows": 500, "heldout_repeated_components": 10,
             "repeated_families": 50}


def build():
    rows, targets, ligands = read_dataset(KIBA)
    components = cdhit_components(targets, "kiba")
    observations, no_cut = build_observations(rows, targets, ligands, components)

    effects: dict[tuple[str, str], dict] = defaultdict(list)
    for item in observations:
        effects[(item["family_key"], item["target"])].append(item)
    target_effects = []
    for (family, target), items in sorted(effects.items()):
        target_effects.append({
            "family": family, "target": target,
            "component": items[0]["component"],
            "delta_y": float(np.median([x["delta_y"] for x in items])),
            "observations": len(items),
        })

    by_family: dict[str, list[dict]] = defaultdict(list)
    for effect in target_effects:
        by_family[effect["family"]].append(effect)

    d_rows = []
    for family, family_effects in sorted(by_family.items()):
        ordered = sorted(family_effects, key=lambda x: x["target"])
        for i, left in enumerate(ordered):
            for right in ordered[i + 1:]:
                d_rows.append({
                    "family": family,
                    "target_left": left["target"],
                    "target_right": right["target"],
                    "component_left": left["component"],
                    "component_right": right["component"],
                    "value": left["delta_y"] - right["delta_y"],
                    "cross_component": left["component"] != right["component"],
                })

    all_components = sorted(set(components.values()))
    rng = np.random.default_rng(SEED)
    order = rng.permutation(len(all_components))
    heldout = tuple(sorted(all_components[int(i)]
                           for i in order[:HELDOUT_COMPONENTS]))
    fit = tuple(sorted(c for c in all_components if c not in set(heldout)))
    fit_set, heldout_set = set(fit), set(heldout)

    fit_rows = [r for r in d_rows
                if r["component_left"] in fit_set
                and r["component_right"] in fit_set]
    fit_cross = [r for r in fit_rows if r["cross_component"]]
    heldout_rows = [r for r in d_rows
                    if r["component_left"] in heldout_set
                    and r["component_right"] in heldout_set]
    fit_families = {e["family"] for e in target_effects
                    if e["component"] in fit_set}
    heldout_repeated = [r for r in heldout_rows if r["family"] in fit_families]
    heldout_cold = [r for r in heldout_rows if r["family"] not in fit_families]

    order = rng.permutation(len(fit_cross))
    holdout_idx = set(order[:max(1, int(len(fit_cross) * FIT_UNSAMPLED_FRACTION))])
    train = [r for i, r in enumerate(fit_cross) if i not in holdout_idx]
    fit_unsampled = [r for i, r in enumerate(fit_cross) if i in holdout_idx]

    admission = {
        "heldout_repeated_rows": len(heldout_repeated),
        "heldout_repeated_components": len({
            c for r in heldout_repeated
            for c in (r["component_left"], r["component_right"])}),
        "repeated_families": len({r["family"] for r in heldout_repeated}),
    }
    gate = {name: {"measured": value, "threshold": ADMISSION[name],
                   "pass": value >= ADMISSION[name]}
            for name, value in admission.items()}

    payload = {
        "schema": "MetaSieve.StageW.W1Data.v1",
        "stage": "stageW_soft_mmp",
        "w1_preregistration_sha256": PREREG_SHA,
        "seed": SEED,
        "split": {
            "fit_components": len(fit),
            "heldout_components": len(heldout),
            "fit": list(fit),
            "heldout": list(heldout),
        },
        "counts": {
            "observations": len(observations),
            "target_effects": len(target_effects),
            "d_rows": len(d_rows),
            "fit_cross_component_D": len(fit_cross),
            "train_rows": len(train),
            "fit_unsampled_rows": len(fit_unsampled),
            "heldout_all_rows": len(heldout_rows),
            "heldout_repeated_rows": len(heldout_repeated),
            "heldout_repeated_components": admission[
                "heldout_repeated_components"],
            "repeated_families": admission["repeated_families"],
            "heldout_cold_rows": len(heldout_cold),
            "heldout_cold_components": len({
                c for r in heldout_cold
                for c in (r["component_left"], r["component_right"])}),
        },
        "admission_gate": gate,
        "all_pass": all(v["pass"] for v in gate.values()),
    }
    return payload, observations, target_effects, train, fit_unsampled, \
        heldout_repeated, heldout_cold


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "W1_DATA.json")
    parser.add_argument("--observations", type=Path,
                        default=HERE / "W1_OBSERVATIONS.jsonl.gz")
    args = parser.parse_args()

    payload, observations, target_effects, train, fit_unsampled, \
        repeated, cold = build()

    with gzip.open(args.observations, "wt", encoding="utf-8") as handle:
        for item in observations:
            handle.write(json.dumps(item, sort_keys=True) + "\n")

    args.output.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "admission": payload["admission_gate"],
        "all_pass": payload["all_pass"],
        "counts": payload["counts"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
