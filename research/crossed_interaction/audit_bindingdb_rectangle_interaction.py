"""Label-only 2x2 rectangle magnitude audit for BindingDB CQ panels.

X0/X1 for the Crossed-SAR Interaction Section route. This script does not train
an encoder. It asks whether complete target x ligand rectangles have substantial
observed-label double-difference magnitude after protein-only and ligand-only
main effects cancel algebraically. Without replicate/noise correction this does
not identify latent non-additivity.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_bindingdb_sardelta_cq_bridge import CORPUS
from research.crossed_interaction.train_cq_observable import OUT as CQ_OUT


OUT = CQ_OUT.parent / "bindingdb_rectangle_interaction_x1"


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def rectangle_value(
        y_ta_la: float, y_ta_lb: float,
        y_tb_la: float, y_tb_lb: float) -> float:
    return float(y_ta_lb - y_ta_la - y_tb_lb + y_tb_la)


def panel_rectangles(panel: dict, cell_by_id: dict[str, dict]) -> list[dict]:
    values = {}
    for cell_id in panel["cell_ids"]:
        cell = cell_by_id[cell_id]
        key = (cell["target_id"], cell["ligand_id"])
        if key in values:
            raise ValueError(f"duplicated target-ligand cell in panel {panel['panel_id']}")
        values[key] = cell
    targets = sorted({target for target, _ in values})
    ligands = sorted({ligand for _, ligand in values})
    rows = []
    for target_a, target_b in combinations(targets, 2):
        common = [
            ligand for ligand in ligands
            if (target_a, ligand) in values and (target_b, ligand) in values
        ]
        for ligand_a, ligand_b in combinations(common, 2):
            caa = values[(target_a, ligand_a)]
            cab = values[(target_a, ligand_b)]
            cba = values[(target_b, ligand_a)]
            cbb = values[(target_b, ligand_b)]
            value = rectangle_value(caa["pK"], cab["pK"], cba["pK"], cbb["pK"])
            rows.append({
                "panel_id": panel["panel_id"],
                "split": panel["split"],
                "dependency_component": panel["dependency_component"],
                "target_a": target_a,
                "target_b": target_b,
                "ligand_a": ligand_a,
                "ligand_b": ligand_b,
                "cell_ta_la": caa["cell_id"],
                "cell_ta_lb": cab["cell_id"],
                "cell_tb_la": cba["cell_id"],
                "cell_tb_lb": cbb["cell_id"],
                "rectangle": value,
                "squared_rectangle": float(value * value),
            })
    return rows


def build_rectangles(corpus: Path) -> tuple[list[dict], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    cell_by_id = {cell["cell_id"]: cell for cell in cells}
    rows = []
    panel_counts = []
    for panel in panels:
        panel_rows = panel_rectangles(panel, cell_by_id)
        if panel_rows:
            panel_counts.append({
                "panel_id": panel["panel_id"],
                "split": panel["split"],
                "dependency_component": panel["dependency_component"],
                "rectangles": len(panel_rows),
            })
            rows.extend(panel_rows)
    metadata = {
        "cells": len(cells),
        "panels": len(panels),
        "rectangle_positive_panels": len(panel_counts),
        "rectangles": len(rows),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return rows, metadata


def split_summary(rows: list[dict], split: str) -> dict:
    selected = [row for row in rows if row["split"] == split]
    if not selected:
        return {
            "rectangles": 0,
            "panels": 0,
            "dependency_components": 0,
            "rms_rectangle_pK": 0.0,
            "mean_abs_rectangle_pK": 0.0,
            "mean_square_rectangle": 0.0,
        }
    values = np.asarray([row["rectangle"] for row in selected], dtype=np.float64)
    return {
        "rectangles": len(selected),
        "panels": len({row["panel_id"] for row in selected}),
        "dependency_components": len({row["dependency_component"] for row in selected}),
        "rms_rectangle_pK": float(np.sqrt(np.mean(values * values))),
        "mean_abs_rectangle_pK": float(np.mean(np.abs(values))),
        "mean_square_rectangle": float(np.mean(values * values)),
        "quantiles_abs_rectangle_pK": {
            "q50": float(np.quantile(np.abs(values), 0.50)),
            "q75": float(np.quantile(np.abs(values), 0.75)),
            "q90": float(np.quantile(np.abs(values), 0.90)),
        },
    }


def component_bootstrap(rows: list[dict], split: str, *, draws: int, seed: int) -> dict:
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        if row["split"] == split:
            grouped[row["dependency_component"]].append(row["squared_rectangle"])
    components = sorted(grouped)
    if len(components) < 2:
        raise ValueError("component bootstrap needs at least two components")
    component_mse = np.asarray([
        float(np.mean(grouped[component])) for component in components
    ], dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = component_mse[
        rng.integers(0, len(component_mse), size=(draws, len(component_mse)))
    ].mean(axis=1)
    rms_samples = np.sqrt(samples)
    return {
        "split": split,
        "components": len(components),
        "component_macro_rms": float(np.sqrt(component_mse.mean())),
        "one_sided_95_lcb_rms": float(np.quantile(rms_samples, 0.05)),
        "rms_lcb_gt_0_10": bool(np.quantile(rms_samples, 0.05) > 0.10),
        "rms_lcb_gt_0_20": bool(np.quantile(rms_samples, 0.05) > 0.20),
        "rms_lcb_gt_0_30": bool(np.quantile(rms_samples, 0.05) > 0.30),
    }


def additive_cancellation_check() -> dict:
    max_abs = 0.0
    for target_a, target_b in combinations(range(4), 2):
        for ligand_a, ligand_b in combinations(range(5), 2):
            y_ta_la = 1.7 * target_a - 0.4 * ligand_a + 3.0
            y_ta_lb = 1.7 * target_a - 0.4 * ligand_b + 3.0
            y_tb_la = 1.7 * target_b - 0.4 * ligand_a + 3.0
            y_tb_lb = 1.7 * target_b - 0.4 * ligand_b + 3.0
            max_abs = max(max_abs, abs(rectangle_value(
                y_ta_la, y_ta_lb, y_tb_la, y_tb_lb)))
    return {
        "synthetic_additive_max_abs_rectangle": float(max_abs),
        "pass": bool(max_abs <= 1e-12),
    }


def run(
        corpus: Path = CORPUS, output: Path = OUT,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    rows, metadata = build_rectangles(corpus)
    train_summary = split_summary(rows, "train")
    development_summary = split_summary(rows, "development")
    development_bootstrap = component_bootstrap(
        rows, "development", draws=bootstrap_draws, seed=seed)
    gates = {
        "x0_development_rectangles_ge_1000": development_summary["rectangles"] >= 1000,
        "x0_development_components_ge_5": development_summary["dependency_components"] >= 5,
        "additive_cancellation_exact": additive_cancellation_check()["pass"],
        "x1_development_rms_lcb_gt_0_10": development_bootstrap["rms_lcb_gt_0_10"],
        "x1_development_rms_lcb_gt_0_20": development_bootstrap["rms_lcb_gt_0_20"],
    }
    verdict = (
        "BINDINGDB_RECTANGLE_INTERACTION_X1_PASS_LABEL_HEADROOM"
        if all(gates.values())
        else "BINDINGDB_RECTANGLE_INTERACTION_X1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBRectangleInteractionX1.v1",
        "hypothesis": (
            "Complete 2x2 target x ligand rectangles have label-side observed "
            "double-difference magnitude after protein-only and ligand-only "
            "effects cancel algebraically."),
        "literature_mechanism": {
            "double_mutant_cycle": (
                "2x2 free-energy cycles isolate non-additive coupling terms"),
            "matched_molecular_pairs": (
                "ligand transformations are interpretable SAR units; here their "
                "effect is contrasted across targets"),
        },
        "corpus": metadata,
        "config": {
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "gate_rms_thresholds_pK": [0.10, 0.20, 0.30],
            "labels_used_for_training": False,
        },
        "train_summary": train_summary,
        "development_summary": development_summary,
        "development_bootstrap": development_bootstrap,
        "additive_cancellation_check": additive_cancellation_check(),
        "gates": gates,
        "development_training_authorized": False,
        "g2_label_side_followup_authorized": verdict.endswith("PASS_LABEL_HEADROOM"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    (output / "rectangles_development.jsonl").write_text(
        "".join(
            json.dumps(row, sort_keys=True) + "\n"
            for row in rows if row["split"] == "development"
        ),
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
