"""Label-only reference predictors on the frozen nested cold-target bank.

These are not models. They bound how much of the frozen-bank MSE is reachable
without any learned representation, so that a learned gain can be attributed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK


def component_target_mean(rows: list[dict], field: str) -> float:
    target_values: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        if not np.isfinite(row[field]):
            continue
        target_values.setdefault((row["component"], row["target"]), []).append(
            row[field])
    component_values: dict[str, list[float]] = {}
    for (component, _), values in target_values.items():
        component_values.setdefault(component, []).append(float(np.mean(values)))
    if not component_values:
        return float("nan")
    return float(np.mean([np.mean(v) for v in component_values.values()]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=20)
    parser.add_argument("--targets-per-component", type=int, default=1)
    parser.add_argument("--draws", type=int, default=1)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    train_values = np.asarray(
        [cell["pK"] for cell in data.cells if cell["split"] == "meta_train"])
    global_mean = float(train_values.mean())
    ligand_sum: dict[str, list[float]] = {}
    for cell in data.cells:
        if cell["split"] == "meta_train":
            ligand_sum.setdefault(cell["ligand_id"], []).append(float(cell["pK"]))
    ligand_prior = {k: float(np.mean(v)) for k, v in ligand_sum.items()}

    banks = data.fixed_nested_episode_banks(
        "meta_test", (0, 1, 2, 3, 5), args.query_size, args.draws,
        args.evaluation_seed, args.targets_per_component)
    report = {}
    coverage = []
    for k, specs in banks.items():
        rows = []
        for spec in specs:
            truth = np.asarray([data.cells[i]["pK"] for i in spec.query])
            ligand_ids = [data.cells[i]["ligand_id"] for i in spec.query]
            known = np.asarray([ligand_prior.get(x, global_mean) for x in ligand_ids])
            seen = np.asarray([x in ligand_prior for x in ligand_ids], dtype=float)
            coverage.append(float(seen.mean()))
            support = (np.asarray([data.cells[i]["pK"] for i in spec.support])
                       if k else None)
            entry = {
                "component": spec.component, "target": spec.target,
                "global_mean": float(((truth - global_mean) ** 2).mean()),
                "ligand_prior": float(((truth - known) ** 2).mean()),
                "oracle_target_mean": float(((truth - truth.mean()) ** 2).mean()),
            }
            if k:
                offset = float((support - np.asarray(
                    [ligand_prior.get(data.cells[i]["ligand_id"], global_mean)
                     for i in spec.support])).mean())
                entry["support_mean"] = float(((truth - support.mean()) ** 2).mean())
                entry["ligand_prior_plus_support_offset"] = float(
                    ((truth - (known + offset)) ** 2).mean())
            else:
                entry["support_mean"] = float("nan")
                entry["ligand_prior_plus_support_offset"] = float("nan")
            rows.append(entry)
        report[str(k)] = {
            field: component_target_mean(rows, field)
            for field in ("global_mean", "ligand_prior", "oracle_target_mean",
                          "support_mean", "ligand_prior_plus_support_offset")
        }
        report[str(k)]["episodes"] = len(rows)
    payload = {
        "schema": "MetaSieve.QPSMPReferenceBaselines.v1",
        "evaluation_seed": args.evaluation_seed,
        "targets_per_component": args.targets_per_component,
        "draws_per_target": args.draws,
        "query_size": args.query_size,
        "meta_train_global_mean_pk": global_mean,
        "meta_test_query_ligand_seen_in_meta_train_fraction": float(
            np.mean(coverage)),
        "reference_mse_pk": report,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
