"""Uniform arm comparison with MSE, concordance index and Spearman.

`train_qpsmp.evaluate` reports squared error only. The Stage 4 admission failure
was a *ranking* failure, so every later decision needs CI and Spearman computed
on the same episodes, for the same arms, from saved checkpoints.
"""
from __future__ import annotations

import argparse
from dataclasses import fields, replace
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, complete_foreign_prediction, counterfactual_label_assignments,
    forward, normalized_episode, resolve_architecture, training_label_scale,
    wrong_protein_prediction,
)


SUPPORT_SIZES = (0, 1, 2, 3, 5)


def component_target_metric(rows: list[dict], metric: str) -> float:
    by_target: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        value = row[metric]
        if value is None or not np.isfinite(value):
            continue
        by_target.setdefault((row["component"], row["target"]), []).append(float(value))
    by_component: dict[str, list[float]] = {}
    for (component, _), values in by_target.items():
        by_component.setdefault(component, []).append(float(np.mean(values)))
    if not by_component:
        return float("nan")
    return float(np.mean([np.mean(v) for v in by_component.values()]))


def evaluate_arm(checkpoint: Path, data: QPSMPData, split: str,
                 targets_per_component: int, device: str | None) -> dict:
    payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
    valid = {field.name for field in fields(TrainConfig)}
    values = {k: v for k, v in payload["config"].items() if k in valid}
    if device is not None:
        values["device"] = device
    config = TrainConfig(**values)
    model = resolve_architecture(config.arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks, adapter_scale=config.adapter_scale,
        use_cartesian=config.use_cartesian)
    model.load_state_dict(payload["model_state"])
    model.to(config.device).eval()
    scale = training_label_scale(data)
    banks = data.fixed_nested_episode_banks(
        split, SUPPORT_SIZES, config.query_size, config.test_draws_per_target,
        config.evaluation_seed, targets_per_component)
    rows = []
    with torch.no_grad():
        for k, specs in banks.items():
            for spec in specs:
                episode = compact_episode(normalized_episode(
                    data.materialize(spec), scale))
                full = forward(model, episode)
                arms = {
                    "full": full.prediction,
                    "zero_shot": forward(model, episode, adapt=False).prediction,
                    "level_only": full.level_baseline,
                    "sar_cut": full.prediction - full.sar_adaptation,
                }
                if k:
                    wrong_labels = counterfactual_label_assignments(
                        full, episode.support_y)[0]
                    arms["permuted_state"] = forward(model, replace(
                        episode, support_y=wrong_labels)).prediction
                    arms["foreign_code_state"] = complete_foreign_prediction(
                        model, data, episode, spec.donor_target, scale)
                    arms["wrong_protein_state"] = wrong_protein_prediction(
                        model, data, episode, spec.donor_target)
                truth = episode.query_y.detach().cpu().numpy()
                for arm, prediction in arms.items():
                    predicted = prediction.detach().cpu().float().numpy()
                    mse = float(scale.squared_error_pk(
                        prediction, episode.query_y).mean())
                    ci, comparable = concordance_index(predicted, truth)
                    rows.append({
                        "component": spec.component, "target": spec.target,
                        "k": k, "arm": arm, "mse_pk": mse,
                        "ci": ci if comparable else None,
                        "spearman": spearman(predicted, truth),
                    })
    summary = {}
    for k in SUPPORT_SIZES:
        entry = {}
        for arm in sorted({row["arm"] for row in rows if row["k"] == k}):
            selected = [r for r in rows if r["k"] == k and r["arm"] == arm]
            entry[arm] = {
                "mse_pk": component_target_metric(selected, "mse_pk"),
                "ci": component_target_metric(selected, "ci"),
                "spearman": component_target_metric(selected, "spearman"),
            }
        entry["episodes"] = len({(r["component"], r["target"]) for r in rows
                                 if r["k"] == k})
        summary[str(k)] = entry
    return {
        "checkpoint": str(checkpoint.resolve()), "arch": config.arch,
        "difference_loss_weight": getattr(config, "difference_loss_weight", 0.0),
        "split": split, "targets_per_component": targets_per_component,
        "trainable_parameters": int(sum(
            p.numel() for p in model.parameters() if p.requires_grad)),
        "metrics": summary,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", action="append", required=True,
                        help="name=path/to/checkpoint.pt (repeatable)")
    parser.add_argument("--split", default="meta_val",
                        choices=("meta_train", "meta_val", "meta_test"))
    parser.add_argument("--targets-per-component", type=int, default=999999)
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    report = {"schema": "MetaSieve.ArmRankingComparison.v1",
              "split": args.split, "arms": {}}
    for item in args.arm:
        name, _, path = item.partition("=")
        report["arms"][name] = evaluate_arm(
            Path(path), data, args.split, args.targets_per_component, args.device)
        entry = report["arms"][name]["metrics"]
        print(f"== {name} ({report['arms'][name]['arch']}, "
              f"diff={report['arms'][name]['difference_loss_weight']})")
        for k in ("0", "1", "2", "3", "5"):
            full = entry[k]["full"]
            level = entry[k].get("level_only", full)
            print(f"   k={k} mse={full['mse_pk']:.4f} ci={full['ci']:.4f} "
                  f"rho={full['spearman']:.4f} | level mse={level['mse_pk']:.4f} "
                  f"ci={level['ci']:.4f}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows_path = args.output.with_suffix(".rows.jsonl")
    with rows_path.open("w", encoding="utf-8") as handle:
        for name, entry in report["arms"].items():
            for row in entry.pop("rows"):
                handle.write(json.dumps({"arm_name": name, **row},
                                        sort_keys=True) + "\n")
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")


if __name__ == "__main__":
    main()
