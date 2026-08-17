"""Evaluate one QPSMP checkpoint on a fixed nested-k episode bank."""
from __future__ import annotations

import argparse
from dataclasses import fields
import json
from pathlib import Path
import sys

import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, evaluate, normalized_episode, resolve_architecture,
    training_label_scale,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--targets-per-component", type=int, default=None,
                        help="override the checkpoint's bank width; the frozen "
                             "protocol bank uses the checkpoint value")
    parser.add_argument("--draws-per-target", type=int, default=None)
    parser.add_argument("--device", default=None,
                        help="override the checkpoint's device, e.g. cpu")
    parser.add_argument("--split", default="meta_val",
                        choices=("meta_train", "meta_val", "meta_test"),
                        help="evaluate on this split; meta_test is sealed and "
                             "requires --include-meta-test (contract 2026-08-16)")
    parser.add_argument("--split-directory", type=Path, default=None,
                        help="governed split directory (double-cold protocol)")
    parser.add_argument("--include-meta-test", action="store_true",
                        help="physically unseal meta_test cells in QPSMPData")
    parser.add_argument("--meta-test-authorization", default=None,
                        help="written reason recorded in the artifact; "
                             "mandatory with --include-meta-test "
                             "(contract 2026-08-16)")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.split == "meta_test" and not args.include_meta_test:
        parser.error("evaluating the sealed meta_test split requires "
                     "--include-meta-test")
    if args.include_meta_test and not (args.meta_test_authorization or "").strip():
        parser.error("--include-meta-test requires a written "
                     "--meta-test-authorization (contract 2026-08-16)")

    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    valid = {field.name for field in fields(TrainConfig)}
    values = {key: value for key, value in payload["config"].items() if key in valid}
    values["evaluation_seed"] = args.evaluation_seed
    if args.device is not None:
        values["device"] = args.device
    config = TrainConfig(**values)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=args.include_meta_test,
                     meta_test_authorization=args.meta_test_authorization)
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
    model.to(config.device)
    scale = training_label_scale(data)
    targets_per_component = (config.eval_targets_per_component
                             if args.targets_per_component is None
                             else args.targets_per_component)
    draws_per_target = (config.test_draws_per_target
                        if args.draws_per_target is None else args.draws_per_target)
    specs = data.fixed_nested_episode_banks(
        args.split, (0, 1, 2, 3, 5), config.query_size,
        draws_per_target, args.evaluation_seed, targets_per_component)
    metrics = {
        str(k): evaluate(
            model, data,
            tuple(compact_episode(normalized_episode(
                data.materialize(spec), scale)) for spec in bank),
            controls=k > 0, label_scale=scale)
        for k, bank in specs.items()
    }
    result = {"schema": "MetaSieve.PairedNestedEvaluation.v1",
              "evaluation_seed": args.evaluation_seed,
              "arch": config.arch,
              "split": args.split,
              "targets_per_component": targets_per_component,
              "draws_per_target": draws_per_target,
              "checkpoint": str(args.checkpoint.resolve()), "test": metrics}
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
