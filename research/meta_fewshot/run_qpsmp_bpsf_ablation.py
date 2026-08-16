"""GPU development ablation for research-only BPSF variants."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys

import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.meta_fewshot.qpsmp_bpsf_v2 import research_model_factory
from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    evaluate, train,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", choices=("baseline", "relevance", "shared"),
                        required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260871)
    parser.add_argument("--stage-a-steps", type=int, default=60)
    parser.add_argument("--stage-b-steps", type=int, default=60)
    parser.add_argument("--episodes-per-step", type=int, default=4)
    parser.add_argument("--train-cache-size", type=int, default=64)
    parser.add_argument("--query-size", type=int, default=8)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.mkdir(parents=True)
    stage_a = TrainConfig(
        seed=args.seed, steps=args.stage_a_steps,
        episodes_per_step=args.episodes_per_step,
        train_cache_size=args.train_cache_size, query_size=args.query_size,
        val_interval=max(10, args.stage_a_steps // 5),
        eval_targets_per_component=2, device=args.device,
        zero_support_only=True, zero_shot_loss_weight=0.0)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    factory = research_model_factory(
        shared_latent=args.variant == "shared",
        relevance_weighted=args.variant == "relevance")
    model, stage_a_diagnostics, scale = train(
        data, stage_a, support_sizes=(1, 2, 3, 5),
        progress_path=args.output / "stage_a_progress.jsonl", model_factory=factory)
    parent = args.output / "stage_a_checkpoint.pt"
    torch.save({"model_state": model.state_dict(), "config": asdict(stage_a)}, parent)
    stage_b = replace(
        stage_a, steps=args.stage_b_steps,
        val_interval=max(10, args.stage_b_steps // 5),
        zero_support_only=False, section_only=True,
        pretrained_checkpoint=str(parent.resolve()))
    model, stage_b_diagnostics, scale = train(
        data, stage_b, support_sizes=(1, 2, 3, 5),
        progress_path=args.output / "stage_b_progress.jsonl", model_factory=factory)
    parent.unlink()
    bank = data.fixed_nested_episode_banks(
        "meta_test", (1, 2, 3, 5), args.query_size, 1, args.seed, 2)
    metrics = {str(k): evaluate(model, data, episodes, True, scale)
               for k, episodes in bank.items()}
    result = {
        "schema": "MetaSieve.QPSMPBPSFV2Research.v1",
        "scope": "consumed-development-research-only",
        "variant": args.variant,
        "stage_a_config": asdict(stage_a), "stage_b_config": asdict(stage_b),
        "training": {"stage_a": stage_a_diagnostics,
                     "stage_b": stage_b_diagnostics},
        "metrics": metrics,
        "promotion_authorized": False,
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
