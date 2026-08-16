"""Stage 0 audit diagnostics for the retained QPSMP baseline.

Read-only instrumentation. It adds gradient-coverage, activation-scale,
parameter, memory, and control diagnostics without changing the active model,
the frozen evaluation bank, or any retained result.
"""
from __future__ import annotations

import argparse
from dataclasses import fields, replace
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, LabelScale,
    TrainConfig, compact_episode, forward, normalized_episode,
    resolve_architecture, training_label_scale,
)


ROOT = Path(__file__).resolve().parents[1]


def build_model(payload: dict, data: QPSMPData):
    valid = {field.name for field in fields(TrainConfig)}
    config = TrainConfig(**{k: v for k, v in payload["config"].items() if k in valid})
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
    return model, config


def learned_kernel_scalars(model: torch.nn.Module) -> dict:
    """Read the few-shot kernel scalars of whichever architecture is loaded."""
    softplus = torch.nn.functional.softplus
    kernel = getattr(getattr(model, "meta", None), "term", None)
    if kernel is not None:
        return {
            "temperature": float(softplus(kernel.log_temperature) + 1.0),
            "local_scale": float(torch.sigmoid(kernel.local_scale_logit)),
            "shrinkage_prior_strength": float(softplus(kernel.log_shrinkage)),
        }
    transport = getattr(model, "transport", None)
    if transport is None:
        return {}
    return {
        "temperature": float(softplus(transport.log_temperature) + 1.0),
        "shrinkage_prior_strength": float(softplus(transport.log_shrinkage)),
    }


def parameter_census(model: torch.nn.Module) -> dict:
    groups: dict[str, int] = {}
    for name, parameter in model.named_parameters():
        top = name.split(".")[0]
        groups[top] = groups.get(top, 0) + parameter.numel()
    return {
        "total": int(sum(groups.values())),
        "trainable": int(sum(p.numel() for p in model.parameters() if p.requires_grad)),
        "by_module": {k: int(v) for k, v in sorted(groups.items())},
    }


def gradient_coverage(model: QPSMPBioModel, data: QPSMPData,
                      label_scale: LabelScale, support_sizes: tuple[int, ...],
                      seed: int, episodes: int) -> dict:
    """Per-k finite-gradient coverage of every trainable tensor."""
    report: dict[str, dict] = {}
    rng = np.random.default_rng(seed)
    for k in support_sizes:
        model.zero_grad(set_to_none=True)
        total = 0.0
        for _ in range(episodes):
            spec = data.draw_episode("meta_train", k, 8, rng)
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            output = forward(model, episode, adapt=k > 0)
            truth = episode.query_y.to(device=output.prediction.device,
                                       dtype=output.prediction.dtype)
            loss = (output.prediction - truth).square().mean()
            loss.backward()
            total += float(loss.detach())
        norms, dead, nonfinite = {}, [], []
        for name, parameter in model.named_parameters():
            if parameter.grad is None:
                dead.append(name)
                norms[name] = 0.0
                continue
            value = float(parameter.grad.detach().float().norm())
            norms[name] = value
            if not np.isfinite(value):
                nonfinite.append(name)
            elif value == 0.0:
                dead.append(name)
        report[str(k)] = {
            "loss_mean": total / episodes,
            "global_grad_norm": float(np.sqrt(sum(
                v ** 2 for v in norms.values() if np.isfinite(v)))),
            "dead_parameter_count": len(dead),
            "dead_parameters": sorted(dead),
            "nonfinite_parameters": sorted(nonfinite),
            "grad_norm_by_module": _by_module(norms),
        }
    model.zero_grad(set_to_none=True)
    return report


def _by_module(norms: dict[str, float]) -> dict[str, float]:
    groups: dict[str, float] = {}
    for name, value in norms.items():
        top = name.split(".")[0]
        groups[top] = groups.get(top, 0.0) + (value ** 2 if np.isfinite(value) else 0.0)
    return {k: float(np.sqrt(v)) for k, v in sorted(groups.items())}


def activation_scales(model: QPSMPBioModel, data: QPSMPData,
                      label_scale: LabelScale, seed: int, episodes: int) -> dict:
    values: dict[str, list[float]] = {}
    rng = np.random.default_rng(seed)
    with torch.no_grad():
        for _ in range(episodes):
            spec = data.draw_episode("meta_train", 5, 8, rng)
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            out = forward(model, episode)
            for name, tensor in (
                    ("zero_shot", out.zero_shot),
                    ("cross_zero_shot", out.cross_zero_shot),
                    ("ligand_only", out.ligand_only),
                    ("additive", out.additive),
                    ("level_adjustment", out.level_adjustment),
                    ("sar_adaptation", out.sar_adaptation),
                    ("adaptation", out.adaptation),
                    ("query_primitive", out.query_basis),
                    ("task_state", out.task_state),
                    ("support_residual", out.support_residual_quotient)):
                if tensor.numel():
                    values.setdefault(name, []).append(
                        float(tensor.detach().float().abs().mean()))
    return {k: {"abs_mean": float(np.mean(v)), "abs_mean_std": float(np.std(v))}
            for k, v in sorted(values.items())}


def protein_sensitivity(model: QPSMPBioModel, data: QPSMPData,
                        label_scale: LabelScale, seed: int, episodes: int) -> dict:
    """How much does the zero-shot endpoint move when the protein changes?"""
    rng = np.random.default_rng(seed)
    correct, wrong, spread, wrong_mse, correct_mse = [], [], [], [], []
    ligand_only_spread = []
    with torch.no_grad():
        for _ in range(episodes):
            spec = data.draw_episode("meta_test", 5, 8, rng)
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            out = forward(model, episode, adapt=False)
            pooled, tokens, mask = data.protein_for_target(spec.donor_target)
            swapped = replace(
                episode, protein_pooled=pooled, protein_tokens=tokens,
                protein_mask=mask,
                protein_chemistry=data.protein_chemistry_for_target(spec.donor_target))
            out_wrong = forward(model, swapped, adapt=False)
            delta = (out.prediction - out_wrong.prediction).abs()
            correct.append(float(out.prediction.float().std()))
            wrong.append(float(out_wrong.prediction.float().std()))
            spread.append(float(delta.float().mean()))
            ligand_only_spread.append(float(
                (out.ligand_only - out_wrong.ligand_only).abs().float().mean()))
            correct_mse.append(float(label_scale.squared_error_pk(
                out.prediction, episode.query_y).mean()))
            wrong_mse.append(float(label_scale.squared_error_pk(
                out_wrong.prediction, episode.query_y).mean()))
    return {
        "episodes": episodes,
        "zero_shot_protein_swap_abs_delta_pk": float(
            np.mean(spread) * label_scale.scale),
        "ligand_only_protein_swap_abs_delta_pk": float(
            np.mean(ligand_only_spread) * label_scale.scale),
        "zero_shot_within_episode_std_pk": float(np.mean(correct) * label_scale.scale),
        "zero_shot_mse_pk": float(np.mean(correct_mse)),
        "wrong_protein_zero_shot_mse_pk": float(np.mean(wrong_mse)),
        "wrong_protein_zero_shot_gap_mse_pk": float(
            np.mean(wrong_mse) - np.mean(correct_mse)),
    }


def target_level_decomposition(data: QPSMPData) -> dict:
    """Variance of pK explained by target mean on the meta-test split."""
    report = {}
    for split in ("meta_train", "meta_test"):
        by_target: dict[str, list[float]] = {}
        for cell in data.cells:
            if cell["split"] == split:
                by_target.setdefault(cell["target_id"], []).append(float(cell["pK"]))
        values = np.concatenate([np.asarray(v) for v in by_target.values()])
        within = float(np.mean([np.var(np.asarray(v)) for v in by_target.values()
                                if len(v) > 1]))
        report[split] = {
            "targets": len(by_target),
            "cells": int(values.size),
            "total_variance": float(values.var()),
            "mean_within_target_variance": within,
            "between_target_variance_share": float(
                1.0 - within / values.var()),
        }
    return report


def oracle_ceilings(data: QPSMPData, seed: int, draws: int) -> dict:
    """Label-only reference ceilings on the frozen nested bank shape."""
    banks = data.fixed_nested_episode_banks(
        "meta_test", (0, 1, 2, 3, 5), 20, draws, seed, 1)
    out = {}
    for k, specs in banks.items():
        rows = []
        for spec in specs:
            query = np.asarray([data.cells[i]["pK"] for i in spec.query])
            if k:
                support = np.asarray([data.cells[i]["pK"] for i in spec.support])
                rows.append({
                    "component": spec.component, "target": spec.target,
                    "support_mean_mse": float(((query - support.mean()) ** 2).mean()),
                    "oracle_target_mean_mse": float(((query - query.mean()) ** 2).mean()),
                })
            else:
                rows.append({
                    "component": spec.component, "target": spec.target,
                    "support_mean_mse": float("nan"),
                    "oracle_target_mean_mse": float(((query - query.mean()) ** 2).mean()),
                })
        entry = {}
        for field in ("support_mean_mse", "oracle_target_mean_mse"):
            finite = [r for r in rows if np.isfinite(r[field])]
            if not finite:
                entry[field] = None
                continue
            by_component: dict[str, list[float]] = {}
            for row in finite:
                by_component.setdefault(row["component"], []).append(row[field])
            entry[field] = float(np.mean(
                [np.mean(v) for v in by_component.values()]))
        entry["episodes"] = len(rows)
        out[str(k)] = entry
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--episodes", type=int, default=8)
    parser.add_argument("--seed", type=int, default=907001)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model, config = build_model(payload, data)
    model.to(args.device)
    label_scale = training_label_scale(data)
    if args.device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats()

    started = time.monotonic()
    result = {
        "schema": "MetaSieve.QPSMPStage0Diagnostics.v1",
        "checkpoint": str(args.checkpoint.resolve()),
        "device": args.device,
        "label_scale": {"mean": label_scale.mean, "scale": label_scale.scale},
        "parameters": parameter_census(model),
        "learned_kernel_scalars": learned_kernel_scalars(model),
        "activation_scales": activation_scales(
            model, data, label_scale, args.seed, args.episodes),
        "protein_sensitivity": protein_sensitivity(
            model, data, label_scale, args.seed + 1, args.episodes),
        "gradient_coverage": gradient_coverage(
            model, data, label_scale, (0, 1, 2, 3, 5), args.seed + 2, 2),
        "label_structure": target_level_decomposition(data),
        "oracle_ceilings": oracle_ceilings(data, config.evaluation_seed, 1),
        "wall_seconds": None,
        "peak_cuda_memory_mb": None,
    }
    result["wall_seconds"] = time.monotonic() - started
    result["peak_cuda_memory_mb"] = (
        torch.cuda.max_memory_allocated() / 2 ** 20
        if args.device.startswith("cuda") else 0.0)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "parameters": result["parameters"]["trainable"],
        "kernel": result["learned_kernel_scalars"],
        "protein_sensitivity": result["protein_sensitivity"],
        "dead_by_k": {k: v["dead_parameter_count"]
                      for k, v in result["gradient_coverage"].items()},
    }, indent=2))


if __name__ == "__main__":
    main()
