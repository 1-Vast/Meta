"""Frozen-checkpoint audit of weighting schemes on the *residual* the model transports.

The transferable-signal audit weights raw support affinities, but the model
transports `r_k = y_k - f0(P, L_k)`. That is a different quantity: `f0` may
already absorb part of the chemical similarity structure, in which case a
Tanimoto weighting of residuals gains less than a Tanimoto weighting of labels.

This isolates the weighting scheme from the residual decomposition by holding
one frozen checkpoint's `f0` fixed and varying only how the residuals are
combined, on a fixed `meta_val` episode bank, with the checkpoint's own learned
shrinkage.
"""
from __future__ import annotations

import argparse
from dataclasses import fields
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, normalized_episode, resolve_architecture, training_label_scale,
)


SUPPORT_SIZES = (1, 2, 3, 5)


def component_mean(rows: list[dict], field: str) -> float:
    by_target: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        value = row.get(field)
        if value is None or not np.isfinite(value):
            continue
        by_target.setdefault((row["component"], row["target"]), []).append(value)
    by_component: dict[str, list[float]] = {}
    for (component, _), values in by_target.items():
        by_component.setdefault(component, []).append(float(np.mean(values)))
    if not by_component:
        return float("nan")
    return float(np.mean([np.mean(v) for v in by_component.values()]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--beta", type=float, default=8.0)
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    valid = {field.name for field in fields(TrainConfig)}
    values = {k: v for k, v in payload["config"].items() if k in valid}
    if args.device is not None:
        values["device"] = args.device
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
    transport = getattr(model, "transport", None)

    banks = data.fixed_nested_episode_banks(
        args.split, (0, *SUPPORT_SIZES), config.query_size,
        config.test_draws_per_target, args.evaluation_seed, None)

    report: dict[str, dict] = {}
    with torch.no_grad():
        for k in SUPPORT_SIZES:
            rows = []
            for spec in banks[k]:
                episode = compact_episode(normalized_episode(
                    data.materialize(spec), scale))
                raw = torch.cat((episode.support_atoms, episode.query_atoms), 0)
                bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
                mask = torch.cat((episode.support_mask, episode.query_mask), 0)
                model_dtype = next(model.parameters()).dtype

                def cast(value: torch.Tensor) -> torch.Tensor:
                    # Protein-bank tensors are stored as float16; `encode` is
                    # called directly here so the usual forward-path cast is
                    # not applied for us.
                    return value.unsqueeze(0).to(config.device, model_dtype)

                endpoint = model.encode(
                    cast(episode.protein_pooled), cast(episode.protein_tokens),
                    cast(episode.protein_mask), cast(raw), cast(bonds),
                    cast(mask), cast(episode.protein_chemistry),
                )[0][0].float().cpu().numpy()
                support_zero, query_zero = endpoint[:k], endpoint[k:]
                support_y = episode.support_y.numpy()
                query_y = episode.query_y.numpy()
                residual = support_y - support_zero            # [K]
                support_fp = episode.support_fingerprint.numpy()
                query_fp = episode.query_fingerprint.numpy()
                inter = query_fp @ support_fp.T
                union = (query_fp.sum(-1)[:, None] + support_fp.sum(-1)[None, :]
                         - inter)
                similarity = inter / np.maximum(union, 1e-9)   # [Q,K]
                shrink = (float(transport.shrinkage(k, torch.zeros(1)))
                          if transport is not None and hasattr(transport, "shrinkage")
                          else k / (k + 2.0))
                logits = args.beta * similarity
                weight = np.exp(logits - logits.max(-1, keepdims=True))
                weight = weight / weight.sum(-1, keepdims=True)
                nearest = similarity.argmax(-1)

                def mse(prediction: np.ndarray) -> float:
                    return float((((prediction - query_y) * scale.scale) ** 2).mean())

                per_support = query_zero[:, None] + shrink * residual[None, :]
                rows.append({
                    "component": spec.component, "target": spec.target,
                    "zero_shot": mse(query_zero),
                    "mean_residual": mse(query_zero + shrink * residual.mean()),
                    "tanimoto_residual": mse(
                        query_zero + shrink * (weight * residual[None, :]).sum(-1)),
                    "nearest_residual": mse(
                        query_zero + shrink * residual[nearest]),
                    "oracle_residual": float(
                        (((per_support - query_y[:, None]) * scale.scale) ** 2
                         ).min(-1).mean()),
                })
            report[str(k)] = {
                field: component_mean(rows, field)
                for field in ("zero_shot", "mean_residual", "tanimoto_residual",
                              "nearest_residual", "oracle_residual")
            }
            report[str(k)]["episodes"] = len(rows)

    out = {
        "schema": "MetaSieve.ResidualTransportAudit.v1",
        "checkpoint": str(args.checkpoint.resolve()), "arch": config.arch,
        "split": args.split, "evaluation_seed": args.evaluation_seed,
        "fingerprint": "production 1024-bit Morgan r=2 from QPSMPData",
        "softmax_beta": args.beta,
        "note": ("oracle_residual picks the best single support per query using "
                 "query labels; diagnostic upper bound only"),
        "estimator_mse_pk": report,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(out["estimator_mse_pk"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
