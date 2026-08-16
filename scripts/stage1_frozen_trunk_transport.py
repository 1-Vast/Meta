"""Stage 1: swap only the support transport on a frozen incumbent trunk.

Stage 6 left one question unresolved: the similarity model beat the incumbent on
`meta_val` and lost on `meta_test`, but those were *independently trained*
trunks, so trunk co-adaptation and mechanism were confounded.

This removes the confound. Each incumbent (`grammar`) checkpoint is frozen and
used as the exact zero-shot trunk, its own shrinkage contract is reused, and
only the support transport is replaced at inference:

* `mean`      - the support-residual mean (the incumbent's own level baseline)
* `tanimoto`  - fixed beta=8 Morgan/Tanimoto residual weighting
* `nearest`   - hard nearest-support residual
* `incumbent` - the checkpoint's own learned transport, for reference

Identical episodes, identical `f0`, identical shrinkage. Any difference is the
transport and nothing else.
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

from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, forward, normalized_episode, resolve_architecture,
    training_label_scale,
)


SUPPORT_SIZES = (1, 2, 3, 5)
TRANSPORTS = ("mean", "tanimoto", "nearest", "incumbent")


def build(checkpoint: Path, data: QPSMPData, device: str | None):
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
    return model, config


def endpoint_for(model, episode, device) -> np.ndarray:
    """`f0` for support ligands followed by query ligands, from the frozen trunk."""
    dtype = next(model.parameters()).dtype
    raw = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)

    def cast(value: torch.Tensor) -> torch.Tensor:
        return value.unsqueeze(0).to(device, dtype)

    return model.encode(
        cast(episode.protein_pooled), cast(episode.protein_tokens),
        cast(episode.protein_mask), cast(raw), cast(bonds), cast(mask),
        cast(episode.protein_chemistry))[0][0].float().cpu().numpy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", action="append", required=True,
                        help="name=path/to/checkpoint.pt (repeatable)")
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--beta", type=float, default=8.0)
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    scale = training_label_scale(data)
    rows: list[dict] = []

    for item in args.arm:
        name, _, path = item.partition("=")
        model, config = build(Path(path), data, args.device)
        transport = getattr(model, "transport", None)
        banks = data.fixed_nested_episode_banks(
            args.split, (0, *SUPPORT_SIZES), config.query_size,
            config.test_draws_per_target, args.evaluation_seed, None)
        with torch.no_grad():
            for k in SUPPORT_SIZES:
                for spec in banks[k]:
                    episode = compact_episode(normalized_episode(
                        data.materialize(spec), scale))
                    endpoint = endpoint_for(model, episode, config.device)
                    support_zero, query_zero = endpoint[:k], endpoint[k:]
                    support_y = episode.support_y.numpy()
                    truth = episode.query_y.numpy()
                    residual = support_y - support_zero
                    support_fp = episode.support_fingerprint.numpy()
                    query_fp = episode.query_fingerprint.numpy()
                    inter = query_fp @ support_fp.T
                    union = (query_fp.sum(-1)[:, None]
                             + support_fp.sum(-1)[None, :] - inter)
                    similarity = inter / np.maximum(union, 1e-9)
                    shrink = float(transport.shrinkage(k, torch.zeros(1)))
                    logits = args.beta * similarity
                    weight = np.exp(logits - logits.max(-1, keepdims=True))
                    weight = weight / weight.sum(-1, keepdims=True)
                    predictions = {
                        "mean": query_zero + shrink * residual.mean(),
                        "tanimoto": query_zero + shrink * (
                            weight * residual[None, :]).sum(-1),
                        "nearest": query_zero + shrink * residual[
                            similarity.argmax(-1)],
                        "incumbent": forward(
                            model, episode).prediction.detach().cpu().float().numpy(),
                    }
                    for transport_name, prediction in predictions.items():
                        ci, comparable = concordance_index(prediction, truth)
                        rows.append({
                            "arm_name": name, "component": spec.component,
                            "target": spec.target, "k": k,
                            "arm": transport_name,
                            "mse_pk": float((((prediction - truth)
                                              * scale.scale) ** 2).mean()),
                            "ci": ci if comparable else None,
                            "spearman": spearman(prediction, truth),
                        })
        del model
        if str(config.device).startswith("cuda"):
            torch.cuda.empty_cache()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows_path = args.output.with_suffix(".rows.jsonl")
    with rows_path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    def component_mean(selected: list[dict], metric: str) -> float:
        by_target: dict[tuple[str, str], list[float]] = {}
        for row in selected:
            value = row.get(metric)
            if value is None or not np.isfinite(value):
                continue
            by_target.setdefault((row["component"], row["target"]), []).append(value)
        by_component: dict[str, list[float]] = {}
        for (component, _), values in by_target.items():
            by_component.setdefault(component, []).append(float(np.mean(values)))
        return (float(np.mean([np.mean(v) for v in by_component.values()]))
                if by_component else float("nan"))

    summary = {}
    for k in SUPPORT_SIZES:
        summary[str(k)] = {}
        for transport_name in TRANSPORTS:
            selected = [r for r in rows if r["k"] == k and r["arm"] == transport_name]
            summary[str(k)][transport_name] = {
                metric: component_mean(selected, metric)
                for metric in ("mse_pk", "ci", "spearman")}
    payload = {"schema": "MetaSieve.Stage1FrozenTrunkTransport.v1",
               "split": args.split, "beta": args.beta,
               "arms": [item.partition("=")[0] for item in args.arm],
               "summary": summary}
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"split={args.split}")
    for k in SUPPORT_SIZES:
        line = "  k=%d " % k
        for transport_name in TRANSPORTS:
            entry = summary[str(k)][transport_name]
            line += "| %s mse=%.4f ci=%.4f rho=%.4f " % (
                transport_name, entry["mse_pk"], entry["ci"], entry["spearman"])
        print(line)


if __name__ == "__main__":
    main()
