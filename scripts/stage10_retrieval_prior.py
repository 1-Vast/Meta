"""Stage C: blend a train-only retrieval prior into the zero-shot endpoint.

Stage A showed the trained endpoint is calibration-dominated and contributes
almost no within-target ranking, while a `meta_train`-only retrieval predictor
beats it by 25.5% at k=0. This tests the smallest deployable consequence:

```text
f0'(q) = (1 - w_q) * f0(q) + w_q * retrieval(q)
r_k    = y_k - f0'(L_k)
f(q)   = f0'(q) + s(n) * sum_k softmax_k(8 * Tanimoto) r_k     (transport unchanged)
```

`w = 0` reproduces the accepted model exactly, which is asserted as a test.
Every retrieval quantity is per-query and uses `meta_train` records only: no
query label, no target identity, no query-set statistic, no training.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
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
    compact_episode, normalized_episode, resolve_architecture, training_label_scale,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)


def tanimoto_rows(query: np.ndarray, bank: np.ndarray) -> np.ndarray:
    inter = query @ bank.T
    union = query.sum(-1)[:, None] + bank.sum(-1)[None, :] - inter
    return inter / np.maximum(union, 1e-9)


def softmax_rows(scores: np.ndarray, beta: float) -> np.ndarray:
    logits = beta * scores
    shifted = np.exp(logits - logits.max(-1, keepdims=True))
    return shifted / shifted.sum(-1, keepdims=True)


def component_target_mean(rows: list[dict], field: str) -> float:
    by_target: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        value = row.get(field)
        if value is None or not np.isfinite(value):
            continue
        by_target.setdefault((row["component"], row["target"]), []).append(float(value))
    by_component: dict[str, list[float]] = {}
    for (component, _), values in by_target.items():
        by_component.setdefault(component, []).append(float(np.mean(values)))
    return (float(np.mean([np.mean(v) for v in by_component.values()]))
            if by_component else float("nan"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", action="append", required=True)
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--retrieval-beta", type=float, default=24.0)
    parser.add_argument("--transport-beta", type=float, default=8.0)
    parser.add_argument("--protein-beta", type=float, default=16.0)
    parser.add_argument("--weights", type=float, nargs="+",
                        default=[0.0, 0.25, 0.5, 0.75, 1.0])
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    scale = training_label_scale(data)
    fingerprints = data.fingerprints

    train_cells = [c for c in data.cells if c["split"] == "meta_train"]
    ligand_values: dict[str, list[float]] = defaultdict(list)
    by_target: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for cell in train_cells:
        ligand_values[cell["ligand_id"]].append(float(cell["pK"]))
        by_target[cell["target_id"]].append((cell["ligand_id"], float(cell["pK"])))
    train_ligands = sorted(ligand_values)
    train_fp = np.stack([fingerprints[k] for k in train_ligands])
    train_mean = np.asarray([float(np.mean(ligand_values[k])) for k in train_ligands])
    train_targets = sorted(by_target)
    pooled = {t: np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
              for t in set(train_targets) | set(data.tasks[args.split])}
    train_pooled = np.stack([pooled[t] for t in train_targets])
    train_pooled /= np.maximum(np.linalg.norm(train_pooled, axis=-1, keepdims=True), 1e-9)

    def retrieval_for(ligand_ids: list[str], target: str) -> dict[str, np.ndarray]:
        """Per-query retrieval predictions in pK. meta_train index only."""
        query_fp = np.stack([fingerprints[k] for k in ligand_ids])
        similarity = tanimoto_rows(query_fp, train_fp)
        ligand = (softmax_rows(similarity, args.retrieval_beta)
                  * train_mean[None, :]).sum(-1)
        vector = pooled[target] / max(float(np.linalg.norm(pooled[target])), 1e-9)
        protein_sim = train_pooled @ vector
        top = np.argsort(-protein_sim)[:16]
        numerator = np.zeros(len(ligand_ids))
        denominator = 1e-9
        for index in top:
            entries = by_target[train_targets[index]]
            bank_fp = np.stack([fingerprints[k] for k, _ in entries])
            bank_y = np.asarray([v for _, v in entries])
            weight = softmax_rows(tanimoto_rows(query_fp, bank_fp),
                                  args.retrieval_beta)
            share = float(np.exp(args.protein_beta *
                                 (protein_sim[index] - protein_sim[top].max())))
            numerator += share * (weight * bank_y[None, :]).sum(-1)
            denominator += share
        dual = numerator / denominator
        return {"ligand": ligand, "dual": dual, "blend": 0.5 * (ligand + dual),
                "novelty": similarity.max(-1)}

    rows: list[dict] = []
    for item in args.arm:
        name, _, path = item.partition("=")
        payload = torch.load(Path(path), map_location="cpu", weights_only=False)
        valid = {f.name for f in fields(TrainConfig)}
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
            adaptive_blocks=config.adaptive_blocks,
            adapter_scale=config.adapter_scale, use_cartesian=config.use_cartesian)
        model.load_state_dict(payload["model_state"])
        model.to(config.device).eval()
        dtype = next(model.parameters()).dtype
        banks = data.fixed_nested_episode_banks(
            args.split, SUPPORT_SIZES, config.query_size,
            config.test_draws_per_target, args.evaluation_seed, None)
        with torch.no_grad():
            for k in SUPPORT_SIZES:
                for spec in banks[k]:
                    episode = compact_episode(normalized_episode(
                        data.materialize(spec), scale))
                    raw = torch.cat((episode.support_atoms, episode.query_atoms), 0)
                    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
                    mask = torch.cat((episode.support_mask, episode.query_mask), 0)

                    def cast(value):
                        return value.unsqueeze(0).to(config.device, dtype)

                    endpoint = model.encode(
                        cast(episode.protein_pooled), cast(episode.protein_tokens),
                        cast(episode.protein_mask), cast(raw), cast(bonds),
                        cast(mask), cast(episode.protein_chemistry),
                    )[0][0].float().cpu().numpy()
                    endpoint = endpoint * scale.scale + scale.mean
                    ligand_ids = [data.cells[i]["ligand_id"]
                                  for i in (*spec.support, *spec.query)]
                    prior = retrieval_for(ligand_ids, spec.target)
                    support_y = (episode.support_y.numpy() * scale.scale
                                 + scale.mean)
                    truth = episode.query_y.numpy() * scale.scale + scale.mean
                    support_fp = episode.support_fingerprint.numpy()
                    query_fp = episode.query_fingerprint.numpy()
                    shrink = float(model.transport.shrinkage(k, torch.zeros(1))) if k else 0.0
                    if k:
                        weight = softmax_rows(
                            tanimoto_rows(query_fp, support_fp), args.transport_beta)
                    for source in ("ligand", "dual", "blend"):
                        for w in args.weights:
                            adjusted = (1.0 - w) * endpoint + w * prior[source]
                            support_zero, query_zero = adjusted[:k], adjusted[k:]
                            prediction = query_zero
                            if k:
                                residual = support_y - support_zero
                                prediction = query_zero + shrink * (
                                    weight * residual[None, :]).sum(-1)
                            ci, comparable = concordance_index(prediction, truth)
                            rows.append({
                                "arm_name": name, "component": spec.component,
                                "target": spec.target, "k": k,
                                "arm": f"{source}_w{w:g}",
                                "mse_pk": float(((prediction - truth) ** 2).mean()),
                                "ci": ci if comparable else None,
                                "spearman": spearman(prediction, truth),
                                "novelty": float(prior["novelty"][k:].mean()),
                            })
        del model
        if str(config.device).startswith("cuda"):
            torch.cuda.empty_cache()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".rows.jsonl").write_text(
        "\n".join(json.dumps(r, sort_keys=True) for r in rows) + "\n",
        encoding="utf-8")
    summary = {}
    for k in SUPPORT_SIZES:
        summary[str(k)] = {}
        for arm in sorted({r["arm"] for r in rows}):
            selected = [r for r in rows if r["k"] == k and r["arm"] == arm]
            summary[str(k)][arm] = {m: component_target_mean(selected, m)
                                    for m in ("mse_pk", "ci", "spearman")}
    args.output.write_text(json.dumps(
        {"schema": "MetaSieve.RetrievalPriorBlend.v1", "split": args.split,
         "retrieval_beta": args.retrieval_beta, "summary": summary},
        indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for k in SUPPORT_SIZES:
        print(f"k={k}")
        for arm in sorted(summary[str(k)]):
            entry = summary[str(k)][arm]
            print("   %-12s mse=%.4f ci=%.4f rho=%.4f" % (
                arm, entry["mse_pk"], entry["ci"], entry["spearman"]))


if __name__ == "__main__":
    main()
