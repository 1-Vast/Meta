"""Stage A: leakage-safe, training-free decomposition of cold-target k=0 error.

Splits each target's squared error exactly into a calibration part and a shape
part, and compares protein-blind, ligand-blind, retrieval and model predictors
on the same episodes. Every retrieval index is built from `meta_train` records
only.
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


ORACLE_ESTIMATORS = ("target_oracle", "model_f0_oracle_level")


def kmer_profile(sequence: str, k: int = 3) -> dict[str, float]:
    counts: dict[str, float] = defaultdict(float)
    for index in range(len(sequence) - k + 1):
        counts[sequence[index:index + k]] += 1.0
    norm = np.sqrt(sum(v * v for v in counts.values())) or 1.0
    return {key: value / norm for key, value in counts.items()}


def sparse_cosine(left: dict[str, float], right: dict[str, float]) -> float:
    if len(left) > len(right):
        left, right = right, left
    return float(sum(value * right.get(key, 0.0) for key, value in left.items()))


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


def load_model(checkpoint: Path, data: QPSMPData, device: str | None):
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", action="append", default=[],
                        help="accepted checkpoint(s) supplying model_f0")
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=20)
    parser.add_argument("--ligand-beta", type=float, default=8.0)
    parser.add_argument("--protein-beta", type=float, default=16.0)
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    fingerprints = data.fingerprints

    # ---------- meta_train-only retrieval indices -------------------------
    train_cells = [cell for cell in data.cells if cell["split"] == "meta_train"]
    global_mean = float(np.mean([cell["pK"] for cell in train_cells]))
    ligand_values: dict[str, list[float]] = defaultdict(list)
    target_values: dict[str, list[float]] = defaultdict(list)
    for cell in train_cells:
        ligand_values[cell["ligand_id"]].append(float(cell["pK"]))
        target_values[cell["target_id"]].append(float(cell["pK"]))
    ligand_prior = {k: float(np.mean(v)) for k, v in ligand_values.items()}
    target_prior = {k: float(np.mean(v)) for k, v in target_values.items()}

    train_ligands = sorted(ligand_prior)
    train_ligand_fp = np.stack([fingerprints[key] for key in train_ligands])
    train_ligand_mean = np.asarray([ligand_prior[key] for key in train_ligands])

    train_targets = sorted(target_prior)
    train_target_mean = np.asarray([target_prior[key] for key in train_targets])
    pooled = {}
    for target in set(train_targets) | set(data.tasks[args.split]):
        pooled[target] = np.asarray(
            data.protein_for_target(target)[0], dtype=np.float32)
    train_pooled = np.stack([pooled[t] for t in train_targets])
    train_pooled = train_pooled / np.maximum(
        np.linalg.norm(train_pooled, axis=-1, keepdims=True), 1e-9)
    train_kmer = [kmer_profile(data._protein_sequences[t]) for t in train_targets]

    # dual index: per train target, its own ligand set
    train_cells_by_target: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for cell in train_cells:
        train_cells_by_target[cell["target_id"]].append(
            (cell["ligand_id"], float(cell["pK"])))

    # ---------- models ----------------------------------------------------
    models = []
    for path in args.checkpoint:
        model, config = load_model(Path(path), data, args.device)
        models.append((Path(path).parent.name, model, config))
    scale = training_label_scale(data)

    banks = data.fixed_nested_episode_banks(
        args.split, (0,), args.query_size, 1, args.evaluation_seed, None)
    rows: list[dict] = []
    sensitivity: list[float] = []

    for spec in banks[0]:
        query_ids = [data.cells[i]["ligand_id"] for i in spec.query]
        truth = np.asarray([data.cells[i]["pK"] for i in spec.query])
        query_fp = np.stack([fingerprints[key] for key in query_ids])
        ligand_sim = tanimoto_rows(query_fp, train_ligand_fp)      # [Q,Ltrain]
        ligand_novelty = float(ligand_sim.max(-1).mean())

        target_vector = pooled[spec.target]
        target_vector = target_vector / max(float(np.linalg.norm(target_vector)), 1e-9)
        protein_sim_esm = train_pooled @ target_vector             # [Ttrain]
        query_kmer = kmer_profile(data._protein_sequences[spec.target])
        protein_sim_kmer = np.asarray(
            [sparse_cosine(query_kmer, profile) for profile in train_kmer])

        predictions: dict[str, np.ndarray] = {
            "global_mean": np.full_like(truth, global_mean),
            "ligand_prior": np.asarray(
                [ligand_prior.get(key, global_mean) for key in query_ids]),
            "ligand_neighbor": (softmax_rows(ligand_sim, args.ligand_beta)
                                * train_ligand_mean[None, :]).sum(-1),
            "protein_neighbor_esm": np.full_like(
                truth, float(softmax_rows(protein_sim_esm[None, :],
                                          args.protein_beta)[0] @ train_target_mean)),
            "protein_neighbor_kmer": np.full_like(
                truth, float(softmax_rows(protein_sim_kmer[None, :],
                                          args.protein_beta)[0] @ train_target_mean)),
            "target_oracle": np.full_like(truth, float(truth.mean())),
        }

        # dual retrieval: top protein neighbours, then ligand-weighted within them
        top = np.argsort(-protein_sim_esm)[:16]
        for beta in (8.0, 16.0, 24.0, 32.0):
            numerator = np.zeros_like(truth)
            denominator = np.zeros_like(truth) + 1e-9
            for index in top:
                entries = train_cells_by_target[train_targets[index]]
                if not entries:
                    continue
                bank_fp = np.stack([fingerprints[key] for key, _ in entries])
                bank_y = np.asarray([value for _, value in entries])
                weight = softmax_rows(tanimoto_rows(query_fp, bank_fp), beta)
                protein_weight = float(np.exp(
                    args.protein_beta * (protein_sim_esm[index]
                                         - protein_sim_esm[top].max())))
                numerator += protein_weight * (weight * bank_y[None, :]).sum(-1)
                denominator += protein_weight
            key = "dual_neighbor" if beta == args.ligand_beta else f"dual_neighbor_b{beta:g}"
            predictions[key] = numerator / denominator
        for beta in (2.0, 4.0, 12.0, 16.0, 24.0, 32.0):
            predictions[f"ligand_neighbor_b{beta:g}"] = (
                softmax_rows(ligand_sim, beta) * train_ligand_mean[None, :]).sum(-1)
        # Fully deployable, train-only, no model and no query label: take the
        # target level from sharp ligand retrieval (best calibration) and the
        # within-target shape from dual retrieval (best shape).
        for level_source in ("ligand_neighbor_b16", "ligand_neighbor_b24"):
            for shape_source in ("dual_neighbor_b16", "dual_neighbor_b24"):
                level, shape = predictions[level_source], predictions[shape_source]
                predictions[f"retrieval[{shape_source}|{level_source}]"] = (
                    shape - shape.mean() + level.mean())

        for name, model, config in models:
            episode = compact_episode(normalized_episode(
                data.materialize(spec), scale))
            with torch.no_grad():
                out = model(
                    episode.protein_pooled, episode.protein_tokens,
                    episode.protein_mask, episode.support_atoms,
                    episode.support_bonds, episode.support_mask, episode.support_y,
                    episode.query_atoms, episode.query_bonds, episode.query_mask,
                    adapt=False, protein_chemistry=episode.protein_chemistry,
                    support_fingerprint=episode.support_fingerprint,
                    query_fingerprint=episode.query_fingerprint)
                f0 = out.prediction.detach().cpu().float().numpy()
                f0 = f0 * scale.scale + scale.mean
                donor = data.protein_for_target(spec.donor_target)
                swapped = model(
                    donor[0], donor[1], donor[2], episode.support_atoms,
                    episode.support_bonds, episode.support_mask, episode.support_y,
                    episode.query_atoms, episode.query_bonds, episode.query_mask,
                    adapt=False,
                    protein_chemistry=data.protein_chemistry_for_target(
                        spec.donor_target),
                    support_fingerprint=episode.support_fingerprint,
                    query_fingerprint=episode.query_fingerprint)
                swapped = swapped.prediction.detach().cpu().float().numpy()
                swapped = swapped * scale.scale + scale.mean
            sensitivity.append(float(np.abs(f0 - swapped).mean()))
            predictions[f"model_f0[{name}]"] = f0
            predictions[f"model_f0_oracle_level[{name}]"] = (
                f0 - f0.mean() + truth.mean())
            # Deployable hybrids: keep the model's target level (the best
            # available calibrator) and take the within-target shape from a
            # train-only retriever. Uses no query label.
            for shape_source in ("ligand_neighbor", "dual_neighbor",
                                 "ligand_neighbor_b16", "dual_neighbor_b16"):
                shape = predictions[shape_source]
                predictions[f"hybrid[{shape_source}|{name}]"] = (
                    shape - shape.mean() + f0.mean())

        for estimator, prediction in predictions.items():
            error = prediction - truth
            ci, comparable = concordance_index(prediction, truth)
            rows.append({
                "component": spec.component, "target": spec.target,
                "estimator": estimator,
                "mse_pk": float((error ** 2).mean()),
                "calibration_pk": float(error.mean() ** 2),
                "shape_pk": float(error.var()),
                "ci": ci if comparable else None,
                "spearman": spearman(prediction, truth),
                "ligand_novelty": ligand_novelty,
                "protein_max_sim_esm": float(protein_sim_esm.max()),
                "protein_max_sim_kmer": float(protein_sim_kmer.max()),
                "assay_size": len(data.tasks[args.split][spec.target]),
                "label_range": float(truth.max() - truth.min()),
            })

    estimators = sorted({row["estimator"] for row in rows})
    summary = {}
    for estimator in estimators:
        selected = [r for r in rows if r["estimator"] == estimator]
        summary[estimator] = {
            field: component_target_mean(selected, field)
            for field in ("mse_pk", "calibration_pk", "shape_pk", "ci", "spearman")}
        summary[estimator]["is_oracle"] = any(
            estimator.startswith(prefix) for prefix in ORACLE_ESTIMATORS)

    payload = {
        "schema": "MetaSieve.K0ErrorDecomposition.v1",
        "split": args.split, "episodes": len(banks[0]),
        "retrieval_index": "meta_train records only",
        "meta_train_global_mean_pk": global_mean,
        "model_protein_swap_abs_delta_pk": (float(np.mean(sensitivity))
                                            if sensitivity else None),
        "summary": summary,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".rows.jsonl").write_text(
        "\n".join(json.dumps(row, sort_keys=True) for row in rows) + "\n",
        encoding="utf-8")
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")

    print(f"{args.split}: {len(banks[0])} episodes, index = meta_train only")
    print("%-38s %8s %8s %8s %7s %7s" % (
        "estimator", "MSE", "calib", "shape", "CI", "rho"))
    for estimator in estimators:
        entry = summary[estimator]
        flag = " (oracle)" if entry["is_oracle"] else ""
        print("%-38s %8.4f %8.4f %8.4f %7.4f %7.4f%s" % (
            estimator, entry["mse_pk"], entry["calibration_pk"],
            entry["shape_pk"], entry["ci"], entry["spearman"], flag))
    if sensitivity:
        print(f"\nmodel protein-swap |delta| = "
              f"{float(np.mean(sensitivity)):.4f} pK")


if __name__ == "__main__":
    main()
