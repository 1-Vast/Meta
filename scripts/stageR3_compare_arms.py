"""Stage R3/R4: evaluate every arm on one identical bank and compare them paired.

Handles both architectures behind one interface so the comparison cannot drift:
the incumbent `similarity_only` trunk and the level-shape factorized model are
scored on exactly the same nested episode bank with exactly the same controls and
the same aggregation. Component-level paired bootstraps resample components, the
unit of independent evidence.

Nothing here reads `meta_test`.
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

from model.level_shape import LevelShapeModel
from model.similarity_grammar import tanimoto
from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData
from scripts.stageR0_retrieval_falsification import (
    component_bootstrap, component_target_mean, tanimoto_rows,
)
from scripts.train_level_shape import (
    RouteConfig, episode_tensors, matched_donors, normalized, protein_inputs,
)
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, resolve_architecture, training_label_scale,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)
ARMS = ("full", "zero_shot", "level_only", "permuted", "wrong_protein",
        "ligand_only")


def load_arm(path: Path, data: QPSMPData, device: str):
    payload = torch.load(path, map_location="cpu", weights_only=False)
    config = payload["config"]
    protein_dim = int(data.protein_bank.manifest["hidden_dim"])
    if "routing" in config:                       # level-shape checkpoint
        valid = {f.name for f in fields(RouteConfig)}
        route = RouteConfig(**{k: v for k, v in config.items() if k in valid})
        model = LevelShapeModel(
            protein_dim=protein_dim, hidden_dim=route.hidden_dim,
            task_dim=route.task_dim, ligand_layers=route.ligand_layers,
            pair_dim=route.pair_dim, pair_heads=route.pair_heads,
            anchors=route.anchors)
        kind = "level_shape"
    else:
        valid = {f.name for f in fields(TrainConfig)}
        train = TrainConfig(**{k: v for k, v in config.items() if k in valid})
        model = resolve_architecture(train.arch)(
            protein_dim=protein_dim, hidden_dim=train.hidden_dim,
            task_dim=train.task_dim, ligand_layers=train.ligand_layers,
            pair_dim=train.pair_dim, pair_blocks=train.pair_blocks,
            pair_latents=train.pair_latents, pair_heads=train.pair_heads,
            pair_chunk_size=train.pair_chunk_size,
            support_hidden_dim=train.support_hidden_dim,
            support_blocks=train.support_blocks, adapter_rank=train.adapter_rank,
            adaptive_blocks=train.adaptive_blocks,
            adapter_scale=train.adapter_scale, use_cartesian=train.use_cartesian)
        kind = train.arch
    model.load_state_dict(payload["model_state"])
    return model.to(device).eval(), kind, int(config.get("seed", -1))


def predictions_for(model, kind: str, data: QPSMPData, spec, episode,
                    donor: str, device: str, dtype) -> dict[str, np.ndarray]:
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    correct = protein_inputs(data, spec.target, device, dtype)
    wrong = protein_inputs(data, donor, device, dtype)
    if kind == "level_shape":
        channels = model.encode_ligand(atoms, bonds, mask)
        endpoint, prior, level, _, _ = model.endpoint_with_channels(
            channels, *correct)
        wrong_endpoint = model.endpoint_with_channels(channels, *wrong)[0]
        ligand_only = prior[:, support:] + level[:, support:]
    else:
        endpoint, ligand_value = model.encode(
            *correct[:3], atoms, bonds, mask, correct[3])[:2]
        wrong_endpoint = model.encode(
            *wrong[:3], atoms, bonds, mask, wrong[3])[0]
        ligand_only = ligand_value[:, support:]
    query_endpoint = endpoint[:, support:]
    query_wrong = wrong_endpoint[:, support:]
    transport = torch.zeros_like(query_endpoint)
    permuted = torch.zeros_like(query_endpoint)
    level_only = torch.zeros_like(query_endpoint)
    if support:
        similarity = tanimoto(
            episode.query_fingerprint.to(device, dtype).unsqueeze(0),
            episode.support_fingerprint.to(device, dtype).unsqueeze(0))
        residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                    - endpoint[:, :support])
        shrink = model.transport.shrinkage(support, residual)
        if kind == "level_shape":
            transport = shrink * model.transport(residual, similarity)[0]
            permuted = shrink * model.transport(
                residual.roll(1, dims=-1), similarity)[0]
        else:
            transport = shrink * model.transport(
                None, None, residual, similarity)[0]
            permuted = shrink * model.transport(
                None, None, residual.roll(1, dims=-1), similarity)[0]
        level_only = shrink * residual.mean(-1, keepdim=True)
    return {
        "full": query_endpoint + transport,
        "zero_shot": query_endpoint,
        "level_only": query_endpoint + level_only,
        "permuted": query_endpoint + permuted,
        "wrong_protein": query_wrong + transport,
        "ligand_only": ligand_only,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", action="append", required=True,
                        help="name=path/to/checkpoint.pt")
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=16)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--reference", default="A0")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory)
    scale = training_label_scale(data)
    # Donors come from the evaluation split, not meta_train, so the control
    # varies protein identity alone and not seen-versus-unseen; the whitening
    # transform is fitted on meta_train only (contract 2026-08-16).
    donors = matched_donors(data, args.split, donor_pool=args.split,
                            whitening_pool="meta_train")
    fingerprints = data.fingerprints
    train_fp = np.stack([fingerprints[c["ligand_id"]].numpy()
                         for c in data.cells if c["split"] == "meta_train"])
    banks = data.fixed_nested_episode_banks(
        args.split, SUPPORT_SIZES, args.query_size, 1, args.evaluation_seed, None)

    rows: list[dict] = []
    for item in args.arm:
        name, _, path = item.partition("=")
        model, kind, seed = load_arm(Path(path), data, args.device)
        dtype = next(model.parameters()).dtype
        with torch.no_grad():
            for k, specs in banks.items():
                for spec in specs:
                    episode = compact_episode(
                        normalized(data.materialize(spec), scale))
                    values = predictions_for(model, kind, data, spec, episode,
                                             donors[spec.target], args.device, dtype)
                    truth = (episode.query_y.numpy() * scale.scale + scale.mean)
                    novelty = tanimoto_rows(
                        np.stack([fingerprints[data.cells[i]["ligand_id"]].numpy()
                                  for i in spec.query]), train_fp).max(-1)
                    row = {"arm": name, "kind": kind, "seed": seed, "k": k,
                           "component": spec.component, "target": spec.target,
                           "mean_novelty": float(novelty.mean()),
                           "low_tier_fraction": float((novelty < 0.4).mean())}
                    for label, tensor in values.items():
                        prediction = (tensor.squeeze(0).float().cpu().numpy()
                                      * scale.scale + scale.mean)
                        row[f"{label}_mse_pk"] = float(
                            ((prediction - truth) ** 2).mean())
                        low = novelty < 0.4
                        row[f"{label}_mse_pk_lt40"] = (
                            float(((prediction[low] - truth[low]) ** 2).mean())
                            if low.any() else None)
                        if label in ("full", "zero_shot"):
                            ci, comparable = concordance_index(prediction, truth)
                            row[f"{label}_ci"] = ci if comparable else None
                            row[f"{label}_spearman"] = spearman(prediction, truth)
                        if label == "full":
                            error = prediction - truth
                            row["calibration_pk"] = float(error.mean() ** 2)
                            row["shape_pk"] = float(error.var())
                            row["query_spread_pk"] = float(prediction.std())
                    rows.append(row)
        del model
        if args.device.startswith("cuda"):
            torch.cuda.empty_cache()

    names = [item.partition("=")[0] for item in args.arm]
    fields_of_interest = [f for f in rows[0]
                          if f.endswith(("_mse_pk", "_mse_pk_lt40", "_ci",
                                         "_spearman", "_pk"))]
    summary: dict[str, dict] = {}
    for name in names:
        summary[name] = {}
        for k in SUPPORT_SIZES:
            selected = [r for r in rows if r["arm"] == name and r["k"] == k]
            summary[name][str(k)] = {
                field: component_target_mean(
                    (r["component"], r["target"], r.get(field))
                    for r in selected)
                for field in fields_of_interest}

    def per_target(name: str, field: str, k: int) -> dict[tuple[str, str], float]:
        """Seeds are averaged inside a target before components are resampled,
        so every interval below is conditional on the trained seeds."""
        collected: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows:
            if row["arm"] == name and row["k"] == k and row.get(field) is not None:
                collected[(row["component"], row["target"])].append(
                    float(row[field]))
        return {key: float(np.mean(values)) for key, values in collected.items()}

    def paired(left: str, right: str, field: str, k: int) -> dict:
        """right minus left, so a positive mean favours `left`."""
        a, b = per_target(left, field, k), per_target(right, field, k)
        values = [(component, target, b[(component, target)] - value)
                  for (component, target), value in a.items()
                  if (component, target) in b]
        return component_bootstrap(values, args.bootstrap_draws, 20260815)

    def internal(name: str, control: str, k: int) -> dict:
        """control minus full, so a positive mean means the control is worse."""
        a = per_target(name, "full_mse_pk", k)
        b = per_target(name, f"{control}_mse_pk", k)
        return component_bootstrap(
            [(component, target, b[(component, target)] - value)
             for (component, target), value in a.items()
             if (component, target) in b],
            args.bootstrap_draws, 20260815)

    contrasts: dict[str, dict] = {}
    for name in names:
        if name == args.reference:
            continue
        contrasts[f"{name}_vs_{args.reference}"] = {
            str(k): {field: paired(name, args.reference, field, k)
                     for field in ("full_mse_pk", "full_mse_pk_lt40",
                                   "zero_shot_mse_pk", "full_ci",
                                   "full_spearman")}
            for k in SUPPORT_SIZES}
    for name in names:
        contrasts[f"{name}_internal"] = {
            str(k): {f"{control}_gap": internal(name, control, k)
                     for control in ("wrong_protein", "permuted",
                                     "ligand_only", "level_only")}
            for k in SUPPORT_SIZES}
    for name in names:
        contrasts[f"{name}_seed_direction"] = {
            str(k): {
                field: [float(np.mean([r[field] for r in rows
                                       if r["arm"] == name and r["k"] == k
                                       and r["seed"] == seed
                                       and r.get(field) is not None]))
                        for seed in sorted({r["seed"] for r in rows
                                            if r["arm"] == name})]
                for field in ("full_mse_pk", "full_ci")}
            for k in SUPPORT_SIZES}

    payload = {
        "schema": "MetaSieve.StageR3ArmComparison.v1",
        "split_directory": str(args.split_directory), "split": args.split,
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "population": {
            "targets": len({r["target"] for r in rows}),
            "components": len({r["component"] for r in rows}),
            "episodes_per_k": len(banks[0])},
        "wrong_protein_donor": (
            f"most similar {args.split} target from a different homology "
            "component under meta_train-fitted whitened ESM"),
        "wrong_protein_donor_pool": args.split,
        "whitening_pool": "meta_train",
        "summary": summary, "contrasts": contrasts,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".rows.jsonl").write_text(
        "\n".join(json.dumps(r, sort_keys=True) for r in rows) + "\n",
        encoding="utf-8")
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")

    print("%s: %d targets, %d components" % (
        args.split, payload["population"]["targets"],
        payload["population"]["components"]))
    for k in SUPPORT_SIZES:
        print(f"\nk={k}")
        print("  %-6s %9s %9s %9s %8s %8s %9s" % (
            "arm", "full", "lt40", "wrongP", "CI", "rho", "ligonly"))
        for name in names:
            entry = summary[name][str(k)]
            print("  %-6s %9.4f %9.4f %9.4f %8.4f %8.4f %9.4f" % (
                name, entry["full_mse_pk"], entry["full_mse_pk_lt40"],
                entry["wrong_protein_mse_pk"], entry["full_ci"],
                entry["full_spearman"], entry["ligand_only_mse_pk"]))


if __name__ == "__main__":
    main()
