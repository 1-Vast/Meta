"""Stage R6 (Stage 2): compare every arm on one identical bank, paired.

Handles three architectures behind one interface so the comparison cannot
drift: the incumbent `similarity_only` trunk, the level-shape factorized
model, and the relative-transport model — all scored on exactly the same
nested episode bank with the same controls and the same aggregation.
Component-level paired bootstraps resample components, the unit of
independent evidence.

Added for Stage R6 (contract 2026-08-16): per-query novelty tiers (max
Morgan Tanimoto to the meta_train ligand block), activity-cliff sign
accuracy, and the `nogate` arm isolating the relative gate. Nothing here
reads `meta_test`.
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
from model.reltransport import RelTransportModel
from model.similarity_grammar import tanimoto
from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData
from scripts.stageR0_retrieval_falsification import (
    component_bootstrap, component_target_mean, tanimoto_rows,
)
from scripts.train_level_shape import (
    RouteConfig, matched_donors, normalized, protein_inputs,
)
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, resolve_architecture, training_label_scale,
)
from scripts.train_reltransport import RelConfig, transport_block
from scripts.qpsmp_data import stable_seed

SUPPORT_SIZES = (0, 1, 2, 3, 5)
ARMS = ("full", "zero_shot", "level_only", "permuted", "wrong_protein",
        "ligand_only", "nogate")


def load_arm(path: Path, data: QPSMPData, device: str):
    payload = torch.load(path, map_location="cpu", weights_only=False)
    config = payload["config"]
    protein_dim = int(data.protein_bank.manifest["hidden_dim"])
    if "identify_weight" in config and "routing" not in config:
        # Grammar-trunk shape-first checkpoint (train_grammar_shape): the
        # exact incumbent architecture with the shape-first training method.
        from model.similarity_grammar import SimilarityGrammarModel
        model = SimilarityGrammarModel(
            protein_dim=protein_dim, hidden_dim=config["hidden_dim"],
            task_dim=config["task_dim"], ligand_layers=config["ligand_layers"],
            pair_dim=config["pair_dim"], pair_latents=config["pair_latents"],
            pair_heads=config["pair_heads"], use_learned_key=False)
        kind = "grammar_shape"
    elif "routing" in config and "gate" in config:      # relative transport
        valid = {f.name for f in fields(RelConfig)}
        rel = RelConfig(**{k: v for k, v in config.items() if k in valid})
        model = RelTransportModel(
            protein_dim=protein_dim, hidden_dim=rel.hidden_dim,
            task_dim=rel.task_dim, ligand_layers=rel.ligand_layers,
            pair_dim=rel.pair_dim, pair_latents=rel.pair_latents,
            pair_heads=rel.pair_heads, anchors=rel.anchors, rank=rel.rank)
        kind = "reltransport"
    elif "routing" in config:                          # level-shape checkpoint
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


def predict_reltransport(model, data, spec, episode, donor, device, dtype):
    support = episode.support_atoms.shape[0]
    atoms = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)
    atoms = atoms.to(device, dtype).unsqueeze(0)
    bonds = bonds.to(device, dtype).unsqueeze(0)
    mask = mask.to(device, dtype).unsqueeze(0)
    correct = protein_inputs(data, spec.target, device, dtype)
    wrong = protein_inputs(data, donor, device, dtype)
    parts = model.forward_parts(
        *correct[:3], atoms, bonds, mask, correct[3])
    donor_parts = model.forward_parts(
        *wrong[:3], atoms, bonds, mask, wrong[3])
    endpoint, prior, level, _, _, embed, u, u_gate = parts
    query_endpoint = endpoint[:, support:]
    wrong_endpoint = donor_parts[0][:, support:]
    ligand_only = prior[:, support:] + level[:, support:]
    transport = torch.zeros_like(query_endpoint)
    permuted = torch.zeros_like(query_endpoint)
    nogate = torch.zeros_like(query_endpoint)
    wrong_ligand = torch.zeros_like(query_endpoint)
    level_only = torch.zeros_like(query_endpoint)
    if support:
        query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
        support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
        similarity = tanimoto(query_fp, support_fp)
        residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                    - endpoint[:, :support])
        shrink = model.transport.shrinkage(support, residual)
        query_embed = embed[:, support:]
        support_embed = embed[:, :support]
        transport, _ = transport_block(
            model, support, query_embed, support_embed, u, u_gate, residual,
            similarity)
        transport = shrink * transport
        nogate, _ = transport_block(
            model, support, query_embed, support_embed, u, u_gate, residual,
            similarity, relative_on=False)
        nogate = shrink * nogate
        permuted_residual = (residual.roll(1, dims=-1) if support > 1
                             else -residual)
        permuted, _ = transport_block(
            model, support, query_embed, support_embed, u, u_gate,
            permuted_residual, similarity)
        permuted = shrink * permuted
        level_only = shrink * residual.mean(-1, keepdim=True)
        # F4 control: replace the support ligand with a cross-component donor
        # ligand, keep the labels, recompute the transport under the correct
        # protein. The prediction must move (the correction depends on the
        # support ligand) — and the arm must be worse than the correct support.
        donor_cells = [int(index) for index in data.tasks[spec.split][spec.donor_target]]
        order = np.random.default_rng(stable_seed(
            "reltransport-donor-ligand", spec.split, spec.target,
            spec.support, spec.query)).permutation(donor_cells)
        donor_indices = tuple(map(int, order[:support]))
        values = [data.ligand_bank.get(data.cells[index]["ligand_id"])
                  for index in donor_indices]
        max_atoms = max(value[0].shape[0] for value in values)
        w_atoms = torch.stack([torch.nn.functional.pad(
            torch.from_numpy(value[0]),
            (0, 0, 0, max_atoms - value[0].shape[0])) for value in values])
        w_bonds = torch.stack([torch.nn.functional.pad(
            torch.from_numpy(value[1]),
            (0, 0, 0, max_atoms - value[1].shape[0],
             0, max_atoms - value[1].shape[0])) for value in values])
        w_mask = torch.stack([torch.nn.functional.pad(
            torch.from_numpy(value[2]),
            (0, max_atoms - value[2].shape[0])) for value in values])
        wrong_parts = model.forward_parts(
            *correct[:3], w_atoms.to(device, dtype).unsqueeze(0),
            w_bonds.to(device, dtype).unsqueeze(0),
            w_mask.to(device, dtype).unsqueeze(0), correct[3])
        wrong_residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                          - wrong_parts[0]).detach()
        wrong_ligand_t, _ = transport_block(
            model, support, query_embed, wrong_parts[5], u, u_gate,
            wrong_residual, similarity)
        wrong_ligand = query_endpoint + shrink * wrong_ligand_t
    return {
        "full": query_endpoint + transport,
        "zero_shot": query_endpoint,
        "level_only": query_endpoint + level_only,
        "permuted": query_endpoint + permuted,
        "wrong_protein": wrong_endpoint + transport,
        "ligand_only": ligand_only,
        "nogate": query_endpoint + nogate,
        "wrong_ligand": wrong_ligand,
    }


def predict_level_shape(model, data, spec, episode, donor, device, dtype):
    support = episode.support_atoms.shape[0]
    atoms = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)
    atoms = atoms.to(device, dtype).unsqueeze(0)
    bonds = bonds.to(device, dtype).unsqueeze(0)
    mask = mask.to(device, dtype).unsqueeze(0)
    correct = protein_inputs(data, spec.target, device, dtype)
    wrong = protein_inputs(data, donor, device, dtype)
    channels = model.encode_ligand(atoms, bonds, mask)
    endpoint, prior, level, _, _ = model.endpoint_with_channels(
        channels, *correct)
    wrong_endpoint = model.endpoint_with_channels(channels, *wrong)[0]
    query_endpoint = endpoint[:, support:]
    query_wrong = wrong_endpoint[:, support:]
    ligand_only = prior[:, support:] + level[:, support:]
    transport = torch.zeros_like(query_endpoint)
    permuted = torch.zeros_like(query_endpoint)
    level_only = torch.zeros_like(query_endpoint)
    if support:
        query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
        support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
        similarity = tanimoto(query_fp, support_fp)
        residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                    - endpoint[:, :support])
        shrink = model.transport.shrinkage(support, residual)
        transport = shrink * model.transport(residual, similarity)[0]
        permuted = shrink * model.transport(
            residual.roll(1, dims=-1), similarity)[0]
        level_only = shrink * residual.mean(-1, keepdim=True)
    return {
        "full": query_endpoint + transport,
        "zero_shot": query_endpoint,
        "level_only": query_endpoint + level_only,
        "permuted": query_endpoint + permuted,
        "wrong_protein": query_wrong + transport,
        "ligand_only": ligand_only,
        "nogate": query_endpoint + transport,   # no gate in this architecture
    }


def predict_grammar(model, data, spec, episode, donor, device, dtype):
    support = episode.support_atoms.shape[0]
    atoms = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)
    atoms = atoms.to(device, dtype).unsqueeze(0)
    bonds = bonds.to(device, dtype).unsqueeze(0)
    mask = mask.to(device, dtype).unsqueeze(0)
    correct = protein_inputs(data, spec.target, device, dtype)
    wrong = protein_inputs(data, donor, device, dtype)
    endpoint, ligand_value = model.encode(
        *correct[:3], atoms, bonds, mask, correct[3])[:2]
    wrong_endpoint = model.encode(
        *wrong[:3], atoms, bonds, mask, wrong[3])[0]
    query_endpoint = endpoint[:, support:]
    query_wrong = wrong_endpoint[:, support:]
    ligand_only = ligand_value[:, support:]
    transport = torch.zeros_like(query_endpoint)
    permuted = torch.zeros_like(query_endpoint)
    level_only = torch.zeros_like(query_endpoint)
    if support:
        query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
        support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
        similarity = tanimoto(query_fp, support_fp)
        residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                    - endpoint[:, :support])
        shrink = model.transport.shrinkage(support, residual)
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
        "nogate": query_endpoint + transport,   # no gate in the incumbent
    }


def sign_accuracy(prediction: np.ndarray, truth: np.ndarray) -> float:
    comparable, correct = 0, 0
    for left in range(len(truth)):
        for right in range(left + 1, len(truth)):
            true_delta = float(truth[left] - truth[right])
            if true_delta == 0:
                continue
            comparable += 1
            correct += float((prediction[left] - prediction[right])
                             * true_delta > 0)
    return float(correct / comparable) if comparable else float("nan")


def cliff_sign_accuracy(prediction: np.ndarray, truth: np.ndarray,
                        similarity: np.ndarray, threshold: float,
                        gap: float) -> float:
    comparable, correct = 0, 0
    for left in range(len(truth)):
        for right in range(left + 1, len(truth)):
            if similarity[left, right] < threshold:
                continue
            true_delta = float(truth[left] - truth[right])
            if abs(true_delta) < gap:
                continue
            comparable += 1
            correct += float((prediction[left] - prediction[right])
                             * true_delta > 0)
    return float(correct / comparable) if comparable else float("nan")


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
                     split_directory=args.split_directory,
                     include_meta_test=False)
    scale = training_label_scale(data)
    # Donors from the evaluation split (identity alone varies); whitening
    # fitted on meta_train only (contract 2026-08-16).
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
                    if kind == "reltransport":
                        values = predict_reltransport(
                            model, data, spec, episode, donors[spec.target],
                            args.device, dtype)
                    elif kind == "level_shape":
                        values = predict_level_shape(
                            model, data, spec, episode, donors[spec.target],
                            args.device, dtype)
                    else:
                        values = predict_grammar(
                            model, data, spec, episode, donors[spec.target],
                            args.device, dtype)
                    truth = (episode.query_y.numpy() * scale.scale + scale.mean)
                    novelty = tanimoto_rows(
                        np.stack([fingerprints[data.cells[i]["ligand_id"]].numpy()
                                  for i in spec.query]), train_fp).max(-1)
                    query_sim = tanimoto_rows(
                        np.stack([fingerprints[data.cells[i]["ligand_id"]].numpy()
                                  for i in spec.query]),
                        np.stack([fingerprints[data.cells[i]["ligand_id"]].numpy()
                                  for i in spec.query]))
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
                            row[f"{label}_sign_accuracy"] = sign_accuracy(
                                prediction, truth)
                            row[f"{label}_cliff_sign_accuracy"] = \
                                cliff_sign_accuracy(
                                    prediction, truth, query_sim,
                                    threshold=0.6, gap=1.0)
                        if label == "full":
                            error = prediction - truth
                            row["calibration_pk"] = float(error.mean() ** 2)
                            row["shape_pk"] = float(error.var())
                            row["query_spread_pk"] = float(prediction.std())
                            if "wrong_ligand" in values:
                                wrong_ligand_pred = (
                                    values["wrong_ligand"].squeeze(0)
                                    .float().cpu().numpy()
                                    * scale.scale + scale.mean)
                                row["support_replacement_abs_delta_pk"] = float(
                                    np.abs(prediction - wrong_ligand_pred).mean())
                    rows.append(row)
        del model
        if args.device.startswith("cuda"):
            torch.cuda.empty_cache()

    names = [item.partition("=")[0] for item in args.arm]
    fields_of_interest = [f for f in rows[0]
                          if f.endswith(("_mse_pk", "_mse_pk_lt40", "_ci",
                                         "_spearman", "_sign_accuracy",
                                         "_cliff_sign_accuracy", "_pk"))]
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
        collected: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows:
            if row["arm"] == name and row["k"] == k and row.get(field) is not None:
                collected[(row["component"], row["target"])].append(
                    float(row[field]))
        return {key: float(np.mean(values)) for key, values in collected.items()}

    def paired(left: str, right: str, field: str, k: int) -> dict:
        a, b = per_target(left, field, k), per_target(right, field, k)
        values = [(component, target, b[(component, target)] - value)
                  for (component, target), value in a.items()
                  if (component, target) in b]
        return component_bootstrap(values, args.bootstrap_draws, 20260816)

    def internal(name: str, control: str, k: int) -> dict:
        a = per_target(name, "full_mse_pk", k)
        b = per_target(name, f"{control}_mse_pk", k)
        if not b:
            return {"absent": True}
        return component_bootstrap(
            [(component, target, b[(component, target)] - value)
             for (component, target), value in a.items()
             if (component, target) in b],
            args.bootstrap_draws, 20260816)

    contrasts: dict[str, dict] = {}
    for name in names:
        if name == args.reference:
            continue
        contrasts[f"{name}_vs_{args.reference}"] = {
            str(k): {field: paired(name, args.reference, field, k)
                     for field in ("full_mse_pk", "full_mse_pk_lt40",
                                   "zero_shot_mse_pk", "full_ci",
                                   "full_spearman", "full_sign_accuracy",
                                   "full_cliff_sign_accuracy",
                                   "zero_shot_cliff_sign_accuracy")}
            for k in SUPPORT_SIZES}
    for name in names:
        contrasts[f"{name}_internal"] = {
            str(k): {f"{control}_gap": internal(name, control, k)
                     for control in ("wrong_protein", "permuted",
                                     "ligand_only", "level_only", "nogate",
                                     "wrong_ligand")}
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
        "schema": "MetaSieve.StageR6ArmComparison.v1",
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
        print("  %-8s %9s %9s %9s %8s %8s %8s %9s" % (
            "arm", "full", "lt40", "wrongP", "CI", "rho", "sign", "cliff"))
        for name in names:
            entry = summary[name][str(k)]
            print("  %-8s %9.4f %9.4f %9.4f %8.4f %8.4f %8.4f %9.4f" % (
                name, entry["full_mse_pk"], entry["full_mse_pk_lt40"],
                entry["wrong_protein_mse_pk"], entry["full_ci"],
                entry["full_spearman"], entry["full_sign_accuracy"],
                entry["full_cliff_sign_accuracy"]))


if __name__ == "__main__":
    main()
