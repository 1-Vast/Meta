"""Shape-first routed training applied to the *incumbent* grammar trunk.

The R7/R8 cycle measured two facts: the incumbent trunk (A0,
`similarity_only`) owns the best k=0 calibration in this project (1.236),
and the shape-first training method owns the first real within-target shape
gain (0.943 -> 0.896) — but the two have never met: every shape-first run so
far used the factorized architecture, whose routed level path calibrates
worse than A0's unified trunk.

This script closes that cell with **zero architecture change**: the exact
`SimilarityGrammarModel` (frozen design, same parameters, same Tanimoto
transport) trained with the shape-first method.

    p_level = ligand_value + protein_value + interaction.detach()
              + transport.detach()
    p_shape = ligand_value + protein_value.detach() + interaction
              + transport
    L = mean(p_level - y)^2                    (level, routed)
      + 1.0 * pairwise_ranking(p_shape, y)     (RankNet, cliff-weighted)
      + 1.0 * var(p_shape - y)                 (shape variance)
      + counterfactual contrasts               (wrong protein shape/level,
                                                permuted labels k>=2,
                                                wrong support ligand k=1)

One backward pass, one optimizer update, single stage. The frozen A0
checkpoints (identical architecture, ordinary training) are the same-arm
ablation, so the training-method contribution is measured one-to-one.
No query-specific transport claim is made here: the transport is the
retained Stage 6/7 Tanimoto baseline, exactly as in A0.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from model.similarity_grammar import SimilarityGrammarModel, tanimoto
from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData, stable_seed
from scripts.train_level_shape import (
    matched_donors, normalized, protein_inputs,
)
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, LabelScale,
    compact_episode, file_sha256, training_label_scale,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)


@dataclass(frozen=True)
class GrammarShapeConfig:
    seed: int = 20260815
    evaluation_seed: int = 73101
    steps: int = 1200
    episodes_per_step: int = 3
    query_size: int = 16
    min_query_size: int = 4
    learning_rate: float = 6e-4
    weight_decay: float = 1e-5
    grad_clip: float = 1.0
    hidden_dim: int = 192
    task_dim: int = 48
    ligand_layers: int = 4
    pair_dim: int = 96
    pair_latents: int = 24
    pair_heads: int = 8
    cliff_pair_weight: float = 4.0
    cliff_tanimoto: float = 0.6
    cliff_gap_pk: float = 1.0
    counterfactual_weight: float = 0.25
    contrast_temperature: float = 0.1
    identify_weight: float = 0.3
    val_interval: int = 200
    lr_warmup_fraction: float = 0.05
    lr_final_fraction: float = 0.1
    amp: bool = True
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


def level_term(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    return (prediction - truth).mean(-1).square().mean()


def shape_variance(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    error = prediction - truth
    return (error - error.mean(-1, keepdim=True)).square().mean()


def pairwise_ranking(prediction: torch.Tensor, truth: torch.Tensor,
                     temperature: float,
                     pair_weight: torch.Tensor | None = None) -> torch.Tensor:
    delta_y = truth.unsqueeze(-1) - truth.unsqueeze(-2)
    delta_p = prediction.unsqueeze(-1) - prediction.unsqueeze(-2)
    comparable = delta_y != 0
    if not bool(comparable.any()):
        return prediction.new_zeros(())
    signed = delta_y.sign() * delta_p / temperature
    loss = F.softplus(-signed[comparable])
    if pair_weight is not None and pair_weight.numel() == loss.numel():
        loss = loss * pair_weight[comparable]
    return loss.mean()


def cliff_pair_weights(similarity: torch.Tensor, truth: torch.Tensor,
                       cliff_tanimoto: float, cliff_gap: float,
                       cliff_weight: float) -> torch.Tensor:
    delta_y = (truth.unsqueeze(-1) - truth.unsqueeze(-2)).abs()
    cliff = (similarity >= cliff_tanimoto) & (delta_y >= cliff_gap)
    return 1.0 + (cliff_weight - 1.0) * cliff.to(similarity.dtype)


def contrast(correct: torch.Tensor, wrong: torch.Tensor,
             temperature: float) -> torch.Tensor:
    logits = -torch.stack((correct, wrong)) / temperature
    return F.cross_entropy(logits.unsqueeze(0),
                           logits.new_zeros(1, dtype=torch.long))


def learning_rate_factor(step: int, config: GrammarShapeConfig) -> float:
    warmup = max(1, int(config.steps * config.lr_warmup_fraction))
    if step <= warmup:
        return step / warmup
    progress = (step - warmup) / max(1, config.steps - warmup)
    final = config.lr_final_fraction
    return final + (1.0 - final) * 0.5 * (1.0 + np.cos(np.pi * min(progress, 1.0)))


def episode_tensors(episode, device: str, dtype: torch.dtype):
    support = episode.support_atoms.shape[0]
    width = max(episode.support_atoms.shape[-2],
                episode.query_atoms.shape[-2], 1)

    def pad_atoms(values):
        return torch.nn.functional.pad(values, (0, 0, 0, width - values.shape[-2]))

    def pad_bonds(values):
        return torch.nn.functional.pad(
            values, (0, 0, 0, width - values.shape[-2],
                     0, width - values.shape[-2]))

    def pad_mask(values):
        return torch.nn.functional.pad(values, (0, width - values.shape[-1]))

    support_atoms = pad_atoms(episode.support_atoms)
    support_bonds = pad_bonds(episode.support_bonds)
    support_mask = pad_mask(episode.support_mask)
    query_atoms = pad_atoms(episode.query_atoms)
    query_bonds = pad_bonds(episode.query_bonds)
    query_mask = pad_mask(episode.query_mask)
    atoms = torch.cat((support_atoms, query_atoms), 0)
    bonds = torch.cat((support_bonds, query_bonds), 0)
    mask = torch.cat((support_mask, query_mask), 0)
    return (support, atoms.to(device, dtype).unsqueeze(0),
            bonds.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0))


def donor_ligand_graphs(data: QPSMPData, spec, device: str, dtype):
    donor_cells = [int(index) for index in data.tasks[spec.split][spec.donor_target]]
    seed = stable_seed("grammar-shape-donor-ligand", spec.split, spec.target,
                       spec.support, spec.query)
    order = np.random.default_rng(seed).permutation(donor_cells)
    indices = tuple(map(int, order[:len(spec.support)]))
    values = [data.ligand_bank.get(data.cells[index]["ligand_id"])
              for index in indices]
    max_atoms = max(value[0].shape[0] for value in values)
    atoms = torch.stack([torch.nn.functional.pad(
        torch.from_numpy(value[0]),
        (0, 0, 0, max_atoms - value[0].shape[0])) for value in values])
    bonds = torch.stack([torch.nn.functional.pad(
        torch.from_numpy(value[1]),
        (0, 0, 0, max_atoms - value[1].shape[0],
         0, max_atoms - value[1].shape[0])) for value in values])
    mask = torch.stack([torch.nn.functional.pad(
        torch.from_numpy(value[2]),
        (0, max_atoms - value[2].shape[0])) for value in values])
    fingerprints = torch.stack([
        data.fingerprints[data.cells[index]["ligand_id"]]
        for index in indices])
    # CPU tensors: the caller's episode stays on CPU and `episode_tensors`
    # moves everything to the device in one place.
    return (atoms.unsqueeze(0), bonds.unsqueeze(0), mask.unsqueeze(0),
            fingerprints.unsqueeze(0))


def forward_episode(model, data, episode, device, dtype, donor):
    """Full forward with the correct and donor proteins."""
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    correct = protein_inputs(data, episode.spec.target, device, dtype)
    wrong = protein_inputs(data, donor, device, dtype)
    query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
    support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
    similarity = tanimoto(query_fp, support_fp) if support else None
    full = model(
        *correct[:3], atoms[:, :support] if support else atoms[:, :0],
        bonds[:, :support] if support else bonds[:, :0],
        mask[:, :support] if support else mask[:, :0],
        episode.support_y.to(device, dtype).unsqueeze(0),
        atoms[:, support:], bonds[:, support:], mask[:, support:],
        protein_chemistry=correct[3],
        support_fingerprint=support_fp, query_fingerprint=query_fp)
    full_wrong = model(
        *wrong[:3], atoms[:, :support] if support else atoms[:, :0],
        bonds[:, :support] if support else bonds[:, :0],
        mask[:, :support] if support else mask[:, :0],
        episode.support_y.to(device, dtype).unsqueeze(0),
        atoms[:, support:], bonds[:, support:], mask[:, support:],
        protein_chemistry=wrong[3],
        support_fingerprint=support_fp, query_fingerprint=query_fp)
    return support, full, full_wrong, similarity


def episode_loss(model, data, episode, donors, config, dtype,
                 scale) -> tuple[torch.Tensor, dict]:
    device = config.device
    support, full, full_wrong, similarity = forward_episode(
        model, data, episode, device, dtype, donors[episode.spec.target])
    query_y = episode.query_y.to(device, dtype).unsqueeze(0)
    ligand_only = full.ligand_only
    protein_value = full.additive - full.ligand_only
    interaction = full.zero_shot - full.additive
    transport = full.adaptation
    wrong_ligand_only = full_wrong.ligand_only
    wrong_protein_value = full_wrong.additive - full_wrong.ligand_only
    wrong_interaction = full_wrong.zero_shot - full_wrong.additive
    wrong_transport = full_wrong.adaptation

    p_level = ligand_only + protein_value + interaction.detach() \
        + transport.detach()
    p_shape = ligand_only + protein_value.detach() + interaction + transport
    loss_level = level_term(p_level, query_y)
    pair_weight = None
    query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
    if query_y.shape[-1] >= 2:
        query_similarity = tanimoto(query_fp, query_fp)
        cliff_gap = config.cliff_gap_pk / scale.scale
        pair_weight = cliff_pair_weights(
            query_similarity, query_y, config.cliff_tanimoto, cliff_gap,
            config.cliff_pair_weight)
    loss_rank = pairwise_ranking(p_shape, query_y, 1.0, pair_weight)
    loss_var = shape_variance(p_shape, query_y)
    # Routing leaves the interaction branch's per-target mean free (null
    # direction of level and shape alike); the label-free pin closes it.
    loss_identify = config.identify_weight * interaction.mean(-1).square().mean()
    loss = loss_level + loss_rank + loss_var + loss_identify \
        + 0.05 * full.support_match_loss

    parts = {"level": float(loss_level.detach()),
             "rank": float(loss_rank.detach()),
             "var": float(loss_var.detach())}

    if config.counterfactual_weight > 0:
        # protein-shape: ligand prior detached, no ligand-only shortcut.
        correct_shape = shape_variance(
            ligand_only.detach() + interaction + transport, query_y)
        wrong_shape = shape_variance(
            wrong_ligand_only.detach() + wrong_interaction + wrong_transport,
            query_y)
        loss_protein_shape = contrast(
            correct_shape, wrong_shape, config.contrast_temperature)
        # protein-level: interaction detached, only protein_value responds.
        correct_level = level_term(
            ligand_only.detach() + protein_value + interaction.detach()
            + transport.detach(), query_y)
        wrong_level = level_term(
            ligand_only.detach() + wrong_protein_value + interaction.detach()
            + transport.detach(), query_y)
        loss_protein_level = contrast(
            correct_level, wrong_level, config.contrast_temperature)
        loss = loss + config.counterfactual_weight * (
            loss_protein_shape + loss_protein_level)
        parts["protein_shape"] = float(loss_protein_shape.detach())
        parts["protein_level"] = float(loss_protein_level.detach())

        frozen = full.zero_shot.detach()
        correct_mse = ((frozen + transport) - query_y).square().mean()
        loss_binding = correct_mse.new_zeros(())
        if support > 1:
            rolled = episode.support_y.roll(1, dims=0).to(device, dtype).unsqueeze(0)
            _, wrong_labels, _, _ = forward_episode(
                model, data,
                replace(episode, support_y=rolled.squeeze(0)),
                device, dtype, donors[episode.spec.target])
            loss_binding = loss_binding + contrast(
                correct_mse,
                ((frozen + wrong_labels.adaptation) - query_y).square().mean(),
                config.contrast_temperature)
        elif support == 1:
            w_atoms, w_bonds, w_mask, w_fp = donor_ligand_graphs(
                data, episode.spec, device, dtype)
            correct = protein_inputs(data, episode.spec.target, device, dtype)
            wrong_episode = replace(
                episode, support_atoms=w_atoms.squeeze(0),
                support_bonds=w_bonds.squeeze(0),
                support_mask=w_mask.squeeze(0),
                support_fingerprint=w_fp.squeeze(0))
            _, wrong_ligand, _, _ = forward_episode(
                model, data, wrong_episode, device, dtype,
                donors[episode.spec.target])
            loss_binding = loss_binding + contrast(
                correct_mse,
                ((frozen + wrong_ligand.adaptation) - query_y).square().mean(),
                config.contrast_temperature)
        loss = loss + config.counterfactual_weight * loss_binding
        parts["binding"] = float(loss_binding.detach())

    parts["prediction_mse"] = float(
        (full.prediction - query_y).square().mean().detach())
    return loss, parts


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


def evaluate(model, data, split, donors, scale, config, draws=1) -> dict:
    device = config.device
    dtype = next(model.parameters()).dtype
    banks = data.fixed_nested_episode_banks(
        split, SUPPORT_SIZES, config.query_size, draws,
        config.evaluation_seed, None)
    rows = []
    model.eval()
    with torch.no_grad():
        for k, specs in banks.items():
            for spec in specs:
                episode = compact_episode(normalized(data.materialize(spec), scale))
                support, full, full_wrong, _ = forward_episode(
                    model, data, episode, device, dtype, donors[spec.target])
                query_y = episode.query_y.to(device, dtype).unsqueeze(0)
                # Label-locked support residuals: r_k = y_k - f0(L_k). The
                # grammar output stores them as (locked - level_adjustment),
                # so the locked residual is the sum of the two.
                locked = (full.support_residual_quotient
                          + full.level_adjustment)
                shrink = (model.transport.shrinkage(support, locked)
                          if support else full.zero_shot.new_zeros(()))
                level_only = (full.zero_shot + shrink * locked.mean(-1, keepdim=True)
                              if support else full.zero_shot)
                wrong = full_wrong.prediction

                def pk(value):
                    return (value.squeeze(0).float().cpu().numpy()
                            * scale.scale + scale.mean)

                truth = pk(query_y)
                predictions = {
                    "full": pk(full.prediction),
                    "zero_shot": pk(full.zero_shot),
                    "level_only": pk(level_only),
                    "wrong_protein": pk(wrong),
                    "ligand_only": pk(full.ligand_only),
                }
                error = predictions["full"] - truth
                row = {"k": k, "component": spec.component,
                       "target": spec.target,
                       "calibration_pk": float(error.mean() ** 2),
                       "shape_pk": float(error.var()),
                       "endpoint_spread_pk": float(
                           predictions["zero_shot"].std())}
                for name, value in predictions.items():
                    row[f"{name}_mse_pk"] = float(((value - truth) ** 2).mean())
                ci, comparable = concordance_index(predictions["full"], truth)
                row["ci"] = ci if comparable else None
                row["spearman"] = spearman(predictions["full"], truth)
                row["sign_accuracy"] = sign_accuracy(
                    predictions["full"], truth)
                zero_ci, zero_comparable = concordance_index(
                    predictions["zero_shot"], truth)
                row["zero_shot_ci"] = zero_ci if zero_comparable else None
                rows.append(row)
    model.train()
    return {"rows": rows}


def component_target_mean(rows: list[dict], field: str, k: int | None = None
                          ) -> float:
    from collections import defaultdict
    by_target: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        if k is not None and row["k"] != k:
            continue
        value = row.get(field)
        if value is not None and np.isfinite(value):
            by_target[(row["component"], row["target"])].append(float(value))
    by_component: dict[str, list[float]] = defaultdict(list)
    for (component, _), values in by_target.items():
        by_component[component].append(float(np.mean(values)))
    if not by_component:
        return float("nan")
    return float(np.mean([np.mean(v) for v in by_component.values()]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=GrammarShapeConfig.seed)
    parser.add_argument("--steps", type=int, default=GrammarShapeConfig.steps)
    parser.add_argument("--episodes-per-step", type=int,
                        default=GrammarShapeConfig.episodes_per_step)
    parser.add_argument("--learning-rate", type=float,
                        default=GrammarShapeConfig.learning_rate)
    parser.add_argument("--cliff-pair-weight", type=float,
                        default=GrammarShapeConfig.cliff_pair_weight)
    parser.add_argument("--counterfactual-weight", type=float,
                        default=GrammarShapeConfig.counterfactual_weight)
    parser.add_argument("--val-interval", type=int,
                        default=GrammarShapeConfig.val_interval)
    parser.add_argument("--device", default=GrammarShapeConfig.device)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)

    config = GrammarShapeConfig(
        seed=args.seed, steps=args.steps,
        episodes_per_step=args.episodes_per_step,
        learning_rate=args.learning_rate,
        cliff_pair_weight=args.cliff_pair_weight,
        counterfactual_weight=args.counterfactual_weight,
        val_interval=args.val_interval, device=args.device)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=False)
    result = train(data, config, args.output)
    model, scale = result.pop("model"), result["scale"]
    donors_eval = result.pop("donors_eval")
    checkpoint_path = args.output / "checkpoint.pt"
    torch.save({"model_state": model.state_dict(), "config": asdict(config),
                "split_directory": str(args.split_directory)},
               checkpoint_path)
    rows = evaluate(model, data, "meta_val", donors_eval, scale, config)["rows"]
    (args.output / "PREDICTIONS_meta_val.jsonl").write_text(
        "\n".join(json.dumps(r, sort_keys=True) for r in rows) + "\n",
        encoding="utf-8")
    fields = [f for f in rows[0] if f.endswith("_mse_pk")] + [
        "ci", "spearman", "sign_accuracy", "zero_shot_ci",
        "calibration_pk", "shape_pk", "endpoint_spread_pk"]
    summary = {str(k): {f: component_target_mean(rows, f, k) for f in fields}
               for k in SUPPORT_SIZES}
    payload = {
        "schema": "MetaSieve.GrammarShapeFirstTraining.v1",
        "config": asdict(config),
        "split_directory": str(args.split_directory),
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "donors": {
            "training_counterfactual_pool": "meta_train",
            "evaluation_wrong_protein_pool": "meta_val",
            "whitening_pool": "meta_train",
            "metric": "esm_whitened (train-fitted)",
        },
        "meta_test": {"included": False, "evaluated": False},
        "checkpoint_sha256": file_sha256(checkpoint_path),
        "training": {k: v for k, v in result.items()
                     if k not in {"scale", "donors_eval"}},
        "label_scale": asdict(result["scale"]),
        "meta_val": summary,
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("\n%-3s %9s %9s %9s %9s %8s %8s" % (
        "k", "full", "zero", "ligonly", "wrongP", "CI", "sign"))
    for k in SUPPORT_SIZES:
        entry = summary[str(k)]
        print("%-3d %9.4f %9.4f %9.4f %9.4f %8.4f %8.4f" % (
            k, entry["full_mse_pk"], entry["zero_shot_mse_pk"],
            entry["ligand_only_mse_pk"], entry["wrong_protein_mse_pk"],
            entry["ci"], entry["sign_accuracy"]))


def train(data: QPSMPData, config: GrammarShapeConfig, output: Path) -> dict:
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = SimilarityGrammarModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers,
        pair_dim=config.pair_dim, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads,
        use_learned_key=False).to(config.device)
    dtype = next(model.parameters()).dtype
    scale = training_label_scale(data)
    donors_train = matched_donors(data, "meta_train", donor_pool="meta_train")
    donors_eval = matched_donors(data, "meta_val", donor_pool="meta_val",
                                 whitening_pool="meta_train")
    parameters = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.AdamW(parameters, lr=config.learning_rate,
                                  weight_decay=config.weight_decay)
    amp_enabled = config.amp and config.device.startswith("cuda")
    scaler = torch.amp.GradScaler("cuda", enabled=amp_enabled)
    best_state, best_value, best_step = None, float("inf"), 0
    trace = []
    started = time.monotonic()
    progress_path = output / "progress.jsonl"
    for step in range(1, config.steps + 1):
        factor = learning_rate_factor(step, config)
        for group in optimizer.param_groups:
            group["lr"] = config.learning_rate * factor
        optimizer.zero_grad(set_to_none=True)
        support_size = SUPPORT_SIZES[(step - 1) % len(SUPPORT_SIZES)]
        episodes = []
        for _ in range(config.episodes_per_step):
            spec = data.draw_episode(
                "meta_train", support_size,
                int(rng.integers(config.min_query_size, config.query_size + 1)),
                rng, min_query_size=config.min_query_size)
            episodes.append(compact_episode(normalized(data.materialize(spec), scale)))
        accumulated = {}
        with torch.autocast(
                device_type="cuda",
                dtype=(torch.bfloat16 if torch.cuda.is_bf16_supported()
                       else torch.float16), enabled=amp_enabled):
            for episode in episodes:
                loss, parts = episode_loss(
                    model, data, episode, donors_train, config, dtype, scale)
                scaler.scale(loss / len(episodes)).backward()
                for key, value in parts.items():
                    accumulated[key] = accumulated.get(key, 0.0) + value / len(episodes)
        scaler.unscale_(optimizer)
        grad_norm = float(torch.nn.utils.clip_grad_norm_(
            model.parameters(), config.grad_clip))
        scaler.step(optimizer)
        scaler.update()
        trace.append({"step": step, "grad_norm": grad_norm, **accumulated})
        if step % config.val_interval == 0 or step == config.steps:
            rows = evaluate(model, data, "meta_val", donors_eval, scale, config)["rows"]
            value = float(np.mean([component_target_mean(rows, "full_mse_pk", k)
                                   for k in SUPPORT_SIZES]))
            record = {"step": step, "val_mean_mse_pk": value,
                      "val_k0_mse_pk": component_target_mean(rows, "full_mse_pk", 0),
                      "val_k0_ci": component_target_mean(rows, "ci", 0),
                      "elapsed_seconds": time.monotonic() - started}
            print(json.dumps(record), flush=True)
            with progress_path.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(record) + "\n")
            if value < best_value:
                best_state = copy.deepcopy(model.state_dict())
                best_value, best_step = value, step
    if best_state is None:
        raise RuntimeError("training produced no validation checkpoint")
    model.load_state_dict(best_state)
    return {"model": model, "scale": scale,
            "donors_eval": donors_eval,
            "best_val_mean_mse_pk": best_value, "best_step": best_step,
            "loss_trace": trace[-50:],
            "trainable_parameters": sum(p.numel() for p in parameters),
            "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                    if config.device.startswith("cuda") else 0.0),
            "wall_seconds": time.monotonic() - started}


if __name__ == "__main__":
    main()
