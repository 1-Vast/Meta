"""Shape-first routed training for the direct-shape model (R13 candidate).

Same training recipe as the R9/R12 C2 configuration — routed level term,
RankNet ranking with cliff weight 2, shape variance 1.5, counterfactual
contrasts, identifiability pin — with one structural difference: the
relative supervision now targets the deployed ordering quantity itself,
`s(e_i) - s(e_j)` for the direct interaction-head shape, instead of a
bilinear potential that only approximates it. The transport is the retained
Tanimoto+key baseline; no query-specific gate exists.
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

from model.shape_direct import ShapeDirectModel
from model.similarity_grammar import tanimoto
from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData, stable_seed
from scripts.train_level_shape import (
    matched_donors, normalized, protein_inputs,
)
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, LabelScale,
    compact_episode, file_sha256, training_label_scale,
)
from scripts.train_reltransport import (
    cliff_pair_weights, contrast, level_term, learning_rate_factor,
    pairwise_ranking, shape_variance,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)


@dataclass(frozen=True)
class ShapeDirectConfig:
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
    anchors: int = 16
    shape_hidden: int = 96
    shape_variance_weight: float = 1.5
    difference_weight: float = 1.0
    cliff_pair_weight: float = 2.0
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


def episode_tensors(episode, device: str, dtype: torch.dtype):
    support = episode.support_atoms.shape[0]
    width = max(episode.support_atoms.shape[-2],
                episode.query_atoms.shape[-2], 1)

    def pad_atoms(values):
        return F.pad(values, (0, 0, 0, width - values.shape[-2]))

    def pad_bonds(values):
        return F.pad(values, (0, 0, 0, width - values.shape[-2],
                              0, width - values.shape[-2]))

    def pad_mask(values):
        return F.pad(values, (0, width - values.shape[-1]))

    atoms = torch.cat((pad_atoms(episode.support_atoms),
                       pad_atoms(episode.query_atoms)), 0)
    bonds = torch.cat((pad_bonds(episode.support_bonds),
                       pad_bonds(episode.query_bonds)), 0)
    mask = torch.cat((pad_mask(episode.support_mask),
                      pad_mask(episode.query_mask)), 0)
    return (support, atoms.to(device, dtype).unsqueeze(0),
            bonds.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0))


def donor_ligand_graphs(data: QPSMPData, spec, support: int):
    donor_cells = [int(index) for index in data.tasks[spec.split][spec.donor_target]]
    seed = stable_seed("shape-direct-donor-ligand", spec.split, spec.target,
                       spec.support, spec.query)
    order = np.random.default_rng(seed).permutation(donor_cells)
    indices = tuple(map(int, order[:support]))
    values = [data.ligand_bank.get(data.cells[index]["ligand_id"])
              for index in indices]
    max_atoms = max(value[0].shape[0] for value in values)
    atoms = torch.stack([F.pad(torch.from_numpy(value[0]),
                               (0, 0, 0, max_atoms - value[0].shape[0]))
                         for value in values])
    bonds = torch.stack([F.pad(torch.from_numpy(value[1]),
                               (0, 0, 0, max_atoms - value[1].shape[0],
                                0, max_atoms - value[1].shape[0]))
                         for value in values])
    mask = torch.stack([F.pad(torch.from_numpy(value[2]),
                              (0, max_atoms - value[2].shape[0]))
                        for value in values])
    fingerprints = torch.stack([
        data.fingerprints[data.cells[index]["ligand_id"]] for index in indices])
    return atoms, bonds, mask, fingerprints


def difference_supervision(shape: torch.Tensor, full_y: torch.Tensor,
                           scale: LabelScale, config: ShapeDirectConfig,
                           full_similarity: torch.Tensor | None
                           ) -> torch.Tensor:
    """s(e_i) - s(e_j) ~ y_i - y_j over every in-target pair, cliff-weighted.

    This is the deployed ordering quantity itself — the exact difference of
    the shape branch — so the supervision is direct by construction.
    """
    delta = shape.unsqueeze(-1) - shape.unsqueeze(-2)
    target = full_y.unsqueeze(-1) - full_y.unsqueeze(-2)
    pair_weight = None
    if full_similarity is not None:
        cliff_gap = config.cliff_gap_pk / scale.scale
        pair_weight = cliff_pair_weights(
            full_similarity, full_y, config.cliff_tanimoto, cliff_gap,
            config.cliff_pair_weight)
    error = (delta - target) ** 2
    mask = torch.ones_like(error) - torch.eye(
        error.shape[-1], device=error.device, dtype=error.dtype)
    if pair_weight is not None:
        error = error * pair_weight
    return (error * mask).sum() / mask.sum().clamp_min(1.0)


def episode_loss(model, data, episode, donors, config, dtype,
                 scale) -> tuple[torch.Tensor, dict]:
    device = config.device
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, episode.spec.target, device, dtype)
    donor_inputs = protein_inputs(data, donors[episode.spec.target], device, dtype)
    parts = model.forward_parts(pooled, tokens, protein_mask, atoms, bonds,
                                mask, chemistry)
    donor_parts = model.forward_parts(
        *donor_inputs[:3], atoms, bonds, mask, donor_inputs[3])
    endpoint, prior, level, shape, _, embed, _ = parts
    query_y = episode.query_y.to(device, dtype).unsqueeze(0)
    support_y = episode.support_y.to(device, dtype).unsqueeze(0)
    full_y = torch.cat((support_y, query_y), -1)
    query_prior, query_level = prior[:, support:], level[:, support:]
    query_shape = shape[:, support:]
    query_endpoint = endpoint[:, support:]
    query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
    support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
    query_similarity = tanimoto(query_fp, query_fp)
    full_similarity = tanimoto(
        torch.cat((support_fp, query_fp), 1),
        torch.cat((support_fp, query_fp), 1))
    transport = torch.zeros_like(query_endpoint)
    if support:
        residual = (support_y - endpoint[:, :support]).detach()
        shrink = model.transport.shrinkage(support, residual)
        transport, _ = model.transport(
            embed[:, :support], embed[:, support:], residual,
            tanimoto(query_fp, support_fp))
        transport = shrink * transport

    p_level = query_prior + query_level + query_shape.detach() \
        + transport.detach()
    p_shape = query_prior + query_level.detach() + query_shape + transport
    loss_level = level_term(p_level, query_y)
    cliff_gap = config.cliff_gap_pk / scale.scale
    pair_weight = cliff_pair_weights(
        query_similarity, query_y, config.cliff_tanimoto, cliff_gap,
        config.cliff_pair_weight)
    loss_rank = pairwise_ranking(p_shape, query_y, 1.0, pair_weight)
    loss_var = config.shape_variance_weight * shape_variance(p_shape, query_y)
    loss_diff = config.difference_weight * difference_supervision(
        shape, full_y, scale, config, full_similarity)
    loss_identify = config.identify_weight * query_shape.mean(-1).square().mean()
    loss = loss_level + loss_rank + loss_var + loss_diff + loss_identify
    parts_names = {"level": float(loss_level.detach()),
                   "rank": float(loss_rank.detach()),
                   "var": float(loss_var.detach()),
                   "diff": float(loss_diff.detach())}

    if config.counterfactual_weight > 0:
        donor_prior, donor_level = donor_parts[1], donor_parts[2]
        donor_shape = donor_parts[3]
        loss_protein_shape = contrast(
            shape_variance(query_prior.detach() + query_shape + transport,
                           query_y),
            shape_variance(donor_prior[:, support:].detach()
                           + donor_shape[:, support:] + transport.detach(),
                           query_y),
            config.contrast_temperature)
        loss_protein_level = contrast(
            level_term(query_prior.detach() + query_level
                       + query_shape.detach() + transport.detach(), query_y),
            level_term(query_prior.detach() + donor_level[:, support:]
                       + query_shape.detach() + transport.detach(), query_y),
            config.contrast_temperature)
        loss = loss + config.counterfactual_weight * (
            loss_protein_shape + loss_protein_level)
        parts_names["protein_shape"] = float(loss_protein_shape.detach())
        parts_names["protein_level"] = float(loss_protein_level.detach())

        frozen = query_endpoint.detach()
        correct_mse = ((frozen + transport) - query_y).square().mean()
        loss_binding = correct_mse.new_zeros(())
        if support > 1:
            rolled = support_y.roll(1, dims=-1)
            wrong_residual = (rolled - endpoint[:, :support]).detach()
            shrink = model.transport.shrinkage(support, wrong_residual)
            wrong_transport, _ = model.transport(
                embed[:, :support], embed[:, support:], wrong_residual,
                tanimoto(query_fp, support_fp))
            wrong_transport = shrink * wrong_transport
            loss_binding = loss_binding + contrast(
                correct_mse,
                ((frozen + wrong_transport) - query_y).square().mean(),
                config.contrast_temperature)
        elif support == 1:
            w_atoms, w_bonds, w_mask, w_fp = donor_ligand_graphs(
                data, episode.spec, support)
            wrong_episode = replace(
                episode, support_atoms=w_atoms, support_bonds=w_bonds,
                support_mask=w_mask, support_fingerprint=w_fp)
            w_support, w_atoms_t, w_bonds_t, w_mask_t = episode_tensors(
                wrong_episode, device, dtype)
            w_parts = model.forward_parts(
                pooled, tokens, protein_mask, w_atoms_t, w_bonds_t, w_mask_t,
                chemistry)
            wrong_residual = (support_y - w_parts[0][:, :1]).detach()
            shrink = model.transport.shrinkage(1, wrong_residual)
            wrong_transport, _ = model.transport(
                w_parts[5][:, :1], embed[:, support:], wrong_residual,
                tanimoto(query_fp, w_fp.to(device, dtype).unsqueeze(0)))
            wrong_transport = shrink * wrong_transport
            loss_binding = loss_binding + contrast(
                correct_mse,
                ((frozen + wrong_transport) - query_y).square().mean(),
                config.contrast_temperature)
        loss = loss + config.counterfactual_weight * loss_binding
        parts_names["binding"] = float(loss_binding.detach())

    parts_names["prediction_mse"] = float(
        (query_endpoint + transport - query_y).square().mean().detach())
    return loss, parts_names


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
                support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
                pooled, tokens, protein_mask, chemistry = protein_inputs(
                    data, spec.target, device, dtype)
                donor_inputs = protein_inputs(data, donors[spec.target], device, dtype)
                parts = model.forward_parts(pooled, tokens, protein_mask, atoms,
                                            bonds, mask, chemistry)
                donor_parts = model.forward_parts(
                    *donor_inputs[:3], atoms, bonds, mask, donor_inputs[3])
                endpoint, prior, level, shape, _, embed, _ = parts
                query_y = episode.query_y.to(device, dtype).unsqueeze(0)
                query_endpoint = endpoint[:, support:]
                wrong_endpoint = donor_parts[0][:, support:]
                ligand_only = prior[:, support:] + level[:, support:]
                transport = torch.zeros_like(query_endpoint)
                permuted = torch.zeros_like(query_endpoint)
                level_shift = torch.zeros_like(query_endpoint)
                if support:
                    residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                                - endpoint[:, :support])
                    shrink = model.transport.shrinkage(support, residual)
                    transport, _ = model.transport(
                        embed[:, :support], embed[:, support:], residual,
                        tanimoto(episode.query_fingerprint.to(device, dtype).unsqueeze(0),
                                 episode.support_fingerprint.to(device, dtype).unsqueeze(0)))
                    transport = shrink * transport
                    permuted_residual = (residual.roll(1, dims=-1) if support > 1
                                         else -residual)
                    permuted, _ = model.transport(
                        embed[:, :support], embed[:, support:],
                        permuted_residual,
                        tanimoto(episode.query_fingerprint.to(device, dtype).unsqueeze(0),
                                 episode.support_fingerprint.to(device, dtype).unsqueeze(0)))
                    permuted = shrink * permuted
                    level_shift = shrink * residual.mean(-1, keepdim=True)

                def pk(value):
                    return (value.squeeze(0).float().cpu().numpy()
                            * scale.scale + scale.mean)

                truth = pk(query_y)
                predictions = {
                    "full": pk(query_endpoint + transport),
                    "zero_shot": pk(query_endpoint),
                    "level_only": pk(query_endpoint + level_shift),
                    "permuted": pk(query_endpoint + permuted),
                    "wrong_protein": pk(wrong_endpoint + transport),
                    "ligand_only": pk(ligand_only),
                }
                error = predictions["full"] - truth
                row = {"k": k, "component": spec.component, "target": spec.target,
                       "calibration_pk": float(error.mean() ** 2),
                       "shape_pk": float(error.var()),
                       "endpoint_spread_pk": float(predictions["zero_shot"].std()),
                       "shape_abs_mean_pk": float(
                           shape[:, support:].squeeze(0).float().cpu().numpy()
                           .__abs__().mean() * scale.scale),
                       "level_pk": float(level[0, 0]) * scale.scale + scale.mean}
                for name, value in predictions.items():
                    row[f"{name}_mse_pk"] = float(((value - truth) ** 2).mean())
                ci, comparable = concordance_index(predictions["full"], truth)
                row["ci"] = ci if comparable else None
                row["spearman"] = spearman(predictions["full"], truth)
                row["sign_accuracy"] = sign_accuracy(predictions["full"], truth)
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


def train(data: QPSMPData, config: ShapeDirectConfig, output: Path) -> dict:
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = ShapeDirectModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_latents=config.pair_latents, pair_heads=config.pair_heads,
        anchors=config.anchors, shape_hidden=config.shape_hidden
    ).to(config.device)
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
    return {"model": model, "scale": scale, "donors_eval": donors_eval,
            "best_val_mean_mse_pk": best_value, "best_step": best_step,
            "loss_trace": trace[-50:],
            "trainable_parameters": sum(p.numel() for p in parameters),
            "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                    if config.device.startswith("cuda") else 0.0),
            "wall_seconds": time.monotonic() - started}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=ShapeDirectConfig.seed)
    parser.add_argument("--steps", type=int, default=ShapeDirectConfig.steps)
    parser.add_argument("--episodes-per-step", type=int,
                        default=ShapeDirectConfig.episodes_per_step)
    parser.add_argument("--learning-rate", type=float,
                        default=ShapeDirectConfig.learning_rate)
    parser.add_argument("--cliff-pair-weight", type=float,
                        default=ShapeDirectConfig.cliff_pair_weight)
    parser.add_argument("--shape-variance-weight", type=float,
                        default=ShapeDirectConfig.shape_variance_weight)
    parser.add_argument("--difference-weight", type=float,
                        default=ShapeDirectConfig.difference_weight)
    parser.add_argument("--val-interval", type=int,
                        default=ShapeDirectConfig.val_interval)
    parser.add_argument("--device", default=ShapeDirectConfig.device)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)
    config = ShapeDirectConfig(
        seed=args.seed, steps=args.steps,
        episodes_per_step=args.episodes_per_step,
        learning_rate=args.learning_rate,
        cliff_pair_weight=args.cliff_pair_weight,
        shape_variance_weight=args.shape_variance_weight,
        difference_weight=args.difference_weight,
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
        "calibration_pk", "shape_pk", "endpoint_spread_pk",
        "shape_abs_mean_pk", "level_pk"]
    summary = {str(k): {f: component_target_mean(rows, f, k) for f in fields}
               for k in SUPPORT_SIZES}
    payload = {
        "schema": "MetaSieve.ShapeDirectTraining.v1",
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


if __name__ == "__main__":
    main()
