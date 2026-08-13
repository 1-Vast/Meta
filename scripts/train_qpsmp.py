"""CPU-budget episodic training for the biological QPSMP meta-learner.

This script is an implementation smoke, not a G2/G3 admission analysis.  Its
reported gate authorizations are deliberately always false.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from model.qpsmp_meta import QPSMPBioModel
from scripts.qpsmp_data import EpisodeBatch, EpisodeSpec, QPSMPData


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT_LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact"
OUT = ROOT / "report/meta_fewshot/qpsmp_meta_smoke"


@dataclass(frozen=True)
class TrainConfig:
    seed: int = 20260812
    support_size: int = 5
    query_size: int = 8
    hidden_dim: int = 32
    task_dim: int = 4
    ligand_layers: int = 1
    steps: int = 20
    episodes_per_step: int = 2
    train_cache_size: int = 64
    eval_batch_size: int = 16
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    val_interval: int = 10
    val_draws_per_target: int = 1
    test_draws_per_target: int = 1
    eval_targets_per_component: int = 1
    grad_clip: float = 5.0
    zero_shot_loss_weight: float = 0.25
    section_mode: str = "support_span"
    interaction_mode: str = "pooled"
    zero_support_only: bool = False
    pretrained_checkpoint: str | None = None
    geometry_checkpoint: str | None = None
    section_only: bool = False
    amp: bool = True
    pair_blocks: int = 2
    pair_latents: int = 8
    pair_heads: int = 4
    pair_chunk_size: int = 32
    episode_cache: str | None = None
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


@dataclass(frozen=True)
class LabelScale:
    mean: float
    scale: float

    def normalize(self, values: torch.Tensor) -> torch.Tensor:
        return (values - self.mean) / self.scale

    def squared_error_pk(self, prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
        truth = truth.to(device=prediction.device, dtype=prediction.dtype)
        return ((prediction - truth) * self.scale).square()


def training_label_scale(data: QPSMPData) -> LabelScale:
    values = np.asarray([cell["pK"] for cell in data.cells
                         if cell["split"] == "meta_train"], dtype=np.float64)
    scale = float(values.std())
    if not np.isfinite(scale) or scale < 1e-6:
        raise ValueError("meta-train labels have invalid scale")
    return LabelScale(float(values.mean()), scale)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_parent_checkpoint(model: QPSMPBioModel, config: TrainConfig) -> str | None:
    if config.pretrained_checkpoint is None:
        return None
    path = Path(config.pretrained_checkpoint)
    payload = torch.load(path, map_location="cpu", weights_only=False)
    parent_config = payload.get("config", {})
    for field in (
            "hidden_dim", "task_dim", "ligand_layers", "interaction_mode",
            "pair_blocks", "pair_latents", "pair_heads"):
        if parent_config.get(field) != getattr(config, field):
            raise ValueError(f"parent checkpoint has incompatible {field}")
    incompatible = model.load_state_dict(payload["model_state"], strict=False)
    allowed_missing = (
        {key for key in incompatible.missing_keys if key.startswith("meta.qp_ams.")}
        if config.section_mode == "qp_ams" else set())
    if set(incompatible.missing_keys) != allowed_missing or incompatible.unexpected_keys:
        raise ValueError(
            "parent checkpoint state is incompatible: "
            f"missing={incompatible.missing_keys}, unexpected={incompatible.unexpected_keys}")
    return file_sha256(path)


def load_geometry_checkpoint(model: QPSMPBioModel, config: TrainConfig) -> str | None:
    if config.geometry_checkpoint is None:
        return None
    if config.interaction_mode != "bpsf":
        raise ValueError("geometry checkpoint requires interaction_mode=bpsf")
    path = Path(config.geometry_checkpoint)
    payload = torch.load(path, map_location="cpu", weights_only=False)
    teacher = payload.get("config", {})
    required = {
        "hidden_dim": config.hidden_dim,
        "section_dim": config.task_dim,
        "pair_blocks": config.pair_blocks,
        "pair_latents": config.pair_latents,
        "pair_heads": config.pair_heads,
        "gine_layers": config.ligand_layers,
    }
    for field, expected in required.items():
        if teacher.get(field) != expected:
            raise ValueError(f"geometry checkpoint has incompatible {field}")
    source = payload["model_state"]
    mappings = (
        ("protein.", "protein_encoder."),
        ("ligand.", "ligand_encoder."),
        ("bridge.trunk.", "pair_section."),
    )
    target = model.state_dict()
    loaded = set()
    with torch.no_grad():
        for source_prefix, target_prefix in mappings:
            for name, value in source.items():
                if not name.startswith(source_prefix):
                    continue
                target_name = target_prefix + name[len(source_prefix):]
                if target_name not in target or target[target_name].shape != value.shape:
                    raise ValueError(f"geometry tensor is incompatible: {name}")
                target[target_name].copy_(value)
                loaded.add(target_name)
    if not loaded or not any(name.startswith("pair_section.") for name in loaded):
        raise ValueError("geometry checkpoint did not contain a pair trunk")
    model.load_state_dict(target)
    return file_sha256(path)


def freeze_for_section_training(model: QPSMPBioModel) -> None:
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    for parameter in model.meta.section_head.parameters():
        parameter.requires_grad_(True)
    if model.meta.section_mode == "ridge":
        model.meta.support_span.ridge_raw.requires_grad_(True)
    elif model.meta.section_mode == "neural":
        for parameter in model.meta.adapter.parameters():
            parameter.requires_grad_(True)
    elif model.meta.section_mode == "qp_ams":
        for parameter in model.meta.qp_ams.parameters():
            parameter.requires_grad_(True)
    elif model.meta.section_mode == "section_former":
        for name, parameter in model.pair_section.latent.named_parameters():
            if not name.startswith("endpoint."):
                parameter.requires_grad_(True)
        for parameter in model.meta.section_former.parameters():
            parameter.requires_grad_(True)


def normalized_episode(episode: EpisodeBatch, scale: LabelScale) -> EpisodeBatch:
    return replace(episode, support_y=scale.normalize(episode.support_y),
                   query_y=scale.normalize(episode.query_y))


def compact_episode(episode: EpisodeBatch) -> EpisodeBatch:
    active = max(int(episode.support_mask.sum(1).max()),
                 int(episode.query_mask.sum(1).max()))
    def atoms(values: torch.Tensor) -> torch.Tensor:
        values = values[:, :active]
        return torch.nn.functional.pad(values, (0, 0, 0, active - values.shape[1])).clone()
    def bonds(values: torch.Tensor) -> torch.Tensor:
        values = values[:, :active, :active]
        missing = active - values.shape[1]
        return torch.nn.functional.pad(values, (0, 0, 0, missing, 0, missing)).clone()
    def mask(values: torch.Tensor) -> torch.Tensor:
        values = values[:, :active]
        return torch.nn.functional.pad(values, (0, active - values.shape[1])).clone()
    return replace(
        episode,
        support_atoms=atoms(episode.support_atoms),
        support_bonds=bonds(episode.support_bonds),
        support_mask=mask(episode.support_mask),
        query_atoms=atoms(episode.query_atoms),
        query_bonds=bonds(episode.query_bonds),
        query_mask=mask(episode.query_mask))


def stack_episodes(episodes: list[EpisodeBatch]) -> EpisodeBatch:
    if not episodes:
        raise ValueError("cannot stack an empty episode list")
    support_size = episodes[0].support_y.numel()
    query_size = episodes[0].query_y.numel()
    if any(item.support_y.numel() != support_size or item.query_y.numel() != query_size
           for item in episodes):
        raise ValueError("batched episodes must have equal support and query sizes")
    atom_count = max(item.support_atoms.shape[1] for item in episodes)

    def pad_graphs(values: torch.Tensor, *, bonds: bool = False,
                   mask: bool = False) -> torch.Tensor:
        missing = atom_count - values.shape[1]
        if missing == 0:
            return values
        if bonds:
            return torch.nn.functional.pad(values, (0, 0, 0, missing, 0, missing))
        if mask:
            return torch.nn.functional.pad(values, (0, missing))
        return torch.nn.functional.pad(values, (0, 0, 0, missing))

    return EpisodeBatch(
        tuple(item.spec for item in episodes),
        torch.stack([item.protein_pooled for item in episodes]),
        torch.stack([item.protein_tokens for item in episodes]),
        torch.stack([item.protein_mask for item in episodes]),
        torch.stack([pad_graphs(item.support_atoms) for item in episodes]),
        torch.stack([pad_graphs(item.support_bonds, bonds=True) for item in episodes]),
        torch.stack([pad_graphs(item.support_mask, mask=True) for item in episodes]),
        torch.stack([item.support_y for item in episodes]),
        torch.stack([pad_graphs(item.query_atoms) for item in episodes]),
        torch.stack([pad_graphs(item.query_bonds, bonds=True) for item in episodes]),
        torch.stack([pad_graphs(item.query_mask, mask=True) for item in episodes]),
        torch.stack([item.query_y for item in episodes]))


def forward(model: QPSMPBioModel, episode: EpisodeBatch, *, adapt: bool = True):
    if episode.support_atoms.ndim == 3:
        episode = compact_episode(episode)
    active = max(int(episode.support_mask.sum(-1).max()),
                 int(episode.query_mask.sum(-1).max()))
    if episode.support_atoms.ndim == 4:
        support_atoms = episode.support_atoms[:, :, :active]
        support_bonds = episode.support_bonds[:, :, :active, :active]
        support_mask = episode.support_mask[:, :, :active]
        query_atoms = episode.query_atoms[:, :, :active]
        query_bonds = episode.query_bonds[:, :, :active, :active]
        query_mask = episode.query_mask[:, :, :active]
    else:
        support_atoms = episode.support_atoms[:, :active]
        support_bonds = episode.support_bonds[:, :active, :active]
        support_mask = episode.support_mask[:, :active]
        query_atoms = episode.query_atoms[:, :active]
        query_bonds = episode.query_bonds[:, :active, :active]
        query_mask = episode.query_mask[:, :active]
    return model(
        episode.protein_pooled, episode.protein_tokens, episode.protein_mask,
        support_atoms, support_bonds, support_mask, episode.support_y,
        query_atoms, query_bonds, query_mask, adapt=adapt)


def replay_state(target_output, state: torch.Tensor) -> torch.Tensor:
    """Replace only the SAR state; retain target zero-shot, level, and basis."""
    return (target_output.level_baseline
            + target_output.sar_scale * (target_output.query_basis @ state))


def component_target_mean(rows: list[dict], field: str) -> float:
    target_values: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        target_values.setdefault((row["component"], row["target"]), []).append(row[field])
    component_values: dict[str, list[float]] = {}
    for (component, _), values in target_values.items():
        component_values.setdefault(component, []).append(float(np.mean(values)))
    return float(np.mean([np.mean(values) for values in component_values.values()]))


def centered_mse_pk(prediction: torch.Tensor, truth: torch.Tensor,
                    label_scale: LabelScale) -> torch.Tensor:
    truth = truth.to(device=prediction.device, dtype=prediction.dtype)
    error = prediction - truth
    error = error - error.mean(dim=-1, keepdim=error.ndim > 1)
    return error.square().mean() * label_scale.scale ** 2


def donor_state(model: QPSMPBioModel, data: QPSMPData, episode: EpisodeBatch,
                donor_target: str, label_scale: LabelScale, *,
                wrong_protein: bool) -> torch.Tensor:
    if wrong_protein:
        pooled, tokens, mask = data.protein_for_target(donor_target)
        support = replace(episode, protein_pooled=pooled,
                          protein_tokens=tokens, protein_mask=mask)
        return forward(model, support).task_state
    donor_indices = data.tasks[episode.spec.split][donor_target]
    count = len(episode.spec.support)
    order = np.random.default_rng(
        sum(episode.spec.support) + len(donor_target)).permutation(donor_indices)
    donor_spec = EpisodeSpec(
        episode.spec.split, data.cells[int(order[0])]["protein_group_40"], donor_target,
        tuple(map(int, order[:count])), (int(order[count]),), episode.spec.target)
    donor_episode = normalized_episode(data.materialize(donor_spec), label_scale)
    return forward(model, donor_episode).task_state


def wrong_protein_zero_shot(model: QPSMPBioModel, data: QPSMPData,
                            episode: EpisodeBatch, donor_target: str) -> torch.Tensor:
    pooled, tokens, mask = data.protein_for_target(donor_target)
    wrong_episode = replace(
        episode, protein_pooled=pooled, protein_tokens=tokens, protein_mask=mask)
    return forward(model, wrong_episode, adapt=False).zero_shot


def evaluate(model: QPSMPBioModel, data: QPSMPData,
             bank: tuple[EpisodeSpec | EpisodeBatch, ...], controls: bool,
             label_scale: LabelScale) -> dict:
    rows = []
    model.eval()
    with torch.no_grad():
        for item in bank:
            if isinstance(item, EpisodeSpec):
                spec = item
                episode = normalized_episode(data.materialize(spec), label_scale)
            else:
                episode = item
                spec = episode.spec
            full = forward(model, episode)
            frozen = forward(model, episode, adapt=False)
            evidence_null = replay_state(full, torch.zeros_like(full.task_state))
            sar_cut = full.prediction - full.sar_adaptation
            level_only = episode.support_y.mean().expand_as(episode.query_y)
            def mse_pk(prediction: torch.Tensor) -> float:
                return float(label_scale.squared_error_pk(prediction, episode.query_y).mean())
            values = {
                "full_mse_pk": mse_pk(full.prediction),
                "zero_shot_mse_pk": mse_pk(frozen.prediction),
                "sar_cut_mse_pk": mse_pk(sar_cut),
                "level_only_mse_pk": mse_pk(level_only),
                "no_interaction_mse_pk": mse_pk(full.additive),
                "ligand_only_mse_pk": mse_pk(full.ligand_only),
                "sar_centered_mse_pk": float(centered_mse_pk(
                    full.zero_shot + full.sar_adaptation,
                    episode.query_y, label_scale)),
                "level_adjustment_abs_mean_pk": float(
                    full.level_adjustment.abs().mean() * label_scale.scale),
                "sar_adaptation_abs_mean_pk": float(
                    full.sar_adaptation.abs().mean() * label_scale.scale),
                "evidence_score_mean": float(full.evidence_score),
                "level_shrinkage": float(full.level_shrinkage),
                "shape_scale": float(full.shape_scale),
                "sar_scale": float(full.sar_scale),
                "cross_zero_shot_abs_mean_pk": float(
                    full.cross_zero_shot.abs().mean() * label_scale.scale),
            }
            if controls:
                permuted = replace(episode, support_y=episode.support_y.roll(1))
                permuted_output = forward(model, permuted)
                foreign = donor_state(
                    model, data, episode, spec.donor_target, label_scale,
                    wrong_protein=False)
                wrong = donor_state(
                    model, data, episode, spec.donor_target, label_scale,
                    wrong_protein=True)
                values.update({
                    "permuted_mse_pk": mse_pk(replay_state(full, permuted_output.task_state)),
                    "foreign_state_mse_pk": mse_pk(replay_state(full, foreign)),
                    "wrong_protein_state_mse_pk": mse_pk(replay_state(full, wrong)),
                    "wrong_protein_zero_shot_mse_pk": mse_pk(
                        wrong_protein_zero_shot(
                            model, data, episode, spec.donor_target)),
                })
            rows.append({"component": spec.component, "target": spec.target, **values})
    metrics = {field: component_target_mean(rows, field)
               for field in rows[0] if field not in {"component", "target"}}
    metrics["level_plus_sar_gain_mse_pk"] = metrics["zero_shot_mse_pk"] - metrics["full_mse_pk"]
    metrics["sar_gain_mse_pk"] = metrics["sar_cut_mse_pk"] - metrics["full_mse_pk"]
    if controls:
        metrics["binding_did_mse_pk"] = metrics["permuted_mse_pk"] - metrics["full_mse_pk"]
        metrics["foreign_state_gap_mse_pk"] = metrics["foreign_state_mse_pk"] - metrics["full_mse_pk"]
        metrics["wrong_protein_gap_mse_pk"] = metrics["wrong_protein_state_mse_pk"] - metrics["full_mse_pk"]
        metrics["wrong_protein_zero_shot_gap_mse_pk"] = (
            metrics["wrong_protein_zero_shot_mse_pk"] - metrics["zero_shot_mse_pk"])
    metrics["weighting"] = "equal_component_then_equal_target_then_equal_draw"
    metrics["episodes"] = len(rows)
    return metrics


def train(data: QPSMPData, config: TrainConfig,
          support_sizes: tuple[int, ...] | None = None,
          progress_path: Path | None = None):
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = QPSMPBioModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers,
        section_mode=config.section_mode,
        interaction_mode=config.interaction_mode,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        dtype=torch.float32)
    model.to(config.device)
    geometry_sha256 = load_geometry_checkpoint(model, config)
    parent_sha256 = load_parent_checkpoint(model, config)
    if geometry_sha256 is not None and parent_sha256 is not None:
        raise ValueError("geometry and parent checkpoints are mutually exclusive initializers")
    if config.section_only:
        if config.section_mode not in {
                "ridge", "neural", "qp_ams", "section_former"} or parent_sha256 is None:
            raise ValueError("section-only training requires a section mode and parent checkpoint")
        freeze_for_section_training(model)
    trainable_parameters = [parameter for parameter in model.parameters()
                            if parameter.requires_grad]
    optimizer = torch.optim.AdamW(trainable_parameters, lr=config.learning_rate,
                                  weight_decay=config.weight_decay)
    amp_enabled = config.amp and config.device.startswith("cuda")
    scaler = torch.amp.GradScaler("cuda", enabled=amp_enabled)
    label_scale = training_label_scale(data)
    train_support_sizes = support_sizes or (config.support_size,)
    if not train_support_sizes or any(k < 1 for k in train_support_sizes):
        raise ValueError("support sizes must be positive")
    cache_contract = {
        "seed": config.seed, "support_sizes": list(train_support_sizes),
        "query_size": config.query_size, "train_cache_size": config.train_cache_size,
        "val_draws_per_target": config.val_draws_per_target,
        "eval_targets_per_component": config.eval_targets_per_component,
        "label_scale": asdict(label_scale),
    }
    cache_path = Path(config.episode_cache) if config.episode_cache else None
    if cache_path is not None and cache_path.exists():
        cached = torch.load(cache_path, map_location="cpu", weights_only=False)
        if cached.get("contract") != cache_contract:
            raise ValueError("episode cache contract does not match training configuration")
        train_cache = cached["train"]
        val_banks = cached["val"]
    else:
        val_specs = {
            k: data.fixed_episode_bank(
                "meta_val", k, config.query_size,
                config.val_draws_per_target, config.seed,
                config.eval_targets_per_component)
            for k in train_support_sizes
        }
        val_banks = {
            k: tuple(compact_episode(normalized_episode(data.materialize(spec), label_scale))
                     for spec in specs)
            for k, specs in val_specs.items()
        }
        train_cache: dict[int, tuple[EpisodeBatch, ...]] = {}
        for k in train_support_sizes:
            episodes = []
            while len(episodes) < config.train_cache_size:
                spec = data.draw_episode("meta_train", k, config.query_size, rng)
                if len(spec.query) != config.query_size:
                    continue
                episodes.append(compact_episode(
                    normalized_episode(data.materialize(spec), label_scale)))
            train_cache[k] = tuple(episodes)
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            torch.save({"schema": "MetaSieve.QPSMPEpisodeCache.v1",
                        "contract": cache_contract, "train": train_cache,
                        "val": val_banks}, cache_path)
    best_state, best_value, best_step = None, float("inf"), 0
    trace = []
    started = time.monotonic()
    for step in range(1, config.steps + 1):
        model.train()
        optimizer.zero_grad()
        support_size = train_support_sizes[(step - 1) % len(train_support_sizes)]
        cache = train_cache[support_size]
        offset = ((step - 1) * config.episodes_per_step) % len(cache)
        selected = [cache[(offset + index) % len(cache)]
                    for index in range(config.episodes_per_step)]
        episode = stack_episodes(selected)
        with torch.autocast(
                device_type="cuda", dtype=torch.float16, enabled=amp_enabled):
            full = forward(model, episode, adapt=not config.zero_support_only)
            endpoint_prediction = full.zero_shot if config.zero_support_only else full.prediction
            query_y = episode.query_y.to(
                device=endpoint_prediction.device, dtype=endpoint_prediction.dtype)
            endpoint_loss = (endpoint_prediction - query_y).square().mean()
            query_residual = full.zero_shot - query_y
            zero_shot_loss = (
                query_residual - query_residual.mean(dim=-1, keepdim=True)
            ).square().mean()
            if config.section_only:
                sar_error = full.zero_shot + full.sar_adaptation - query_y
                sar_error = sar_error - sar_error.mean(dim=-1, keepdim=True)
                loss = sar_error.square().mean()
            else:
                loss = endpoint_loss + config.zero_shot_loss_weight * zero_shot_loss
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        scaler.step(optimizer)
        scaler.update()
        trace.append(float(loss.detach()))
        if step % config.val_interval == 0 or step == config.steps:
            selection_field = (
                "sar_centered_mse_pk" if config.section_only else
                "zero_shot_mse_pk" if config.zero_support_only else "full_mse_pk")
            values = [evaluate(
                model, data, bank, controls=False,
                label_scale=label_scale)[selection_field] for bank in val_banks.values()]
            value = float(np.mean(values))
            progress = {
                "step": step, "validation_mse_pk": value,
                "best_validation_mse_pk": min(best_value, value),
                "elapsed_seconds": time.monotonic() - started,
                "device": config.device,
            }
            line = json.dumps(progress)
            print(line, flush=True)
            if progress_path is not None:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(line + "\n")
            if value < best_value:
                best_state, best_value, best_step = copy.deepcopy(model.state_dict()), value, step
    if best_state is None:
        raise RuntimeError("training produced no validation checkpoint")
    model.load_state_dict(best_state)
    return model, {"best_val_component_target_mse_pk": best_value,
                   "best_step": best_step, "loss_trace": trace,
                   "loss_trace_units": "standardized_squared_error",
                   "label_scale": asdict(label_scale),
                   "support_sizes": list(train_support_sizes),
                   "validation_episode_bank_sizes": {
                       str(k): len(bank) for k, bank in val_banks.items()
                   },
                   "train_episode_cache_sizes": {
                       str(k): len(bank) for k, bank in train_cache.items()
                   },
                   "parent_checkpoint_sha256": parent_sha256,
                   "geometry_checkpoint_sha256": geometry_sha256,
                   "amp_enabled": amp_enabled,
                   "episode_cache": str(cache_path) if cache_path else None,
                   "peak_cuda_memory_mb": (
                       torch.cuda.max_memory_allocated() / 2 ** 20
                       if config.device.startswith("cuda") else 0.0),
                   "trainable_parameter_names": [name for name, parameter in
                                                   model.named_parameters()
                                                   if parameter.requires_grad]}, label_scale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=TrainConfig.seed)
    parser.add_argument("--support-size", type=int, default=TrainConfig.support_size)
    parser.add_argument("--query-size", type=int, default=TrainConfig.query_size)
    parser.add_argument("--hidden-dim", type=int, default=TrainConfig.hidden_dim)
    parser.add_argument("--task-dim", type=int, default=TrainConfig.task_dim)
    parser.add_argument("--steps", type=int, default=TrainConfig.steps)
    parser.add_argument("--learning-rate", type=float, default=TrainConfig.learning_rate)
    parser.add_argument("--episodes-per-step", type=int, default=TrainConfig.episodes_per_step)
    parser.add_argument("--train-cache-size", type=int, default=TrainConfig.train_cache_size)
    parser.add_argument("--val-interval", type=int, default=TrainConfig.val_interval)
    parser.add_argument("--val-draws-per-target", type=int, default=TrainConfig.val_draws_per_target)
    parser.add_argument("--test-draws-per-target", type=int, default=TrainConfig.test_draws_per_target)
    parser.add_argument("--eval-targets-per-component", type=int,
                        default=TrainConfig.eval_targets_per_component)
    parser.add_argument("--zero-shot-loss-weight", type=float,
                        default=TrainConfig.zero_shot_loss_weight)
    parser.add_argument("--section-mode", choices=(
        "support_span", "ridge", "neural", "qp_ams", "section_former"),
                        default=TrainConfig.section_mode)
    parser.add_argument("--interaction-mode", choices=("pooled", "atom_residue", "bpsf"),
                        default=TrainConfig.interaction_mode)
    parser.add_argument("--zero-support-only", action="store_true")
    parser.add_argument("--pretrained-checkpoint", type=Path)
    parser.add_argument("--geometry-checkpoint", type=Path)
    parser.add_argument("--section-only", action="store_true")
    parser.add_argument("--no-amp", action="store_true")
    parser.add_argument("--pair-blocks", type=int, default=TrainConfig.pair_blocks)
    parser.add_argument("--pair-latents", type=int, default=TrainConfig.pair_latents)
    parser.add_argument("--pair-heads", type=int, default=TrainConfig.pair_heads)
    parser.add_argument("--pair-chunk-size", type=int, default=TrainConfig.pair_chunk_size)
    parser.add_argument("--episode-cache", type=Path)
    parser.add_argument("--device", default=TrainConfig.device)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)
    config = TrainConfig(
        seed=args.seed, support_size=args.support_size, query_size=args.query_size,
        hidden_dim=args.hidden_dim, task_dim=args.task_dim,
        steps=args.steps, episodes_per_step=args.episodes_per_step,
        learning_rate=args.learning_rate,
        train_cache_size=args.train_cache_size,
        val_interval=args.val_interval,
        val_draws_per_target=args.val_draws_per_target,
        test_draws_per_target=args.test_draws_per_target,
        eval_targets_per_component=args.eval_targets_per_component,
        zero_shot_loss_weight=args.zero_shot_loss_weight,
        section_mode=args.section_mode,
        interaction_mode=args.interaction_mode,
        zero_support_only=args.zero_support_only,
        pretrained_checkpoint=(str(args.pretrained_checkpoint.resolve())
                               if args.pretrained_checkpoint else None),
        geometry_checkpoint=(str(args.geometry_checkpoint.resolve())
                             if args.geometry_checkpoint else None),
        section_only=args.section_only,
        amp=not args.no_amp,
        pair_blocks=args.pair_blocks, pair_latents=args.pair_latents,
        pair_heads=args.pair_heads, pair_chunk_size=args.pair_chunk_size,
        episode_cache=(str(args.episode_cache.resolve()) if args.episode_cache else None),
        device=args.device)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    model, training, label_scale = train(
        data, config, progress_path=args.output / "progress.jsonl")
    test_bank = data.fixed_episode_bank(
        "meta_test", config.support_size, config.query_size,
        config.test_draws_per_target, config.seed, config.eval_targets_per_component)
    metrics = evaluate(model, data, test_bank, controls=True, label_scale=label_scale)
    checkpoint = args.output / "checkpoint.pt"
    torch.save({"model_state": model.state_dict(), "config": asdict(config)}, checkpoint)
    result = {
        "schema": "MetaSieve.QPSMPMetaSmoke.v1",
        "scope": "implementation_smoke_only",
        "data": {"corpus": str(CORPUS), "protein_bank_records": len(data.protein_bank),
                 "ligand_bank_records": len(data.ligand_bank)},
        "config": asdict(config), "training": training, "test": metrics,
        "evaluation_population": "fixed hash-selected targets within every eligible component",
        "controls": {"foreign": "only task_state is replaced; target level and query channels stay fixed",
                     "wrong_protein": "only task_state inferred under the donor protein is replaced",
                     "evidence_null": "SAR state is zero; zero-shot and level calibration stay fixed",
                     "sar_cut": "the SAR term alone is removed; zero-shot and level stay fixed",
                     "permuted": "support labels are cyclically permuted"},
        "gate_authorization": {"G2": False, "G3a": False, "G3b": False},
        "authorization_reason": "A training smoke cannot authorize preregistered inferential gates.",
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
