"""CPU-budget episodic training for the biological QPSMP meta-learner.

This script is an implementation smoke, not a G2/G3 admission analysis.  Its
reported gate authorizations are deliberately always false.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
import sys

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
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    val_interval: int = 10
    val_draws_per_target: int = 1
    test_draws_per_target: int = 1
    eval_targets_per_component: int = 1
    grad_clip: float = 5.0
    zero_shot_loss_weight: float = 0.25
    section_mode: str = "support_span"
    interaction_mode: str = "atom_residue"
    zero_support_only: bool = False
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


def normalized_episode(episode: EpisodeBatch, scale: LabelScale) -> EpisodeBatch:
    return replace(episode, support_y=scale.normalize(episode.support_y),
                   query_y=scale.normalize(episode.query_y))


def forward(model: QPSMPBioModel, episode: EpisodeBatch, *, adapt: bool = True):
    active = max(int(episode.support_mask.sum(1).max()),
                 int(episode.query_mask.sum(1).max()))
    return model(
        episode.protein_pooled, episode.protein_tokens, episode.protein_mask,
        episode.support_atoms[:, :active],
        episode.support_bonds[:, :active, :active],
        episode.support_mask[:, :active], episode.support_y,
        episode.query_atoms[:, :active],
        episode.query_bonds[:, :active, :active],
        episode.query_mask[:, :active], adapt=adapt)


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
             bank: tuple[EpisodeSpec, ...], controls: bool,
             label_scale: LabelScale) -> dict:
    rows = []
    model.eval()
    with torch.no_grad():
        for spec in bank:
            episode = normalized_episode(data.materialize(spec), label_scale)
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
          support_sizes: tuple[int, ...] | None = None):
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = QPSMPBioModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers,
        section_mode=config.section_mode,
        interaction_mode=config.interaction_mode, dtype=torch.float32)
    model.to(config.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate,
                                  weight_decay=config.weight_decay)
    label_scale = training_label_scale(data)
    train_support_sizes = support_sizes or (config.support_size,)
    if not train_support_sizes or any(k < 1 for k in train_support_sizes):
        raise ValueError("support sizes must be positive")
    val_banks = {
        k: data.fixed_episode_bank(
            "meta_val", k, config.query_size,
            config.val_draws_per_target, config.seed,
            config.eval_targets_per_component)
        for k in train_support_sizes
    }
    best_state, best_value, best_step = None, float("inf"), 0
    trace = []
    for step in range(1, config.steps + 1):
        model.train()
        optimizer.zero_grad()
        losses = []
        for episode_index in range(config.episodes_per_step):
            support_size = train_support_sizes[
                ((step - 1) * config.episodes_per_step + episode_index)
                % len(train_support_sizes)
            ]
            spec = data.draw_episode(
                "meta_train", support_size, config.query_size, rng)
            episode = normalized_episode(data.materialize(spec), label_scale)
            full = forward(model, episode, adapt=not config.zero_support_only)
            endpoint_prediction = full.zero_shot if config.zero_support_only else full.prediction
            query_y = episode.query_y.to(
                device=endpoint_prediction.device, dtype=endpoint_prediction.dtype)
            endpoint_loss = (endpoint_prediction - query_y).square().mean()
            if episode.query_y.numel() > 1:
                query_residual = full.zero_shot - query_y
                zero_shot_loss = (
                    query_residual - query_residual.mean()
                ).square().mean()
            else:
                zero_shot_loss = endpoint_loss.new_zeros(())
            losses.append(endpoint_loss + config.zero_shot_loss_weight * zero_shot_loss)
        loss = torch.stack(losses).mean()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        optimizer.step()
        trace.append(float(loss.detach()))
        if step % config.val_interval == 0 or step == config.steps:
            selection_field = "zero_shot_mse_pk" if config.zero_support_only else "full_mse_pk"
            values = [evaluate(
                model, data, bank, controls=False,
                label_scale=label_scale)[selection_field] for bank in val_banks.values()]
            value = float(np.mean(values))
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
                   }}, label_scale


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=TrainConfig.seed)
    parser.add_argument("--support-size", type=int, default=TrainConfig.support_size)
    parser.add_argument("--query-size", type=int, default=TrainConfig.query_size)
    parser.add_argument("--steps", type=int, default=TrainConfig.steps)
    parser.add_argument("--episodes-per-step", type=int, default=TrainConfig.episodes_per_step)
    parser.add_argument("--val-interval", type=int, default=TrainConfig.val_interval)
    parser.add_argument("--val-draws-per-target", type=int, default=TrainConfig.val_draws_per_target)
    parser.add_argument("--test-draws-per-target", type=int, default=TrainConfig.test_draws_per_target)
    parser.add_argument("--eval-targets-per-component", type=int,
                        default=TrainConfig.eval_targets_per_component)
    parser.add_argument("--zero-shot-loss-weight", type=float,
                        default=TrainConfig.zero_shot_loss_weight)
    parser.add_argument("--section-mode", choices=("support_span", "ridge", "neural"),
                        default=TrainConfig.section_mode)
    parser.add_argument("--interaction-mode", choices=("pooled", "atom_residue"),
                        default=TrainConfig.interaction_mode)
    parser.add_argument("--zero-support-only", action="store_true")
    parser.add_argument("--device", default=TrainConfig.device)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    config = TrainConfig(
        seed=args.seed, support_size=args.support_size, query_size=args.query_size,
        steps=args.steps, episodes_per_step=args.episodes_per_step,
        val_interval=args.val_interval,
        val_draws_per_target=args.val_draws_per_target,
        test_draws_per_target=args.test_draws_per_target,
        eval_targets_per_component=args.eval_targets_per_component,
        zero_shot_loss_weight=args.zero_shot_loss_weight,
        section_mode=args.section_mode,
        interaction_mode=args.interaction_mode,
        zero_support_only=args.zero_support_only,
        device=args.device)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK)
    model, training, label_scale = train(data, config)
    test_bank = data.fixed_episode_bank(
        "meta_test", config.support_size, config.query_size,
        config.test_draws_per_target, config.seed, config.eval_targets_per_component)
    metrics = evaluate(model, data, test_bank, controls=True, label_scale=label_scale)
    args.output.mkdir(parents=True, exist_ok=False)
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
