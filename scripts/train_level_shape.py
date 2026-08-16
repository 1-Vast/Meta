"""Core Innovation B: counterfactual level-shape gradient-routed episodic training.

The optimization principle
--------------------------
For any episode the squared error decomposes **exactly**:

    MSE = mean(p - y)^2 + var(p - y)
          \___________/   \________/
           level            shape

Stage 9 measured that split on cold targets: 59% level, 41% shape, with the
trained endpoint contributing essentially nothing to shape. Ordinary training
lets one scalar head buy almost all of its squared error from the level term,
and a level shift is constant across the queries of a target, so it can never
change ordering. Four architectural interventions failed against that gradient.

This method matches the **algebraic** decomposition of the objective to the
**architectural** decomposition of the model, and routes each term only to the
components that own it:

    p_level = ligand_prior + target_level + centered.detach()
    p_shape = ligand_prior + target_level.detach() + centered
    L_route = mean(p_level - y)^2 + var(p_shape - y)

`L_route` is numerically identical to `MSE(p)` — this is not an added auxiliary
loss and it does not reweight anything. Only the gradient paths differ:

* the level term cannot train the interaction branch, so calibration pressure
  can no longer collapse it into a per-target constant;
* the shape term reaches `centered_interaction` and `ligand_prior`;
* `target_level` is constant across queries, so `var` is analytically
  independent of it — detaching it there is exact rather than approximate.

Counterfactual supervision, all in the same step
------------------------------------------------
Three contrasts, each routed so that it cannot be satisfied by the wrong module:

* **protein-shape**: correct protein must beat a *similarity-matched* wrong
  protein on the shape term, with `ligand_prior` detached — so a ligand-only
  shortcut cannot satisfy it, structurally rather than by penalty;
* **protein-level**: correct protein must beat the same donor on the level term,
  with the interaction detached — so only `target_level` can respond;
* **support-binding**: correct support labels must beat a permutation, with the
  endpoint detached — so only the transport can respond, and since permutation
  leaves `mean(residual)` invariant this isolates the query-specific channel.

Wrong-protein donors are the most similar training target from a **different**
homology component under Stage R2's `esm_whitened` metric — the hardest
available control, not a random cross-component protein.

One backward pass, one optimizer step, single stage. No warmup, no alternating
phase, no separately trained adapter, no inner loop, no test-time gradient.
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

from model.level_shape import LevelShapeModel
from model.similarity_grammar import tanimoto
from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import EpisodeBatch, EpisodeSpec, QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, LabelScale,
    compact_episode, file_sha256, training_label_scale,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)
MODULE_GROUPS = ("ligand_encoder", "channels", "protein_encoder", "level_head",
                 "ligand_head", "interaction", "anchor", "transport")


@dataclass(frozen=True)
class RouteConfig:
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
    pair_heads: int = 4
    anchors: int = 32
    ranking_loss_weight: float = 0.5
    ranking_temperature: float = 1.0
    counterfactual_weight: float = 0.25
    contrast_temperature: float = 0.1
    routing: bool = True
    counterfactual: bool = True
    val_interval: int = 100
    diagnostic_interval: int = 200
    lr_warmup_fraction: float = 0.05
    lr_final_fraction: float = 0.1
    amp: bool = True
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


# ---------------------------------------------------------------- objectives


def level_term(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    return (prediction - truth).mean(-1).square().mean()


def shape_term(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    error = prediction - truth
    return (error - error.mean(-1, keepdim=True)).square().mean()


def pairwise_ranking(prediction: torch.Tensor, truth: torch.Tensor,
                     temperature: float) -> torch.Tensor:
    delta_y = truth.unsqueeze(-1) - truth.unsqueeze(-2)
    delta_p = prediction.unsqueeze(-1) - prediction.unsqueeze(-2)
    comparable = delta_y != 0
    if not bool(comparable.any()):
        return prediction.new_zeros(())
    signed = delta_y.sign() * delta_p / temperature
    return F.softplus(-signed[comparable]).mean()


def contrast(correct: torch.Tensor, wrong: torch.Tensor,
             temperature: float) -> torch.Tensor:
    """Prefer the correct configuration over its counterfactual."""
    logits = -torch.stack((correct, wrong)) / temperature
    return F.cross_entropy(logits.unsqueeze(0),
                           logits.new_zeros(1, dtype=torch.long))


def learning_rate_factor(step: int, config: RouteConfig) -> float:
    warmup = max(1, int(config.steps * config.lr_warmup_fraction))
    if step <= warmup:
        return step / warmup
    progress = (step - warmup) / max(1, config.steps - warmup)
    final = config.lr_final_fraction
    return final + (1.0 - final) * 0.5 * (1.0 + np.cos(np.pi * min(progress, 1.0)))


# ---------------------------------------------------------------- data helpers


def normalized(episode: EpisodeBatch, scale: LabelScale) -> EpisodeBatch:
    return replace(episode, support_y=scale.normalize(episode.support_y),
                   query_y=scale.normalize(episode.query_y))


def matched_donors(data: QPSMPData, split: str,
                   donor_pool: str = "meta_train",
                   whitening_pool: str = "meta_train") -> dict[str, str]:
    """Most similar cross-component donor for every target in `split`.

    Stage R2 selected `esm_whitened` for explicit protein similarity: raw pooled
    cosine occupies a band of width 0.21, so a "most similar" donor chosen with
    it is barely distinguishable from a random one. Whitening widens the usable
    spread fourteen-fold, which makes this control genuinely hard.

    Two pool parameters, both part of the frozen evaluation contract
    (2026-08-16):

    * `donor_pool` supplies the donor **candidates**. During training both the
      correct and the donor protein are `meta_train` targets, so that contrast
      is clean with a `meta_train` pool. At **evaluation** a `meta_train` donor
      would confound wrong identity with seen-versus-unseen: the model has
      fitted the donor's level and not the query target's, so a wrong protein
      can win for a reason that has nothing to do with specificity. Evaluation
      therefore draws candidates from the evaluation split itself, where both
      proteins are equally unseen.
    * `whitening_pool` supplies the mean and covariance of the whitening
      transform. The contract requires `meta_train` here, **always** — the
      evaluation population must never contribute its own whitening statistics.
    """
    whitening_targets = sorted(data.tasks[whitening_pool])
    donor_targets = sorted(data.tasks[donor_pool])
    component = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    needed = set(whitening_targets) | set(donor_targets) | set(data.tasks[split])
    pooled = {t: np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
              for t in needed}
    matrix = np.stack([pooled[t] for t in whitening_targets])
    center = matrix.mean(0, keepdims=True)
    deviation = matrix - center
    covariance = deviation.T @ deviation / max(len(matrix) - 1, 1)
    values, vectors = np.linalg.eigh(covariance.astype(np.float64))
    whiten = ((vectors / np.sqrt(np.maximum(values, 1e-3)))
              @ vectors.T).astype(np.float32)

    def transform(vector: np.ndarray) -> np.ndarray:
        out = (vector - center[0]) @ whiten.T
        return out / max(float(np.linalg.norm(out)), 1e-9)

    bank = np.stack([transform(pooled[t]) for t in donor_targets])
    donors: dict[str, str] = {}
    for target in sorted(set(data.tasks[split])):
        similarity = bank @ transform(pooled[target])
        for index in np.argsort(-similarity):
            candidate = donor_targets[int(index)]
            if component[candidate] != component[target]:
                donors[target] = candidate
                break
        if target not in donors:
            raise ValueError(f"no cross-component donor exists for {target}")
    return donors


def protein_inputs(data: QPSMPData, target: str, device: str,
                   dtype: torch.dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return (pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0))


def episode_tensors(episode: EpisodeBatch, device: str, dtype: torch.dtype):
    support = episode.support_atoms.shape[0]
    atoms = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)
    return (support, atoms.to(device, dtype).unsqueeze(0),
            bonds.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0))


def episode_loss(model: LevelShapeModel, data: QPSMPData, episode: EpisodeBatch,
                 donors: dict[str, str], config: RouteConfig,
                 dtype: torch.dtype) -> tuple[torch.Tensor, dict]:
    device = config.device
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, episode.spec.target, device, dtype)
    channels = model.encode_ligand(atoms, bonds, mask)
    endpoint, prior, level, centered, _ = model.endpoint_with_channels(
        channels, pooled, tokens, protein_mask, chemistry)
    query_y = episode.query_y.to(device, dtype).unsqueeze(0)
    support_y = episode.support_y.to(device, dtype).unsqueeze(0)

    def split(value):
        return value[:, :support], value[:, support:]

    support_endpoint, query_endpoint = split(endpoint)
    _, query_prior = split(prior)
    _, query_level = split(level)
    _, query_centered = split(centered)

    transport = torch.zeros_like(query_endpoint)
    permuted_transport = None
    if support:
        similarity = tanimoto(
            episode.query_fingerprint.to(device, dtype).unsqueeze(0),
            episode.support_fingerprint.to(device, dtype).unsqueeze(0))
        residual = (support_y - support_endpoint).detach()
        shrink = model.transport.shrinkage(support, residual)
        transport = shrink * model.transport(residual, similarity)[0]
        if support > 1:
            rolled = residual.roll(1, dims=-1)
            permuted_transport = shrink * model.transport(rolled, similarity)[0]

    prediction = query_endpoint + transport
    if config.routing:
        p_level = query_prior + query_level + query_centered.detach() \
            + transport.detach()
        p_shape = query_prior + query_level.detach() + query_centered + transport
    else:
        p_level = p_shape = prediction
    loss_level = level_term(p_level, query_y)
    loss_shape = shape_term(p_shape, query_y)
    loss_rank = pairwise_ranking(p_shape, query_y, config.ranking_temperature)
    loss = loss_level + loss_shape + config.ranking_loss_weight * loss_rank

    parts = {"level": float(loss_level.detach()),
             "shape": float(loss_shape.detach()),
             "rank": float(loss_rank.detach())}

    if config.routing:
        # Routing creates one unidentifiable direction: a constant added to
        # `centered` within a target is detached out of the level term and
        # invisible to the shape term, which is a variance. Left free it drifts
        # and `target_level` silently compensates, which destroys the meaning of
        # the factorization even though the sum still fits.
        #
        # This term pins it. It contains **no label**, so it cannot improve fit
        # or trade off against anything; it is an identifiability constraint on
        # the anchor reference set, not an objective. It is switched on exactly
        # when routing is, because routing is what creates the freedom.
        loss_identify = query_centered.mean(-1).square().mean()
        loss = loss + loss_identify
        parts["identify"] = float(loss_identify.detach())
    if config.counterfactual:
        donor = donors[episode.spec.target]
        donor_inputs = protein_inputs(data, donor, device, dtype)
        _, wrong_prior, wrong_level, wrong_centered, _ = \
            model.endpoint_with_channels(channels, *donor_inputs)
        _, wrong_query_level = split(wrong_level)
        _, wrong_query_centered = split(wrong_centered)
        del wrong_prior                       # protein-blind by construction

        # shape contrast: ligand_prior detached, so no ligand-only shortcut
        correct_shape = shape_term(
            query_prior.detach() + query_centered, query_y)
        wrong_shape = shape_term(
            query_prior.detach() + wrong_query_centered, query_y)
        loss_protein_shape = contrast(
            correct_shape, wrong_shape, config.contrast_temperature)

        # level contrast: interaction detached, so only target_level responds
        correct_level = level_term(
            query_prior.detach() + query_level + query_centered.detach(), query_y)
        wrong_level_error = level_term(
            query_prior.detach() + wrong_query_level + query_centered.detach(),
            query_y)
        loss_protein_level = contrast(
            correct_level, wrong_level_error, config.contrast_temperature)

        loss = loss + config.counterfactual_weight * (
            loss_protein_shape + loss_protein_level)
        parts["protein_shape"] = float(loss_protein_shape.detach())
        parts["protein_level"] = float(loss_protein_level.detach())

        if permuted_transport is not None:
            # endpoint detached, so only the transport can satisfy this
            frozen = query_endpoint.detach()
            loss_binding = contrast(
                ((frozen + transport) - query_y).square().mean(),
                ((frozen + permuted_transport) - query_y).square().mean(),
                config.contrast_temperature)
            loss = loss + config.counterfactual_weight * loss_binding
            parts["binding"] = float(loss_binding.detach())

    parts["prediction_mse"] = float(
        (prediction - query_y).square().mean().detach())
    return loss, parts


# ---------------------------------------------------------------- diagnostics


def module_of(name: str) -> str:
    for group in MODULE_GROUPS:
        if name.startswith(group):
            return group
    return "other"


def gradient_report(model: LevelShapeModel, data: QPSMPData,
                    episode: EpisodeBatch, donors: dict[str, str],
                    config: RouteConfig, dtype: torch.dtype) -> dict:
    """Per-objective, per-module gradient norms and pairwise cosines."""
    device = config.device
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, episode.spec.target, device, dtype)
    channels = model.encode_ligand(atoms, bonds, mask)
    endpoint, prior, level, centered, _ = model.endpoint_with_channels(
        channels, pooled, tokens, protein_mask, chemistry)
    query_y = episode.query_y.to(device, dtype).unsqueeze(0)
    query_prior, query_level = prior[:, support:], level[:, support:]
    query_centered = centered[:, support:]
    del endpoint

    objectives = {}
    if config.routing:
        objectives["level"] = level_term(
            query_prior + query_level + query_centered.detach(), query_y)
        objectives["shape"] = shape_term(
            query_prior + query_level.detach() + query_centered, query_y)
    else:
        joint = query_prior + query_level + query_centered
        objectives["level"] = level_term(joint, query_y)
        objectives["shape"] = shape_term(joint, query_y)

    names = [n for n, p in model.named_parameters() if p.requires_grad]
    parameters = [p for _, p in model.named_parameters() if p.requires_grad]
    flat: dict[str, dict[str, np.ndarray]] = {}
    for key, value in objectives.items():
        grads = torch.autograd.grad(value, parameters, retain_graph=True,
                                    allow_unused=True)
        collected: dict[str, list[np.ndarray]] = {g: [] for g in MODULE_GROUPS}
        collected["other"] = []
        for name, parameter, grad in zip(names, parameters, grads):
            # An unused parameter must contribute a zero vector of its own
            # length, not a placeholder: the per-objective vectors are compared
            # elementwise by cosine and have to stay aligned.
            vector = (np.zeros(parameter.numel(), dtype=np.float64) if grad is None
                      else grad.detach().float().reshape(-1).cpu().numpy())
            collected[module_of(name)].append(vector.astype(np.float64))
        flat[key] = {group: (np.concatenate(vectors) if vectors
                             else np.zeros(1)) for group, vectors in
                     collected.items()}

    report: dict[str, dict] = {"norms": {}, "cosine_level_shape": {}}
    for group in (*MODULE_GROUPS, "other"):
        for key in objectives:
            report["norms"].setdefault(key, {})[group] = float(
                np.linalg.norm(flat[key][group]))
        a, b = flat["level"][group], flat["shape"][group]
        denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
        report["cosine_level_shape"][group] = (
            float(a @ b / denominator) if denominator > 1e-12 else None)
    return report


def aggregate_gradient_reports(reports: list[dict]) -> dict:
    """Mean over episodes of one step's per-episode gradient geometry.

    Contract 2026-08-16: no single-episode snapshot may be reported as a
    mechanism statement. Every saved diagnostic step records the mean over all
    of its episodes plus the per-episode records. The module groups are read
    from the reports, so the helper is shared across architectures.
    """
    groups = sorted({group for report in reports
                     for key in report["norms"]
                     for group in report["norms"][key]})
    aggregate: dict[str, dict] = {"norms": {"level": {}, "shape": {}},
                                  "cosine_level_shape": {}}
    for key in ("level", "shape"):
        for group in groups:
            aggregate["norms"][key][group] = float(np.mean([
                report["norms"][key][group] for report in reports]))
    for group in groups:
        cosines = [report["cosine_level_shape"][group] for report in reports
                   if report["cosine_level_shape"][group] is not None]
        aggregate["cosine_level_shape"][group] = (
            float(np.mean(cosines)) if cosines else None)
    return aggregate


def aggregate_gradient_diagnostics(diagnostics: list[dict]) -> dict:
    """Summarize gradient geometry across all recorded steps of one seed.

    Reports the mean and conflict frequency (fraction of episode snapshots with
    a negative level-shape cosine) per module group, so reports can cite
    step-aggregated numbers instead of single-step snapshots. Seed aggregation
    across runs happens at the stage level over these summaries.
    """
    groups = sorted({group
                     for step in diagnostics
                     for group in step["mean_over_episodes"]
                     ["cosine_level_shape"]})
    summary: dict[str, dict] = {}
    for group in groups:
        cosines = [step["mean_over_episodes"]["cosine_level_shape"][group]
                   for step in diagnostics
                   if step["mean_over_episodes"]
                   ["cosine_level_shape"][group] is not None]
        per_episode = [entry["cosine_level_shape"][group]
                       for step in diagnostics for entry in step["per_episode"]
                       if entry["cosine_level_shape"][group] is not None]
        norms = {"level": [], "shape": []}
        for step in diagnostics:
            for key in ("level", "shape"):
                norms[key].append(
                    step["mean_over_episodes"]["norms"][key][group])
        summary[group] = {
            "steps_recorded": len(diagnostics),
            "episode_snapshots": len(per_episode),
            "conflict_frequency": (
                float(np.mean([value < 0 for value in per_episode]))
                if per_episode else None),
            "cosine_mean": (float(np.mean(cosines)) if cosines else None),
            "cosine_min": (float(np.min(cosines)) if cosines else None),
            "level_norm_mean": float(np.mean(norms["level"])),
            "shape_norm_mean": float(np.mean(norms["shape"])),
        }
    return summary


def gradient_coverage(model: LevelShapeModel, data: QPSMPData, rng,
                      scale: LabelScale, config: RouteConfig,
                      dtype: torch.dtype) -> dict:
    """Zero-gradient census of the final checkpoint on a fresh k=5 episode."""
    model.train()
    spec = data.draw_episode("meta_train", 5,
                             int(rng.integers(config.min_query_size,
                                              config.query_size + 1)),
                             rng, min_query_size=config.min_query_size)
    episode = compact_episode(normalized(data.materialize(spec), scale))
    support, atoms, bonds, mask = episode_tensors(episode, config.device, dtype)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, spec.target, config.device, dtype)
    channels = model.encode_ligand(atoms, bonds, mask)
    endpoint, prior, level, centered, _ = model.endpoint_with_channels(
        channels, pooled, tokens, protein_mask, chemistry)
    query_y = episode.query_y.to(config.device, dtype).unsqueeze(0)
    prediction = endpoint[:, support:]
    if support:
        similarity = tanimoto(
            episode.query_fingerprint.to(config.device, dtype).unsqueeze(0),
            episode.support_fingerprint.to(config.device, dtype).unsqueeze(0))
        residual = (episode.support_y.to(config.device, dtype).unsqueeze(0)
                    - endpoint[:, :support]).detach()
        prediction = prediction + model.transport.shrinkage(support, residual) \
            * model.transport(residual, similarity)[0]
    prediction.square().mean().backward()
    zero = [name for name, parameter in model.named_parameters()
            if parameter.requires_grad
            and (parameter.grad is None or not parameter.grad.abs().sum())]
    model.zero_grad(set_to_none=True)
    total = sum(parameter.numel() for parameter in model.parameters()
                if parameter.requires_grad)
    zero_params = sum(model.get_parameter(name).numel() for name in zero)
    return {"nonzero_fraction": float(1.0 - zero_params / max(total, 1)),
            "zero_gradient_parameters": zero}


# ---------------------------------------------------------------- evaluation


def evaluate(model: LevelShapeModel, data: QPSMPData, split: str,
             donors: dict[str, str], scale: LabelScale, config: RouteConfig,
             draws: int = 1) -> dict:
    device = config.device
    dtype = next(model.parameters()).dtype
    banks = data.fixed_nested_episode_banks(
        split, SUPPORT_SIZES, config.query_size, draws,
        config.evaluation_seed, None)
    rows: list[dict] = []
    model.eval()
    with torch.no_grad():
        for k, specs in banks.items():
            for spec in specs:
                episode = compact_episode(normalized(data.materialize(spec), scale))
                support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
                channels = model.encode_ligand(atoms, bonds, mask)
                endpoint, prior, level, centered, _ = model.endpoint_with_channels(
                    channels, *protein_inputs(data, spec.target, device, dtype))
                wrong_endpoint = model.endpoint_with_channels(
                    channels, *protein_inputs(
                        data, donors[spec.target], device, dtype))[0]
                query_y = episode.query_y.to(device, dtype).unsqueeze(0)
                query_endpoint = endpoint[:, support:]
                query_wrong = wrong_endpoint[:, support:]
                ligand_only = (prior[:, support:]
                               + level[:, support:])   # interaction removed
                transport = torch.zeros_like(query_endpoint)
                permuted = torch.zeros_like(query_endpoint)
                level_shift = torch.zeros_like(query_endpoint)
                if support:
                    similarity = tanimoto(
                        episode.query_fingerprint.to(device, dtype).unsqueeze(0),
                        episode.support_fingerprint.to(device, dtype).unsqueeze(0))
                    residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                                - endpoint[:, :support])
                    shrink = model.transport.shrinkage(support, residual)
                    transport = shrink * model.transport(residual, similarity)[0]
                    permuted = shrink * model.transport(
                        residual.roll(1, dims=-1), similarity)[0]
                    # Level-only: one constant per target, built from the
                    # support residuals alone — never from the query panel.
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
                    "wrong_protein": pk(query_wrong + transport),
                    # the interaction branch removed: level + ligand tendency
                    "ligand_only": pk(ligand_only),
                }
                error = predictions["full"] - truth
                row = {"k": k, "component": spec.component, "target": spec.target,
                       "calibration_pk": float(error.mean() ** 2),
                       "shape_pk": float(error.var()),
                       "centered_spread_pk": float(
                           centered[:, support:].squeeze(0).float().cpu().numpy().std()
                           * scale.scale),
                       "level_pk": float(level[0, 0]) * scale.scale + scale.mean,
                       # Activation statistics (contract 2026-08-16).
                       "centered_mean_pk": float(
                           centered[:, support:].squeeze(0).float().cpu().numpy().mean()
                           * scale.scale),
                       "centered_abs_mean_pk": float(
                           centered[:, support:].squeeze(0).float().cpu().numpy()
                           .__abs__().mean() * scale.scale),
                       "transport_abs_mean_pk": float(
                           transport.squeeze(0).float().cpu().numpy()
                           .__abs__().mean() * scale.scale),
                       "endpoint_spread_pk": float(
                           predictions["zero_shot"].std())}
                for name, value in predictions.items():
                    row[f"{name}_mse_pk"] = float(((value - truth) ** 2).mean())
                ci, comparable = concordance_index(predictions["full"], truth)
                row["ci"] = ci if comparable else None
                row["spearman"] = spearman(predictions["full"], truth)
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


# ---------------------------------------------------------------- training


def train(data: QPSMPData, config: RouteConfig, output: Path) -> dict:
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = LevelShapeModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_heads=config.pair_heads, anchors=config.anchors).to(config.device)
    dtype = next(model.parameters()).dtype
    scale = training_label_scale(data)
    # Two donor maps, per the 2026-08-16 contract: training counterfactuals use
    # meta_train donors (both proteins seen, contrast varies identity alone);
    # evaluation wrong-protein arms use same-split donors (both proteins
    # equally unseen). Whitening is fitted on meta_train in both maps.
    donors_train = matched_donors(data, "meta_train", donor_pool="meta_train")
    donors_eval = matched_donors(data, "meta_val", donor_pool="meta_val",
                                 whitening_pool="meta_train")
    parameters = [p for p in model.parameters() if p.requires_grad]
    optimizer = torch.optim.AdamW(parameters, lr=config.learning_rate,
                                  weight_decay=config.weight_decay)
    amp_enabled = config.amp and config.device.startswith("cuda")
    scaler = torch.amp.GradScaler("cuda", enabled=amp_enabled)

    best_state, best_value, best_step = None, float("inf"), 0
    trace, diagnostics = [], []
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
                    model, data, episode, donors_train, config, dtype)
                scaler.scale(loss / len(episodes)).backward()
                for key, value in parts.items():
                    accumulated[key] = accumulated.get(key, 0.0) + value / len(episodes)
        scaler.unscale_(optimizer)
        grad_norm = float(torch.nn.utils.clip_grad_norm_(
            model.parameters(), config.grad_clip))
        scaler.step(optimizer)
        scaler.update()
        trace.append({"step": step, "grad_norm": grad_norm, **accumulated})

        if step % config.diagnostic_interval == 0 or step == config.steps:
            # Gradient geometry is recorded for *every* episode of the step and
            # averaged, so no report may present a single-episode snapshot as a
            # general claim (contract 2026-08-16). Seed and step aggregation
            # happens at the stage level over these records.
            per_episode = [gradient_report(
                model, data, episode, donors_train, config, dtype)
                for episode in episodes]
            diagnostics.append({
                "step": step, "episode_count": len(per_episode),
                "mean_over_episodes": aggregate_gradient_reports(per_episode),
                "per_episode": per_episode})
        if step % config.val_interval == 0 or step == config.steps:
            rows = evaluate(model, data, "meta_val", donors_eval, scale,
                            config)["rows"]
            value = float(np.mean([component_target_mean(rows, "full_mse_pk", k)
                                   for k in SUPPORT_SIZES]))
            record = {"step": step, "val_mean_mse_pk": value,
                      "val_k0_mse_pk": component_target_mean(rows, "full_mse_pk", 0),
                      "val_k0_ci": component_target_mean(rows, "zero_shot_ci", 0),
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
    coverage = gradient_coverage(model, data, rng, scale, config, dtype)
    return {"model": model, "scale": scale,
            "donors_train": donors_train, "donors_eval": donors_eval,
            "best_val_mean_mse_pk": best_value, "best_step": best_step,
            "loss_trace": trace[-50:], "gradient_diagnostics": diagnostics,
            "gradient_summary": aggregate_gradient_diagnostics(diagnostics),
            "gradient_coverage": coverage,
            "trainable_parameters": sum(p.numel() for p in parameters),
            "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                    if config.device.startswith("cuda") else 0.0),
            "wall_seconds": time.monotonic() - started}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=RouteConfig.seed)
    parser.add_argument("--steps", type=int, default=RouteConfig.steps)
    parser.add_argument("--episodes-per-step", type=int,
                        default=RouteConfig.episodes_per_step)
    parser.add_argument("--learning-rate", type=float,
                        default=RouteConfig.learning_rate)
    parser.add_argument("--anchors", type=int, default=RouteConfig.anchors)
    parser.add_argument("--no-routing", action="store_true")
    parser.add_argument("--no-counterfactual", action="store_true")
    parser.add_argument("--val-interval", type=int, default=RouteConfig.val_interval)
    parser.add_argument("--device", default=RouteConfig.device)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)

    config = RouteConfig(
        seed=args.seed, steps=args.steps,
        episodes_per_step=args.episodes_per_step,
        learning_rate=args.learning_rate, anchors=args.anchors,
        routing=not args.no_routing, counterfactual=not args.no_counterfactual,
        val_interval=args.val_interval, device=args.device)
    # The sealed confirmation split is dropped physically: this script cannot
    # read meta_test (contract 2026-08-16).
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=False)
    result = train(data, config, args.output)
    model, scale = result.pop("model"), result["scale"]
    donors_train = result.pop("donors_train")
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
        "ci", "spearman", "zero_shot_ci", "calibration_pk", "shape_pk",
        "centered_spread_pk", "centered_mean_pk", "centered_abs_mean_pk",
        "transport_abs_mean_pk", "endpoint_spread_pk", "level_pk"]
    summary = {str(k): {f: component_target_mean(rows, f, k) for f in fields}
               for k in SUPPORT_SIZES}
    payload = {
        "schema": "MetaSieve.LevelShapeRoutedTraining.v2",
        "config": asdict(config),
        "split_directory": str(args.split_directory),
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "arm": {"routing": config.routing, "counterfactual": config.counterfactual},
        "donors": {
            "training_counterfactual_pool": "meta_train",
            "evaluation_wrong_protein_pool": "meta_val",
            "whitening_pool": "meta_train",
            "metric": "esm_whitened (train-fitted)",
            "criterion": "most similar target from a different homology "
                         "component",
        },
        "meta_test": {
            "included": False,
            "evaluated": False,
            "seal": "physical: QPSMPData include_meta_test=False",
        },
        "checkpoint_sha256": file_sha256(checkpoint_path),
        "training": {k: v for k, v in result.items()
                     if k not in {"scale", "donors_train", "donors_eval"}},
        "label_scale": asdict(result["scale"]),
        "meta_val": summary,
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("\n%-3s %9s %9s %9s %9s %8s %8s" % (
        "k", "full", "zero", "ligonly", "wrongP", "CI", "rho"))
    for k in SUPPORT_SIZES:
        entry = summary[str(k)]
        print("%-3d %9.4f %9.4f %9.4f %9.4f %8.4f %8.4f" % (
            k, entry["full_mse_pk"], entry["zero_shot_mse_pk"],
            entry["ligand_only_mse_pk"], entry["wrong_protein_mse_pk"],
            entry["ci"], entry["spearman"]))


if __name__ == "__main__":
    main()
