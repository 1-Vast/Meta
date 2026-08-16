"""Core Innovation B: counterfactual gradient-routed shape-first training.

The squared error decomposes exactly into a level (calibration) term and a
shape term; ordinary training lets one scalar head buy almost all of its error
from the level, and Stage R3/R4 measured the two objectives pulling the shared
trunk in opposite directions (gradient cosine -0.334). Routing removed the
conflict and bought calibration — and left shape flat, because the shape
objective was still a variance, whose optimum a weak interaction branch can
partially ignore.

This method changes what the shape signal *is*, and routes every term only to
the modules that own it:

* **shape-first**: the shape objective is within-target pairwise ranking
  (RankNet-style logistic, ActFound/PBCNet-style shift invariance) plus direct
  relative supervision on the antisymmetric potential
  `delta(P, i, j) ~ y_i - y_j` over every in-target pair (support-support,
  support-query, query-query). A per-target constant has zero gradient under
  both, so the interaction trunk is forced to produce ligand-dependent
  ordering — exactly the quantity the zero-shot endpoint has been missing;
* **cliff-aware pair weighting**: in-target pairs with Tanimoto >= 0.6 and a
  >= 1.0 pK label gap — activity cliffs, the failure mode Stage R0 localized —
  are upweighted in both the ranking and the relative supervision, so small
  chemical edits that reverse activity dominate the shape gradient;
* **routing**: the level term cannot train the interaction trunk; the shape
  terms cannot train the level head (a constant is analytically invisible to
  them, so detaching it there is exact rather than approximate);
* **three counterfactual contrasts, each routed so it cannot be satisfied by
  the wrong module**: correct protein must beat a similarity-matched wrong
  protein on the shape objective (ligand prior detached) and on the level
  objective (interaction detached); correct support binding must beat a
  permuted / magnitude-matched-wrong label assignment *and* a wrong support
  ligand (endpoint detached), so only the transport can respond.

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

from model.reltransport import RelTransportModel
from model.similarity_grammar import tanimoto
from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import EpisodeBatch, EpisodeSpec, QPSMPData, stable_seed
from scripts.train_level_shape import (
    aggregate_gradient_diagnostics, aggregate_gradient_reports,
    matched_donors, normalized, protein_inputs,
)
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, LabelScale,
    compact_episode, file_sha256, training_label_scale,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)
MODULE_GROUPS = ("ligand_encoder", "protein_encoder", "grammar", "embed",
                 "relative", "ligand_head", "level_head", "anchor", "transport")


@dataclass(frozen=True)
class RelConfig:
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
    rank: int = 96
    ranking_loss_weight: float = 1.0      # primary shape signal
    ranking_loss: str = "ranknet"         # "ranknet" | "margin"
    ranking_margin: float = 0.1           # hinge margin (normalized pK)
    shape_variance_weight: float = 1.0
    relative_loss_weight: float = 0.5     # delta ~ y_i - y_j supervision
    cliff_pair_weight: float = 4.0
    cliff_tanimoto: float = 0.6
    cliff_gap_pk: float = 1.0
    counterfactual_weight: float = 0.25
    contrast_temperature: float = 0.1
    identify_weight: float = 0.3        # label-free pin on the shape mean
    routing: bool = True
    counterfactual: bool = True
    gate: bool = True                    # False = Stage 2 ablation A3
    ordinary: bool = False               # True = Stage 2 arm A1: the R3R4
                                         # "ordinary" recipe (MSE + 0.5 ranking),
                                         # no relative supervision, no routing,
                                         # no counterfactual
    val_interval: int = 100
    diagnostic_interval: int = 200
    lr_warmup_fraction: float = 0.05
    lr_final_fraction: float = 0.1
    amp: bool = True
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


# ---------------------------------------------------------------- objectives


def level_term(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    return (prediction - truth).mean(-1).square().mean()


def shape_variance(prediction: torch.Tensor, truth: torch.Tensor) -> torch.Tensor:
    error = prediction - truth
    return (error - error.mean(-1, keepdim=True)).square().mean()


def pairwise_ranking(prediction: torch.Tensor, truth: torch.Tensor,
                     temperature: float,
                     pair_weight: torch.Tensor | None = None,
                     mode: str = "ranknet",
                     margin: float = 0.1) -> torch.Tensor:
    """Pairwise ranking loss with optional per-pair weights.

    `mode == "ranknet"`: RankNet logistic (softplus) loss.
    `mode == "margin"`: margin-ranking (hinge) loss
    `max(0, margin - sign(dy) * dp)` — it keeps pushing correctly ordered
    pairs until their predicted margin exceeds `margin`, directly countering
    the margin compression the R9 pair audit measured (C1 mean |margin|
    0.097 vs A0 0.121). Literature: margin-based ranking in metric learning
    (Hadsell et al. 2006) and the pairwise-loss/global-epistasis analysis of
    [Diaz et al., arXiv:2305.03136], which frames pairwise losses over
    biological fitness/affinity exactly as a level/shape decomposition.
    """
    delta_y = truth.unsqueeze(-1) - truth.unsqueeze(-2)
    delta_p = prediction.unsqueeze(-1) - prediction.unsqueeze(-2)
    comparable = delta_y != 0
    if not bool(comparable.any()):
        return prediction.new_zeros(())
    signed = delta_y.sign() * delta_p
    if mode == "margin":
        loss = torch.relu(margin - signed[comparable])
    else:
        loss = F.softplus(-signed[comparable] / temperature)
    if pair_weight is not None and pair_weight.numel() == loss.numel():
        loss = loss * pair_weight[comparable]
    return loss.mean()


def contrast(correct: torch.Tensor, wrong: torch.Tensor,
             temperature: float) -> torch.Tensor:
    """Prefer the correct configuration over its counterfactual."""
    logits = -torch.stack((correct, wrong)) / temperature
    return F.cross_entropy(logits.unsqueeze(0),
                           logits.new_zeros(1, dtype=torch.long))


def cliff_pair_weights(similarity: torch.Tensor, truth: torch.Tensor,
                       cliff_tanimoto: float, cliff_gap: float,
                       cliff_weight: float) -> torch.Tensor:
    """[B,Q,Q] weights, >1 exactly on activity-cliff pairs."""
    delta_y = (truth.unsqueeze(-1) - truth.unsqueeze(-2)).abs()
    cliff = (similarity >= cliff_tanimoto) & (delta_y >= cliff_gap)
    return 1.0 + (cliff_weight - 1.0) * cliff.to(similarity.dtype)


def learning_rate_factor(step: int, config: RelConfig) -> float:
    warmup = max(1, int(config.steps * config.lr_warmup_fraction))
    if step <= warmup:
        return step / warmup
    progress = (step - warmup) / max(1, config.steps - warmup)
    final = config.lr_final_fraction
    return final + (1.0 - final) * 0.5 * (1.0 + np.cos(np.pi * min(progress, 1.0)))


# ---------------------------------------------------------------- data helpers


def episode_tensors(episode: EpisodeBatch, device: str, dtype: torch.dtype):
    support = episode.support_atoms.shape[0]
    atoms = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)
    return (support, atoms.to(device, dtype).unsqueeze(0),
            bonds.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0))


def donor_ligand_indices(data: QPSMPData, spec: EpisodeSpec,
                         count: int) -> tuple[int, ...]:
    """Deterministic wrong-support-ligand draw from a cross-component target."""
    donor_cells = [int(index) for index in data.tasks[spec.split][spec.donor_target]]
    seed = stable_seed("reltransport-donor-ligand", spec.split, spec.target,
                       spec.support, spec.query)
    order = np.random.default_rng(seed).permutation(donor_cells)
    return tuple(map(int, order[:count]))


def donor_ligand_graphs(data: QPSMPData, spec: EpisodeSpec,
                        device: str, dtype: torch.dtype) -> torch.Tensor:
    indices = donor_ligand_indices(data, spec, len(spec.support))
    values = [data.ligand_bank.get(data.cells[index]["ligand_id"])
              for index in indices]
    max_atoms = max(value[0].shape[0] for value in values)
    atoms = torch.stack([torch.nn.functional.pad(
        torch.from_numpy(value[0]), (0, 0, 0, max_atoms - value[0].shape[0]))
        for value in values])
    bonds = torch.stack([torch.nn.functional.pad(
        torch.from_numpy(value[1]),
        (0, 0, 0, max_atoms - value[1].shape[0],
         0, max_atoms - value[1].shape[0])) for value in values])
    mask = torch.stack([torch.nn.functional.pad(
        torch.from_numpy(value[2]),
        (0, max_atoms - value[2].shape[0])) for value in values])
    return (atoms.to(device, dtype).unsqueeze(0),
            bonds.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0))


# ---------------------------------------------------------------- per-episode loss


def transport_block(model: RelTransportModel, support: int, query_embed,
                    support_embed, u, u_gate, residual, similarity, *,
                    relative_on: bool = True) -> tuple:
    """Query-specific residual transport for one (query, support, protein).

    rho is the linear zero-initialised transferability gate; `relative_on`
    False pins rho == 1 (Stage 2 ablation A3).
    """
    delta = model.relative.delta_matrix(u, query_embed, support_embed)
    if relative_on:
        rho = model.gate.gate_matrix(u_gate, query_embed, support_embed)
    else:
        rho = torch.ones_like(delta)
    shrink = model.transport.shrinkage(support, residual)
    transport, weight = model.transport(
        support_embed, query_embed, residual, similarity, rho)
    return shrink * transport, delta


def shape_objective(model: RelTransportModel, config: RelConfig,
                    prediction: torch.Tensor, truth: torch.Tensor,
                    scale: LabelScale,
                    similarity: torch.Tensor | None = None) -> torch.Tensor:
    """Ranking-primary shape objective with cliff-aware pair weighting."""
    pair_weight = None
    if similarity is not None:
        cliff_gap = config.cliff_gap_pk / scale.scale
        pair_weight = cliff_pair_weights(
            similarity, truth, config.cliff_tanimoto, cliff_gap,
            config.cliff_pair_weight)
    rank = pairwise_ranking(prediction, truth, 1.0, pair_weight,
                            mode=config.ranking_loss,
                            margin=config.ranking_margin)
    variance = shape_variance(prediction, truth)
    return config.ranking_loss_weight * rank \
        + config.shape_variance_weight * variance


def relative_supervision(model: RelTransportModel, config: RelConfig,
                         parts, full_y: torch.Tensor, scale: LabelScale,
                         full_similarity: torch.Tensor | None = None
                         ) -> torch.Tensor:
    """delta(P, i, j) ~ y_i - y_j over every in-target pair, cliff-weighted."""
    _, _, _, _, _, embed, u, _ = parts
    delta = model.relative.delta_matrix(u, embed, embed)
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


def episode_loss(model: RelTransportModel, data: QPSMPData,
                 episode: EpisodeBatch, donors: dict[str, str],
                 config: RelConfig, dtype: torch.dtype,
                 scale: LabelScale) -> tuple[torch.Tensor, dict]:
    device = config.device
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, episode.spec.target, device, dtype)
    donor = donors[episode.spec.target]
    donor_inputs = protein_inputs(data, donor, device, dtype)
    parts = model.forward_parts(pooled, tokens, protein_mask, atoms, bonds,
                                mask, chemistry)
    donor_parts = model.forward_parts(
        *donor_inputs[:3], atoms, bonds, mask, donor_inputs[3])
    endpoint, prior, level, shape, _, embed, u, _ = parts
    query_y = episode.query_y.to(device, dtype).unsqueeze(0)
    support_y = episode.support_y.to(device, dtype).unsqueeze(0)
    full_y = torch.cat((support_y, query_y), -1)

    def split(value):
        return value[:, :support], value[:, support:]

    _, query_endpoint = split(endpoint)
    _, query_prior = split(prior)
    _, query_level = split(level)
    _, query_shape = split(shape)
    _, query_embed = split(embed)

    query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
    support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
    query_similarity = tanimoto(query_fp, query_fp)
    full_similarity = tanimoto(
        torch.cat((support_fp, query_fp), 1),
        torch.cat((support_fp, query_fp), 1))

    transport = torch.zeros_like(query_endpoint)
    permuted_transport = None
    wrong_ligand_transport = None
    wrong_protein_transport = None
    if support:
        support_embed = embed[:, :support]
        similarity = tanimoto(query_fp, support_fp)
        residual = (support_y - endpoint[:, :support]).detach()
        u_gate = parts[7]
        transport, _ = transport_block(
            model, support, query_embed, support_embed, u, u_gate, residual,
            similarity, relative_on=config.gate)
        if config.counterfactual:
            # k>=2: permuted labels (mean(r) invariant -> isolates the
            # query-specific channel). k=1 has no permutation; the
            # magnitude-matched label flip is deliberately NOT trained: under
            # a query-specific gate it teaches per-episode error ratios and
            # destabilizes the transport (measured on the synthetic gates).
            # The k=1 wrong-label behaviour is evaluated, not trained.
            if support > 1:
                permuted_transport, _ = transport_block(
                    model, support, query_embed, support_embed, u, u_gate,
                    residual.roll(1, dims=-1), similarity,
                    relative_on=config.gate)
            donor_embed = donor_parts[5]
            wrong_protein_transport, _ = transport_block(
                model, support, query_embed, donor_embed[:, :support],
                donor_parts[6], donor_parts[7], residual, similarity,
                relative_on=config.gate)
            wrong_atoms, wrong_bonds, wrong_mask = donor_ligand_graphs(
                data, episode.spec, device, dtype)
            wrong_parts = model.forward_parts(
                pooled, tokens, protein_mask, wrong_atoms, wrong_bonds,
                wrong_mask, chemistry)
            wrong_residual = (support_y - wrong_parts[0]).detach()
            wrong_ligand_transport, _ = transport_block(
                model, support, query_embed, wrong_parts[5], u, u_gate,
                wrong_residual, similarity, relative_on=config.gate)

    prediction = query_endpoint + transport
    if config.ordinary:
        # Stage 2 arm A1: the conventional recipe applied to the same
        # architecture — squared error (level + shape) with a 0.5-weighted
        # pairwise ranking auxiliary on the full prediction. No relative
        # supervision, no routing, no counterfactuals. This is the identical
        # "ordinary training" used as the B1 arm of the Stage R3/R4 ladder.
        loss = level_term(prediction, query_y) \
            + shape_variance(prediction, query_y) \
            + 0.5 * pairwise_ranking(prediction, query_y, 1.0)
        return loss, {"level": float(
            level_term(prediction, query_y).detach()),
            "shape": float(shape_variance(prediction, query_y).detach()),
            "prediction_mse": float(
                (prediction - query_y).square().mean().detach())}
    if config.routing:
        p_level = query_prior + query_level + query_shape.detach() \
            + transport.detach()
        p_shape = query_prior + query_level.detach() + query_shape + transport
    else:
        p_level = p_shape = prediction
    loss_level = level_term(p_level, query_y)
    loss_shape = shape_objective(model, config, p_shape, query_y, scale,
                                 query_similarity)
    loss_rel = config.relative_loss_weight * relative_supervision(
        model, config, parts, full_y, scale, full_similarity)
    loss = loss_level + loss_shape + loss_rel

    parts_names = {"level": float(loss_level.detach()),
                   "shape": float(loss_shape.detach()),
                   "rel": float(loss_rel.detach())}

    if config.routing:
        # Routing leaves one identifiability freedom: the shape branch's
        # per-target mean is detached out of the level term and invisible to
        # variance/ranking (shift-invariant), so the anchors can drift
        # outward while `target_level` silently compensates — the Stage R3/R4
        # predrift defect, reproduced by the R6 screening seed 1 (ligand_only
        # exploded to ~100 MSE^2). This label-free pin closes the freedom; it
        # contains no label, so it cannot buy fit. A small weight suffices to
        # stop the null-direction walk without fighting calibration.
        loss_identify = config.identify_weight * query_shape.mean(-1).square().mean()
        loss = loss + loss_identify
        parts_names["identify"] = float(loss_identify.detach())

    if config.counterfactual:
        # protein-shape: ligand prior detached, so no ligand-only shortcut.
        wrong_prior = donor_parts[1][:, support:]
        wrong_shape = donor_parts[3][:, support:]
        correct_p_shape = query_prior.detach() + query_shape + transport
        wrong_p_shape = wrong_prior.detach() + wrong_shape \
            + (wrong_protein_transport if wrong_protein_transport is not None
               else transport.detach())
        loss_protein_shape = contrast(
            shape_objective(model, config, correct_p_shape, query_y, scale,
                            query_similarity),
            shape_objective(model, config, wrong_p_shape, query_y, scale,
                            query_similarity),
            config.contrast_temperature)
        # protein-level: interaction detached, so only target_level responds.
        wrong_level = donor_parts[2][:, support:]
        loss_protein_level = contrast(
            level_term(query_prior.detach() + query_level
                       + query_shape.detach() + transport.detach(), query_y),
            level_term(query_prior.detach() + wrong_level
                       + query_shape.detach() + transport.detach(), query_y),
            config.contrast_temperature)
        loss = loss + config.counterfactual_weight * (
            loss_protein_shape + loss_protein_level)
        parts_names["protein_shape"] = float(loss_protein_shape.detach())
        parts_names["protein_level"] = float(loss_protein_level.detach())

        # support binding: endpoint detached, so only the transport responds.
        # The contrast target is the full squared error (deployment metric):
        # the shape-first principle governs the interaction trunk elsewhere,
        # but support utility is decided by absolute prediction error. At
        # k>=2 a label permutation leaves mean(r) unchanged, so this contrast
        # isolates the query-specific channel exactly as in Stage R3/R4.
        frozen = query_endpoint.detach()
        correct_mse = ((frozen + transport) - query_y).square().mean()
        loss_binding = correct_mse.new_zeros(())
        if permuted_transport is not None:
            wrong_mse = ((frozen + permuted_transport) - query_y).square().mean()
            loss_binding = loss_binding + contrast(
                correct_mse, wrong_mse, config.contrast_temperature)
        if wrong_ligand_transport is not None:
            wrong_ligand_mse = ((frozen + wrong_ligand_transport)
                                - query_y).square().mean()
            loss_binding = loss_binding + contrast(
                correct_mse, wrong_ligand_mse, config.contrast_temperature)
        loss = loss + config.counterfactual_weight * loss_binding
        parts_names["binding"] = float(loss_binding.detach())

    parts_names["prediction_mse"] = float(
        (prediction - query_y).square().mean().detach())
    return loss, parts_names


# ---------------------------------------------------------------- diagnostics


def module_of(name: str) -> str:
    for group in MODULE_GROUPS:
        if name.startswith(group):
            return group
    return "other"


def gradient_report(model: RelTransportModel, data: QPSMPData,
                    episode: EpisodeBatch, donors: dict[str, str],
                    config: RelConfig, dtype: torch.dtype,
                    scale: LabelScale) -> dict:
    """Per-objective, per-module gradient norms and pairwise cosines."""
    device = config.device
    support, atoms, bonds, mask = episode_tensors(episode, device, dtype)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, episode.spec.target, device, dtype)
    parts = model.forward_parts(pooled, tokens, protein_mask, atoms, bonds,
                                mask, chemistry)
    endpoint, prior, level, shape, _, _, _, _ = parts
    query_y = episode.query_y.to(device, dtype).unsqueeze(0)
    query_prior, query_level = prior[:, support:], level[:, support:]
    query_shape = shape[:, support:]
    del endpoint

    objectives = {}
    if config.routing and not config.ordinary:
        objectives["level"] = level_term(
            query_prior + query_level + query_shape.detach(), query_y)
        objectives["shape"] = shape_objective(
            model, config, query_prior + query_level.detach() + query_shape,
            query_y, scale)
    else:
        joint = query_prior + query_level + query_shape
        objectives["level"] = level_term(joint, query_y)
        if config.ordinary:
            objectives["shape"] = shape_variance(joint, query_y) \
                + 0.5 * pairwise_ranking(joint, query_y, 1.0)
        else:
            objectives["shape"] = shape_objective(
                model, config, joint, query_y, scale)

    names = [n for n, p in model.named_parameters() if p.requires_grad]
    parameters = [p for _, p in model.named_parameters() if p.requires_grad]
    flat: dict[str, dict[str, np.ndarray]] = {}
    for key, value in objectives.items():
        grads = torch.autograd.grad(value, parameters, retain_graph=True,
                                    allow_unused=True)
        collected: dict[str, list[np.ndarray]] = {g: [] for g in MODULE_GROUPS}
        collected["other"] = []
        for name, parameter, grad in zip(names, parameters, grads):
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


def gradient_coverage(model: RelTransportModel, data: QPSMPData, rng,
                      scale: LabelScale, config: RelConfig,
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
    parts = model.forward_parts(pooled, tokens, protein_mask, atoms, bonds,
                                mask, chemistry)
    query_y = episode.query_y.to(config.device, dtype).unsqueeze(0)
    support_y = episode.support_y.to(config.device, dtype).unsqueeze(0)
    _, query_endpoint = torch.split(parts[0], (support, parts[0].shape[1] - support), 1)
    query_fp = episode.query_fingerprint.to(config.device, dtype).unsqueeze(0)
    support_fp = episode.support_fingerprint.to(config.device, dtype).unsqueeze(0)
    similarity = tanimoto(query_fp, support_fp)
    embed = parts[5]
    support_embed, query_embed = torch.split(embed, (support, embed.shape[1] - support), 1)
    residual = (support_y - parts[0][:, :support]).detach()
    delta = model.relative.delta_matrix(parts[6], query_embed, support_embed)
    rho = model.gate.gate_matrix(parts[7], query_embed, support_embed)
    shrink = model.transport.shrinkage(support, residual)
    transport, _ = model.transport(support_embed, query_embed, residual,
                                   similarity, rho)
    prediction = query_endpoint + shrink * transport
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


def sign_accuracy(prediction: np.ndarray, truth: np.ndarray) -> float:
    comparable, correct = 0, 0
    for left in range(len(truth)):
        for right in range(left + 1, len(truth)):
            true_delta = float(truth[left] - truth[right])
            if true_delta == 0:
                continue
            comparable += 1
            pred_delta = float(prediction[left] - prediction[right])
            correct += float(pred_delta * true_delta > 0)
    return float(correct / comparable) if comparable else float("nan")


def evaluate(model: RelTransportModel, data: QPSMPData, split: str,
             donors: dict[str, str], scale: LabelScale, config: RelConfig,
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
                pooled, tokens, protein_mask, chemistry = protein_inputs(
                    data, spec.target, device, dtype)
                donor_inputs = protein_inputs(data, donors[spec.target],
                                              device, dtype)
                parts = model.forward_parts(
                    pooled, tokens, protein_mask, atoms, bonds, mask, chemistry)
                donor_parts = model.forward_parts(
                    *donor_inputs[:3], atoms, bonds, mask, donor_inputs[3])
                endpoint, prior, level, shape, _, embed, u, _ = parts
                query_y = episode.query_y.to(device, dtype).unsqueeze(0)
                query_endpoint = endpoint[:, support:]
                wrong_endpoint = donor_parts[0][:, support:]
                ligand_only = prior[:, support:] + level[:, support:]
                transport = torch.zeros_like(query_endpoint)
                permuted = torch.zeros_like(query_endpoint)
                nogate = torch.zeros_like(query_endpoint)
                level_shift = torch.zeros_like(query_endpoint)
                if support:
                    query_fp = episode.query_fingerprint.to(device, dtype).unsqueeze(0)
                    support_fp = episode.support_fingerprint.to(device, dtype).unsqueeze(0)
                    similarity = tanimoto(query_fp, support_fp)
                    residual = (episode.support_y.to(device, dtype).unsqueeze(0)
                                - endpoint[:, :support])
                    shrink = model.transport.shrinkage(support, residual)
                    support_embed = embed[:, :support]
                    query_embed = embed[:, support:]
                    u_gate = parts[7]
                    transport, _ = transport_block(
                        model, support, query_embed, support_embed, u, u_gate,
                        residual, similarity, relative_on=config.gate)
                    transport = shrink * transport
                    nogate, _ = transport_block(
                        model, support, query_embed, support_embed, u, u_gate,
                        residual, similarity, relative_on=False)
                    nogate = shrink * nogate
                    permuted_residual = (residual.roll(1, dims=-1) if support > 1
                                         else -residual)
                    permuted, _ = transport_block(
                        model, support, query_embed, support_embed, u, u_gate,
                        permuted_residual, similarity)
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
                    "nogate": pk(query_endpoint + nogate),
                }
                error = predictions["full"] - truth
                row = {"k": k, "component": spec.component, "target": spec.target,
                       "calibration_pk": float(error.mean() ** 2),
                       "shape_pk": float(error.var()),
                       "endpoint_spread_pk": float(predictions["zero_shot"].std()),
                       "shape_abs_mean_pk": float(
                           shape[:, support:].squeeze(0).float().cpu().numpy()
                           .__abs__().mean() * scale.scale),
                       "transport_abs_mean_pk": float(
                           transport.squeeze(0).float().cpu().numpy()
                           .__abs__().mean() * scale.scale),
                       "level_pk": float(level[0, 0]) * scale.scale + scale.mean}
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
                row["zero_shot_sign_accuracy"] = sign_accuracy(
                    predictions["zero_shot"], truth)
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


def train(data: QPSMPData, config: RelConfig, output: Path) -> dict:
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = RelTransportModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_latents=config.pair_latents, pair_heads=config.pair_heads,
        anchors=config.anchors, rank=config.rank).to(config.device)
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

        if step % config.diagnostic_interval == 0 or step == config.steps:
            per_episode = [gradient_report(
                model, data, episode, donors_train, config, dtype, scale)
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
                      "val_k0_sign": component_target_mean(
                          rows, "zero_shot_sign_accuracy", 0),
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
    parser.add_argument("--seed", type=int, default=RelConfig.seed)
    parser.add_argument("--steps", type=int, default=RelConfig.steps)
    parser.add_argument("--episodes-per-step", type=int,
                        default=RelConfig.episodes_per_step)
    parser.add_argument("--learning-rate", type=float,
                        default=RelConfig.learning_rate)
    parser.add_argument("--anchors", type=int, default=RelConfig.anchors)
    parser.add_argument("--rank", type=int, default=RelConfig.rank)
    parser.add_argument("--ranking-loss-weight", type=float,
                        default=RelConfig.ranking_loss_weight)
    parser.add_argument("--ranking-loss", default=RelConfig.ranking_loss,
                        choices=("ranknet", "margin"))
    parser.add_argument("--ranking-margin", type=float,
                        default=RelConfig.ranking_margin)
    parser.add_argument("--shape-variance-weight", type=float,
                        default=RelConfig.shape_variance_weight)
    parser.add_argument("--relative-loss-weight", type=float,
                        default=RelConfig.relative_loss_weight)
    parser.add_argument("--cliff-pair-weight", type=float,
                        default=RelConfig.cliff_pair_weight)
    parser.add_argument("--counterfactual-weight", type=float,
                        default=RelConfig.counterfactual_weight)
    parser.add_argument("--no-routing", action="store_true")
    parser.add_argument("--no-counterfactual", action="store_true")
    parser.add_argument("--no-gate", action="store_true",
                        help="disable the delta terms in the transport "
                             "correction (Stage 2 ablation A3)")
    parser.add_argument("--ordinary", action="store_true",
                        help="conventional training of the same architecture "
                             "(Stage 2 arm A1): MSE + 0.5 ranking, no relative "
                             "supervision, no routing, no counterfactuals")
    parser.add_argument("--val-interval", type=int, default=RelConfig.val_interval)
    parser.add_argument("--device", default=RelConfig.device)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)

    config = RelConfig(
        seed=args.seed, steps=args.steps,
        episodes_per_step=args.episodes_per_step,
        learning_rate=args.learning_rate, anchors=args.anchors,
        rank=args.rank, ranking_loss_weight=args.ranking_loss_weight,
        ranking_loss=args.ranking_loss, ranking_margin=args.ranking_margin,
        shape_variance_weight=args.shape_variance_weight,
        relative_loss_weight=args.relative_loss_weight,
        cliff_pair_weight=args.cliff_pair_weight,
        counterfactual_weight=args.counterfactual_weight,
        routing=not args.no_routing, counterfactual=not args.no_counterfactual,
        gate=not args.no_gate,
        ordinary=args.ordinary,
        val_interval=args.val_interval, device=args.device)
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
        "ci", "spearman", "sign_accuracy", "zero_shot_ci",
        "zero_shot_sign_accuracy", "calibration_pk", "shape_pk",
        "endpoint_spread_pk", "shape_abs_mean_pk", "transport_abs_mean_pk",
        "level_pk"]
    summary = {str(k): {f: component_target_mean(rows, f, k) for f in fields}
               for k in SUPPORT_SIZES}
    payload = {
        "schema": "MetaSieve.RelTransportShapeFirstTraining.v1",
        "config": asdict(config),
        "split_directory": str(args.split_directory),
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "arm": {"routing": config.routing, "counterfactual": config.counterfactual,
                "gate": config.gate, "ordinary": config.ordinary},
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
    print("\n%-3s %9s %9s %9s %9s %8s %8s %8s" % (
        "k", "full", "zero", "ligonly", "wrongP", "CI", "rho", "sign"))
    for k in SUPPORT_SIZES:
        entry = summary[str(k)]
        print("%-3d %9.4f %9.4f %9.4f %9.4f %8.4f %8.4f %8.4f" % (
            k, entry["full_mse_pk"], entry["zero_shot_mse_pk"],
            entry["ligand_only_mse_pk"], entry["wrong_protein_mse_pk"],
            entry["ci"], entry["spearman"], entry["sign_accuracy"]))


if __name__ == "__main__":
    main()
