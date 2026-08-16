"""GPU episodic training for the active QPSMP-BPSF meta-learner.

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
import torch.nn.functional as F

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from model.interaction_grammar import InteractionGrammarModel
from model.qpsmp_meta import QPSMPBioModel
from model.relative_grammar import RelativeGrammarModel
from model.locality_grammar import LocalityGrammarModel
from model.similarity_grammar import SimilarityGrammarModel
from scripts.qpsmp_data import EpisodeBatch, EpisodeSpec, QPSMPData


def _relative_no_reliability(*args, **kwargs):
    return RelativeGrammarModel(*args, use_reliability=False, **kwargs)


def _similarity_only(*args, **kwargs):
    return SimilarityGrammarModel(*args, use_learned_key=False, **kwargs)


def _locality(*args, **kwargs):
    return LocalityGrammarModel(*args, use_learned_key=False, **kwargs)


ARCHITECTURES = {
    "bpsf": QPSMPBioModel,
    "grammar": InteractionGrammarModel,
    "relative": RelativeGrammarModel,
    "relative_noreliability": _relative_no_reliability,
    "similarity": SimilarityGrammarModel,
    "similarity_only": _similarity_only,
    "locality": _locality,
}
RELATIVE_ARCHITECTURES = {"relative", "relative_noreliability"}
FINGERPRINT_ARCHITECTURES = {"similarity", "similarity_only", "locality"}


def resolve_architecture(name: str):
    if name not in ARCHITECTURES:
        raise ValueError(f"unknown architecture: {name}")
    return ARCHITECTURES[name]


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
COMPACT_LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact"
OUT = ROOT / "report/meta_fewshot/qpsmp_meta_smoke"


@dataclass(frozen=True)
class TrainConfig:
    seed: int = 20260812
    evaluation_seed: int = 73101
    support_size: int = 5
    query_size: int = 20
    min_query_size: int = 4
    hidden_dim: int = 192
    task_dim: int = 48
    adapter_rank: int = 4
    adaptive_blocks: int = 2
    adapter_scale: float = 0.25
    ligand_layers: int = 4
    steps: int = 80
    episodes_per_step: int = 3
    learning_rate: float = 3e-4
    backbone_lr_scale: float = 0.25
    weight_decay: float = 1e-5
    val_interval: int = 10
    val_draws_per_target: int = 1
    test_draws_per_target: int = 1
    eval_targets_per_component: int = 1
    # Checkpoint selection needs more episodes than the 6-episode protocol test
    # bank; a separate width keeps the frozen test bank untouched.
    val_targets_per_component: int | None = None
    grad_clip: float = 1.0
    zero_shot_loss_weight: float = 1.0
    ranking_loss_weight: float = 0.5
    shape_loss_weight: float = 0.5
    support_match_loss_weight: float = 0.05
    binding_loss_weight: float = 0.25
    binding_temperature: float = 0.1
    protein_contrast_loss_weight: float = 0.5
    ranking_temperature: float = 1.0
    # Form of the within-target ranking term. "ranknet" is the incumbent
    # pairwise softplus; "listce" is the R14 regression-compatible listwise
    # cross-entropy, whose optimum coincides with the regression optimum so it
    # contributes exactly zero gradient at `prediction == truth`. See
    # report/meta_fewshot/stageR14_diagnostics_20260816/.
    ranking_loss_form: str = "ranknet"
    representation_warmup_fraction: float = 0.0
    admission_binding_margin_pk: float = 0.01
    zero_support_only: bool = False
    pretrained_checkpoint: str | None = None
    geometry_checkpoint: str | None = None
    amp: bool = True
    pair_dim: int = 96
    pair_blocks: int = 4
    pair_latents: int = 24
    pair_heads: int = 8
    pair_chunk_size: int = 12
    support_hidden_dim: int = 192
    support_blocks: int = 2
    use_cartesian: bool = False
    episode_cache: str | None = None
    # Optional governed split that replaces the corpus `split` field. The
    # default keeps the retained CD-HIT40 protocol exactly as it was.
    split_directory: str | None = None
    arch: str = "bpsf"
    difference_loss_weight: float = 0.0
    lr_schedule: str = "constant"
    lr_warmup_fraction: float = 0.05
    lr_final_fraction: float = 0.1
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
            "hidden_dim", "task_dim", "ligand_layers", "pair_dim",
            "pair_blocks", "pair_latents", "pair_heads",
            "support_hidden_dim", "support_blocks", "adapter_rank",
            "adaptive_blocks", "adapter_scale"):
        if parent_config.get(field) != getattr(config, field):
            raise ValueError(f"parent checkpoint has incompatible {field}")
    model.load_state_dict(payload["model_state"], strict=True)
    return file_sha256(path)


def load_geometry_checkpoint(model: QPSMPBioModel, config: TrainConfig) -> str | None:
    if config.geometry_checkpoint is None:
        return None
    path = Path(config.geometry_checkpoint)
    payload = torch.load(path, map_location="cpu", weights_only=False)
    teacher = payload.get("config", {})
    required = {
        "hidden_dim": config.hidden_dim,
        "section_dim": config.task_dim,
        "pair_dim": config.pair_dim,
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
                    # Retained mechanism slots add outputs but do not alter the
                    # coordinate-free prefix learned by the geometry teacher.
                    if target_name.startswith("pair_section.latent"):
                        continue
                    raise ValueError(f"geometry tensor is incompatible: {name}")
                target[target_name].copy_(value)
                loaded.add(target_name)
    if not loaded or not any(name.startswith("pair_section.") for name in loaded):
        raise ValueError("geometry checkpoint did not contain a pair trunk")
    model.load_state_dict(target)
    return file_sha256(path)


def learning_rate_factor(step: int, config: TrainConfig) -> float:
    """Multiplier applied to every parameter group at `step` (1-indexed)."""
    if config.lr_schedule == "constant":
        return 1.0
    if config.lr_schedule != "cosine":
        raise ValueError(f"unknown learning-rate schedule: {config.lr_schedule}")
    warmup = max(1, int(config.steps * config.lr_warmup_fraction))
    if step <= warmup:
        return step / warmup
    progress = (step - warmup) / max(1, config.steps - warmup)
    final = config.lr_final_fraction
    return final + (1.0 - final) * 0.5 * (1.0 + np.cos(np.pi * min(progress, 1.0)))


def normalized_episode(episode: EpisodeBatch, scale: LabelScale) -> EpisodeBatch:
    return replace(episode, support_y=scale.normalize(episode.support_y),
                   query_y=scale.normalize(episode.query_y))


def compact_episode(episode: EpisodeBatch) -> EpisodeBatch:
    if episode.support_coordinates is not None or episode.query_coordinates is not None:
        if episode.support_coordinates is None or episode.query_coordinates is None:
            raise ValueError("support/query Cartesian coordinates must be paired")
        # Coordinate node order is tied to the original padded atom width.
        # Keep that width rather than silently invalidating the packed edges.
        return episode
    support_active = (int(episode.support_mask.sum(1).max())
                      if episode.support_mask.shape[0] else 0)
    active = max(support_active, int(episode.query_mask.sum(1).max()))
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


def forward(model: QPSMPBioModel, episode: EpisodeBatch, *, adapt: bool = True,
            task_state_override: torch.Tensor | None = None):
    if episode.support_atoms.ndim == 3:
        episode = compact_episode(episode)
    has_geometry = (episode.support_coordinates is not None
                    or episode.query_coordinates is not None)
    if has_geometry:
        if episode.support_coordinates is None or episode.query_coordinates is None:
            raise ValueError("support/query Cartesian coordinates must be paired")
        active = episode.query_atoms.shape[-2]
        if episode.support_atoms.shape[-2] != active:
            raise ValueError("Cartesian episodes require a common padded atom width")
    else:
        support_active = (int(episode.support_mask.sum(-1).max())
                          if episode.support_mask.numel() else 0)
        active = max(support_active, int(episode.query_mask.sum(-1).max()))
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
    if not has_geometry:
        def pad_atoms(values: torch.Tensor) -> torch.Tensor:
            return F.pad(values, (0, 0, 0, active - values.shape[-2]))
        def pad_bonds(values: torch.Tensor) -> torch.Tensor:
            missing = active - values.shape[-2]
            return F.pad(values, (0, 0, 0, missing, 0, missing))
        def pad_mask(values: torch.Tensor) -> torch.Tensor:
            return F.pad(values, (0, active - values.shape[-1]))
        support_atoms, query_atoms = pad_atoms(support_atoms), pad_atoms(query_atoms)
        support_bonds, query_bonds = pad_bonds(support_bonds), pad_bonds(query_bonds)
        support_mask, query_mask = pad_mask(support_mask), pad_mask(query_mask)
    geometry_coordinates = None
    if has_geometry:
        pair_axis = 1 if episode.support_coordinates.ndim == 4 else 0
        geometry_coordinates = torch.cat(
            (episode.support_coordinates, episode.query_coordinates), pair_axis)
    extra = {}
    if getattr(model, "transport", None).__class__.__name__ == "SimilarityTransport":
        extra = {"support_fingerprint": episode.support_fingerprint,
                 "query_fingerprint": episode.query_fingerprint}
    return model(
        episode.protein_pooled, episode.protein_tokens, episode.protein_mask,
        support_atoms, support_bonds, support_mask, episode.support_y,
        query_atoms, query_bonds, query_mask, adapt=adapt,
        protein_chemistry=episode.protein_chemistry, **extra,
        task_state_override=task_state_override,
        geometry_coordinates=geometry_coordinates,
        geometry_edge_index=episode.geometry_edge_index,
        geometry_available=episode.geometry_available,
        geometry_common_frame=episode.geometry_common_frame)


def centered_task_error(prediction: torch.Tensor,
                        truth: torch.Tensor) -> torch.Tensor:
    prediction = prediction - prediction.mean(-1, keepdim=True)
    truth = truth - truth.mean(-1, keepdim=True)
    return (prediction - truth).square().mean()


def pairwise_ranking_loss(prediction: torch.Tensor, truth: torch.Tensor,
                          temperature: float) -> torch.Tensor:
    delta_y = truth.unsqueeze(-1) - truth.unsqueeze(-2)
    delta_p = prediction.unsqueeze(-1) - prediction.unsqueeze(-2)
    comparable = delta_y != 0
    if not bool(comparable.any()):
        return prediction.new_zeros(())
    signed = delta_y.sign() * delta_p / temperature
    return torch.nn.functional.softplus(-signed[comparable]).mean()


#: Fixed shift below the pK label range, used by the listwise transform. It is
#: not a swept hyperparameter: any constant below the corpus minimum gives the
#: same stationary point, and 2.0 pK sits under every measured Ki in the bank.
LISTCE_SHIFT_PK = 2.0


def regression_compatible_ranking_loss(prediction: torch.Tensor,
                                       truth: torch.Tensor,
                                       label_scale: "LabelScale") -> torch.Tensor:
    """Within-target listwise cross-entropy that shares the regression optimum.

    The incumbent `pairwise_ranking_loss` is minimised by *any* scale-free
    monotone arrangement, so its optimum differs from the regression optimum
    and it keeps pushing for wider margins even at `prediction == truth`. The
    R14 Phase 2 audit measured the cost of that residual gradient: within-target
    Pearson `r` at k=0 is lower than the incumbent's in 8 of 8 arms.

    Following the regression-compatible construction of Bai et al.
    (arXiv:2211.01494), the normalizer uses the model's own value link instead
    of `exp`. With a fixed shift `m` below the label range, `T(x) = x - m` and
    weights `w = y - m`:

        L = -(1/Σ_j w_j) Σ_i w_i · log[ T(p_i) / Σ_j T(p_j) ]
        ∂L/∂p_k ∝ w_k/T(p_k) - Σ_j w_j / Σ_j T(p_j)

    which vanishes for every `k` exactly when `T(p) ∝ w`. The regression
    optimum `p = y` is such a point, so the two terms share a minimum rather
    than competing for it, and the regression term alone pins the amplitude.

    Inputs are in normalized label units; the transform is applied in pK,
    which is an affine reparameterization and therefore leaves the stationary
    point unchanged.
    """
    truth = truth.to(device=prediction.device, dtype=prediction.dtype)
    truth_pk = truth * label_scale.scale + label_scale.mean
    prediction_pk = prediction * label_scale.scale + label_scale.mean
    weights = (truth_pk - LISTCE_SHIFT_PK).clamp_min(1e-6)
    transformed = (prediction_pk - LISTCE_SHIFT_PK).clamp_min(1e-6)
    total_weight = weights.sum(-1, keepdim=True)
    if not bool((total_weight > 0).all()):
        return prediction.new_zeros(())
    normalized_log = transformed.log() - transformed.sum(-1, keepdim=True).log()
    return -((weights * normalized_log).sum(-1) / total_weight.squeeze(-1)).mean()


def ranking_term(prediction: torch.Tensor, truth: torch.Tensor,
                 config: "TrainConfig", label_scale: "LabelScale") -> torch.Tensor:
    """Dispatch the configured within-target ranking form. One changed variable."""
    if config.ranking_loss_form == "listce":
        return regression_compatible_ranking_loss(prediction, truth, label_scale)
    if config.ranking_loss_form == "ranknet":
        return pairwise_ranking_loss(prediction, truth, config.ranking_temperature)
    raise ValueError(f"unknown ranking_loss_form {config.ranking_loss_form!r}; "
                     "expected 'ranknet' or 'listce'")


def relative_difference_loss(output, query_y: torch.Tensor) -> torch.Tensor:
    """Per-(query, support) supervision of the signed difference operator.

    The prediction is exact when `delta(k->q) == rho_q - r_k`, where
    `rho_q = y_q - f0(q)` and `r_k = y_k - f0(L_k)`.  Supervising every entry of
    the difference tensor rather than only its weighted sum is what gives the
    relative operator a direct training signal.  Query labels are used as a loss
    target only; they are never a model input, and both sides of the target are
    detached so the operator stays label locked.
    """
    delta = output.support_evidence
    if delta.ndim < 2 or delta.shape[-1] == 0:
        return delta.new_zeros(())
    residual = output.support_residual_quotient + output.level_adjustment
    query_y = query_y.to(device=delta.device, dtype=delta.dtype)
    query_residual = query_y - output.zero_shot
    target = (query_residual.unsqueeze(-1) - residual.unsqueeze(-2)).detach()
    if target.shape != delta.shape:
        raise ValueError(
            f"difference target {tuple(target.shape)} does not match "
            f"operator {tuple(delta.shape)}")
    return (delta - target).square().mean()


def binding_contrastive_loss(errors: list[torch.Tensor], temperature: float) -> torch.Tensor:
    """Prefer the true ligand-label binding over matched counterfactuals."""
    if len(errors) < 2:
        raise ValueError("binding contrast requires a counterfactual error")
    logits = -torch.stack(errors) / temperature
    return F.cross_entropy(logits.unsqueeze(0), logits.new_zeros(1, dtype=torch.long))


def matched_wrong_labels(output, support_y: torch.Tensor) -> torch.Tensor:
    """Flip the zero-shot residual magnitude for a one-support counterfactual."""
    if support_y.shape[-1] != 1:
        raise ValueError("magnitude-matched labels are defined only for k=1")
    raw_residual = output.support_residual_quotient + output.level_adjustment[..., :1]
    support_y = support_y.to(device=raw_residual.device, dtype=raw_residual.dtype)
    return support_y - 2.0 * raw_residual.detach()


def counterfactual_label_assignments(output, support_y: torch.Tensor) -> list[torch.Tensor]:
    """Enumerate non-identity cyclic bindings; k=1 uses an equal-magnitude flip."""
    support_count = support_y.shape[-1]
    if support_count == 0:
        return []
    if support_count == 1:
        return [matched_wrong_labels(output, support_y)]
    return [support_y.roll(shift, dims=-1)
            for shift in range(1, support_count)]


def batch_counterfactual_episode(
        episode: EpisodeBatch, assignments: list[torch.Tensor]) -> EpisodeBatch:
    """Stack label assignments on a new batch axis for one model forward."""
    if not assignments:
        raise ValueError("counterfactual batch cannot be empty")
    count = len(assignments)
    def repeat(value):
        if not isinstance(value, torch.Tensor):
            return value
        return value.unsqueeze(0).expand(count, *value.shape)
    if episode.support_coordinates is not None:
        raise ValueError("batched counterfactual geometry is not implemented")
    return replace(
        episode,
        protein_pooled=repeat(episode.protein_pooled),
        protein_tokens=repeat(episode.protein_tokens),
        protein_mask=repeat(episode.protein_mask),
        protein_chemistry=repeat(episode.protein_chemistry),
        support_fingerprint=repeat(episode.support_fingerprint),
        query_fingerprint=repeat(episode.query_fingerprint),
        support_atoms=repeat(episode.support_atoms),
        support_bonds=repeat(episode.support_bonds),
        support_mask=repeat(episode.support_mask),
        support_y=torch.stack(assignments),
        query_atoms=repeat(episode.query_atoms),
        query_bonds=repeat(episode.query_bonds),
        query_mask=repeat(episode.query_mask),
        query_y=repeat(episode.query_y))


def admission_score(metrics: dict[str, float], binding_margin_pk: float) -> float:
    """Balance zero-shot quality with useful, label-sensitive adaptation."""
    full = metrics["full_mse_pk"]
    level = metrics["sar_cut_mse_pk"]
    binding_gap = metrics["permuted_mse_pk"] - full
    return (full + 0.25 * metrics["zero_shot_mse_pk"]
            + max(0.0, full - level)
            + max(0.0, binding_margin_pk - binding_gap))


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


def donor_episode(model: QPSMPBioModel, data: QPSMPData,
                  episode: EpisodeBatch, donor_target: str,
                  label_scale: LabelScale) -> EpisodeBatch:
    """Materialize a complete cross-component foreign support episode."""
    donor_indices = data.tasks[episode.spec.split][donor_target]
    count = len(episode.spec.support)
    order = np.random.default_rng(
        sum(episode.spec.support) + len(donor_target)).permutation(donor_indices)
    if len(order) < count + 1:
        raise ValueError("foreign target cannot supply a complete support episode")
    donor_spec = EpisodeSpec(
        episode.spec.split, data.cells[int(order[0])]["protein_group_40"],
        donor_target, tuple(map(int, order[:count])), (int(order[count]),),
        episode.spec.target)
    return normalized_episode(data.materialize(donor_spec), label_scale)


def complete_foreign_prediction(model: QPSMPBioModel, data: QPSMPData,
                                recipient: EpisodeBatch, donor_target: str,
                                label_scale: LabelScale) -> torch.Tensor:
    """Replace support ligand/label pairs under the recipient protein."""
    donor = donor_episode(model, data, recipient, donor_target, label_scale)
    hybrid = replace(
        donor, protein_pooled=recipient.protein_pooled,
        protein_tokens=recipient.protein_tokens,
        protein_mask=recipient.protein_mask,
        protein_chemistry=recipient.protein_chemistry,
        query_atoms=recipient.query_atoms, query_bonds=recipient.query_bonds,
        query_mask=recipient.query_mask, query_y=recipient.query_y,
        query_fingerprint=recipient.query_fingerprint)
    return forward(model, compact_episode(hybrid)).prediction


def wrong_protein_prediction(model: QPSMPBioModel, data: QPSMPData,
                             episode: EpisodeBatch, donor_target: str,
                             *, adapt: bool = True) -> torch.Tensor:
    """Replace the protein throughout the full recipient episode."""
    pooled, tokens, mask = data.protein_for_target(donor_target)
    wrong_episode = replace(
        episode, protein_pooled=pooled, protein_tokens=tokens, protein_mask=mask,
        protein_chemistry=data.protein_chemistry_for_target(donor_target))
    return forward(model, wrong_episode, adapt=adapt).prediction


def wrong_protein_zero_shot(model: QPSMPBioModel, data: QPSMPData,
                            episode: EpisodeBatch, donor_target: str) -> torch.Tensor:
    return wrong_protein_prediction(
        model, data, episode, donor_target, adapt=False)


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
            sar_cut = full.prediction - full.sar_adaptation
            level_only = full.level_baseline
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
                # Degeneracy diagnostic: a zero-shot endpoint that does not
                # vary across the queries of one episode is a constant, and
                # every apparent few-shot gain is then pure level calibration.
                "zero_shot_query_spread_pk": float(
                    frozen.prediction.std(-1, unbiased=False).mean()
                    * label_scale.scale),
                "query_truth_spread_pk": float(
                    episode.query_y.float().std(-1, unbiased=False).mean()
                    * label_scale.scale),
            }
            if controls:
                assignments = counterfactual_label_assignments(
                    full, episode.support_y)
                wrong_episode = batch_counterfactual_episode(episode, assignments)
                wrong_output = forward(model, wrong_episode)
                counterfactual_mse = label_scale.squared_error_pk(
                    wrong_output.prediction, wrong_episode.query_y).mean(-1)
                values.update({
                    "permuted_mse_pk": float(counterfactual_mse.mean()),
                    "foreign_code_state_mse_pk": mse_pk(
                        complete_foreign_prediction(
                            model, data, episode, spec.donor_target,
                            label_scale)),
                    "wrong_protein_state_mse_pk": mse_pk(
                        wrong_protein_prediction(
                            model, data, episode, spec.donor_target)),
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
        metrics["foreign_code_state_gap_mse_pk"] = (
            metrics["foreign_code_state_mse_pk"] - metrics["full_mse_pk"])
        metrics["wrong_protein_gap_mse_pk"] = metrics["wrong_protein_state_mse_pk"] - metrics["full_mse_pk"]
        metrics["wrong_protein_zero_shot_gap_mse_pk"] = (
            metrics["wrong_protein_zero_shot_mse_pk"] - metrics["zero_shot_mse_pk"])
    metrics["weighting"] = "equal_component_then_equal_target_then_equal_draw"
    metrics["episodes"] = len(rows)
    return metrics


def train(data: QPSMPData, config: TrainConfig,
          support_sizes: tuple[int, ...] | None = None,
          progress_path: Path | None = None,
          model_factory=None):
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model_factory = model_factory or resolve_architecture(config.arch)
    model = model_factory(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers,
        pair_dim=config.pair_dim, pair_blocks=config.pair_blocks,
        pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks,
        adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks,
        adapter_scale=config.adapter_scale,
        use_cartesian=config.use_cartesian,
        dtype=torch.float32)
    model.to(config.device)
    geometry_sha256 = load_geometry_checkpoint(model, config)
    parent_sha256 = load_parent_checkpoint(model, config)
    if geometry_sha256 is not None and parent_sha256 is not None:
        raise ValueError("geometry and parent checkpoints are mutually exclusive initializers")
    trainable_parameters = [parameter for parameter in model.parameters()
                            if parameter.requires_grad]
    fast_prefixes = ("meta.term.", "transport.")
    if config.difference_loss_weight and config.arch not in RELATIVE_ARCHITECTURES:
        raise ValueError(
            "difference_loss_weight requires a relative-transport architecture")
    fast_parameters = [parameter for name, parameter in model.named_parameters()
                       if parameter.requires_grad
                       and name.startswith(fast_prefixes)]
    fast_ids = {id(parameter) for parameter in fast_parameters}
    slow_parameters = [parameter for parameter in trainable_parameters
                       if id(parameter) not in fast_ids]
    optimizer = torch.optim.AdamW([
        {"params": slow_parameters,
         "lr": config.learning_rate * config.backbone_lr_scale},
        {"params": fast_parameters, "lr": config.learning_rate},
    ], weight_decay=config.weight_decay)
    base_learning_rates = [group["lr"] for group in optimizer.param_groups]
    amp_enabled = config.amp and config.device.startswith("cuda")
    scaler = torch.amp.GradScaler("cuda", enabled=amp_enabled)
    label_scale = training_label_scale(data)
    train_support_sizes = support_sizes or (0, 1, 2, 3, 5)
    if not train_support_sizes or any(k < 0 for k in train_support_sizes):
        raise ValueError("support sizes cannot be negative")
    if config.min_query_size < 1 or config.query_size < config.min_query_size:
        raise ValueError("query-size range is invalid")
    val_targets_per_component = (config.eval_targets_per_component
                                 if config.val_targets_per_component is None
                                 else config.val_targets_per_component)
    cache_contract = {
        "schema": "MetaSieve.QPSMPValidationCache.v2",
        "seed": config.evaluation_seed,
        "support_sizes": list(train_support_sizes),
        "query_size": config.query_size,
        "val_draws_per_target": config.val_draws_per_target,
        "eval_targets_per_component": val_targets_per_component,
        "label_scale": asdict(label_scale),
    }
    cache_path = Path(config.episode_cache) if config.episode_cache else None
    if cache_path is not None and cache_path.exists():
        cached = torch.load(cache_path, map_location="cpu", weights_only=False)
        if cached.get("contract") != cache_contract:
            raise ValueError("episode cache contract does not match training configuration")
        val_banks = cached["val"]
    else:
        val_specs = data.fixed_nested_episode_banks(
            "meta_val", tuple(train_support_sizes), config.query_size,
            config.val_draws_per_target, config.evaluation_seed,
            val_targets_per_component)
        val_banks = {
            k: tuple(compact_episode(normalized_episode(data.materialize(spec), label_scale))
                     for spec in specs)
            for k, specs in val_specs.items()
        }
        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            torch.save({"schema": "MetaSieve.QPSMPValidationCache.v2",
                        "contract": cache_contract, "val": val_banks}, cache_path)

    diversity = {"components": set(), "targets": set(), "support_cells": set(),
                 "query_cells": set(), "ligands": set(), "episodes": 0}
    best_state, best_value, best_step = None, float("inf"), 0
    trace = []
    started = time.monotonic()
    for step in range(1, config.steps + 1):
        model.train()
        factor = learning_rate_factor(step, config)
        for group, base in zip(optimizer.param_groups, base_learning_rates):
            group["lr"] = base * factor
        optimizer.zero_grad(set_to_none=True)
        support_size = train_support_sizes[(step - 1) % len(train_support_sizes)]
        selected = []
        for _ in range(config.episodes_per_step):
            requested_query = int(rng.integers(
                config.min_query_size, config.query_size + 1))
            spec = data.draw_episode(
                "meta_train", support_size, requested_query, rng)
            selected.append(compact_episode(
                normalized_episode(data.materialize(spec), label_scale)))
        for item in selected:
            spec = item.spec
            diversity["episodes"] += 1
            diversity["components"].add(spec.component)
            diversity["targets"].add(spec.target)
            diversity["support_cells"].update(spec.support)
            diversity["query_cells"].update(spec.query)
            diversity["ligands"].update(
                data.cells[index]["ligand_id"]
                for index in (*spec.support, *spec.query))
        warmup_steps = max(
            0, int(config.steps * config.representation_warmup_fraction))
        phase_a = step <= warmup_steps
        loss_value = 0.0
        with torch.autocast(
                device_type="cuda",
                dtype=(torch.bfloat16 if torch.cuda.is_bf16_supported()
                       else torch.float16), enabled=amp_enabled):
            for episode in selected:
                adapt = not (config.zero_support_only or phase_a)
                full = forward(model, episode, adapt=adapt)
                query_y = episode.query_y.to(
                    device=full.prediction.device, dtype=full.prediction.dtype)
                endpoint_prediction = full.prediction if adapt else full.zero_shot
                loss_full = F.smooth_l1_loss(endpoint_prediction, query_y)
                loss_zero = F.smooth_l1_loss(full.zero_shot, query_y)
                loss_rank = ranking_term(
                    full.prediction, query_y, config, label_scale)
                correct_error = centered_task_error(full.prediction, query_y)
                loss_binding = loss_full.new_zeros(())
                loss_protein = loss_full.new_zeros(())
                if adapt and support_size > 0:
                    binding_errors = [(full.prediction - query_y).square().mean()]
                    assignments = counterfactual_label_assignments(
                        full, episode.support_y)
                    wrong_episode = batch_counterfactual_episode(
                        episode, assignments)
                    wrong = forward(model, wrong_episode)
                    wrong_truth = wrong_episode.query_y.to(
                        device=wrong.prediction.device, dtype=wrong.prediction.dtype)
                    binding_errors.extend(
                        (wrong.prediction - wrong_truth).square().mean(-1).unbind())
                    loss_binding = binding_contrastive_loss(
                        binding_errors, config.binding_temperature)
                    wrong_zero = wrong_protein_zero_shot(
                        model, data, episode, episode.spec.donor_target)
                    correct_zero_error = (full.zero_shot - query_y).square().mean()
                    wrong_zero_error = (wrong_zero - query_y).square().mean()
                    loss_protein = binding_contrastive_loss(
                        [correct_zero_error, wrong_zero_error],
                        config.binding_temperature)
                if phase_a:
                    episode_loss = (
                        loss_zero
                        + config.ranking_loss_weight * ranking_term(
                            full.zero_shot, query_y, config, label_scale)
                        + config.support_match_loss_weight * full.support_match_loss)
                else:
                    episode_loss = (
                        loss_full + config.zero_shot_loss_weight * loss_zero
                        + config.ranking_loss_weight * loss_rank
                        + config.shape_loss_weight * correct_error
                        + config.support_match_loss_weight * full.support_match_loss
                        + config.binding_loss_weight * loss_binding
                        + config.protein_contrast_loss_weight * loss_protein)
                    if config.difference_loss_weight and adapt and support_size > 0:
                        episode_loss = episode_loss + (
                            config.difference_loss_weight
                            * relative_difference_loss(full, query_y))
                loss_value += float(episode_loss.detach())
                scaler.scale(episode_loss / len(selected)).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        scaler.step(optimizer)
        scaler.update()
        trace.append(loss_value / len(selected))
        if step % config.val_interval == 0 or step == config.steps:
            validation_metrics = {k: evaluate(
                model, data, bank,
                controls=not config.zero_support_only and k > 0,
                label_scale=label_scale) for k, bank in val_banks.items()}
            values = ([item["zero_shot_mse_pk"]
                       for item in validation_metrics.values()]
                      if config.zero_support_only else
                      [(admission_score(item, config.admission_binding_margin_pk)
                        if k > 0 else item["full_mse_pk"])
                       for k, item in validation_metrics.items()])
            value = float(np.mean(values))
            progress = {
                "step": step, "phase": "representation" if phase_a else "meta",
                "validation_admission_score": value,
                "best_validation_admission_score": min(best_value, value),
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
    return model, {"best_val_admission_score": best_value,
                   "best_step": best_step, "loss_trace": trace,
                   "loss_trace_units": "standardized_squared_error",
                   "label_scale": asdict(label_scale),
                   "support_sizes": list(train_support_sizes),
                   "validation_episode_bank_sizes": {
                       str(k): len(bank) for k, bank in val_banks.items()
                   },
                   "val_targets_per_component": val_targets_per_component,
                   "training_sampler": "online_fresh_component_uniform",
                   "training_query_size_range": [
                       config.min_query_size, config.query_size],
                   "training_diversity": {
                       key: (len(value) if isinstance(value, set) else value)
                       for key, value in diversity.items()
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
    parser.add_argument("--evaluation-seed", type=int,
                        default=TrainConfig.evaluation_seed)
    parser.add_argument("--support-size", type=int, default=TrainConfig.support_size)
    parser.add_argument("--query-size", type=int, default=TrainConfig.query_size)
    parser.add_argument("--min-query-size", type=int,
                        default=TrainConfig.min_query_size)
    parser.add_argument("--hidden-dim", type=int, default=TrainConfig.hidden_dim)
    parser.add_argument("--task-dim", type=int, default=TrainConfig.task_dim)
    parser.add_argument("--adapter-rank", type=int, default=TrainConfig.adapter_rank)
    parser.add_argument("--adaptive-blocks", type=int, default=TrainConfig.adaptive_blocks)
    parser.add_argument("--adapter-scale", type=float, default=TrainConfig.adapter_scale)
    parser.add_argument("--ligand-layers", type=int, default=TrainConfig.ligand_layers)
    parser.add_argument("--steps", type=int, default=TrainConfig.steps)
    parser.add_argument("--learning-rate", type=float, default=TrainConfig.learning_rate)
    parser.add_argument("--backbone-lr-scale", type=float,
                        default=TrainConfig.backbone_lr_scale)
    parser.add_argument("--episodes-per-step", type=int, default=TrainConfig.episodes_per_step)
    parser.add_argument("--val-interval", type=int, default=TrainConfig.val_interval)
    parser.add_argument("--val-draws-per-target", type=int, default=TrainConfig.val_draws_per_target)
    parser.add_argument("--test-draws-per-target", type=int, default=TrainConfig.test_draws_per_target)
    parser.add_argument("--eval-targets-per-component", type=int,
                        default=TrainConfig.eval_targets_per_component)
    parser.add_argument("--val-targets-per-component", type=int, default=None)
    parser.add_argument("--zero-shot-loss-weight", type=float,
                        default=TrainConfig.zero_shot_loss_weight)
    parser.add_argument("--ranking-loss-weight", type=float,
                        default=TrainConfig.ranking_loss_weight)
    parser.add_argument("--ranking-loss-form", choices=("ranknet", "listce"),
                        default=TrainConfig.ranking_loss_form,
                        help="form of the within-target ranking term; 'listce' "
                             "is the R14 regression-compatible construction "
                             "whose optimum coincides with the regression "
                             "optimum (default: ranknet, the incumbent)")
    parser.add_argument("--shape-loss-weight", type=float,
                        default=TrainConfig.shape_loss_weight)
    parser.add_argument("--support-match-loss-weight", type=float,
                        default=TrainConfig.support_match_loss_weight)
    parser.add_argument("--binding-loss-weight", type=float,
                        default=TrainConfig.binding_loss_weight)
    parser.add_argument("--binding-temperature", type=float,
                        default=TrainConfig.binding_temperature)
    parser.add_argument("--protein-contrast-loss-weight", type=float,
                        default=TrainConfig.protein_contrast_loss_weight)
    parser.add_argument("--representation-warmup-fraction", type=float,
                        default=TrainConfig.representation_warmup_fraction)
    parser.add_argument("--zero-support-only", action="store_true")
    parser.add_argument("--pretrained-checkpoint", type=Path)
    parser.add_argument("--geometry-checkpoint", type=Path)
    parser.add_argument("--no-amp", action="store_true")
    parser.add_argument("--pair-dim", type=int, default=TrainConfig.pair_dim)
    parser.add_argument("--pair-blocks", type=int, default=TrainConfig.pair_blocks)
    parser.add_argument("--pair-latents", type=int, default=TrainConfig.pair_latents)
    parser.add_argument("--pair-heads", type=int, default=TrainConfig.pair_heads)
    parser.add_argument("--pair-chunk-size", type=int, default=TrainConfig.pair_chunk_size)
    parser.add_argument("--support-hidden-dim", type=int,
                        default=TrainConfig.support_hidden_dim)
    parser.add_argument("--support-blocks", type=int, default=TrainConfig.support_blocks)
    parser.add_argument("--use-cartesian", action="store_true",
                        help="enable optional common-frame Cartesian inputs; the active BindingDB banks provide none")
    parser.add_argument("--episode-cache", type=Path)
    parser.add_argument("--split-directory", type=Path, default=None)
    parser.add_argument("--include-meta-test", action="store_true",
                        help="physically unseal the meta_test cells in QPSMPData "
                             "(contract 2026-08-16: off by default)")
    parser.add_argument("--eval-meta-test", action="store_true",
                        help="evaluate meta_test after training; requires "
                             "--include-meta-test (off by default)")
    parser.add_argument("--arch", default=TrainConfig.arch,
                        choices=sorted(ARCHITECTURES))
    parser.add_argument("--difference-loss-weight", type=float,
                        default=TrainConfig.difference_loss_weight)
    parser.add_argument("--lr-schedule", default=TrainConfig.lr_schedule,
                        choices=("constant", "cosine"))
    parser.add_argument("--lr-warmup-fraction", type=float,
                        default=TrainConfig.lr_warmup_fraction)
    parser.add_argument("--lr-final-fraction", type=float,
                        default=TrainConfig.lr_final_fraction)
    parser.add_argument("--device", default=TrainConfig.device)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)
    config = TrainConfig(
        seed=args.seed, evaluation_seed=args.evaluation_seed,
        support_size=args.support_size, query_size=args.query_size,
        min_query_size=args.min_query_size,
        hidden_dim=args.hidden_dim, task_dim=args.task_dim,
        adapter_rank=args.adapter_rank, adaptive_blocks=args.adaptive_blocks,
        adapter_scale=args.adapter_scale,
        ligand_layers=args.ligand_layers,
        steps=args.steps, episodes_per_step=args.episodes_per_step,
        learning_rate=args.learning_rate,
        backbone_lr_scale=args.backbone_lr_scale,
        val_interval=args.val_interval,
        val_draws_per_target=args.val_draws_per_target,
        test_draws_per_target=args.test_draws_per_target,
        eval_targets_per_component=args.eval_targets_per_component,
        val_targets_per_component=args.val_targets_per_component,
        zero_shot_loss_weight=args.zero_shot_loss_weight,
        ranking_loss_weight=args.ranking_loss_weight,
        ranking_loss_form=args.ranking_loss_form,
        shape_loss_weight=args.shape_loss_weight,
        support_match_loss_weight=args.support_match_loss_weight,
        binding_loss_weight=args.binding_loss_weight,
        binding_temperature=args.binding_temperature,
        protein_contrast_loss_weight=args.protein_contrast_loss_weight,
        representation_warmup_fraction=args.representation_warmup_fraction,
        zero_support_only=args.zero_support_only,
        pretrained_checkpoint=(str(args.pretrained_checkpoint.resolve())
                               if args.pretrained_checkpoint else None),
        geometry_checkpoint=(str(args.geometry_checkpoint.resolve())
                             if args.geometry_checkpoint else None),
        amp=not args.no_amp,
        pair_dim=args.pair_dim, pair_blocks=args.pair_blocks,
        pair_latents=args.pair_latents,
        pair_heads=args.pair_heads, pair_chunk_size=args.pair_chunk_size,
        support_hidden_dim=args.support_hidden_dim,
        support_blocks=args.support_blocks,
        use_cartesian=args.use_cartesian,
        episode_cache=(str(args.episode_cache.resolve()) if args.episode_cache else None),
        split_directory=(str(args.split_directory.resolve())
                         if args.split_directory else None),
        arch=args.arch, difference_loss_weight=args.difference_loss_weight,
        lr_schedule=args.lr_schedule,
        lr_warmup_fraction=args.lr_warmup_fraction,
        lr_final_fraction=args.lr_final_fraction, device=args.device)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=args.include_meta_test)
    if args.eval_meta_test and not args.include_meta_test:
        parser.error("--eval-meta-test requires --include-meta-test")
    model, training, label_scale = train(
        data, config, progress_path=args.output / "progress.jsonl")
    checkpoint = args.output / "checkpoint.pt"
    torch.save({"model_state": model.state_dict(), "config": asdict(config)}, checkpoint)
    metrics = {}
    if args.eval_meta_test:
        test_specs = data.fixed_nested_episode_banks(
            "meta_test", (0, 1, 2, 3, 5), config.query_size,
            config.test_draws_per_target, config.evaluation_seed,
            config.eval_targets_per_component)
        metrics = {
            str(k): evaluate(
                model, data,
                tuple(compact_episode(normalized_episode(
                    data.materialize(spec), label_scale)) for spec in specs),
                controls=k > 0, label_scale=label_scale)
            for k, specs in test_specs.items()
        }
    result = {
        "schema": "MetaSieve.LCIPFELMTTrainingRun.v2",
        "scope": "implementation_smoke_only",
        "data": {"corpus": str(CORPUS), "protein_bank_records": len(data.protein_bank),
                 "ligand_bank_records": len(data.ligand_bank)},
        "config": asdict(config), "training": training,
        "meta_test": {
            "included": bool(args.include_meta_test),
            "evaluated": bool(args.eval_meta_test),
            "seal": "physical: QPSMPData include_meta_test flag "
                    "(contract 2026-08-16)",
        },
        "checkpoint_sha256": file_sha256(checkpoint),
        "test": metrics,
        "evaluation_population": "fixed hash-selected targets within every eligible component",
        "controls": {"foreign_code": "donor support ligand/label pairs replace recipient support and are recomputed under the fixed recipient protein and query ligands",
                     "wrong_protein": "ELMT evidence inferred under a donor protein replaces recipient evidence",
                     "sar_cut": "ELMT primitive composition is removed while zero-shot and scalar level calibration remain",
                     "permuted": "support labels are cyclically permuted; k=1 uses an equal-magnitude residual flip"},
        "gate_authorization": {"G2": False, "G3a": False, "G3b": False},
        "authorization_reason": "A training smoke cannot authorize preregistered inferential gates.",
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
