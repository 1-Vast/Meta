"""Stage A trainer: one code path for A0, A1 and A2.

`--inner-steps 0` must reproduce the accepted baseline bitwise. That is what
makes `A0` matched by construction rather than by inspection, and
`test_zero_inner_steps_reproduces_the_production_loss` pins it. Every other
difference between the arms is a named flag.

Arms:

* `A0` — `--inner-steps 0`, uniform sampling. The accepted recipe.
* `A1` — `--inner-steps 1`, uniform sampling.
* `A2` — `--inner-steps 1`, `--task-selection`.

Adaptation runs in float32 outside autocast: an inner gradient on 97 parameters
under float16 can underflow to zero, which would make the whole experiment
measure "no adaptation" while appearing to run correctly.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from model.similarity_grammar import tanimoto                        # noqa: E402
from scripts.qpsmp_data import QPSMPData                             # noqa: E402
from scripts.train_qpsmp import (                                    # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    admission_score, binding_contrastive_loss, centered_task_error,
    compact_episode, evaluate, learning_rate_factor, normalized_episode,
    ranking_term, resolve_architecture, training_label_scale,
)
from tools.research.stageA_innerloop.inner_loop import (             # noqa: E402
    AdaptationConfig, adapt, base_weights, encode_task, readout,
    support_query_gradient_cosine,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"


@dataclass(frozen=True)
class MetaConfig:
    """The base recipe plus exactly the fields this experiment adds."""
    base: TrainConfig = field(default_factory=lambda: TrainConfig(
        arch="similarity_only", steps=1200, seed=20260815,
        split_directory=str(SPLIT)))
    inner_steps: int = 0
    inner_lr: float = 0.1
    first_order: bool = True
    task_selection: bool = False
    candidate_multiplier: int = 3
    post_adaptation_loss_weight: float = 1.0

    def adaptation(self) -> AdaptationConfig:
        return AdaptationConfig(inner_steps=self.inner_steps,
                                inner_lr=self.inner_lr,
                                first_order=self.first_order)


def episode_tensors(model, episode, device, dtype):
    """Move one materialized episode onto the device, batched."""
    def batched(value):
        return value.unsqueeze(0).to(device, dtype)
    return {
        "protein_pooled": batched(episode.protein_pooled),
        "protein_tokens": batched(episode.protein_tokens),
        "protein_mask": batched(episode.protein_mask),
        "protein_chemistry": batched(episode.protein_chemistry),
        "support_atoms": batched(episode.support_atoms),
        "support_bonds": batched(episode.support_bonds),
        "support_mask": batched(episode.support_mask),
        "query_atoms": batched(episode.query_atoms),
        "query_bonds": batched(episode.query_bonds),
        "query_mask": batched(episode.query_mask),
        "support_y": batched(episode.support_y),
        "query_y": batched(episode.query_y),
        "support_fingerprint": batched(episode.support_fingerprint),
        "query_fingerprint": batched(episode.query_fingerprint),
    }


def align_atoms(parts: dict) -> dict:
    """Pad support and query to a common atom width, as `forward` does."""
    support_active = (int(parts["support_mask"].sum(-1).max())
                      if parts["support_mask"].numel() else 0)
    active = max(support_active, int(parts["query_mask"].sum(-1).max()))
    out = dict(parts)
    for side in ("support", "query"):
        atoms, bonds, mask = (parts[f"{side}_atoms"], parts[f"{side}_bonds"],
                              parts[f"{side}_mask"])
        atoms, bonds, mask = atoms[:, :, :active], bonds[:, :, :active, :active], mask[:, :, :active]
        missing = active - atoms.shape[-2]
        if missing > 0:
            atoms = F.pad(atoms, (0, 0, 0, missing))
            bonds = F.pad(bonds, (0, 0, 0, missing, 0, missing))
            mask = F.pad(mask, (0, missing))
        out[f"{side}_atoms"], out[f"{side}_bonds"], out[f"{side}_mask"] = atoms, bonds, mask
    return out


def encode_parts(model, parts: dict):
    """The expensive half, exposed so counterfactuals can reuse it.

    A support-label counterfactual changes only the inner loop and the
    transport residual — never the encoder — so re-encoding for it would be
    pure waste and would also risk the control differing from the real arm in
    something other than the labels.
    """
    return encode_task(
        model, parts["protein_pooled"], parts["protein_tokens"],
        parts["protein_mask"], parts["support_atoms"], parts["support_bonds"],
        parts["support_mask"], parts["query_atoms"], parts["query_bonds"],
        parts["query_mask"], parts["protein_chemistry"])


def predict(model, parts: dict, config: AdaptationConfig, *,
            support_y_override: torch.Tensor | None = None,
            keep: str = "both", return_parts: bool = False,
            task=None):
    """Encode once, adapt on support, read out on query, then transport.

    `support_y_override` is how the counterfactual controls are run: permuted
    or matched-wrong support labels flow into both the inner loop and the
    transport, exactly as the correct labels do, so the control differs from
    the real arm in the labels alone.
    """
    from tools.research.stageA_innerloop.inner_loop import partial_weights

    if task is None:
        task = encode_parts(model, parts)
    support_y = (parts["support_y"] if support_y_override is None
                 else support_y_override)
    weights = base_weights(model)

    pre_support = readout(model, task.support_hidden, task.support_additive,
                          task.support_occupancy, weights)
    pre_query = readout(model, task.query_hidden, task.query_additive,
                        task.query_occupancy, weights)

    fast, inner_trace = adapt(model, task, support_y, config, weights)
    if keep != "both":
        fast = partial_weights(model, fast, keep)
    post_support = readout(model, task.support_hidden, task.support_additive,
                           task.support_occupancy, fast)
    post_query = readout(model, task.query_hidden, task.query_additive,
                         task.query_occupancy, fast)

    support_count = task.support_hidden.shape[1]
    if support_count == 0:
        prediction = post_query
        transport = torch.zeros_like(post_query)
    else:
        similarity = tanimoto(parts["query_fingerprint"],
                              parts["support_fingerprint"])
        locked = (support_y - post_support).detach()
        shrink = model.transport.shrinkage(support_count, locked)
        transport, _ = model.transport(
            task.support_embed, task.query_embed, locked, similarity)
        transport = shrink * transport
        prediction = post_query + transport
    output = {"prediction": prediction, "zero_shot": post_query,
              "pre_adaptation_query": pre_query,
              "pre_adaptation_support": pre_support,
              "post_adaptation_support": post_support,
              "transport": transport, "inner_trace": inner_trace,
              "fast_weights": fast, "task": task,
              "query_occupancy": task.query_occupancy}
    if return_parts:
        output["task"] = task
    return output


def counterfactual_supports(support_y: torch.Tensor,
                            support_prediction: torch.Tensor) -> list[torch.Tensor]:
    """Wrong-label bindings, defined exactly as the accepted recipe defines them.

    `scripts/train_qpsmp.py::counterfactual_label_assignments` enumerates
    non-identity cyclic rolls for k > 1, and for k = 1 uses an equal-magnitude
    residual flip — because a one-element permutation is the identity and would
    make the control a copy of the real arm rather than a counterfactual.

    Reimplemented here rather than imported because the production helper reads
    a `QPSMPMetaOutput`, which the adapted path does not build; the arithmetic
    is identical and a test pins it against the original.
    """
    count = support_y.shape[-1]
    if count == 0:
        return []
    if count == 1:
        raw_residual = (support_y - support_prediction).detach()
        return [support_y - 2.0 * raw_residual]
    return [support_y.roll(shift, dims=-1) for shift in range(1, count)]


def auxiliary_losses(model, parts: dict, task, output, config, base,
                     label_scale):
    """The accepted recipe's three auxiliary terms, identical in every arm.

    `A0` is supposed to *be* the accepted baseline, so dropping these would
    quietly make it a different recipe and turn any A1 gain into "the inner
    loop recovers what we deleted". They are therefore computed for all three
    arms, and computed the same way, so the inner loop remains the only
    difference between arms.

    Two of the three are free or nearly free on the cached encoding:

    * `support_match` is the contact-dictionary regularizer, already a function
      of `query_occupancy`;
    * `binding` needs only the readout and transport re-run under
      counterfactual support labels, because a label change never reaches the
      encoder.

    The protein contrast genuinely needs a second encoder pass on the donor
    protein, and is the only added forward.
    """
    support_count = task.support_hidden.shape[1]
    query_y = parts["query_y"]
    support_match = model.dictionary_regularizer(task.query_occupancy)
    binding = query_y.new_zeros(())
    protein = query_y.new_zeros(())
    if support_count > 0:
        errors = [(output["prediction"] - query_y).square().mean()]
        for assignment in counterfactual_supports(
                parts["support_y"], output["post_adaptation_support"]):
            other = predict(model, parts, config, task=task,
                            support_y_override=assignment)
            errors.append((other["prediction"] - query_y).square().mean())
        binding = binding_contrastive_loss(errors, base.binding_temperature)

        donor = dict(parts)
        donor.update(_protein_inputs_cache[parts["_donor_key"]])
        donor_task = encode_parts(model, donor)
        donor_zero = readout(model, donor_task.query_hidden,
                             donor_task.query_additive,
                             donor_task.query_occupancy, base_weights(model))
        correct_error = (output["pre_adaptation_query"] - query_y).square().mean()
        wrong_error = (donor_zero - query_y).square().mean()
        protein = binding_contrastive_loss([correct_error, wrong_error],
                                           base.binding_temperature)
    return support_match, binding, protein


# Donor protein tensors, cached per target so the protein-contrast term costs
# one encoder pass rather than one bank lookup plus a pass.
_protein_inputs_cache: dict[str, dict] = {}


def cache_donor(data, target: str, device, dtype) -> str:
    if target not in _protein_inputs_cache:
        pooled, tokens, mask = data.protein_for_target(target)
        chemistry = data.protein_chemistry_for_target(target)
        _protein_inputs_cache[target] = {
            "protein_pooled": pooled.unsqueeze(0).to(device, dtype),
            "protein_tokens": tokens.unsqueeze(0).to(device, dtype),
            "protein_mask": mask.unsqueeze(0).to(device, dtype),
            "protein_chemistry": chemistry.unsqueeze(0).to(device, dtype)}
    return target


def task_value(model, parts: dict, config: AdaptationConfig) -> dict:
    """A2's score for one candidate meta_train task. No gradient is kept.

    Both terms are computed at the shared initialization on a candidate drawn
    from `meta_train`. The query label is read here to score the task for
    sampling; it never reaches a model input or an inference path.
    """
    task = encode_task(
        model, parts["protein_pooled"], parts["protein_tokens"],
        parts["protein_mask"], parts["support_atoms"], parts["support_bonds"],
        parts["support_mask"], parts["query_atoms"], parts["query_bonds"],
        parts["query_mask"], parts["protein_chemistry"])
    weights = base_weights(model)
    fast, _ = adapt(model, task, parts["support_y"], config, weights)
    with torch.no_grad():
        post = readout(model, task.query_hidden, task.query_additive,
                       task.query_occupancy, fast)
        loss = float(F.mse_loss(post, parts["query_y"]))
    cosine = support_query_gradient_cosine(
        model, task, parts["support_y"], parts["query_y"])
    return {"post_adaptation_query_loss": loss, "gradient_cosine": cosine}


def standardize(values: list[float]) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    spread = float(array.std())
    if not np.isfinite(spread) or spread < 1e-12:
        return np.zeros_like(array)
    return (array - array.mean()) / spread


def select_tasks(scores: list[dict], keep: int) -> tuple[list[int], dict]:
    """Prefer low post-adaptation query loss and high gradient agreement."""
    value = (standardize([-s["post_adaptation_query_loss"] for s in scores])
             + standardize([s["gradient_cosine"] for s in scores]))
    order = np.argsort(-value)[:keep]
    weight = np.exp(value - value.max())
    weight = weight / weight.sum()
    effective = float(1.0 / np.square(weight).sum())
    return [int(i) for i in order], {
        "value_mean": float(value.mean()), "value_std": float(value.std()),
        "effective_tasks": effective, "candidates": len(scores),
        "selected_cosine_mean": float(np.mean(
            [scores[int(i)]["gradient_cosine"] for i in order])),
        "all_cosine_mean": float(np.mean(
            [s["gradient_cosine"] for s in scores])),
    }


def train(data: QPSMPData, meta: MetaConfig, progress_path: Path | None = None):
    config = meta.base
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = resolve_architecture(config.arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks,
        adapter_scale=config.adapter_scale, use_cartesian=config.use_cartesian,
        dtype=torch.float32).to(config.device)

    trainable = [p for p in model.parameters() if p.requires_grad]
    fast_prefixes = ("meta.term.", "transport.")
    fast_parameters = [p for n, p in model.named_parameters()
                       if p.requires_grad and n.startswith(fast_prefixes)]
    fast_ids = {id(p) for p in fast_parameters}
    slow_parameters = [p for p in trainable if id(p) not in fast_ids]
    optimizer = torch.optim.AdamW([
        {"params": slow_parameters,
         "lr": config.learning_rate * config.backbone_lr_scale},
        {"params": fast_parameters, "lr": config.learning_rate},
    ], weight_decay=config.weight_decay)
    base_lrs = [group["lr"] for group in optimizer.param_groups]
    label_scale = training_label_scale(data)
    support_sizes = (0, 1, 2, 3, 5)
    adaptation = meta.adaptation()

    val_specs = data.fixed_nested_episode_banks(
        "meta_val", support_sizes, config.query_size,
        config.val_draws_per_target, config.evaluation_seed,
        config.eval_targets_per_component)
    val_banks = {k: tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in specs)
        for k, specs in val_specs.items()}

    device, dtype = config.device, torch.float32
    best_state, best_value, best_step = None, float("inf"), 0
    trace, selection_log, cosine_log = [], [], []
    forwards = 0
    started = time.monotonic()

    for step in range(1, config.steps + 1):
        model.train()
        factor = learning_rate_factor(step, config)
        for group, base in zip(optimizer.param_groups, base_lrs):
            group["lr"] = base * factor
        optimizer.zero_grad(set_to_none=True)
        support_size = support_sizes[(step - 1) % len(support_sizes)]

        wanted = config.episodes_per_step
        draw_count = (wanted * meta.candidate_multiplier
                      if meta.task_selection else wanted)
        candidates = []
        for _ in range(draw_count):
            requested = int(rng.integers(config.min_query_size,
                                         config.query_size + 1))
            spec = data.draw_episode("meta_train", support_size, requested, rng)
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            parts = align_atoms(episode_tensors(model, episode, device, dtype))
            # The donor is drawn by `draw_episode` from a *different* homology
            # component of meta_train, so the protein contrast never sees a
            # meta_val protein.
            parts["_donor_key"] = cache_donor(data, spec.donor_target,
                                              device, dtype)
            candidates.append(parts)

        if meta.task_selection and support_size > 0:
            scores = [task_value(model, parts, adaptation)
                      for parts in candidates]
            forwards += len(candidates)
            picked, summary = select_tasks(scores, wanted)
            summary["step"] = step
            selection_log.append(summary)
            cosine_log.extend(s["gradient_cosine"] for s in scores)
            selected = [candidates[i] for i in picked]
        else:
            selected = candidates[:wanted]

        loss_value = 0.0
        for parts in selected:
            forwards += 1
            task = encode_parts(model, parts)
            output = predict(model, parts, adaptation, task=task)
            query_y = parts["query_y"]
            loss_post = F.smooth_l1_loss(output["prediction"], query_y)
            loss_pre = F.smooth_l1_loss(output["pre_adaptation_query"], query_y)
            loss_rank = ranking_term(output["prediction"], query_y, config,
                                     label_scale)
            shape = centered_task_error(output["prediction"], query_y)
            support_match, binding, protein = auxiliary_losses(
                model, parts, task, output, adaptation, config, label_scale)
            if parts["support_y"].shape[-1] > 0:
                forwards += 1                      # the donor encoder pass
            episode_loss = (
                meta.post_adaptation_loss_weight * loss_post
                + config.zero_shot_loss_weight * loss_pre
                + config.ranking_loss_weight * loss_rank
                + config.shape_loss_weight * shape
                + config.support_match_loss_weight * support_match
                + config.binding_loss_weight * binding
                + config.protein_contrast_loss_weight * protein)
            loss_value += float(episode_loss.detach())
            (episode_loss / len(selected)).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        optimizer.step()
        trace.append(loss_value / len(selected))

        if step % config.val_interval == 0 or step == config.steps:
            metrics = {k: evaluate(model, data, bank, controls=k > 0,
                                   label_scale=label_scale)
                       for k, bank in val_banks.items()}
            values = [(admission_score(item, config.admission_binding_margin_pk)
                       if k > 0 else item["full_mse_pk"])
                      for k, item in metrics.items()]
            value = float(np.mean(values))
            progress = {"step": step, "validation_admission_score": value,
                        "best": min(best_value, value),
                        "elapsed_seconds": time.monotonic() - started}
            print(json.dumps(progress), flush=True)
            if progress_path is not None:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(progress) + "\n")
            if value < best_value:
                best_state = copy.deepcopy(model.state_dict())
                best_value, best_step = value, step

    if best_state is None:
        raise RuntimeError("training produced no validation checkpoint")
    model.load_state_dict(best_state)
    return model, {
        "best_val_admission_score": best_value, "best_step": best_step,
        "loss_trace": trace, "label_scale": asdict(label_scale),
        "adaptation": adaptation.to_dict(),
        "task_selection": meta.task_selection,
        "candidate_multiplier": meta.candidate_multiplier,
        "encoder_forward_passes": forwards,
        "optimization_steps": config.steps,
        "selection_log": selection_log,
        "gradient_cosine_samples": cosine_log,
        "wall_time_seconds": time.monotonic() - started,
        "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                if device.startswith("cuda") else 0.0),
        "trainable_parameters": int(sum(p.numel() for p in trainable)),
    }, label_scale


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--inner-steps", type=int, default=0)
    parser.add_argument("--inner-lr", type=float, default=0.1)
    parser.add_argument("--task-selection", action="store_true")
    parser.add_argument("--candidate-multiplier", type=int, default=3)
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing output directory")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    meta = MetaConfig(
        base=TrainConfig(arch="similarity_only", steps=arguments.steps,
                         seed=arguments.seed, split_directory=str(SPLIT),
                         device=arguments.device, amp=False),
        inner_steps=arguments.inner_steps, inner_lr=arguments.inner_lr,
        task_selection=arguments.task_selection,
        candidate_multiplier=arguments.candidate_multiplier)
    # A second run writing into a live output directory silently interleaves
    # two models' progress and checkpoints. That happened once here, when a
    # killed shell wrapper walked on to the next arm while a fresh launch was
    # already using the same path, so the guard is a fail-closed default.
    progress_path = arguments.output / "progress.jsonl"
    if progress_path.exists() and not arguments.force:
        raise SystemExit(
            f"{progress_path} already exists — another run may be writing "
            f"here. Remove the directory or pass --force.")
    arguments.output.mkdir(parents=True, exist_ok=True)
    model, report, label_scale = train(data, meta, progress_path=progress_path)
    torch.save({"model_state": model.state_dict(),
                "config": {**asdict(meta.base), "arm": arguments.arm,
                           **meta.adaptation().to_dict(),
                           "task_selection": meta.task_selection},
                "adaptation": meta.adaptation().to_dict(),
                "label_scale": asdict(label_scale)},
               arguments.output / "checkpoint.pt")
    (arguments.output / "RESULT.json").write_text(
        json.dumps({"arm": arguments.arm, "report": report,
                    "config": {**asdict(meta.base), **meta.adaptation().to_dict(),
                               "task_selection": meta.task_selection},
                    "meta_test": data.seal_record()}, indent=1),
        encoding="utf-8")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
