"""Stage B trainer: four arms, and checkpoint selection that never reads meta_val.

Stage A selected checkpoints on the same `meta_val` population it then reported,
which made every reported figure an optimistic development estimate. Here
`meta_train`'s homology components are partitioned once, deterministically, into
a **fit** set and an **internal validation** set. Training episodes come only
from fit components; checkpoint selection reads only internal-validation
episodes; `meta_val` is read exactly once, after the candidate is frozen.

The partition is by component, never by target: two targets of one homology
component are not independent evidence about an unseen component, so a
target-level split would leak.

Arms are `T` / `M` / `H` / `C` (see `arms.py`). Everything else — seed,
architecture, capacity, optimizer, schedule, episode banks, losses, query panels
— is identical across arms.
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

from scripts.qpsmp_data import EpisodeSpec, QPSMPData, stable_seed    # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    centered_task_error, compact_episode, learning_rate_factor,
    normalized_episode, ranking_term, resolve_architecture,
    training_label_scale,
)
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)
from tools.research.stageB_complementary.arms import (                # noqa: E402
    InnerStepSizes, MODES, StageBAdaptation, predict,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"

# Frozen before any Stage B training.
INTERNAL_VAL_FRACTION = 0.12
PARTITION_SEED = 20260818
INTERNAL_BANK_SEED = 20260819
SUPPORT_SIZES = (0, 1, 2, 3, 5)


@dataclass(frozen=True)
class StageBConfig:
    base: TrainConfig = field(default_factory=lambda: TrainConfig(
        arch="similarity_only", steps=1200, seed=20260815,
        split_directory=str(SPLIT), amp=False))
    mode: str = "T"
    inner_steps: int = 1
    inner_lr: float = 0.1
    learned_step: bool = False
    max_step: float = 0.5
    post_adaptation_loss_weight: float = 1.0
    # "internal" is the leak-free default and the only setting any reported
    # arm uses. "meta_val" exists solely for the leakage diagnostic: run with
    # the same fit components, it isolates how much of Stage A's advantage came
    # from selecting checkpoints on the population it then reported.
    selection: str = "internal"

    def adaptation(self) -> StageBAdaptation:
        return StageBAdaptation(
            mode=self.mode, inner_steps=self.inner_steps,
            inner_lr=self.inner_lr, learned_step=self.learned_step,
            max_step=self.max_step)


def partition_components(data: QPSMPData) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Split meta_train components into fit / internal-validation, by component."""
    components = sorted(data.components["meta_train"])
    order = np.random.default_rng(PARTITION_SEED).permutation(len(components))
    held = max(2, int(round(INTERNAL_VAL_FRACTION * len(components))))
    internal = tuple(sorted(components[int(i)] for i in order[:held]))
    fit = tuple(sorted(components[int(i)] for i in order[held:]))
    return fit, internal


def eligible_targets(data: QPSMPData, components: tuple[str, ...],
                     needed: int) -> dict[str, tuple[str, ...]]:
    out: dict[str, tuple[str, ...]] = {}
    for component in components:
        targets = tuple(
            target for target in data.components["meta_train"][component]
            if data._unique_ligand_count(data.tasks["meta_train"][target]) >= needed)
        if targets:
            out[component] = targets
    return out


def draw_fit_episode(data: QPSMPData, components: tuple[str, ...],
                     support_size: int, query_size: int,
                     rng: np.random.Generator) -> EpisodeSpec:
    """`draw_episode` restricted to the fit components.

    Component-uniform then target-uniform, matching the incumbent sampler, so
    the only difference from production sampling is the withheld internal-
    validation components.
    """
    pool = eligible_targets(data, components, support_size + 1)
    if not pool:
        raise ValueError("no fit task can provide the requested episode")
    keys = sorted(pool)
    component = keys[int(rng.integers(len(keys)))]
    targets = pool[component]
    target = targets[int(rng.integers(len(targets)))]
    order = data._unique_ligand_order(data.tasks["meta_train"][target], rng)
    support = order[:support_size]
    query = order[support_size:support_size + min(query_size,
                                                  len(order) - support_size)]
    donors = [key for key in keys if key != component]
    donor_component = donors[int(rng.integers(len(donors)))]
    donor = pool[donor_component][int(rng.integers(len(pool[donor_component])))]
    return EpisodeSpec("meta_train", component, target, tuple(map(int, support)),
                       tuple(map(int, query)), donor)


def internal_validation_bank(data: QPSMPData, components: tuple[str, ...],
                             label_scale, query_size: int = 16,
                             draws: int = 1) -> dict[int, tuple]:
    """A fixed nested bank over the withheld meta_train components."""
    max_support = max(SUPPORT_SIZES)
    pool = eligible_targets(data, components, max_support + 1)
    keys = sorted(pool)
    if len(keys) < 2:
        raise ValueError("internal validation needs at least two components")
    banks: dict[int, list] = {size: [] for size in SUPPORT_SIZES}
    for index, component in enumerate(keys):
        donor_component = keys[(index + 1) % len(keys)]
        for target in pool[component]:
            for draw in range(draws):
                rng = np.random.default_rng(stable_seed(
                    "stageb-internal", INTERNAL_BANK_SEED, target, draw))
                order = data._unique_ligand_order(
                    data.tasks["meta_train"][target], rng)
                query = tuple(map(int, order[
                    max_support:max_support + min(
                        query_size, len(order) - max_support)]))
                if len(query) < 2:
                    continue
                donor_targets = pool[donor_component]
                donor = donor_targets[int(rng.integers(len(donor_targets)))]
                for size in SUPPORT_SIZES:
                    banks[size].append(EpisodeSpec(
                        "meta_train", component, target,
                        tuple(map(int, order[:size])), query, donor))
    return {size: tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in specs)
        for size, specs in banks.items()}


def internal_score(model, banks, adaptation, steps, device, label_scale) -> float:
    """Mean pK MSE over the internal-validation bank. meta_val is never read."""
    model.eval()
    total, count = 0.0, 0
    for size, bank in banks.items():
        for episode in bank:
            parts = align_atoms(episode_tensors(
                model, episode, device, torch.float32))
            with torch.no_grad():
                task = encode_parts(model, parts)
                output = predict(model, parts, task, adaptation, steps)
                error = float(F.mse_loss(output["prediction"],
                                         parts["query_y"]))
            total += error * label_scale.scale ** 2
            count += 1
    model.train()
    return total / max(count, 1)


def train(data: QPSMPData, config: StageBConfig, progress_path: Path | None = None):
    base = config.base
    torch.manual_seed(base.seed)
    rng = np.random.default_rng(base.seed)
    model = resolve_architecture(base.arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=base.hidden_dim, task_dim=base.task_dim,
        ligand_layers=base.ligand_layers, pair_dim=base.pair_dim,
        pair_blocks=base.pair_blocks, pair_latents=base.pair_latents,
        pair_heads=base.pair_heads, pair_chunk_size=base.pair_chunk_size,
        support_hidden_dim=base.support_hidden_dim,
        support_blocks=base.support_blocks, adapter_rank=base.adapter_rank,
        adaptive_blocks=base.adaptive_blocks, adapter_scale=base.adapter_scale,
        dtype=torch.float32).to(base.device)

    adaptation = config.adaptation()
    steps_module = None
    extra_parameters: list = []
    if config.learned_step and config.mode != "T":
        steps_module = InnerStepSizes(config.inner_lr, config.max_step).to(base.device)
        extra_parameters = list(steps_module.parameters())

    trainable = [p for p in model.parameters() if p.requires_grad]
    fast_prefixes = ("meta.term.", "transport.")
    fast_parameters = [p for n, p in model.named_parameters()
                       if p.requires_grad and n.startswith(fast_prefixes)]
    fast_ids = {id(p) for p in fast_parameters}
    slow_parameters = [p for p in trainable if id(p) not in fast_ids]
    groups = [
        {"params": slow_parameters,
         "lr": base.learning_rate * base.backbone_lr_scale},
        {"params": fast_parameters, "lr": base.learning_rate},
    ]
    if extra_parameters:
        groups.append({"params": extra_parameters, "lr": base.learning_rate})
    optimizer = torch.optim.AdamW(groups, weight_decay=base.weight_decay)
    base_lrs = [group["lr"] for group in optimizer.param_groups]
    label_scale = training_label_scale(data)

    fit_components, internal_components = partition_components(data)
    if config.selection == "internal":
        banks = internal_validation_bank(data, internal_components, label_scale)
        selection_note = "meta_train internal-validation components only"
    elif config.selection == "meta_val":
        # LEAKAGE DIAGNOSTIC ONLY. Reproduces Stage A's selection rule on the
        # same fit components, so the contrast against the `internal` arm
        # isolates the selection effect from the smaller training set.
        specs = data.fixed_nested_episode_banks(
            "meta_val", SUPPORT_SIZES, 16, 1, base.evaluation_seed,
            base.eval_targets_per_component)
        banks = {k: tuple(compact_episode(normalized_episode(
            data.materialize(spec), label_scale)) for spec in items)
            for k, items in specs.items()}
        selection_note = ("meta_val (LEAKAGE DIAGNOSTIC — not an admissible "
                          "arm; reproduces the Stage A rule)")
    else:
        raise ValueError(f"unknown selection rule {config.selection!r}")
    print(f"fit components {len(fit_components)}, "
          f"internal-validation components {len(internal_components)}, "
          f"selection={config.selection}, "
          f"bank {sum(len(b) for b in banks.values())} episodes")

    device, dtype = base.device, torch.float32
    best_state, best_steps, best_value, best_step = None, None, float("inf"), 0
    trace, forwards = [], 0
    started = time.monotonic()

    for step in range(1, base.steps + 1):
        model.train()
        factor = learning_rate_factor(step, base)
        for group, start in zip(optimizer.param_groups, base_lrs):
            group["lr"] = start * factor
        optimizer.zero_grad(set_to_none=True)
        support_size = SUPPORT_SIZES[(step - 1) % len(SUPPORT_SIZES)]

        selected = []
        for _ in range(base.episodes_per_step):
            requested = int(rng.integers(base.min_query_size, base.query_size + 1))
            spec = draw_fit_episode(data, fit_components, support_size,
                                    requested, rng)
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            selected.append(align_atoms(
                episode_tensors(model, episode, device, dtype)))

        loss_value = 0.0
        for parts in selected:
            forwards += 1
            task = encode_parts(model, parts)
            output = predict(model, parts, task, adaptation, steps_module)
            query_y = parts["query_y"]
            # Post-adaptation trains fast few-shot; pre-adaptation protects k=0.
            loss_post = F.smooth_l1_loss(output["prediction"], query_y)
            loss_pre = F.smooth_l1_loss(output["zero_shot"], query_y)
            loss_rank = ranking_term(output["prediction"], query_y, base,
                                     label_scale)
            shape = centered_task_error(output["prediction"], query_y)
            support_match = model.dictionary_regularizer(task.query_occupancy)
            episode_loss = (
                config.post_adaptation_loss_weight * loss_post
                + base.zero_shot_loss_weight * loss_pre
                + base.ranking_loss_weight * loss_rank
                + base.shape_loss_weight * shape
                + base.support_match_loss_weight * support_match)
            loss_value += float(episode_loss.detach())
            (episode_loss / len(selected)).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), base.grad_clip)
        optimizer.step()
        trace.append(loss_value / len(selected))

        if step % base.val_interval == 0 or step == base.steps:
            value = internal_score(model, banks, adaptation, steps_module,
                                   device, label_scale)
            progress = {"step": step, "internal_val_mse_pk": value,
                        "best": min(best_value, value),
                        "elapsed_seconds": time.monotonic() - started}
            if steps_module is not None:
                progress["weight_step"] = float(steps_module.weight_step())
                progress["bias_step"] = float(steps_module.bias_step())
            print(json.dumps(progress), flush=True)
            if progress_path is not None:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(progress) + "\n")
            if value < best_value:
                best_state = copy.deepcopy(model.state_dict())
                best_steps = (copy.deepcopy(steps_module.state_dict())
                              if steps_module is not None else None)
                best_value, best_step = value, step

    if best_state is None:
        raise RuntimeError("training produced no internal-validation checkpoint")
    model.load_state_dict(best_state)
    if steps_module is not None and best_steps is not None:
        steps_module.load_state_dict(best_steps)
    report = {
        "mode": config.mode,
        "best_internal_val_mse_pk": best_value, "best_step": best_step,
        "loss_trace": trace, "label_scale": asdict(label_scale),
        "adaptation": adaptation.to_dict(),
        "checkpoint_selection": selection_note,
        "selection_rule": config.selection,
        "fit_components": len(fit_components),
        "internal_val_components": len(internal_components),
        "internal_val_component_names": list(internal_components),
        "encoder_forward_passes": forwards,
        "optimization_steps": base.steps,
        "wall_time_seconds": time.monotonic() - started,
        "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                if device.startswith("cuda") else 0.0),
        "trainable_parameters": int(sum(p.numel() for p in trainable)),
    }
    if steps_module is not None:
        report["learned_steps"] = {
            "weight": float(steps_module.weight_step()),
            "bias": float(steps_module.bias_step()),
            "max_step": config.max_step,
            "initial": config.inner_lr}
    return model, steps_module, report, label_scale


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=MODES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--inner-steps", type=int, default=1)
    parser.add_argument("--inner-lr", type=float, default=0.1)
    parser.add_argument("--learned-step", action="store_true")
    parser.add_argument("--max-step", type=float, default=0.5)
    parser.add_argument("--selection", default="internal",
                        choices=("internal", "meta_val"),
                        help="meta_val is a leakage diagnostic only")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    progress_path = arguments.output / "progress.jsonl"
    if progress_path.exists() and not arguments.force:
        raise SystemExit(f"{progress_path} exists; pass --force to overwrite")

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    config = StageBConfig(
        base=TrainConfig(arch="similarity_only", steps=arguments.steps,
                         seed=arguments.seed, split_directory=str(SPLIT),
                         device=arguments.device, amp=False),
        mode=arguments.mode, inner_steps=arguments.inner_steps,
        inner_lr=arguments.inner_lr, learned_step=arguments.learned_step,
        max_step=arguments.max_step, selection=arguments.selection)
    arguments.output.mkdir(parents=True, exist_ok=True)
    model, steps_module, report, label_scale = train(
        data, config, progress_path=progress_path)
    payload = {"model_state": model.state_dict(),
               "config": {**asdict(config.base), **config.adaptation().to_dict()},
               "adaptation": config.adaptation().to_dict(),
               "label_scale": asdict(label_scale),
               "stageb_version": 1}
    if steps_module is not None:
        payload["inner_step_state"] = steps_module.state_dict()
    torch.save(payload, arguments.output / "checkpoint.pt")
    (arguments.output / "RESULT.json").write_text(
        json.dumps({"mode": arguments.mode, "report": report,
                    "config": {**asdict(config.base),
                               **config.adaptation().to_dict()},
                    "meta_test": data.seal_record()}, indent=1),
        encoding="utf-8")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
