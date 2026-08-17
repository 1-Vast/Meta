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

from scripts import internal_validation                              # noqa: E402
# Re-exported for the stage lineage: Stages D/E, F, I, J, K, L and Q all
# import `draw_fit_episode`, `internal_validation_bank` and
# `partition_components` from this module.
from scripts.internal_validation import (                            # noqa: E402,F401
    draw_fit_episode, eligible_targets, internal_validation_specs,
    partition_components,
)
from scripts.qpsmp_data import QPSMPData                              # noqa: E402
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

# The rule below was promoted verbatim into `scripts/internal_validation.py` on
# 2026-08-18 and is now the maintained trainer's default. Stage B imports it
# from there so the two cannot drift; the constants and bank seeds are
# unchanged, so a rebuilt Stage B bank is bit-identical to the recorded one.
INTERNAL_VAL_FRACTION = internal_validation.INTERNAL_VAL_FRACTION
PARTITION_SEED = internal_validation.PARTITION_SEED
INTERNAL_BANK_SEED = internal_validation.INTERNAL_BANK_SEED
SUPPORT_SIZES = internal_validation.SUPPORT_SIZES


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


def internal_validation_bank(data: QPSMPData, components: tuple[str, ...],
                             label_scale, query_size: int = 16,
                             draws: int = 1) -> dict[int, tuple]:
    """Materialize the promoted internal-validation specs for Stage B's arms."""
    specs = internal_validation_specs(data, components, query_size, draws,
                                      SUPPORT_SIZES)
    return {size: tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in group)
        for size, group in specs.items()}


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
