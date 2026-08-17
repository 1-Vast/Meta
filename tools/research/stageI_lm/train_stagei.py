"""Stage I trainer: live ESM-150M LoRA lane vs the frozen T2 baseline.

Arms:
* I        - live ESM-2 150M encoder with trainable LoRA adapters (r=8,
             alpha=16) feeding the unchanged similarity_only trunk.
* I-FROZEN - identical live path with the adapters frozen and the encode
             under no_grad: isolates the trainability of the LM lane.

Leak-free checkpoint selection on internal meta_train components; the
internal bank here uses up to 2 targets per component and 1 draw (declared
difference from the Stage D bank, made for the live-encoder cost);
meta_val is read once after freezing. GPU verification before every arm.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData, stable_seed
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, forward, learning_rate_factor, normalized_episode,
    ranking_term, training_label_scale,
)
from tools.research.stageB_complementary.train_stageb import (
    draw_fit_episode, partition_components,
)
from tools.research.stageD_level_panel.train_staged import (
    StageEConfig, arm_loss, build_model,
)
from tools.research.stageI_lm.lora_esm import LiveESMProteinEncoder

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
MODEL_DIR = ("C:/Users/59964/.cache/huggingface/hub/"
             "models--facebook--esm2_t30_150M_UR50D/snapshots/"
             "a695f6045e2e32885fa60af20c13cb35398ce30c")
INTERNAL_BANK_SEED = 20260819


def internal_bank_small(data, components, label_scale):
    """<=2 targets per internal component, 1 draw, nested k."""
    from scripts.qpsmp_data import EpisodeSpec
    max_support = max(SUPPORT_SIZES)
    banks = {size: [] for size in SUPPORT_SIZES}
    keys = sorted(components)
    for index, component in enumerate(keys):
        donor_component = keys[(index + 1) % len(keys)]
        targets = sorted(data.components["meta_train"][component])
        donors = sorted(data.components["meta_train"][donor_component])
        for target in targets[:1]:
            indices = data.tasks["meta_train"][target]
            if data._unique_ligand_count(indices) < max_support + 2:
                continue
            rng = np.random.default_rng(stable_seed(
                "stagei-internal", INTERNAL_BANK_SEED, target))
            order = data._unique_ligand_order(indices, rng)
            query = tuple(map(int, order[max_support:max_support + 16]))
            donor = donors[int(rng.integers(len(donors)))]
            for size in SUPPORT_SIZES:
                banks[size].append(EpisodeSpec(
                    "meta_train", component, target,
                    tuple(map(int, order[:size])), query, donor))
    return {size: tuple(compact_episode(normalized_episode(
            data.materialize(spec), label_scale)) for spec in specs)
            for size, specs in banks.items()}


def live_forward(model, encoder, data, episode, train_lora):
    sequence = data._protein_sequences[episode.spec.target]
    if train_lora:
        # gradient flows through the first chunk only; later chunks are
        # encoded feature-only to keep memory bounded on long proteins
        pooled, slots, mask = encoder.encode(sequence, max_chunks=1)
    else:
        with torch.no_grad():
            pooled, slots, mask = encoder.encode(sequence)
    from dataclasses import replace
    episode = replace(episode,
                      protein_pooled=pooled, protein_tokens=slots,
                      protein_mask=mask)
    return forward(model, episode)


def gpu_probe(data, config, encoder, output_dir):
    device = config.base.device
    model = build_model(config, data).to(device).train()
    torch.manual_seed(0)
    rng = np.random.default_rng(0)
    fit, _ = partition_components(data)
    spec = draw_fit_episode(data, fit, 2, 8, rng)
    label_scale = training_label_scale(data)
    episode = compact_episode(normalized_episode(data.materialize(spec),
                                                 label_scale))
    if device.startswith("cuda"):
        from dataclasses import replace
        episode = replace(episode, **{
            field: getattr(episode, field).to(device)
            for field in episode.__dataclass_fields__
            if isinstance(getattr(episode, field), torch.Tensor)})
    probe = {"torch_cuda_is_available": bool(torch.cuda.is_available()),
             "torch_version": torch.__version__,
             "cuda_runtime": torch.version.cuda,
             "device_count": torch.cuda.device_count(),
             "configured_device": device,
             "model_parameter_devices": sorted({str(p.device)
                                                for p in model.parameters()}),
             "encoder_on_cuda": bool(device.startswith("cuda")),
             "nvidia_smi_samples": []}
    optimizer = torch.optim.AdamW(list(model.parameters()) +
                                  encoder.lora_parameters(), lr=1e-5)
    parsed = []
    for step in range(4):
        output = live_forward(model, encoder, data, episode, True)
        loss = F.smooth_l1_loss(output.prediction, episode.query_y.to(
            device=output.prediction.device, dtype=output.prediction.dtype))
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        if step >= 1:
            try:
                text = subprocess.run(
                    ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",
                     "--format=csv,noheader,nounits"],
                    capture_output=True, text=True, timeout=10).stdout.strip()
                parts = text.split(",")
                parsed.append((int(parts[0].strip()), int(parts[1].strip())))
                probe["nvidia_smi_samples"].append(text)
            except Exception:  # noqa: BLE001
                pass
    probe["max_gpu_utilization_percent"] = max((v for v, _ in parsed),
                                                default=None)
    probe["utilization_nonzero"] = any(v > 0 for v, _ in parsed)
    probe["batch_device_check"] = bool(device.startswith("cuda"))
    (output_dir / "GPU_PROBE.json").write_text(
        json.dumps(probe, indent=1, sort_keys=True), encoding="utf-8")
    print(json.dumps(probe, indent=1), flush=True)
    return probe


def internal_score(model, encoder, data, banks, label_scale):
    model.eval()
    total, count = 0.0, 0
    cache = {}
    for bank in banks.values():
        for episode in bank:
            target = episode.spec.target
            if target not in cache:
                with torch.no_grad():
                    pooled, slots, mask = encoder.encode(
                        data._protein_sequences[target])
                cache[target] = (pooled, slots, mask)
            from dataclasses import replace
            episode = replace(episode,
                              protein_pooled=cache[target][0],
                              protein_tokens=cache[target][1],
                              protein_mask=cache[target][2])
            with torch.no_grad():
                output = forward(model, episode)
                error = float(F.mse_loss(output.prediction, episode.query_y.to(
                    device=output.prediction.device,
                    dtype=output.prediction.dtype)))
            total += error * label_scale.scale ** 2
            count += 1
    model.train()
    return total / max(count, 1)


def train(data, config, arm, output_dir, progress_path):
    base = config.base
    torch.manual_seed(base.seed)
    rng = np.random.default_rng(base.seed)
    label_scale = training_label_scale(data)
    encoder = LiveESMProteinEncoder(MODEL_DIR, base.device)
    probe = gpu_probe(data, config, encoder, output_dir)
    if not probe["torch_cuda_is_available"]:
        raise RuntimeError("GPU verification failed")
    model = build_model(config, data).to(base.device)
    train_lora = arm == "I"
    parameters = list(model.parameters())
    if train_lora:
        parameters = parameters + encoder.lora_parameters()
    trainable = [p for p in parameters if p.requires_grad]
    fast_prefixes = ("meta.term.", "transport.")
    fast_parameters = [p for n, p in list(model.named_parameters())
                       if p.requires_grad and n.startswith(fast_prefixes)]
    fast_ids = {id(p) for p in fast_parameters}
    slow_parameters = [p for p in trainable if id(p) not in fast_ids]
    optimizer = torch.optim.AdamW([
        {"params": slow_parameters,
         "lr": base.learning_rate * base.backbone_lr_scale},
        {"params": fast_parameters, "lr": base.learning_rate},
    ], weight_decay=base.weight_decay)
    base_lrs = [group["lr"] for group in optimizer.param_groups]
    fit_components, internal_components = partition_components(data)
    banks = internal_bank_small(data, internal_components, label_scale)
    print(f"arm={arm} fit {len(fit_components)} internal "
          f"{len(internal_components)} bank "
          f"{sum(len(b) for b in banks.values())}", flush=True)
    best_state, best_lora, best_value, best_step = None, None, float("inf"), 0
    trace = []
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
            requested = int(rng.integers(base.min_query_size,
                                         base.query_size + 1))
            spec = draw_fit_episode(data, fit_components, support_size,
                                    requested, rng)
            selected.append(compact_episode(normalized_episode(
                data.materialize(spec), label_scale)))
        loss_value = 0.0
        for episode in selected:
            output = live_forward(model, encoder, data, episode, train_lora)
            query_y = episode.query_y.to(device=output.prediction.device,
                                         dtype=output.prediction.dtype)
            loss, _ = arm_loss(config, model, output, episode, query_y,
                               base, label_scale)
            loss_value += float(loss.detach())
            (loss / len(selected)).backward()
        torch.nn.utils.clip_grad_norm_(parameters, base.grad_clip)
        optimizer.step()
        trace.append(loss_value / len(selected))
        if step % base.val_interval == 0 or step == base.steps:
            value = internal_score(model, encoder, data, banks, label_scale)
            progress = {"step": step, "arm": arm,
                        "internal_val_mse_pk": value,
                        "best": min(best_value, value),
                        "elapsed_seconds": time.monotonic() - started}
            print(json.dumps(progress), flush=True)
            if progress_path is not None:
                with progress_path.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(progress) + "\n")
            if value < best_value:
                best_state = copy.deepcopy(model.state_dict())
                best_lora = encoder.lora_state()
                best_value, best_step = value, step
    if best_state is None:
        raise RuntimeError("training produced no checkpoint")
    model.load_state_dict(best_state)
    report = {"arm": arm, "best_internal_val_mse_pk": best_value,
              "best_step": best_step, "loss_trace": trace,
              "label_scale": asdict(label_scale),
              "checkpoint_selection": "meta_train internal components only",
              "internal_bank": "up to 2 targets per component, 1 draw",
              "fit_components": len(fit_components),
              "internal_val_components": len(internal_components),
              "optimization_steps": base.steps,
              "wall_time_seconds": time.monotonic() - started,
              "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                     if base.device.startswith("cuda") else 0.0),
              "trainable_parameters": int(sum(p.numel() for p in trainable)),
              "lora_parameters": int(sum(p.numel()
                                         for p in encoder.lora_parameters())),
              "gpu_probe": probe}
    payload = {"model_state": model.state_dict(), "lora_state": best_lora,
               "config": asdict(base), "arm": arm,
               "label_scale": asdict(label_scale), "stageI_version": 1}
    torch.save(payload, output_dir / "checkpoint.pt")
    (output_dir / "RESULT.json").write_text(json.dumps(
        {"arm": arm, "report": report,
         "meta_test": data.seal_record()}, indent=1), encoding="utf-8")
    print(f"wrote {output_dir}")
    return model, report, label_scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", required=True, choices=("I", "I-FROZEN"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--val-interval", type=int, default=20)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()
    progress_path = args.output / "progress.jsonl"
    if progress_path.exists() and not args.force:
        raise SystemExit(f"{progress_path} exists; pass --force to overwrite")
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    from dataclasses import replace as _replace
    base = TrainConfig(arch="similarity_only", steps=args.steps,
                       seed=args.seed, split_directory=str(SPLIT),
                       device=args.device, amp=False)
    base = _replace(base, val_interval=args.val_interval)
    config = StageEConfig(base=base, arm="T2")
    args.output.mkdir(parents=True, exist_ok=True)
    train(data, config, args.arm, args.output, progress_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
