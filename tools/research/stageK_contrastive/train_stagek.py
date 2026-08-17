"""Stage K trainer: contrastive coembedding branch vs frozen T2."""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict
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

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    centered_task_error, compact_episode, forward, learning_rate_factor,
    normalized_episode, ranking_term, training_label_scale,
)
from tools.research.stageB_complementary.train_stageb import (
    draw_fit_episode, internal_validation_bank, partition_components,
)
from tools.research.stageK_contrastive.model import CoembedModel

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
CONTRAST_WEIGHT = 0.5
TEMPERATURE = 0.1


def build_model(config, data):
    base = config.base
    return CoembedModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=base.hidden_dim, task_dim=base.task_dim,
        ligand_layers=base.ligand_layers, pair_dim=base.pair_dim,
        pair_blocks=base.pair_blocks, pair_latents=base.pair_latents,
        pair_heads=base.pair_heads, pair_chunk_size=base.pair_chunk_size,
        support_hidden_dim=base.support_hidden_dim,
        support_blocks=base.support_blocks, adapter_rank=base.adapter_rank,
        adaptive_blocks=base.adaptive_blocks,
        adapter_scale=base.adapter_scale,
        coembed_dim=128, dtype=torch.float32)


def contrastive_loss(model, temperature):
    zp, zl = model._stash
    zp_n = F.normalize(zp, dim=-1)
    zl_n = F.normalize(zl, dim=-1)
    logits = torch.einsum("bd,bqd->bq", zp_n, zl_n) / temperature
    labels = torch.zeros(zl.shape[0], dtype=torch.long, device=zl.device)
    return F.cross_entropy(logits, labels)


def regression_alignment(model):
    zp, zl = model._stash
    zp_n = F.normalize(zp, dim=-1)
    zl_n = F.normalize(zl, dim=-1)
    sim = torch.einsum("bd,bqd->bq", zp_n, zl_n)
    q = sim.shape[1]
    pos_loss = F.mse_loss(sim[:, 0], torch.ones_like(sim[:, 0]))
    if q > 1:
        neg_loss = F.mse_loss(sim[:, 1:], torch.zeros_like(sim[:, 1:]))
    else:
        neg_loss = pos_loss.new_zeros(())
    return pos_loss + neg_loss


def gpu_probe(data, config, model, output_dir):
    device = config.base.device
    model = model.to(device).train()
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
             "nvidia_smi_samples": []}
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-5)
    parsed = []
    for step in range(4):
        output = forward(model, episode)
        loss = F.smooth_l1_loss(output.prediction, episode.query_y.to(
            device=output.prediction.device, dtype=output.prediction.dtype))
        loss = loss + contrastive_loss(model, TEMPERATURE)
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


def internal_score(model, banks, label_scale):
    model.eval()
    total, count = 0.0, 0
    with torch.no_grad():
        for bank in banks.values():
            for episode in bank:
                output = forward(model, episode)
                error = float(F.mse_loss(output.prediction, episode.query_y.to(
                    device=output.prediction.device,
                    dtype=output.prediction.dtype)))
                total += error * label_scale.scale ** 2
                count += 1
    model.train()
    return total / max(count, 1)


def train(data, config, output_dir, progress_path):
    base = config.base
    torch.manual_seed(base.seed)
    rng = np.random.default_rng(base.seed)
    label_scale = training_label_scale(data)
    model = build_model(config, data)
    probe = gpu_probe(data, config, model, output_dir)
    if not probe["torch_cuda_is_available"] or not probe["batch_device_check"]:
        raise RuntimeError("GPU verification failed")
    model = model.to(base.device)
    trainable = [p for p in model.parameters() if p.requires_grad]
    fast_prefixes = ("meta.term.", "transport.")
    fast_parameters = [p for n, p in model.named_parameters()
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
    banks = internal_validation_bank(data, internal_components, label_scale)
    print("arm=" + config.arm + " fit " + str(len(fit_components))
          + " internal " + str(len(internal_components)), flush=True)
    best_state, best_value, best_step = None, float("inf"), 0
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
            output = forward(model, episode)
            query_y = episode.query_y.to(device=output.prediction.device,
                                         dtype=output.prediction.dtype)
            loss_post = F.smooth_l1_loss(output.prediction, query_y)
            loss_pre = F.smooth_l1_loss(output.zero_shot, query_y)
            loss_rank = ranking_term(output.prediction, query_y, base,
                                     label_scale)
            loss_centered = centered_task_error(output.prediction, query_y)
            loss = (loss_post + base.zero_shot_loss_weight * loss_pre
                    + base.ranking_loss_weight * loss_rank
                    + base.shape_loss_weight * loss_centered
                    + base.support_match_loss_weight
                    * output.support_match_loss)
            if config.arm == "K":
                loss = loss + CONTRAST_WEIGHT * contrastive_loss(
                    model, TEMPERATURE)
            else:
                loss = loss + CONTRAST_WEIGHT * regression_alignment(model)
            loss_value += float(loss.detach())
            (loss / len(selected)).backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), base.grad_clip)
        optimizer.step()
        trace.append(loss_value / len(selected))
        if step % base.val_interval == 0 or step == base.steps:
            value = internal_score(model, banks, label_scale)
            progress = {"step": step, "arm": config.arm,
                        "internal_val_mse_pk": value,
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
        raise RuntimeError("training produced no checkpoint")
    model.load_state_dict(best_state)
    report = {"arm": config.arm, "best_internal_val_mse_pk": best_value,
              "best_step": best_step, "loss_trace": trace,
              "label_scale": asdict(label_scale),
              "checkpoint_selection": "meta_train internal components only",
              "fit_components": len(fit_components),
              "internal_val_components": len(internal_components),
              "optimization_steps": base.steps,
              "wall_time_seconds": time.monotonic() - started,
              "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                     if base.device.startswith("cuda") else 0.0),
              "trainable_parameters": int(sum(p.numel() for p in trainable)),
              "gpu_probe": probe}
    payload = {"model_state": model.state_dict(), "config": asdict(base),
               "arm": config.arm, "label_scale": asdict(label_scale),
               "stageK_version": 1}
    torch.save(payload, output_dir / "checkpoint.pt")
    (output_dir / "RESULT.json").write_text(json.dumps(
        {"arm": config.arm, "report": report,
         "meta_test": data.seal_record()}, indent=1), encoding="utf-8")
    print(f"wrote {output_dir}")
    return model, report, label_scale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", required=True, choices=("K", "K-REG"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--steps", type=int, default=1200)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()
    progress_path = args.output / "progress.jsonl"
    if progress_path.exists() and not args.force:
        raise SystemExit(f"{progress_path} exists; pass --force to overwrite")
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    config = TrainConfig(arch="similarity_only", steps=args.steps,
                         seed=args.seed, split_directory=str(SPLIT),
                         device=args.device, amp=False)
    from tools.research.stageD_level_panel.train_staged import StageEConfig
    args.output.mkdir(parents=True, exist_ok=True)
    train(data, StageEConfig(base=config, arm=args.arm), args.output,
          progress_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
