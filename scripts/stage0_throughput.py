"""Per-architecture throughput, memory and utilization benchmark.

Measures the real episodic training step (forward, counterfactual forward,
wrong-protein forward, backward) so that budget planning uses the cost the
trainer actually pays.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    batch_counterfactual_episode, binding_contrastive_loss, compact_episode,
    counterfactual_label_assignments, forward, normalized_episode,
    resolve_architecture, training_label_scale, wrong_protein_zero_shot,
)


def benchmark(arch: str, data: QPSMPData, config: TrainConfig, steps: int,
              support_sizes: tuple[int, ...]) -> dict:
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    model = resolve_architecture(arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks, adapter_scale=config.adapter_scale,
        use_cartesian=False).to(config.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate)
    label_scale = training_label_scale(data)
    amp = config.amp and config.device.startswith("cuda")
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    if config.device.startswith("cuda"):
        torch.cuda.reset_peak_memory_stats()
        torch.cuda.synchronize()
    episodes = 0
    started = time.monotonic()
    for step in range(steps):
        support_size = support_sizes[step % len(support_sizes)]
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast("cuda", dtype=torch.bfloat16, enabled=amp):
            for _ in range(config.episodes_per_step):
                spec = data.draw_episode(
                    "meta_train", support_size,
                    int(rng.integers(config.min_query_size, config.query_size + 1)),
                    rng)
                episode = compact_episode(normalized_episode(
                    data.materialize(spec), label_scale))
                out = forward(model, episode, adapt=support_size > 0)
                truth = episode.query_y.to(device=out.prediction.device,
                                           dtype=out.prediction.dtype)
                loss = (F.smooth_l1_loss(out.prediction, truth)
                        + F.smooth_l1_loss(out.zero_shot, truth)
                        + 0.05 * out.support_match_loss)
                if support_size > 0:
                    assignments = counterfactual_label_assignments(
                        out, episode.support_y)
                    wrong_episode = batch_counterfactual_episode(episode, assignments)
                    wrong = forward(model, wrong_episode)
                    wrong_truth = wrong_episode.query_y.to(
                        device=wrong.prediction.device, dtype=wrong.prediction.dtype)
                    errors = [(out.prediction - truth).square().mean()]
                    errors.extend((wrong.prediction - wrong_truth
                                   ).square().mean(-1).unbind())
                    loss = loss + binding_contrastive_loss(errors, 0.1)
                    wrong_zero = wrong_protein_zero_shot(
                        model, data, episode, episode.spec.donor_target)
                    loss = loss + binding_contrastive_loss([
                        (out.zero_shot - truth).square().mean(),
                        (wrong_zero - truth).square().mean()], 0.1)
                scaler.scale(loss / config.episodes_per_step).backward()
                episodes += 1
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        scaler.step(optimizer)
        scaler.update()
    if config.device.startswith("cuda"):
        torch.cuda.synchronize()
    elapsed = time.monotonic() - started
    return {
        "arch": arch,
        "trainable_parameters": int(sum(
            p.numel() for p in model.parameters() if p.requires_grad)),
        "steps": steps,
        "episodes": episodes,
        "seconds": elapsed,
        "seconds_per_step": elapsed / steps,
        "seconds_per_episode": elapsed / episodes,
        "peak_cuda_memory_mb": (torch.cuda.max_memory_allocated() / 2 ** 20
                                if config.device.startswith("cuda") else 0.0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--device", default=TrainConfig.device)
    parser.add_argument("--archs", nargs="+", default=["bpsf", "grammar"])
    args = parser.parse_args()
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    config = TrainConfig(
        seed=20260812, hidden_dim=192, task_dim=48, ligand_layers=4,
        pair_dim=96, pair_blocks=4, pair_latents=24, pair_heads=8,
        pair_chunk_size=8, query_size=20, min_query_size=12,
        episodes_per_step=4, device=args.device)
    rows = []
    for arch in args.archs:
        rows.append(benchmark(arch, data, config, args.steps, (0, 1, 2, 3, 5)))
        print(json.dumps(rows[-1], indent=2))
        if args.device.startswith("cuda"):
            torch.cuda.empty_cache()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(
        {"schema": "MetaSieve.QPSMPThroughput.v1", "device": args.device,
         "benchmarks": rows}, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
