"""Reproducible CUDA throughput benchmark for the V1 training paths."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time

import torch

from model.metasieve_v1 import MetaSieveV1
from research.meta_fewshot.train_main_v1 import (
    V1TrainConfig,
    load_data,
    train_one,
)


def benchmark(*, steps: int, repeats: int, val_interval: int,
              val_draws: int) -> dict:
    if not torch.cuda.is_available():
        raise RuntimeError("V1 GPU benchmark requires CUDA")
    cells, tensors, _, normalization = load_data("cuda")
    config = replace(
        V1TrainConfig(), steps=steps, val_interval=val_interval,
        val_draws=val_draws)
    torch.cuda.synchronize()
    arms = {
        "uniform_clean": ("uniform", False),
        "ats_clean": ("ats", False),
    }
    measurements = {}
    for arm, (schedule, noise) in arms.items():
        values = []
        for repeat in range(repeats):
            seed = 20260831
            torch.manual_seed(seed)
            model = MetaSieveV1(288, config.section_dim, config.ridge).cuda()
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.synchronize()
            started = time.perf_counter()
            train_one(
                model, cells, tensors, family="correct", schedule=schedule,
                add_support_noise=noise, seed=seed, config=config,
                y_scale=normalization["y_scale"])
            torch.cuda.synchronize()
            elapsed = time.perf_counter() - started
            values.append({
                "repeat": repeat,
                "elapsed_seconds": elapsed,
                "milliseconds_per_step_end_to_end": (
                    1000.0 * elapsed / steps),
                "peak_allocated_mib": (
                    torch.cuda.max_memory_allocated() / (1024 ** 2)),
                "peak_reserved_mib": (
                    torch.cuda.max_memory_reserved() / (1024 ** 2)),
            })
        measurements[arm] = values
    return {
        "schema": "MetaSieve.V1GpuBenchmark.v1",
        "device": torch.cuda.get_device_name(),
        "torch": torch.__version__,
        "config": asdict(config),
        "measurements": measurements,
        "steady_measurement": {
            arm: values[-1] for arm, values in measurements.items()},
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--val-interval", type=int)
    parser.add_argument("--val-draws", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    val_interval = args.val_interval or args.steps
    if min(args.steps, args.repeats, val_interval, args.val_draws) < 1:
        raise ValueError("benchmark counts must be positive")
    result = benchmark(
        steps=args.steps, repeats=args.repeats,
        val_interval=val_interval, val_draws=args.val_draws)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
