"""Verify batched V1 inference against the frozen scalar CUDA artifact."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import torch

from model.metasieve_v1 import MetaSieveV1
from research.meta_fewshot.train_main_v1 import (
    V1TrainConfig,
    build_tasks,
    load_data,
    predict,
)


ROOT = Path(__file__).resolve().parents[2]
FORMAL = ROOT / "report/meta_fewshot/main_v1"


def verify(*, device: str = "cuda", tolerance: float = 2e-5) -> dict:
    result = json.loads(
        (FORMAL / "MAIN_V1_COLD_TARGET_RESULT.json").read_text(encoding="utf-8"))
    config = V1TrainConfig(**result["config"])
    cells, tensors, _, normalization = load_data(device)
    seed = int(result["seeds"][0])
    arm = "uniform_clean_correct"
    checkpoint = torch.load(
        FORMAL / "checkpoints" / f"uniform_clean_seed{seed}.pt",
        map_location=device, weights_only=False)
    model = MetaSieveV1(288, config.section_dim, config.ridge).to(device)
    model.load_state_dict(checkpoint["model_state"])
    tasks = build_tasks(cells, "meta_test", max(config.support_sizes), config.min_query)
    current = []
    for k in config.support_sizes:
        current += predict(
            model, cells, tensors, tasks, family="correct", seed=seed, k=k,
            draws=config.test_draws, max_query=config.test_max_query, arm=arm,
            y_scale=normalization["y_scale"],
            max_support_k=max(config.support_sizes))

    frozen = []
    prediction_path = ROOT / result["prediction_artifact"]["path"]
    with gzip.open(prediction_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["seed"] == seed and row["arm"] == arm:
                frozen.append(row)
    if len(current) != len(frozen):
        raise AssertionError("batched and scalar prediction row counts differ")
    differences = []
    for new, old in zip(current, frozen):
        new_value = new.pop("prediction_standardized")
        old_value = old.pop("prediction_standardized")
        if new != old:
            raise AssertionError("batched inference changed episode identity/order")
        differences.append(abs(new_value - old_value))
    maximum = float(np.max(differences)) if differences else 0.0
    if maximum > tolerance:
        raise AssertionError(
            f"batched inference exceeded tolerance: {maximum} > {tolerance}")
    return {
        "schema": "MetaSieve.V1VectorizationVerification.v1",
        "device": device,
        "rows": len(current),
        "episode_metadata_exact": True,
        "maximum_absolute_prediction_difference": maximum,
        "tolerance": tolerance,
        "pass": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--tolerance", type=float, default=2e-5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify(device=args.device, tolerance=args.tolerance)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
