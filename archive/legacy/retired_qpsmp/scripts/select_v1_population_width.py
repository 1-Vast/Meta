"""Select the narrow ligand-population width on meta-validation only."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import numpy as np
import torch

from model.metasieve_v1 import MetaSieveV1
from research.meta_fewshot.train_main_v0 import CORPUS, FEATURES, load_data, sha256
from research.meta_fewshot.train_main_v1 import SEEDS, V1TrainConfig, train_one


ROOT = Path(__file__).resolve().parents[2]
WIDTHS = (32, 64, 128, 256)


def choose_population_width(rows: list[dict]) -> dict:
    grouped = {}
    for row in rows:
        grouped.setdefault(row["population_hidden_dim"], []).append(
            row["best_combined_val_score"])
    if not grouped:
        raise ValueError("population selection has no validation rows")
    width = min(grouped, key=lambda value: (float(np.mean(grouped[value])), value))
    return {"population_hidden_dim": width,
            "mean_best_combined_val_score": float(np.mean(grouped[width]))}


def run(output: str | Path, *, device: str = "cuda") -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("V1 population selection is registered for CUDA")
    output = Path(output)
    if output.exists():
        raise FileExistsError(f"selection report already exists: {output}")
    cells, tensors, _, normalization = load_data(device)
    rows = []
    started = time.perf_counter()
    for width in WIDTHS:
        config = V1TrainConfig(
            support_only_section=True, population_hidden_dim=width,
            pair_hidden_dim=0, section_dim=4, ridge=1.0)
        for seed in SEEDS:
            torch.manual_seed(seed)
            model = MetaSieveV1(
                288, 4, 1.0, support_only_section=True,
                population_hidden_dim=width).to(device)
            diagnostics, _ = train_one(
                model, cells, tensors, family="correct", schedule="uniform",
                add_support_noise=False, seed=seed, config=config,
                y_scale=normalization["y_scale"])
            rows.append({
                "population_hidden_dim": width,
                "seed": seed,
                "best_step": diagnostics["best_step"],
                "best_combined_val_score": diagnostics[
                    "best_combined_val_score"],
                "clean_val_mse_by_k": diagnostics["clean_val_mse_by_k"],
                "noisy_support_val_mse_by_k": diagnostics[
                    "noisy_support_val_mse_by_k"],
            })
    report = {
        "schema": "MetaSieve.V1PopulationWidthSelection.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "device": device,
        "selection_metric": "three-seed mean best_combined_val_score",
        "meta_test_predictions_read": 0,
        "support_only_section": True,
        "section_dim": 4,
        "ridge": 1.0,
        "pair_hidden_dim": 0,
        "widths": list(WIDTHS),
        "seeds": list(SEEDS),
        "corpus_sha256": sha256(CORPUS / "manifest.json"),
        "features_sha256": sha256(FEATURES),
        "runner_sha256": sha256(Path(__file__)),
        "model_sha256": sha256(ROOT / "model/metasieve_v1.py"),
        "elapsed_seconds": time.perf_counter() - started,
        "selected": choose_population_width(rows),
        "rows": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    print(json.dumps(run(args.output, device=args.device), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
