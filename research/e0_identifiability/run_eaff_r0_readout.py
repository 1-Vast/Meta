"""Run the E-AFF-R0 readout diagnosis.

R0 asks what the registered affinity readout can and cannot see. It uses the
repository's own `concordance` implementation and the real H0C task/component
partition. It reads no affinity label: parts 1 and 2 are algebraic and
simulated, and part 3 reuses per-task metrics H0C already published.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from research.e0_identifiability.metrics import concordance
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-R0_READOUT_DIAGNOSIS"
TEST_LIGANDS_PER_TASK = 20
VARIANCE_RATIOS = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
SEED = 20260808


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def component_macro(per_task: list[dict], key: str) -> float:
    """Reproduce the H0C aggregation: task mean inside component, then over components."""
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in per_task:
        grouped[row["closure_component_id"]].append(row[key])
    return float(np.mean([float(np.mean(values)) for values in grouped.values()]))


def invariance_check(task_sizes: list[int], rng: np.random.Generator) -> dict:
    """Exact invariance of within-task concordance under per-task affine maps.

    Within one task the metric compares only the signs of pairwise differences,
    so adding any constant or multiplying by any positive constant - the exact
    algebraic form of a task-level affinity location and scale - cannot change
    it. This is checked numerically against the repository implementation.
    """
    deviations = {"prediction_shift": 0.0, "prediction_scale": 0.0,
                  "label_shift": 0.0, "label_scale": 0.0}
    for size in task_sizes:
        labels = rng.normal(size=size)
        predictions = rng.normal(size=size)
        base = concordance(labels, predictions)
        shift = float(rng.normal() * 10.0)
        scale = float(abs(rng.normal()) + 0.1)
        deviations["prediction_shift"] = max(
            deviations["prediction_shift"], abs(concordance(labels, predictions + shift) - base))
        deviations["prediction_scale"] = max(
            deviations["prediction_scale"], abs(concordance(labels, predictions * scale) - base))
        deviations["label_shift"] = max(
            deviations["label_shift"], abs(concordance(labels + shift, predictions) - base))
        deviations["label_scale"] = max(
            deviations["label_scale"], abs(concordance(labels * scale, predictions) - base))
    return {
        "max_absolute_deviation": {name: float(value) for name, value in deviations.items()},
        "exactly_invariant": all(value == 0.0 for value in deviations.values()),
        "meaning": (
            "within-task concordance is exactly invariant to task-level affinity "
            "location and scale, so any protein contribution expressed as a task "
            "level is assigned zero credit"),
    }


def location_credit(tasks: list[dict], rng: np.random.Generator) -> list[dict]:
    """Credit assigned to a perfect task-level predictor, by readout.

    Generates y[t,i] = level[t] + within[t,i] and scores oracles that see only
    the level, only the within-task variation, or both. The level channel is the
    one a correct protein sets for a chemical series.
    """
    rows = []
    for ratio in VARIANCE_RATIOS:
        levels = rng.normal(scale=ratio, size=len(tasks))
        per_task = []
        squared: dict[str, list[float]] = defaultdict(list)
        for index, task in enumerate(tasks):
            within = rng.normal(scale=1.0, size=TEST_LIGANDS_PER_TASK)
            observed = levels[index] + within
            predictors = {
                "level_oracle": np.full(TEST_LIGANDS_PER_TASK, levels[index]),
                "within_oracle": within,
                "full_oracle": observed,
                "global_mean": np.zeros(TEST_LIGANDS_PER_TASK),
            }
            per_task.append({
                "closure_component_id": task["closure_component_id"],
                **{name: concordance(observed, value) for name, value in predictors.items()},
            })
            for name, value in predictors.items():
                squared[name].extend(np.square(observed - value).tolist())
        rows.append({
            "between_task_sd": ratio,
            "within_task_sd": 1.0,
            "location_share_of_variance": round(ratio ** 2 / (1.0 + ratio ** 2), 4),
            "within_task_concordance": {
                name: round(component_macro(per_task, name), 5)
                for name in ("level_oracle", "within_oracle", "full_oracle", "global_mean")},
            "location_sensitive_rmse": {
                name: round(float(np.sqrt(np.mean(values))), 5)
                for name, values in squared.items()},
        })
    return rows


def published_geometry_effect(per_task: list[dict]) -> dict:
    """What the frozen geometry actually did to H0C's published per-task scores."""
    local = np.asarray([row["local_ligand"] for row in per_task])
    correct = np.asarray([row["correct"] for row in per_task])
    deranged = np.asarray([row["deranged"] for row in per_task])
    contrasts = {
        "correct_minus_local": correct - local,
        "deranged_minus_local": deranged - local,
        "correct_minus_deranged": correct - deranged,
    }
    return {
        "tasks": len(per_task),
        "contrasts": {
            name: {
                "mean": round(float(values.mean()), 5),
                "median": round(float(np.median(values)), 5),
                "min": round(float(values.min()), 5),
                "max": round(float(values.max()), 5),
                "tasks_unchanged": int(np.sum(np.abs(values) < 1e-12)),
            }
            for name, values in contrasts.items()},
        "meaning": (
            "the geometry term moved almost every task, so the null result is not "
            "an inert feature; its sign carried no consistent affinity direction"),
    }


def run(args) -> dict:
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output exists: {output}")
    h0c_root = Path(args.h0c)
    per_task = list(_read_jsonl(h0c_root / "task_metrics.jsonl"))
    rng = np.random.default_rng(SEED)

    report = {
        "schema": "MetaSieve.EAffR0.v1",
        "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "affinity_labels_read": False,
        "davis_label_reads": 0,
        "recipient_label_reads": 0,
        "training_performed": False,
        "readout_under_test": "within_task_concordance_macro_averaged_over_closure_components",
        "theory_metric": "hausdorff_W1_between_emitted_law_classes",
        "tasks": len(per_task),
        "test_ligands_per_task": TEST_LIGANDS_PER_TASK,
        "seed": SEED,
        "invariance": invariance_check([TEST_LIGANDS_PER_TASK] * len(per_task), rng),
        "location_credit": location_credit(per_task, rng),
        "published_geometry_effect": published_geometry_effect(per_task),
        "interpretation_limits": [
            "R0 diagnoses the readout, not the biology",
            "it does not show that a location channel carries protein-specific affinity",
            "it shows the registered readout would assign that channel zero credit if it did",
            "the simulated generative model is an illustration, not a fit to source data",
        ],
    }
    verdicts = []
    if report["invariance"]["exactly_invariant"]:
        verdicts.append("READOUT_BLIND_TO_TASK_LEVEL_AFFINITY_LOCATION")
    if all(abs(row["within_task_concordance"]["level_oracle"] - 0.5) < 1e-12
           for row in report["location_credit"]):
        verdicts.append("PERFECT_LEVEL_PREDICTOR_SCORES_CHANCE_AT_EVERY_VARIANCE_SHARE")
    report["verdict"] = "|".join(verdicts) if verdicts else "READOUT_DIAGNOSIS_INCONCLUSIVE"

    output.mkdir(parents=True)
    _write_json(output / "report.json", report)
    manifest = {
        "stage": STAGE,
        "affinity_labels_read": False,
        "inputs": {
            "h0c_task_metrics": sha256_file(h0c_root / "task_metrics.jsonl"),
            "h0c_report": sha256_file(h0c_root / "report.json"),
            "metrics_module": sha256_file(Path(__file__).with_name("metrics.py")),
            "registration": sha256_file(Path(__file__).with_name("EAFF_R0_REGISTRATION.md")),
        },
        "outputs": {"report.json": sha256_file(output / "report.json")},
        "label_reads": {"affinity_values": 0, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h0c", default="research/e0_identifiability/artifacts/eaff_h0c_v1_run2")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_r0_v1")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
