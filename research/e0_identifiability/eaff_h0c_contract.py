"""Pure contracts for the E-AFF-H0C interaction-residual diagnostic."""
from __future__ import annotations

from collections import defaultdict

import numpy as np


def centered_interaction(features: np.ndarray, epsilon: float = 1e-12) -> np.ndarray:
    """Remove chemistry/radial marginals and normalize by positive total mass."""
    value = np.asarray(features, dtype=np.float64)
    if value.shape[-3:] != (8, 6, 6):
        raise ValueError("interaction residual requires [...,8,6,6] features")
    chemistry = value.sum(axis=-1, keepdims=True)
    radial = value.sum(axis=(-3, -2), keepdims=True)
    total = value.sum(axis=(-3, -2, -1), keepdims=True)
    if np.any(total <= epsilon):
        raise ValueError("interaction residual requires positive total mass")
    return (value - chemistry * radial / total) / total


def component_summary(per_task: list[dict]) -> dict[str, float]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in per_task:
        grouped[row["closure_component_id"]].append(row)
    keys = ("global_ligand", "local_ligand", "correct", "deranged")
    components = [{key: float(np.mean([row[key] for row in rows])) for key in keys}
                  for rows in grouped.values()]
    means = {key: float(np.mean([row[key] for row in components])) for key in keys}
    return {
        **means,
        "correct_minus_local_ligand": means["correct"] - means["local_ligand"],
        "correct_minus_deranged": means["correct"] - means["deranged"],
        "components": len(components),
    }


def component_bootstrap(per_task: list[dict], seed: int = 29,
                        draws: int = 2000) -> dict[str, list[float]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in per_task:
        grouped[row["closure_component_id"]].append(row)
    matrix = []
    for rows in grouped.values():
        correct = float(np.mean([row["correct"] for row in rows]))
        matrix.append([
            correct - float(np.mean([row["local_ligand"] for row in rows])),
            correct - float(np.mean([row["deranged"] for row in rows])),
        ])
    values = np.asarray(matrix, dtype=np.float64)
    generator = np.random.default_rng(seed)
    samples = np.empty((draws, 2), dtype=np.float64)
    for draw in range(draws):
        samples[draw] = values[generator.integers(0, len(values), len(values))].mean(0)
    names = ("correct_minus_local_ligand", "correct_minus_deranged")
    return {name: [float(np.quantile(samples[:, index], 0.025)),
                   float(np.quantile(samples[:, index], 0.975))]
            for index, name in enumerate(names)}
