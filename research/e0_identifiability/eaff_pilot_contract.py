"""Pure contracts for the E-AFF-P0 fixed-radial affinity pilot."""
from __future__ import annotations

from collections import defaultdict

import numpy as np

from research.e0_identifiability.metrics import concordance


def assert_paffinity_direction(values: np.ndarray, molar_values: np.ndarray) -> None:
    """Fail if pK and molar affinity do not use stronger-is-larger semantics."""
    p = np.asarray(values, dtype=np.float64)
    molar = np.asarray(molar_values, dtype=np.float64)
    if len(p) < 2 or np.any(~np.isfinite(p)) or np.any(molar <= 0):
        raise ValueError("affinity direction audit needs finite positive values")
    order = np.argsort(molar)
    if np.any(np.diff(p[order]) > 1e-10):
        raise ValueError("pAffinity must decrease as molar K increases")


def coupling_null(features: np.ndarray, epsilon: float = 1e-12) -> np.ndarray:
    """Preserve chemistry-pair and radial marginals while deleting their coupling."""
    value = np.asarray(features, dtype=np.float64)
    if value.shape[-3:] != (8, 6, 6):
        raise ValueError("coupling null requires [...,8,6,6] features")
    chemistry = value.sum(axis=-1, keepdims=True)
    radial = value.sum(axis=(-3, -2), keepdims=True)
    total = value.sum(axis=(-3, -2, -1), keepdims=True)
    if np.any(np.abs(total) <= epsilon):
        raise ValueError("coupling null is undefined for zero-mass features")
    return chemistry * radial / total


def task_pair_differences(features: np.ndarray, targets: np.ndarray,
                          task_ids: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return all within-task differences with equal total weight per task."""
    x = np.asarray(features, dtype=np.float64)
    y = np.asarray(targets, dtype=np.float64)
    tasks = np.asarray(task_ids)
    differences, labels, weights = [], [], []
    for task in sorted(set(tasks.tolist())):
        indices = np.flatnonzero(tasks == task)
        left, right = np.triu_indices(len(indices), 1)
        if not len(left):
            continue
        differences.append(x[indices[left]] - x[indices[right]])
        labels.append(y[indices[left]] - y[indices[right]])
        weights.append(np.full(len(left), 1.0 / len(left), dtype=np.float64))
    if not differences:
        raise ValueError("no within-task pairs")
    return np.concatenate(differences), np.concatenate(labels), np.concatenate(weights)


def fit_pair_ridge(features: np.ndarray, targets: np.ndarray, task_ids: np.ndarray,
                   alpha: float) -> tuple[np.ndarray, dict]:
    """Fit a deterministic task-balanced Ridge direction in the original basis."""
    from sklearn.linear_model import Ridge

    dx, dy, weights = task_pair_differences(features, targets, task_ids)
    scale = dx.std(axis=0)
    active = scale > 1e-10
    safe = np.where(active, scale, 1.0)
    model = Ridge(alpha=alpha, fit_intercept=False, solver="svd")
    model.fit(dx[:, active] / safe[active], dy, sample_weight=weights)
    direction = np.zeros(features.shape[1], dtype=np.float64)
    direction[active] = model.coef_ / safe[active]
    diagnostics = {
        "alpha": float(alpha),
        "pairs": int(len(dy)),
        "tasks": int(len(set(np.asarray(task_ids).tolist()))),
        "active_dimensions": int(active.sum()),
        "weighted_pair_rmse": float(np.sqrt(np.average(
            np.square(dx @ direction - dy), weights=weights))),
    }
    return direction, diagnostics


def task_metrics(rows: list[dict], labels: np.ndarray,
                 predictions: dict[str, np.ndarray]) -> list[dict]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[row["task_id"]].append(index)
    result = []
    for task, indices in sorted(grouped.items()):
        first = rows[indices[0]]
        result.append({
            "task_id": task,
            "closure_component_id": first["closure_component_id"],
            "endpoint_family": first["endpoint_family"],
            **{name: concordance(labels[indices], value[indices])
               for name, value in predictions.items()},
        })
    return result


def component_macro_contrasts(per_task: list[dict]) -> dict[str, float]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in per_task:
        grouped[row["closure_component_id"]].append(row)
    component = []
    for rows in grouped.values():
        component.append({key: float(np.mean([row[key] for row in rows]))
                          for key in ("ligand", "correct", "deranged", "null")})
    means = {key: float(np.mean([row[key] for row in component]))
             for key in ("ligand", "correct", "deranged", "null")}
    return {
        **means,
        "correct_minus_ligand": means["correct"] - means["ligand"],
        "correct_minus_deranged": means["correct"] - means["deranged"],
        "correct_minus_null": means["correct"] - means["null"],
        "components": len(component),
    }


def component_bootstrap(per_task: list[dict], seed: int = 17,
                        draws: int = 2000) -> dict[str, list[float]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in per_task:
        grouped[row["closure_component_id"]].append(row)
    components = []
    for rows in grouped.values():
        values = {key: float(np.mean([row[key] for row in rows]))
                  for key in ("ligand", "correct", "deranged", "null")}
        components.append([
            values["correct"] - values["ligand"],
            values["correct"] - values["deranged"],
            values["correct"] - values["null"],
        ])
    matrix = np.asarray(components, dtype=np.float64)
    generator = np.random.default_rng(seed)
    samples = np.empty((draws, 3), dtype=np.float64)
    for draw in range(draws):
        samples[draw] = matrix[generator.integers(0, len(matrix), len(matrix))].mean(0)
    names = ("correct_minus_ligand", "correct_minus_deranged", "correct_minus_null")
    return {name: [float(np.quantile(samples[:, index], 0.025)),
                   float(np.quantile(samples[:, index], 0.975))]
            for index, name in enumerate(names)}
