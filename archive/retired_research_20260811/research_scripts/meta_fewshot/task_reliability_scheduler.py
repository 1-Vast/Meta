"""Cross-fitted task-scheduler audit primitives.

The scorer is deliberately offline: destructive partner/support controls are
accepted only by the evaluation helpers and can never become fit features.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib

import numpy as np
from scipy import stats


DEFAULT_ALPHAS = (1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0)


@dataclass(frozen=True)
class RidgeFold:
    test_indices: np.ndarray
    feature_mean: np.ndarray
    feature_scale: np.ndarray
    coefficient: np.ndarray
    intercept: float
    alpha: float

    def predict(self, features: np.ndarray) -> np.ndarray:
        standardized = (features - self.feature_mean) / self.feature_scale
        return standardized @ self.coefficient + self.intercept


def stable_fold_assignments(groups: np.ndarray, folds: int = 5) -> np.ndarray:
    """Assign whole groups to deterministic, approximately balanced folds."""
    if folds < 2:
        raise ValueError("cross-fitting requires at least two folds")
    unique, counts = np.unique(groups.astype(str), return_counts=True)
    if len(unique) < folds:
        raise ValueError("fewer groups than requested folds")
    order = sorted(
        zip(unique, counts),
        key=lambda item: (-int(item[1]), hashlib.sha256(item[0].encode()).hexdigest()),
    )
    loads = np.zeros(folds, dtype=np.int64)
    mapping: dict[str, int] = {}
    for group, count in order:
        fold = int(np.argmin(loads))
        mapping[group] = fold
        loads[fold] += int(count)
    return np.asarray([mapping[str(group)] for group in groups], dtype=np.int64)


def _fit_ridge(features: np.ndarray, target: np.ndarray,
               alpha: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    mean = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale < 1e-8] = 1.0
    x = (features - mean) / scale
    intercept = float(target.mean())
    centered = target - intercept
    coefficient = np.linalg.solve(
        x.T @ x + alpha * np.eye(x.shape[1]), x.T @ centered)
    return mean, scale, coefficient, intercept


def _group_macro_mse(error: np.ndarray, groups: np.ndarray) -> float:
    return float(np.mean([
        np.mean(error[groups == group] ** 2) for group in np.unique(groups)
    ]))


def _select_alpha(features: np.ndarray, target: np.ndarray, groups: np.ndarray,
                  alphas: tuple[float, ...]) -> float:
    inner_folds = min(3, len(np.unique(groups)))
    assignments = stable_fold_assignments(groups, inner_folds)
    losses = []
    for alpha in alphas:
        prediction = np.empty_like(target, dtype=np.float64)
        for fold in range(inner_folds):
            train = assignments != fold
            held = ~train
            mean, scale, coefficient, intercept = _fit_ridge(
                features[train], target[train], alpha)
            prediction[held] = ((features[held] - mean) / scale) @ coefficient + intercept
        losses.append(_group_macro_mse(prediction - target, groups))
    return float(alphas[int(np.argmin(losses))])


def cross_fitted_ridge(
        features: np.ndarray, target: np.ndarray, groups: np.ndarray, *,
        folds: int = 5,
        alphas: tuple[float, ...] = DEFAULT_ALPHAS) -> tuple[np.ndarray, list[RidgeFold]]:
    """Nested group-cross-fitted ridge predictions and reusable held-fold fits."""
    features = np.asarray(features, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    groups = np.asarray(groups).astype(str)
    if features.ndim != 2 or target.shape != (len(features),):
        raise ValueError("invalid scorer feature or target shape")
    if not np.isfinite(features).all() or not np.isfinite(target).all():
        raise ValueError("scheduler audit arrays must be finite")
    assignments = stable_fold_assignments(groups, folds)
    prediction = np.empty_like(target)
    fitted: list[RidgeFold] = []
    for fold in range(folds):
        train = assignments != fold
        held = ~train
        alpha = _select_alpha(features[train], target[train], groups[train], alphas)
        mean, scale, coefficient, intercept = _fit_ridge(
            features[train], target[train], alpha)
        prediction[held] = ((features[held] - mean) / scale) @ coefficient + intercept
        fitted.append(RidgeFold(
            test_indices=np.flatnonzero(held), feature_mean=mean,
            feature_scale=scale, coefficient=coefficient, intercept=intercept,
            alpha=alpha))
    return prediction, fitted


def apply_cross_fitted(folds: list[RidgeFold], features: np.ndarray) -> np.ndarray:
    prediction = np.empty(len(features), dtype=np.float64)
    covered = np.zeros(len(features), dtype=bool)
    for fitted in folds:
        prediction[fitted.test_indices] = fitted.predict(features[fitted.test_indices])
        covered[fitted.test_indices] = True
    if not covered.all():
        raise ValueError("cross-fitted models do not cover every row")
    return prediction


def permute_informative_rows(
        features: np.ndarray, task_size: np.ndarray, *, seed: int,
        informative_columns: tuple[int, ...] = (0, 1)) -> np.ndarray:
    """Equal-capacity null: permute statistics within task-size quartiles."""
    result = np.asarray(features, dtype=np.float64).copy()
    size = np.asarray(task_size, dtype=np.float64)
    edges = np.unique(np.quantile(size, (0.0, 0.25, 0.5, 0.75, 1.0)))
    strata = (np.digitize(size, edges[1:-1], right=True)
              if len(edges) > 1 else np.zeros(len(size), dtype=np.int64))
    rng = np.random.default_rng(seed)
    for stratum in np.unique(strata):
        rows = np.flatnonzero(strata == stratum)
        donor = rng.permutation(rows)
        result[np.ix_(rows, informative_columns)] = features[
            np.ix_(donor, informative_columns)]
    return result


def component_macro(values: np.ndarray, groups: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unique = np.unique(groups.astype(str))
    macro = np.asarray([
        np.mean(values[groups.astype(str) == group], axis=0) for group in unique
    ])
    return unique, macro


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    value = stats.spearmanr(left, right).statistic
    return float(value) if np.isfinite(value) else 0.0


def residualize(values: np.ndarray, nuisance: np.ndarray) -> np.ndarray:
    nuisance = np.asarray(nuisance, dtype=np.float64)
    design = np.column_stack((np.ones(len(nuisance)), nuisance))
    return values - design @ np.linalg.lstsq(design, values, rcond=None)[0]


def component_bootstrap(
        score: np.ndarray, utility: np.ndarray, null_score: np.ndarray,
        transformed_scores: dict[str, np.ndarray],
        transformed_utilities: dict[str, np.ndarray], nuisance: np.ndarray,
        groups: np.ndarray, *, seed: int, draws: int = 9999) -> dict:
    """Whole-component bootstrap for the frozen scheduler Gate contrasts."""
    unique, joined = component_macro(np.column_stack((
        score, utility, null_score,
        *[transformed_scores[name] for name in sorted(transformed_scores)],
        *[transformed_utilities[name] for name in sorted(transformed_utilities)],
        nuisance,
    )), groups)
    names = sorted(transformed_scores)
    utility_names = sorted(transformed_utilities)
    score_m = joined[:, 0]
    utility_m = joined[:, 1]
    null_m = joined[:, 2]
    offset = 3
    transformed_score_m = {
        name: joined[:, offset + index] for index, name in enumerate(names)}
    offset += len(names)
    transformed_utility_m = {
        name: joined[:, offset + index]
        for index, name in enumerate(utility_names)}
    offset += len(utility_names)
    nuisance_m = joined[:, offset:]

    clean_corr = spearman(score_m, utility_m)
    null_corr = spearman(null_m, utility_m)
    clean_residual = residualize(score_m, nuisance_m)
    utility_residual = residualize(utility_m, nuisance_m)
    null_residual = residualize(null_m, nuisance_m)
    observed = {
        "clean_correlation": clean_corr,
        "matched_null_correlation": null_corr,
        "matched_null_advantage": clean_corr - null_corr,
        "nuisance_residual_advantage": (
            spearman(clean_residual, utility_residual)
            - spearman(null_residual, utility_residual)),
        "score_correlation_losses": {
            name: clean_corr - spearman(values, utility_m)
            for name, values in transformed_score_m.items()},
        "utility_deltas": {
            name: float(np.mean(utility_m - values))
            for name, values in transformed_utility_m.items()},
    }
    rng = np.random.default_rng(seed)
    samples: dict[str, list[float]] = {
        "clean_correlation": [], "matched_null_advantage": [],
        "nuisance_residual_advantage": []}
    for name in names:
        samples[f"score_loss:{name}"] = []
    for name in utility_names:
        samples[f"utility_delta:{name}"] = []
    for _ in range(draws):
        index = rng.integers(0, len(unique), size=len(unique))
        s, u, n = score_m[index], utility_m[index], null_m[index]
        clean = spearman(s, u)
        samples["clean_correlation"].append(clean)
        samples["matched_null_advantage"].append(clean - spearman(n, u))
        nr = nuisance_m[index]
        samples["nuisance_residual_advantage"].append(
            spearman(residualize(s, nr), residualize(u, nr))
            - spearman(residualize(n, nr), residualize(u, nr)))
        for name, values in transformed_score_m.items():
            samples[f"score_loss:{name}"].append(clean - spearman(values[index], u))
        for name, values in transformed_utility_m.items():
            samples[f"utility_delta:{name}"].append(float(np.mean(u - values[index])))
    lower = {name: float(np.quantile(values, 0.05))
             for name, values in samples.items()}
    return {"components": len(unique), "observed": observed, "lcb95": lower}
