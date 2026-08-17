"""Metrics and component-paired bootstrap for signed within-target gaps.

Every reported statistic is computed from **prediction rows** (one row per
evaluation pair), so a contrast between two arms can always be recomputed on
the same resample of components.  That is what "component-paired" means here:
the bootstrap draws a set of components once and both arms are re-scored on
exactly that draw, which removes the between-component variance that dominates
this corpus.

Targets and components are carried as integer codes so the bootstrap can run
2,000 draws without a Python loop over rows.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


class Rows:
    """Prediction rows for one arm on one evaluation bank."""

    __slots__ = ("pair_id", "delta_y", "delta_hat", "target", "component",
                 "target_codes", "component_codes")

    def __init__(self, pair_id: list[str], delta_y, delta_hat,
                 target: list[str], component: list[str]) -> None:
        self.pair_id = list(pair_id)
        self.delta_y = np.asarray(delta_y, dtype=np.float64)
        self.delta_hat = np.asarray(delta_hat, dtype=np.float64)
        self.target = list(target)
        self.component = list(component)
        target_index = {key: i for i, key in enumerate(sorted(set(target)))}
        component_index = {key: i for i, key in enumerate(sorted(set(component)))}
        self.target_codes = np.asarray([target_index[t] for t in target],
                                       dtype=np.int64)
        self.component_codes = np.asarray([component_index[c] for c in component],
                                          dtype=np.int64)

    def __len__(self) -> int:
        return int(self.delta_y.size)

    def select(self, mask) -> "Rows":
        mask = np.asarray(mask)
        keep = np.flatnonzero(mask) if mask.dtype == bool else mask
        return Rows([self.pair_id[int(i)] for i in keep], self.delta_y[keep],
                    self.delta_hat[keep], [self.target[int(i)] for i in keep],
                    [self.component[int(i)] for i in keep])

    def as_dict(self) -> dict:
        return {
            "pair_id": self.pair_id,
            "delta_y": self.delta_y.tolist(),
            "delta_hat": self.delta_hat.tolist(),
            "target": self.target,
            "component": self.component,
        }


def _pearson(truth: np.ndarray, prediction: np.ndarray) -> float:
    if truth.size < 3 or np.std(truth) < 1e-12 or np.std(prediction) < 1e-12:
        return float("nan")
    return float(np.corrcoef(truth, prediction)[0, 1])


def _spearman(truth: np.ndarray, prediction: np.ndarray) -> float:
    if truth.size < 3:
        return float("nan")
    value = stats.spearmanr(truth, prediction).statistic
    return float(value) if np.isfinite(value) else float("nan")


def _concordance(truth: np.ndarray, prediction: np.ndarray) -> float:
    """CI over the evaluation rows, via Kendall tau-b.

    For continuous targets the concordance index is `(tau + 1) / 2`, and tau-b
    handles the ties that a raw O(n^2) loop would have to special-case.
    """
    if truth.size < 3:
        return float("nan")
    tau = stats.kendalltau(truth, prediction).statistic
    return float((tau + 1.0) / 2.0) if np.isfinite(tau) else float("nan")


def _sign_accuracy(truth: np.ndarray, prediction: np.ndarray) -> float:
    usable = np.abs(truth) > 0
    if not usable.any():
        return float("nan")
    return float(np.mean(np.sign(truth[usable]) == np.sign(prediction[usable])))


def _equal_component_target_mean_mse(error: np.ndarray, targets: np.ndarray,
                                     components: np.ndarray) -> float:
    """Mean over components of the mean over that component's target MSEs.

    Targets differ by two orders of magnitude in pair count and components
    differ in target count, so a plain pooled MSE is effectively a report about
    the three largest targets.  This statistic gives every component one vote.
    """
    if error.size == 0:
        return float("nan")
    pair_codes, _ = _pair_codes(targets, components)
    sums = np.bincount(pair_codes, weights=error)
    counts = np.bincount(pair_codes)
    target_mse = sums / np.maximum(counts, 1)
    owner = _first_owner(pair_codes, components, sums.size)
    component_sums = np.bincount(owner, weights=target_mse)
    component_counts = np.bincount(owner)
    keep = component_counts > 0
    return float(np.mean(component_sums[keep] / component_counts[keep]))


def _pair_codes(targets: np.ndarray, components: np.ndarray):
    combined = components.astype(np.int64) * (int(targets.max()) + 1) + targets
    unique, codes = np.unique(combined, return_inverse=True)
    return codes, unique


def _first_owner(pair_codes: np.ndarray, components: np.ndarray, size: int):
    owner = np.zeros(size, dtype=np.int64)
    owner[pair_codes] = components
    return owner


def metrics_of(rows: Rows) -> dict:
    truth, prediction = rows.delta_y, rows.delta_hat
    error = (prediction - truth) ** 2
    return {
        "n": int(truth.size),
        "mse": float(np.mean(error)) if truth.size else float("nan"),
        "equal_component_target_mean_mse": _equal_component_target_mean_mse(
            error, rows.target_codes, rows.component_codes),
        "pearson": _pearson(truth, prediction),
        "spearman": _spearman(truth, prediction),
        "ci": _concordance(truth, prediction),
        "sign_accuracy": _sign_accuracy(truth, prediction),
    }


STATISTICS = ("mse", "equal_component_target_mean_mse", "pearson", "spearman",
              "ci", "sign_accuracy")
# For these, lower is better; the gate wording is stated in terms of gains.
LOWER_IS_BETTER = ("mse", "equal_component_target_mean_mse")


def _component_groups(rows: Rows) -> list[np.ndarray]:
    order = np.argsort(rows.component_codes, kind="stable")
    codes = rows.component_codes[order]
    boundaries = np.flatnonzero(np.diff(codes)) + 1
    return np.split(order, boundaries)


def _metrics_from_arrays(truth, prediction, targets, components) -> dict:
    error = (prediction - truth) ** 2
    return {
        "mse": float(np.mean(error)) if truth.size else float("nan"),
        "equal_component_target_mean_mse": _equal_component_target_mean_mse(
            error, targets, components),
        "pearson": _pearson(truth, prediction),
        "spearman": _spearman(truth, prediction),
        "ci": _concordance(truth, prediction),
        "sign_accuracy": _sign_accuracy(truth, prediction),
    }


def component_paired_bootstrap(left: Rows, right: Rows, draws: int = 2000,
                               seed: int = 20260819) -> dict:
    """Bootstrap the `left - right` contrast by resampling components.

    `left` and `right` must be prediction rows for the **same** evaluation bank
    in the same order.  The check below is not decoration: an earlier stage in
    this repository compared a condition against itself for weeks because the
    pairing was assumed rather than verified.
    """
    if left.pair_id != right.pair_id:
        raise ValueError("paired bootstrap requires identical evaluation rows")
    groups = _component_groups(left)
    rng = np.random.default_rng(seed)
    samples = {name: np.empty(draws, dtype=np.float64) for name in STATISTICS}
    for draw in range(draws):
        picked = rng.integers(0, len(groups), size=len(groups))
        chosen = [groups[int(i)] for i in picked]
        positions = np.concatenate(chosen)
        # A component drawn twice must count twice in the equal-component mean,
        # so each drawn copy gets its own replicate label instead of collapsing
        # back onto the original component code.
        components = np.repeat(np.arange(len(chosen), dtype=np.int64),
                               [group.size for group in chosen])
        targets = left.target_codes[positions]
        left_metrics = _metrics_from_arrays(
            left.delta_y[positions], left.delta_hat[positions], targets, components)
        right_metrics = _metrics_from_arrays(
            right.delta_y[positions], right.delta_hat[positions], targets, components)
        for name in STATISTICS:
            samples[name][draw] = left_metrics[name] - right_metrics[name]
    point_left, point_right = metrics_of(left), metrics_of(right)
    out: dict = {"components": len(groups), "n": len(left), "draws": draws}
    for name in STATISTICS:
        values = samples[name][np.isfinite(samples[name])]
        lo, hi = ((float(np.quantile(values, 0.025)),
                   float(np.quantile(values, 0.975))) if values.size
                  else (float("nan"), float("nan")))
        out[name] = {
            "left": point_left[name],
            "right": point_right[name],
            "delta": point_left[name] - point_right[name],
            "lo": lo,
            "hi": hi,
            "resolved": bool(np.isfinite(lo) and np.isfinite(hi)
                             and (lo > 0.0 or hi < 0.0)),
            "effective_draws": int(values.size),
        }
    return out
