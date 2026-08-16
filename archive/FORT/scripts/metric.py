"""Component-aware evaluation utilities for the shared SIMA-DTA episode registry."""

from __future__ import annotations

import math
from typing import Mapping, Sequence

import torch

from .contract import Episode


def _device() -> torch.device:
    if not torch.cuda.is_available():
        raise RuntimeError("SIMA-DTA numerical evaluation requires CUDA")
    return torch.device("cuda")


def _rank(values: torch.Tensor) -> torch.Tensor:
    order = torch.argsort(values, stable=True)
    sortedvalues = values[order]
    _, inverse, counts = torch.unique_consecutive(
        sortedvalues, return_inverse=True, return_counts=True
    )
    ends = torch.cumsum(counts, dim=0).to(torch.float64)
    starts = ends - counts
    groupranks = (starts + ends - 1.0) / 2.0
    ranks = torch.empty_like(values, dtype=torch.float64)
    ranks[order] = groupranks[inverse]
    return ranks


def _spearman(actual: torch.Tensor, predicted: torch.Tensor) -> float:
    if (
        len(actual) < 2
        or torch.amax(actual) - torch.amin(actual) == 0
        or torch.amax(predicted) - torch.amin(predicted) == 0
    ):
        return float("nan")
    left = _rank(actual)
    right = _rank(predicted)
    return float(torch.corrcoef(torch.stack((left, right)))[0, 1].item())


def _pairmetrics(actual: torch.Tensor, predicted: torch.Tensor) -> tuple[float, float]:
    left, right = torch.triu_indices(len(actual), len(actual), offset=1, device=actual.device)
    actualdiff = actual[left] - actual[right]
    predicteddiff = predicted[left] - predicted[right]
    usable = actualdiff != 0
    if not bool(usable.any()):
        return float("nan"), float("nan")
    product = actualdiff[usable] * predicteddiff[usable]
    credit = (product > 0).to(torch.float64) + 0.5 * (predicteddiff[usable] == 0)
    accuracy = credit.mean()
    concordance = credit.mean()
    return float(accuracy.item()), float(concordance.item())


def _uncertainty(
    actual: torch.Tensor,
    predicted: torch.Tensor,
    variance: torch.Tensor,
) -> dict[str, float]:
    variance = variance.clamp_min(1e-8)
    error = actual - predicted
    standard = torch.sqrt(variance)
    result = {
        "nll": float(
            (0.5 * (torch.log(2.0 * math.pi * variance) + error.square() / variance))
            .mean()
            .item()
        )
    }
    for label, z in (("50", 0.67448975), ("80", 1.28155157), ("95", 1.95996398)):
        result[f"coverage_{label}"] = float((error.abs() <= z * standard).float().mean().item())
        result[f"width_{label}"] = float((2.0 * z * standard).mean().item())
    centeredprediction = predicted - predicted.mean()
    denominator = centeredprediction.square().sum()
    result["calibration_slope"] = (
        float((centeredprediction * (actual - actual.mean())).sum().div(denominator).item())
        if float(denominator) > 0.0
        else float("nan")
    )
    order = torch.argsort(variance)
    squarederror = error[order].square()
    prefixrmse = torch.sqrt(torch.cumsum(squarederror, dim=0) / torch.arange(
        1, len(squarederror) + 1, device=actual.device, dtype=actual.dtype
    ))
    result["risk_coverage_auc"] = float(prefixrmse.mean().item())
    return result


def evaluateprotocol(
    *,
    predictions: Sequence[float],
    labels: Sequence[float],
    episodes: Sequence[Episode],
    prediction_indices: Sequence[int],
    component_by_target: Mapping[str, str],
    variances: Sequence[float] | None = None,
) -> dict[str, float]:
    """Report target-macro metrics; callers supply only already-authorized labels."""

    if not (len(predictions) == len(labels) == len(prediction_indices)):
        raise ValueError("predictions, labels, and indices must have equal length")
    if variances is not None and len(variances) != len(predictions):
        raise ValueError("variances must match predictions")
    index_to_position = {index: position for position, index in enumerate(prediction_indices)}
    target_scores: list[dict[str, float]] = []
    for episode in episodes:
        positions = [index_to_position[index] for index in episode.query_indices if index in index_to_position]
        if not positions:
            continue
        actual = torch.tensor([labels[position] for position in positions], dtype=torch.float64, device=_device())
        predicted = torch.tensor([predictions[position] for position in positions], dtype=torch.float64, device=_device())
        rmse = float(torch.sqrt(torch.mean((actual - predicted) ** 2)).item())
        mae = float(torch.mean(torch.abs(actual - predicted)).item())
        pairwise, concordance = _pairmetrics(actual, predicted)
        scores = {
            "target_macro_rmse": rmse,
            "target_macro_mae": mae,
            "within_target_spearman": _spearman(actual, predicted),
            "pairwise_accuracy": pairwise,
            "concordance_index": concordance,
        }
        if variances is not None:
            variance = torch.tensor(
                [variances[position] for position in positions],
                dtype=torch.float64,
                device=_device(),
            )
            if not bool(torch.isfinite(variance).all()) or bool((variance < 0).any()):
                raise ValueError("predictive variances must be finite and non-negative")
            scores.update(_uncertainty(actual, predicted, variance))
        target_scores.append(scores)
    if not target_scores:
        raise ValueError("no episode query predictions were supplied")
    keys = target_scores[0].keys()
    result = {
        key: float(
            torch.tensor([score[key] for score in target_scores], device=_device()).nanmean().item()
        )
        for key in keys
    }
    result["independent_components"] = float(
        len({component_by_target[episode.target_key] for episode in episodes})
    )
    return result


def pairedcomponents(
    *,
    predictions: Mapping[str, Sequence[float]],
    labels: Sequence[float],
    episodes: Sequence[Episode],
    prediction_indices: Sequence[int],
    component_by_target: Mapping[str, str],
    reference: str,
    replicates: int = 2000,
    seed: int = 1729,
) -> dict[str, dict[str, dict[str, float | list[float]]]]:
    """Paired homology-component bootstrap for the primary mean/ranking metrics."""

    if reference not in predictions:
        raise ValueError("reference arm is missing")
    if any(len(values) != len(labels) for values in predictions.values()):
        raise ValueError("all prediction arms must align with labels")
    positions = {index: position for position, index in enumerate(prediction_indices)}
    metricnames = ("rmse_gain", "mae_gain", "spearman_gain", "pairwise_gain")
    bycomponent: dict[str, dict[str, list[torch.Tensor]]] = {
        name: {} for name in predictions if name != reference
    }
    for episode in episodes:
        selected = [positions[index] for index in episode.query_indices if index in positions]
        if not selected:
            continue
        actual = torch.tensor(
            [labels[position] for position in selected], dtype=torch.float64, device=_device()
        )
        referencevalues = torch.tensor(
            [predictions[reference][position] for position in selected],
            dtype=torch.float64,
            device=_device(),
        )
        referencermse = torch.sqrt(torch.mean((actual - referencevalues).square()))
        referencemae = torch.mean(torch.abs(actual - referencevalues))
        referencespearman = torch.tensor(
            _spearman(actual, referencevalues), dtype=torch.float64, device=_device()
        )
        referencepairwise = torch.tensor(
            _pairmetrics(actual, referencevalues)[0], dtype=torch.float64, device=_device()
        )
        component = component_by_target[episode.target_key]
        for name, values in predictions.items():
            if name == reference:
                continue
            control = torch.tensor(
                [values[position] for position in selected],
                dtype=torch.float64,
                device=_device(),
            )
            delta = torch.stack(
                (
                    torch.sqrt(torch.mean((actual - control).square())) - referencermse,
                    torch.mean(torch.abs(actual - control)) - referencemae,
                    referencespearman
                    - torch.tensor(
                        _spearman(actual, control), dtype=torch.float64, device=_device()
                    ),
                    referencepairwise
                    - torch.tensor(
                        _pairmetrics(actual, control)[0], dtype=torch.float64, device=_device()
                    ),
                )
            )
            bycomponent[name].setdefault(component, []).append(delta)

    output: dict[str, dict[str, dict[str, float | list[float]]]] = {}
    for name, componentvalues in bycomponent.items():
        matrix = torch.stack(
            [torch.stack(values).mean(dim=0) for values in componentvalues.values()]
        )
        output[name] = {}
        for column, metric in enumerate(metricnames):
            values = matrix[:, column]
            values = values[torch.isfinite(values)]
            if len(values) == 0:
                missing = float("nan")
                output[name][metric] = {
                    "mean": missing,
                    "ci95": [missing, missing],
                    "probability_positive": missing,
                    "largest_positive_component_share": missing,
                    "components": 0.0,
                }
                continue
            generator = torch.Generator(device=_device()).manual_seed(seed + column)
            samples = torch.randint(
                len(values),
                (replicates, len(values)),
                generator=generator,
                device=_device(),
            )
            bootstrap = values[samples].mean(dim=1)
            positive = values.clamp_min(0)
            positivesum = positive.sum()
            dominance = (
                positive.max() / positivesum if float(positivesum) > 0.0 else values.new_tensor(float("nan"))
            )
            output[name][metric] = {
                "mean": float(values.mean().item()),
                "ci95": [
                    float(torch.quantile(bootstrap, 0.025).item()),
                    float(torch.quantile(bootstrap, 0.975).item()),
                ],
                "probability_positive": float((bootstrap > 0).float().mean().item()),
                "largest_positive_component_share": float(dominance.item()),
                "components": float(len(values)),
            }
    return output
