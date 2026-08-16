"""Leakage-safe residual bilinear interaction probe.

This is a diagnostic candidate, not a replacement for the registered model.
It asks the smallest useful question after the identifiability audit: does a
fit-only protein vector explain ligand-pair residuals on unseen homology
components?  The learned operator is exactly antisymmetric in the ligand pair
and is compared with protein-free and homology-matched wrong-protein controls.

No checkpoint, development label, structure, target id, or task scheduler is
used.  All fitted statistics and the bilinear matrix are learned from the
TRAIN fit components only; CUDA is mandatory for numerical paths.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from research.shared.priorgate import limitrows, splitcomponents
from scripts.train import loadlabels


ROOT = Path("dataset/public/chembl_37/processed/dualcold")
REPORT = Path("reports/active/residual_bilinear_probe_pKi_seed1729.json")


@dataclass(frozen=True)
class ProbeEpisode:
    target: str
    component: str
    support: tuple[int, ...]
    query: tuple[int, ...]


def _require_cuda() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("residual bilinear probe requires CUDA in the drug environment")


def _seed_everything(seed: int) -> np.random.Generator:
    torch.manual_seed(seed)
    np.random.seed(seed)
    return np.random.default_rng(seed)


def _fit_feature_projection(
    feature: np.ndarray,
    fit_sources: np.ndarray,
    needed_sources: np.ndarray,
    dimension: int,
    seed: int,
) -> tuple[dict[int, torch.Tensor], dict[str, list[float]]]:
    """Fit normalization on fit rows and project only the authorized rows."""

    fit_values = torch.as_tensor(feature[fit_sources], device="cuda", dtype=torch.float32)
    center = fit_values.mean(dim=0)
    scale = fit_values.std(dim=0).clamp_min(1e-6)
    generator = np.random.default_rng(seed + 211)
    matrix = torch.as_tensor(
        generator.normal(size=(feature.shape[1], dimension)).astype(np.float32)
        / math.sqrt(dimension),
        device="cuda",
    )
    output: dict[int, torch.Tensor] = {}
    unique = np.unique(needed_sources.astype(np.int64))
    for start in range(0, len(unique), 8192):
        sources = unique[start : start + 8192]
        values = torch.as_tensor(feature[sources], device="cuda", dtype=torch.float32)
        projected = ((values - center) / scale) @ matrix
        for source, value in zip(sources.tolist(), projected):
            output[int(source)] = value
    return output, {
        "feature_center": center.detach().cpu().tolist(),
        "feature_scale": scale.detach().cpu().tolist(),
    }


def _fit_base(
    feature: np.ndarray,
    fit: pd.DataFrame,
    gate: pd.DataFrame,
    ridge: float,
) -> tuple[dict[int, torch.Tensor], dict[str, object]]:
    """Fit the ligand-only B0 on fit rows and predict fit/gate rows."""

    if ridge <= 0:
        raise ValueError("ridge must be positive")
    fit_sources = fit.source_row.to_numpy(dtype=np.int64)
    values = torch.as_tensor(feature[fit_sources], device="cuda", dtype=torch.float32)
    center = values.mean(dim=0)
    scale = values.std(dim=0).clamp_min(1e-6)
    standardized = (values - center) / scale
    labels = torch.as_tensor(fit.affinity.to_numpy(dtype=np.float32), device="cuda")
    label_center = labels.mean()
    eye = torch.eye(standardized.shape[1], device="cuda", dtype=standardized.dtype)
    weight = torch.linalg.solve(
        standardized.T @ standardized + ridge * eye,
        standardized.T @ (labels - label_center),
    )
    output: dict[int, torch.Tensor] = {}
    all_sources = np.unique(
        np.concatenate(
            [fit.source_row.to_numpy(dtype=np.int64), gate.source_row.to_numpy(dtype=np.int64)]
        )
    )
    for start in range(0, len(all_sources), 8192):
        sources = all_sources[start : start + 8192]
        batch = torch.as_tensor(feature[sources], device="cuda", dtype=torch.float32)
        prediction = ((batch - center) / scale) @ weight + label_center
        for source, value in zip(sources.tolist(), prediction):
            output[int(source)] = value
    residual = ((standardized @ weight + label_center) - labels).square().mean()
    return output, {
        "ridge": ridge,
        "fit_residual_variance": float(residual.detach().cpu()),
        "feature_center": center.detach().cpu().tolist(),
        "feature_scale": scale.detach().cpu().tolist(),
    }


def _strict_episodes(
    frame: pd.DataFrame,
    support_size: int,
    query_cap: int,
) -> tuple[list[ProbeEpisode], int]:
    """Create target-balanced episodes with scaffold/conn/doc/assay closure."""

    episodes: list[ProbeEpisode] = []
    skipped = 0
    for target, group in frame.groupby("target", sort=True):
        ordered = group.sort_values(
            ["scaffold", "conn", "docs", "assays", "source_row"]
        )
        support: list[int] = []
        used = {column: set() for column in ("scaffold", "conn", "docs", "assays")}
        for row in ordered.itertuples(index=False):
            if any(str(getattr(row, column)) in used[column] for column in used):
                continue
            support.append(int(row.source_row))
            for column in used:
                used[column].add(str(getattr(row, column)))
            if len(support) == support_size:
                break
        if len(support) < support_size:
            skipped += 1
            continue
        query: list[int] = []
        for row in ordered.itertuples(index=False):
            source = int(row.source_row)
            if source in support:
                continue
            if any(str(getattr(row, column)) in used[column] for column in used):
                continue
            query.append(source)
            if len(query) == query_cap:
                break
        if not query:
            skipped += 1
            continue
        episodes.append(
            ProbeEpisode(
                target=str(target),
                component=str(group.hcluster.iloc[0]),
                support=tuple(support),
                query=tuple(query),
            )
        )
    return episodes, skipped


def _pair_rows(
    frame: pd.DataFrame,
    pairs_per_target: int,
    seed: int,
) -> list[tuple[str, int, int]]:
    rng = np.random.default_rng(seed + 307)
    pairs: list[tuple[str, int, int]] = []
    for target, group in frame.groupby("target", sort=True):
        rows = list(group.itertuples(index=False))
        candidates = [
            (left, right)
            for position, left in enumerate(rows)
            for right in rows[position + 1 :]
            if left.scaffold != right.scaffold and left.conn != right.conn
        ]
        if not candidates:
            continue
        count = min(pairs_per_target, len(candidates))
        selected = rng.choice(len(candidates), size=count, replace=False)
        pairs.extend((str(target), candidates[int(index)][0].source_row, candidates[int(index)][1].source_row) for index in selected)
    if not pairs:
        raise RuntimeError("fit role has no valid ligand pairs")
    return pairs


def _fit_bilinear(
    pairs: list[tuple[str, int, int]],
    residual: dict[int, torch.Tensor],
    ligand: dict[int, torch.Tensor],
    protein: dict[str, torch.Tensor],
    rank: int,
    ridge: float,
) -> torch.Tensor:
    rows: list[torch.Tensor] = []
    labels: list[torch.Tensor] = []
    for target, left, right in pairs:
        delta = ligand[right] - ligand[left]
        design = torch.outer(protein[target], delta).reshape(-1)
        rows.append(design)
        labels.append(residual[right] - residual[left])
    design = torch.stack(rows)
    target = torch.stack(labels)
    eye = torch.eye(rank * rank, device="cuda", dtype=design.dtype)
    coefficients = torch.linalg.solve(
        design.T @ design + ridge * eye,
        design.T @ target,
    )
    return coefficients.reshape(rank, rank)


def _rank(values: torch.Tensor) -> torch.Tensor:
    order = torch.argsort(values, stable=True)
    result = torch.empty_like(values, dtype=torch.float64)
    result[order] = torch.arange(len(values), device=values.device, dtype=torch.float64)
    return result


def _target_metrics(actual: torch.Tensor, prediction: torch.Tensor) -> dict[str, float]:
    error = actual - prediction
    actual_rank = _rank(actual)
    prediction_rank = _rank(prediction)
    spearman = float(torch.corrcoef(torch.stack((actual_rank, prediction_rank)))[0, 1].cpu()) if len(actual) > 1 and actual.unique().numel() > 1 and prediction.unique().numel() > 1 else float("nan")
    left, right = torch.triu_indices(len(actual), len(actual), offset=1, device=actual.device)
    actual_delta = actual[left] - actual[right]
    prediction_delta = prediction[left] - prediction[right]
    usable = actual_delta != 0
    pairwise = float(((actual_delta[usable] * prediction_delta[usable]) > 0).float().mean().cpu()) if bool(usable.any()) else float("nan")
    return {
        "rmse": float(torch.sqrt(error.square().mean()).cpu()),
        "mae": float(error.abs().mean().cpu()),
        "spearman": spearman,
        "pairwise": pairwise,
    }


def _bootstrap_gain(
    arm: dict[str, dict[str, float]],
    reference: dict[str, dict[str, float]],
    components: dict[str, str],
    metric: str,
    higher_is_better: bool,
    seed: int,
    draws: int = 2000,
) -> dict[str, float | list[float] | None]:
    by_component: dict[str, list[float]] = {}
    for target, values in arm.items():
        if target not in reference or not np.isfinite(values[metric]) or not np.isfinite(reference[target][metric]):
            continue
        raw = values[metric] - reference[target][metric] if higher_is_better else reference[target][metric] - values[metric]
        by_component.setdefault(components[target], []).append(raw)
    means = np.asarray([np.mean(values) for values in by_component.values()], dtype=np.float64)
    if len(means) == 0:
        return {"mean": None, "ci95": [None, None], "components": 0}
    rng = np.random.default_rng(seed)
    sample = means[rng.integers(0, len(means), size=(draws, len(means)))].mean(axis=1)
    return {
        "mean": float(means.mean()),
        "ci95": [float(np.quantile(sample, 0.025)), float(np.quantile(sample, 0.975))],
        "components": int(len(means)),
    }


def _evaluate(
    episodes: list[ProbeEpisode],
    frame: pd.DataFrame,
    base: dict[int, torch.Tensor],
    ligand: dict[int, torch.Tensor],
    protein: dict[str, torch.Tensor],
    bilinear: torch.Tensor,
    matched: dict[str, str | None],
    random_wrong: dict[str, str],
) -> tuple[dict[str, dict[str, dict[str, float]]], dict[str, dict[str, object]], dict[str, int]]:
    rows = frame.set_index("source_row")
    arm_names = ("b0", "calibration", "proteinfree", "correct", "matched_wrong", "random_wrong")
    by_arm: dict[str, dict[str, dict[str, float]]] = {name: {} for name in arm_names}
    valid: dict[str, int] = {name: 0 for name in arm_names}
    for episode in episodes:
        support = list(episode.support)
        query = list(episode.query)
        support_y = torch.as_tensor(rows.loc[support, "affinity"].to_numpy(dtype=np.float32), device="cuda")
        query_y = torch.as_tensor(rows.loc[query, "affinity"].to_numpy(dtype=np.float32), device="cuda")
        support_base = torch.stack([base[source] for source in support])
        query_base = torch.stack([base[source] for source in query])
        residual = support_y - support_base
        query_ligand = torch.stack([ligand[source] for source in query])
        support_ligand = torch.stack([ligand[source] for source in support])
        arms: dict[str, torch.Tensor] = {
            "b0": query_base,
            "calibration": query_base + residual.mean(),
            "proteinfree": query_base + residual.mean(),
        }
        for name, target in (
            ("correct", episode.target),
            ("matched_wrong", matched.get(episode.target)),
            ("random_wrong", random_wrong.get(episode.target)),
        ):
            if target is None or target not in protein:
                continue
            coefficient = protein[target] @ bilinear
            deltas = (
                (query_ligand[:, None, :] - support_ligand[None, :, :])
                * coefficient[None, None, :]
            ).sum(dim=-1)
            arms[name] = query_base + (residual[None, :] + deltas).mean(dim=1)
        for name, prediction in arms.items():
            by_arm[name][episode.target] = _target_metrics(query_y, prediction)
            valid[name] += 1
    components = {episode.target: episode.component for episode in episodes}
    summary: dict[str, dict[str, object]] = {}
    for name, targets in by_arm.items():
        if not targets:
            summary[name] = {"targets": 0}
            continue
        summary[name] = {
            "targets": len(targets),
            "target_macro_rmse": float(np.nanmean([value["rmse"] for value in targets.values()])),
            "target_macro_mae": float(np.nanmean([value["mae"] for value in targets.values()])),
            "within_target_spearman": float(np.nanmean([value["spearman"] for value in targets.values()])),
            "pairwise_accuracy": float(np.nanmean([value["pairwise"] for value in targets.values()])),
        }
    reference = by_arm["proteinfree"]
    gains: dict[str, dict[str, object]] = {}
    for name in ("correct", "matched_wrong", "random_wrong"):
        gains[name] = {
            metric: _bootstrap_gain(
                by_arm[name], reference, components, metric,
                higher_is_better=metric in {"spearman", "pairwise"}, seed=1729 + index,
            )
            for index, metric in enumerate(("rmse", "mae", "spearman", "pairwise"))
        }
    gains["correct_vs_matched_wrong"] = {
        metric: _bootstrap_gain(
            by_arm["correct"], by_arm["matched_wrong"], components, metric,
            higher_is_better=metric in {"spearman", "pairwise"}, seed=1801 + index,
        )
        for index, metric in enumerate(("rmse", "mae", "spearman", "pairwise"))
    }
    return by_arm, {"metrics": summary, "gains_vs_proteinfree": gains}, valid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--endpoint", choices=("pKi", "pKd"), default="pKi")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--holdout", type=float, default=0.2)
    parser.add_argument("--support", type=int, default=5)
    parser.add_argument("--querycap", type=int, default=64)
    parser.add_argument("--pairs-per-target", type=int, default=32)
    parser.add_argument("--fit-targets", type=int, default=0)
    parser.add_argument("--gate-targets", type=int, default=0)
    parser.add_argument("--ligand-dim", type=int, default=16)
    parser.add_argument("--protein-dim", type=int, default=16)
    parser.add_argument("--ridge", type=float, default=10.0)
    parser.add_argument("--out", type=Path, default=REPORT)
    args = parser.parse_args()
    _require_cuda()
    if args.support <= 0 or args.querycap <= 0 or args.pairs_per_target <= 0:
        raise ValueError("support, querycap, and pairs-per-target must be positive")
    if args.ligand_dim <= 0 or args.protein_dim <= 0:
        raise ValueError("projection dimensions must be positive")
    if args.ligand_dim != args.protein_dim:
        raise ValueError("the diagnostic uses a square bilinear operator")
    started = time.perf_counter()
    rng = _seed_everything(args.seed)
    full = loadlabels(args.root, args.endpoint)
    train = full.loc[full.dual_cold_split == "train"].copy()
    fit, gate = splitcomponents(train, args.holdout, args.seed)
    if args.fit_targets:
        fit = limitrows(fit, args.fit_targets, 0, args.seed + 401)
    if args.gate_targets:
        gate = limitrows(gate, args.gate_targets, 0, args.seed + 402)
    fit = fit.reset_index(drop=True)
    gate = gate.reset_index(drop=True)
    feature = np.load(args.root / "ligand_features.npz", allow_pickle=False)["feat"]
    fit_sources = fit.source_row.to_numpy(dtype=np.int64)
    needed_sources = np.concatenate(
        [fit_sources, gate.source_row.to_numpy(dtype=np.int64)]
    )
    ligand, ligand_norm = _fit_feature_projection(
        feature, fit_sources, needed_sources, args.ligand_dim, args.seed
    )
    base, base_info = _fit_base(feature, fit, gate, args.ridge)
    fit_targets = sorted(fit.target.astype(str).unique())
    # Gate targets are unseen in the fit role but their label-free protein vectors
    # are valid inputs.  Refit normalization is forbidden, so use fit moments.
    payload = np.load(args.root / "target_esm2.npz", allow_pickle=False)
    keys = [str(value) for value in payload["keys"]]
    available = {key: index for index, key in enumerate(keys)}
    fit_raw = torch.as_tensor(
        payload["pooled"][[available[target] for target in fit_targets]],
        device="cuda",
        dtype=torch.float32,
    )
    pcenter = fit_raw.mean(dim=0)
    pscale = fit_raw.std(dim=0).clamp_min(1e-6)
    matrix = torch.as_tensor(
        np.random.default_rng(args.seed + 101).normal(
            size=(payload["pooled"].shape[1], args.protein_dim)
        ).astype(np.float32)
        / math.sqrt(args.protein_dim),
        device="cuda",
    )
    all_raw = torch.as_tensor(payload["pooled"], device="cuda", dtype=torch.float32)
    projected = ((all_raw - pcenter) / pscale) @ matrix
    protein = {key: projected[index] for index, key in enumerate(keys)}
    residual = {
        int(source): torch.as_tensor(float(row.affinity), device="cuda") - base[int(source)]
        for row in fit.itertuples(index=False)
        for source in [row.source_row]
    }
    pairs = _pair_rows(fit, args.pairs_per_target, args.seed)
    bilinear = _fit_bilinear(
        pairs, residual, ligand, protein, args.ligand_dim, args.ridge
    )
    episodes, skipped = _strict_episodes(gate, args.support, args.querycap)
    if len(episodes) < 2:
        raise RuntimeError("strict document/assay closure left too few gate episodes")
    gate_targets = sorted({episode.target for episode in episodes})
    vectors = {target: protein[target] for target in gate_targets}
    component_by_target = gate.drop_duplicates("target").set_index("target").hcluster.astype(str).to_dict()
    matched: dict[str, str | None] = {}
    random_wrong: dict[str, str] = {}
    for target in gate_targets:
        same = [candidate for candidate in gate_targets if candidate != target and component_by_target[candidate] == component_by_target[target]]
        matched[target] = min(same, key=lambda candidate: float(torch.linalg.vector_norm(vectors[target] - vectors[candidate]).detach().cpu())) if same else None
        other = [candidate for candidate in gate_targets if component_by_target[candidate] != component_by_target[target]]
        if other:
            random_wrong[target] = other[int(rng.integers(len(other)))]
    by_arm, evaluation, valid = _evaluate(
        episodes, gate, base, ligand, protein, bilinear, matched, random_wrong
    )
    matched_count = sum(value is not None for value in matched.values())
    def positive_ci(gains: dict[str, object]) -> bool:
        lowers = [gains[metric]["ci95"][0] for metric in ("rmse", "mae", "spearman", "pairwise")]
        return bool(lowers) and all(value is not None and value > 0 for value in lowers)

    proteinfree_gain = evaluation["gains_vs_proteinfree"]["correct"]
    matched_gain = evaluation["gains_vs_proteinfree"]["correct_vs_matched_wrong"]
    result = {
        "schema": "idg-rbp-residual-bilinear-probe-v1",
        "protocol": "TRAIN-only fit-component residual antisymmetric bilinear diagnostic",
        "status": "diagnostic_only_not_promoted",
        "endpoint": args.endpoint,
        "seed": args.seed,
        "support": args.support,
        "fit_targets": int(fit.target.nunique()),
        "fit_target_cap": args.fit_targets or None,
        "fit_components": int(fit.hcluster.nunique()),
        "gate_targets": int(gate.target.nunique()),
        "gate_target_cap": args.gate_targets or None,
        "gate_components": int(gate.hcluster.nunique()),
        "strict_gate_episodes": len(episodes),
        "skipped_gate_targets": skipped,
        "matched_wrong_targets": matched_count,
        "random_wrong_targets": len(random_wrong),
        "pairs": len(pairs),
        "ligand_dim": args.ligand_dim,
        "protein_dim": args.protein_dim,
        "ridge": args.ridge,
        "closure": ["scaffold", "ligand_connectivity", "document", "assay"],
        "checkpoint": None,
        "feature_projection": "fit-only normalization plus seeded Gaussian projection",
        "ligand_normalization": ligand_norm,
        "base": base_info,
        "bilinear_parameter_count": int(args.ligand_dim * args.protein_dim),
        "bilinear_frobenius_norm": float(torch.linalg.vector_norm(bilinear).detach().cpu()),
        "valid_episode_counts": valid,
        "evaluation": evaluation,
        "training_seconds": time.perf_counter() - started,
        "gpu": torch.cuda.get_device_name(0),
        "peak_torch_memory_mib": float(torch.cuda.max_memory_allocated() / 2**20),
        "admission": {
            "correct_beats_proteinfree_rmse_ci": bool(
                proteinfree_gain["rmse"]["ci95"][0] is not None
                and proteinfree_gain["rmse"]["ci95"][0] > 0
            ),
            "correct_beats_proteinfree_all_metrics": positive_ci(proteinfree_gain),
            "correct_beats_matched_wrong": positive_ci(matched_gain),
            "requires_matched_wrong_positive_ci": True,
            "promote": False,
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=True), encoding="utf-8")
    print(json.dumps({"evaluation": evaluation, "out": str(args.out)}, indent=2), flush=True)


if __name__ == "__main__":
    main()
