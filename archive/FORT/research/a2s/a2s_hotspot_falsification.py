"""Source-only falsification of heterogeneous coordinate sparsity.

This is a representation diagnostic, not an adaptation mechanism.  Every
basis statistic is fitted on the source ``fit`` role.  Target heads are then
estimated and scored on scaffold-disjoint splits of completely component-
held-out ``probe`` targets.  ``locked`` and recipient roles are never loaded.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch

from research.a2s.a2s_mode_gates import fit_head, source_heads
from research.a2s.a2s_mode_generalization import build_subspace
from research.a2s.a2s_trace import DEFAULT_LOCK, DEFAULT_OOF, DEVICE, Substrate, load_substrate
from research.a2s.a2s_trace_stratum import metric_loss, paired_bootstrap


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_hotspot_falsification_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_hotspot_falsification_records_2026-08-02.parquet"

SPLIT_SEEDS = (1729, 1731, 1733)
ROTATION_SEEDS = (2719, 2729, 2741)
S_SWEEP = (1, 2, 3, 5, 8)
MIN_TARGET_ROWS = 40
MIN_HEAD_TRAIN = 20
MIN_EVAL_ROWS = 8
BOOTSTRAP_DRAWS = 2000
RETENTION_TARGET = 0.60


@dataclass(frozen=True)
class Split:
    target: str
    component: str
    seed: int
    train_rows: np.ndarray
    eval_rows: np.ndarray


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def standardise(values: np.ndarray, fit_mask: np.ndarray) -> np.ndarray:
    mean = values[fit_mask].mean(axis=0, keepdims=True)
    scale = values[fit_mask].std(axis=0, keepdims=True)
    scale[scale < 1e-6] = 1.0
    return ((values - mean) / scale).astype(np.float64)


def pca_basis(values: torch.Tensor, fit_mask: torch.Tensor, dimension: int) -> np.ndarray:
    fit = values[fit_mask]
    centre = fit.mean(dim=0, keepdim=True)
    _, _, right = torch.svd_lowrank(fit - centre, q=dimension + 8, niter=4)
    scores = (values - centre) @ right[:, :dimension]
    fit_scores = scores[fit_mask]
    scores = (scores - fit_scores.mean(0, keepdim=True)) / fit_scores.std(0, keepdim=True).clamp(min=1e-6)
    return scores.cpu().numpy().astype(np.float64)


def random_projection_basis(
    values: torch.Tensor, fit_mask: torch.Tensor, dimension: int, seed: int
) -> np.ndarray:
    generator = torch.Generator(device="cpu").manual_seed(seed)
    projection = torch.randn(values.shape[1], dimension, generator=generator, dtype=torch.float32)
    projection, _ = torch.linalg.qr(projection, mode="reduced")
    projection = projection.to(values.device)
    scores = values @ projection
    fit_scores = scores[fit_mask]
    scores = (scores - fit_scores.mean(0, keepdim=True)) / fit_scores.std(0, keepdim=True).clamp(min=1e-6)
    return scores.cpu().numpy().astype(np.float64)


def pharmacophore_counts(smiles: Iterable[str]) -> np.ndarray:
    """Eight RDKit feature-family counts, cached per unique molecule in memory."""

    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures

    families = (
        "Donor",
        "Acceptor",
        "Aromatic",
        "Hydrophobe",
        "LumpedHydrophobe",
        "PosIonizable",
        "NegIonizable",
        "ZnBinder",
    )
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    cache: dict[str, np.ndarray] = {}
    rows: list[np.ndarray] = []
    for value in smiles:
        key = str(value)
        if key not in cache:
            molecule = Chem.MolFromSmiles(key)
            counts = np.zeros(len(families), dtype=np.float64)
            if molecule is not None:
                lookup = {name: index for index, name in enumerate(families)}
                for feature in factory.GetFeaturesForMol(molecule):
                    index = lookup.get(feature.GetFamily())
                    if index is not None:
                        counts[index] += 1.0
            cache[key] = counts
        rows.append(cache[key])
    return np.vstack(rows)


def build_bases(substrate: Substrate) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    fit_mask_np = (substrate.labeled.role == "fit").to_numpy()
    fit_mask = torch.as_tensor(fit_mask_np, device=DEVICE)
    desc = substrate.desc.cpu().numpy().astype(np.float64)
    morgan26 = pca_basis(substrate.bits, fit_mask, 26)
    morgan16 = morgan26[:, :16]
    original = np.concatenate((desc, morgan16), axis=1)
    pharm = pharmacophore_counts(substrate.labeled.conn.astype(str))
    pharm_desc = standardise(np.concatenate((desc, pharm), axis=1), fit_mask_np)
    bases = {
        "original26": original,
        "descriptors10": standardise(desc, fit_mask_np),
        "morgan_pca26": morgan26,
        "pharm_desc18": pharm_desc,
    }
    for seed in ROTATION_SEEDS:
        rng = np.random.default_rng(seed)
        rotation, _ = np.linalg.qr(rng.normal(size=(original.shape[1], original.shape[1])))
        bases[f"original_rot{seed}"] = original @ rotation
    for seed in ROTATION_SEEDS:
        bases[f"morgan_rp26_{seed}"] = random_projection_basis(substrate.bits, fit_mask, 26, seed)
    metadata = {
        "roles_used": ["fit", "probe"],
        "locked_requested": False,
        "basis_dimensions": {name: int(values.shape[1]) for name, values in bases.items()},
        "statistics_fitted_on": "fit role only",
        "pharmacophore_families": 8,
    }
    return bases, metadata


def target_splits(substrate: Substrate, seed: int) -> list[Split]:
    splits: list[Split] = []
    frame = substrate.labeled.loc[substrate.labeled.role == "probe"]
    for target, group in frame.groupby("target", sort=True):
        if len(group) < MIN_TARGET_ROWS:
            continue
        rows = group.index.to_numpy()
        scaffolds = group.scaffold.astype(str).to_numpy()
        unique = np.asarray(sorted(set(scaffolds)), dtype=object)
        digest = int(sha256(f"{seed}:{target}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(digest)
        held: set[str] = set()
        for index in rng.permutation(len(unique)):
            held.add(str(unique[index]))
            mask = np.isin(scaffolds, list(held))
            if mask.sum() >= max(MIN_EVAL_ROWS, int(0.3 * len(rows))):
                break
        mask = np.isin(scaffolds, list(held))
        if mask.sum() < MIN_EVAL_ROWS or (~mask).sum() < MIN_HEAD_TRAIN:
            continue
        splits.append(
            Split(
                target=str(target),
                component=str(group.component.iloc[0]),
                seed=seed,
                train_rows=rows[~mask],
                eval_rows=rows[mask],
            )
        )
    return splits


def top_coordinates(head: np.ndarray, size: int) -> np.ndarray:
    return np.argsort(-np.abs(head))[: min(size, len(head))]


def truncate_coordinates(head: np.ndarray, size: int) -> np.ndarray:
    result = np.zeros_like(head)
    selected = top_coordinates(head, size)
    result[selected] = head[selected]
    return result


def evaluate_basis(
    substrate: Substrate, basis_name: str, basis: np.ndarray, splits: list[Split]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    heads, _, sigma, _ = source_heads(substrate, basis)
    subspace = build_subspace(heads, sigma)
    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    affinity = substrate.affinity.cpu().numpy().astype(np.float64)
    base_all = substrate.base.cpu().numpy().astype(np.float64)
    records: list[dict[str, object]] = []
    supports: list[dict[str, object]] = []
    for split in splits:
        train_design = basis[split.train_rows]
        eval_design = basis[split.eval_rows]
        labels = affinity[split.eval_rows]
        base = base_all[split.eval_rows]
        if float(np.std(labels)) < 1e-9:
            continue
        head, level = fit_head(train_design, residual[split.train_rows])
        base_ci = float(metric_loss(labels, base)["ci"])
        full_ci = float(metric_loss(labels, base + eval_design @ head + level)["ci"])
        centred = head - subspace.mean_head
        rng = np.random.default_rng(
            int(sha256(f"{basis_name}:{split.seed}:{split.target}".encode()).hexdigest()[:8], 16)
        )
        for size in sorted(set((*S_SWEEP, basis.shape[1]))):
            size = min(size, basis.shape[1])
            coord = truncate_coordinates(head, size)
            directions = subspace.directions[:, :size]
            rank = subspace.mean_head + directions @ (directions.T @ centred)
            random_predictions = []
            for _ in range(8):
                selected = rng.choice(basis.shape[1], size=size, replace=False)
                random_head = np.zeros_like(head)
                random_head[selected] = head[selected]
                random_predictions.append(eval_design @ random_head)
            predictions = {
                "coord": base + eval_design @ coord + level,
                "rank": base + eval_design @ rank + level,
                "random": base + np.mean(random_predictions, axis=0) + level,
            }
            for method, prediction in predictions.items():
                records.append(
                    {
                        "basis": basis_name,
                        "target": split.target,
                        "component": split.component,
                        "seed": split.seed,
                        "size": size,
                        "method": method,
                        "base_ci": base_ci,
                        "full_ci": full_ci,
                        "ci": float(metric_loss(labels, prediction)["ci"]),
                    }
                )
        for size in (2, 3, min(8, basis.shape[1])):
            supports.append(
                {
                    "basis": basis_name,
                    "target": split.target,
                    "component": split.component,
                    "seed": split.seed,
                    "size": size,
                    "coordinates": tuple(sorted(int(value) for value in top_coordinates(head, size))),
                    "mass": float(np.abs(head)[top_coordinates(head, size)].sum() / np.abs(head).sum()),
                }
            )
    return pd.DataFrame.from_records(records), pd.DataFrame.from_records(supports)


def component_bootstrap(frame: pd.DataFrame, value: str) -> dict[str, float]:
    target_mean = (
        frame.groupby(["component", "target"], as_index=False)[value]
        .mean()
    )
    return paired_bootstrap(target_mean, value, draws=BOOTSTRAP_DRAWS)


def support_stability(frame: pd.DataFrame) -> dict[str, float]:
    values: list[float] = []
    for _, group in frame.groupby(["basis", "target", "size"], sort=True):
        sets = [set(value) for value in group.coordinates]
        for left in range(len(sets)):
            for right in range(left + 1, len(sets)):
                union = sets[left] | sets[right]
                values.append(len(sets[left] & sets[right]) / len(union) if union else 1.0)
    return {
        "mean": float(np.mean(values)) if values else float("nan"),
        "median": float(np.median(values)) if values else float("nan"),
        "pairs": int(len(values)),
    }


def summarise(records: pd.DataFrame, supports: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for basis_name, basis_frame in records.groupby("basis", sort=True):
        cell: dict[str, object] = {"sizes": {}}
        full_frame = basis_frame.drop_duplicates(["target", "seed"])[
            ["component", "target", "seed", "base_ci", "full_ci"]
        ].copy()
        full_frame["gain"] = full_frame.full_ci - full_frame.base_ci
        full = component_bootstrap(full_frame, "gain")
        cell["full_gain"] = full
        effective: int | None = None
        for size, size_frame in basis_frame.groupby("size", sort=True):
            methods: dict[str, object] = {}
            for method, method_frame in size_frame.groupby("method", sort=True):
                working = method_frame.copy()
                working["gain"] = working.ci - working.base_ci
                result = component_bootstrap(working, "gain")
                result["retained_fraction"] = (
                    float(result["mean"] / full["mean"]) if abs(full["mean"]) > 1e-9 else None
                )
                methods[method] = result
            cell["sizes"][f"s{int(size)}"] = methods
            coord_fraction = methods.get("coord", {}).get("retained_fraction")
            if effective is None and coord_fraction is not None and coord_fraction >= RETENTION_TARGET:
                effective = int(size)
        cell["effective_s60"] = effective
        cell["stability"] = {
            f"s{int(size)}": support_stability(group)
            for size, group in supports.loc[supports.basis == basis_name].groupby("size", sort=True)
        }
        cell["mean_top3_mass"] = float(
            supports.loc[
                (supports.basis == basis_name) & (supports["size"] == 3), "mass"
            ].mean()
        )
        summary[basis_name] = cell
    return summary


def decide(summary: dict[str, object]) -> dict[str, object]:
    original = summary["original26"]
    coord8 = original["sizes"]["s8"]["coord"]
    rank8 = original["sizes"]["s8"]["rank"]
    rotated = [value for name, value in summary.items() if name.startswith("original_rot")]
    rotation_fraction = float(
        np.mean([cell["sizes"]["s8"]["coord"]["retained_fraction"] for cell in rotated])
    )
    reproducible = bool(
        coord8["lower95"] > 0.0
        and coord8["mean"] > rank8["mean"]
        and coord8["retained_fraction"] > rotation_fraction
    )
    compact = bool((original.get("effective_s60") or 999) <= 5)
    stable = bool(original["stability"]["s3"]["mean"] >= 0.5)
    return {
        "H1_reproducible_coordinate_advantage": reproducible,
        "H1_compact_enough_for_k5": compact,
        "H1_support_stable": stable,
        "original_s8_retained_fraction": coord8["retained_fraction"],
        "mean_rotated_s8_retained_fraction": rotation_fraction,
        "original_effective_s60": original.get("effective_s60"),
        "original_s3_jaccard": original["stability"]["s3"],
        "verdict": (
            "HETEROGENEOUS_SPARSITY_SUPPORTED_BUT_NOT_K5_IDENTIFIABLE"
            if reproducible and not (compact and stable)
            else "HETEROGENEOUS_SPARSITY_ADMITTED"
            if reproducible
            else "HETEROGENEOUS_SPARSITY_NOT_REPRODUCED"
        ),
    }


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    substrate, context = load_substrate(lock_path, oof_cache)
    if set(substrate.labeled.role.unique()) - {"fit", "probe"}:
        raise AssertionError("an unauthorized role entered the hotspot gate")
    bases, basis_metadata = build_bases(substrate)
    splits = [split for seed in SPLIT_SEEDS for split in target_splits(substrate, seed)]
    all_records: list[pd.DataFrame] = []
    all_supports: list[pd.DataFrame] = []
    for name, basis in bases.items():
        records, supports = evaluate_basis(substrate, name, basis, splits)
        all_records.append(records)
        all_supports.append(supports)
    records = pd.concat(all_records, ignore_index=True)
    supports = pd.concat(all_supports, ignore_index=True)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    summary = summarise(records, supports)
    result = {
        "schema": "a2s-hotspot-falsification-v1",
        "status": "SOURCE_ONLY",
        "protocol": {
            "split_seeds": list(SPLIT_SEEDS),
            "rotation_seeds": list(ROTATION_SEEDS),
            "s_sweep": list(S_SWEEP),
            "retention_target": RETENTION_TARGET,
            "basis": basis_metadata,
            "aggregation": "split-seed mean within target, then component bootstrap",
        },
        "data": {
            "targets": int(records.target.nunique()),
            "components": int(records.component.nunique()),
            "roles_opened": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
            "source_context": context,
        },
        "summary": summary,
        "decision": decide(summary),
        "records": str(records_path.relative_to(ROOT)).replace("\\", "/"),
    }
    payload = canonical(result)
    result["content_sha256"] = sha256(payload.encode()).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    args = parser.parse_args()
    result = run(args.lock, args.output, args.records, args.oof_cache)
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
