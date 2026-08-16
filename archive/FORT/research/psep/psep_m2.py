"""M2: can k support labels identify the object, and is it low-rank?

M4 established that a target-specific adaptation object exists and transfers
across documents: at d=266 the unit's own head beats a target-agnostic head of
identical capacity by +0.019 and a wrong-target head by +0.021, under
simultaneous scaffold/document/assay separation.  Those heads were fitted on the
unit's *entire* training set -- on average well over a hundred labels.  Few-shot
adaptation gets k in {1,3,5}.

Two quantities decide whether a few-shot mechanism is possible at all.

**A. The label-budget curve.**  The target-specific gain as a function of k.
The object needs d~266 dimensions to be visible, and k labels can pin down at
most k directions, so the curve is the direct empirical statement of the
rank/sample-complexity obstruction.  Support rows are drawn from the unit's
training documents only, so the separation guarantee is untouched.  Many draws
per unit, and the ridge is swept and taken at its best value for each k --
optimistic, which is the right direction when the risk is a false negative.

**B. The shared-basis falsifier.**  If the per-target heads of *different*
targets lie near a common m-dimensional subspace, then k labels need only locate
a coordinate in that subspace rather than a full 266-dimensional vector, and
few-shot adaptation becomes well-posed.  This is the one lever the A2S programme
left untried, and it registered its own falsifier: the retained fraction of the
gain under a rank-m projection must rise well above the -6 % measured earlier.

The subspace is estimated by principal components of the source heads with the
evaluated component *held out*, so a component never contributes to the basis
used to compress its own head.

Reads the `discover` role only.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits, pair_concordance, paired_bootstrap
from research.psep.psep_m0 import rich_basis
from research.psep.psep_m4 import CAPACITY, GLOBAL_RIDGE, RIDGE, fit_global_heads

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_m2_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_m2_records_2026-08-02.parquet"

K_SWEEP = (1, 3, 5, 10, 20, 50, 100, 200)
K_RIDGES = (10.0, 100.0, 1000.0, 10000.0)
DRAWS = 16
RANKS = (1, 2, 4, 8, 16, 32, 64, 128, 266)


def ridge_fit(design: np.ndarray, target: np.ndarray, ridges: tuple[float, ...]) -> np.ndarray:
    """Ridge weights for several penalties.  Uses the dual form when n < d, which
    is the regime every small-k fit lives in."""

    centre = design.mean(axis=0)
    centred = design - centre
    shifted = target - target.mean()
    rows, columns = centred.shape
    if rows <= columns:
        gram = centred @ centred.T
        values, vectors = np.linalg.eigh(gram)
        projected = vectors.T @ shifted
        return np.stack([
            centred.T @ (vectors @ (projected / (values + ridge))) for ridge in ridges
        ])
    gram = centred.T @ centred
    values, vectors = np.linalg.eigh(gram)
    projected = vectors.T @ (centred.T @ shifted)
    return np.stack([vectors @ (projected / (values + ridge)) for ridge in ridges])


# --------------------------------------------------------------------------- #
# A -- label-budget curve
# --------------------------------------------------------------------------- #


def label_curve(basis: np.ndarray, substrate, splits: list) -> pd.DataFrame:
    documents = substrate.rows.docs.astype(str).to_numpy()
    design = basis[:, :CAPACITY]
    globals_ = fit_global_heads(basis, substrate)
    folds = substrate.rows.oof_fold.to_numpy()

    records: list[dict[str, object]] = []
    for number, split in enumerate(splits):
        train, evaluation = split.train, split.evaluation
        evaluation_docs = documents[evaluation]
        label = substrate.affinity[evaluation]
        base = substrate.base[evaluation]
        reference = pair_concordance(label, base, evaluation_docs)
        if not np.isfinite(reference["ci_within"]):
            continue
        global_weight = globals_[int(folds[evaluation[0]])]
        global_prediction = design[evaluation] @ global_weight
        global_score = pair_concordance(label, base + global_prediction, evaluation_docs)["ci_within"]

        rng = np.random.default_rng(int(sha256(f"{SEED}:m2:{split.unit}".encode()).hexdigest()[:8], 16))
        row: dict[str, object] = {
            "unit": split.unit, "component": split.component, "endpoint": split.endpoint,
            "n_train": int(len(train)), "n_eval": int(len(evaluation)),
            "base__ci_within": reference["ci_within"], "global__ci_within": global_score,
        }
        for k in K_SWEEP:
            if len(train) < k:
                continue
            scores = {ridge: [] for ridge in K_RIDGES}
            draws = DRAWS if k <= 50 else max(4, DRAWS // 4)
            for _ in range(draws):
                support = rng.choice(train, size=k, replace=False)
                if k == 1:
                    # A single label carries no gradient information: the head is
                    # exactly zero after centring.  Recorded, not skipped.
                    for ridge in K_RIDGES:
                        scores[ridge].append(global_score)
                    continue
                weights = ridge_fit(design[support], substrate.residual[support], K_RIDGES)
                for ridge, weight in zip(K_RIDGES, weights):
                    value = pair_concordance(
                        label, base + design[evaluation] @ weight, evaluation_docs
                    )["ci_within"]
                    scores[ridge].append(value)
            for ridge in K_RIDGES:
                row[f"k{k}__r{ridge:g}"] = float(np.nanmean(scores[ridge])) if scores[ridge] else np.nan
        records.append(row)
        if number % 100 == 0:
            print(f"  curve {number}/{len(splits)}", flush=True)
    return pd.DataFrame.from_records(records)


# --------------------------------------------------------------------------- #
# B -- shared low-rank basis
# --------------------------------------------------------------------------- #


def rank_probe(basis: np.ndarray, substrate, splits: list) -> tuple[pd.DataFrame, dict[str, object]]:
    documents = substrate.rows.docs.astype(str).to_numpy()
    design = basis[:, :CAPACITY]
    globals_ = fit_global_heads(basis, substrate)
    folds = substrate.rows.oof_fold.to_numpy()

    heads, keep = [], []
    for split in splits:
        if len(split.train) < 20:
            continue
        heads.append(ridge_fit(design[split.train], substrate.residual[split.train], (RIDGE,))[0])
        keep.append(split)
    matrix = np.stack(heads)
    centre = matrix.mean(axis=0, keepdims=True)
    spectrum = np.linalg.svd(matrix - centre, compute_uv=False)
    variance = spectrum ** 2
    spectrum_stats = {
        "heads": int(len(matrix)),
        "dimension": int(matrix.shape[1]),
        "variance_top1": float(variance[0] / variance.sum()),
        "variance_top2": float(variance[:2].sum() / variance.sum()),
        "variance_top3": float(variance[:3].sum() / variance.sum()),
        "variance_top16": float(variance[:16].sum() / variance.sum()),
        "variance_top64": float(variance[:64].sum() / variance.sum()),
        "participation_ratio": float(variance.sum() ** 2 / (variance ** 2).sum()),
    }

    components = np.asarray([split.component for split in keep])
    records: list[dict[str, object]] = []
    for position, split in enumerate(keep):
        evaluation = split.evaluation
        evaluation_docs = documents[evaluation]
        label = substrate.affinity[evaluation]
        base = substrate.base[evaluation]
        reference = pair_concordance(label, base, evaluation_docs)
        if not np.isfinite(reference["ci_within"]):
            continue
        global_weight = globals_[int(folds[evaluation[0]])]
        global_score = pair_concordance(
            label, base + design[evaluation] @ global_weight, evaluation_docs
        )["ci_within"]

        # Basis from every *other* homology component.
        others = matrix[components != split.component]
        other_centre = others.mean(axis=0, keepdims=True)
        _, _, directions = np.linalg.svd(others - other_centre, full_matrices=False)
        own = heads[position]
        row: dict[str, object] = {
            "unit": split.unit, "component": split.component, "endpoint": split.endpoint,
            "base__ci_within": reference["ci_within"], "global__ci_within": global_score,
            "full__ci_within": pair_concordance(
                label, base + design[evaluation] @ own, evaluation_docs
            )["ci_within"],
        }
        for rank in RANKS:
            if rank > directions.shape[0]:
                continue
            projector = directions[:rank]
            projected = other_centre[0] + projector.T @ (projector @ (own - other_centre[0]))
            row[f"rank{rank}__ci_within"] = pair_concordance(
                label, base + design[evaluation] @ projected, evaluation_docs
            )["ci_within"]
        records.append(row)
        if position % 100 == 0:
            print(f"  rank {position}/{len(keep)}", flush=True)
    return pd.DataFrame.from_records(records), spectrum_stats


# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #


def summarise_curve(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for k in K_SWEEP:
        best = None
        per_ridge = {}
        for ridge in K_RIDGES:
            column = f"k{k}__r{ridge:g}"
            if column not in records.columns:
                continue
            frame = records[["component", "unit", column, "global__ci_within", "base__ci_within"]].copy()
            frame["gain"] = frame[column] - frame["global__ci_within"]
            cell = paired_bootstrap(frame, "gain")
            per_ridge[f"r{ridge:g}"] = cell
            if best is None or (np.isfinite(cell["mean"]) and cell["mean"] > best["mean"]):
                best = {**cell, "ridge": ridge}
        if best is not None:
            summary[f"k{k}"] = {"best": best, "by_ridge": per_ridge}
    return summary


def summarise_rank(records: pd.DataFrame) -> dict[str, object]:
    frame = records.copy()
    frame["full_gain"] = frame["full__ci_within"] - frame["global__ci_within"]
    full = paired_bootstrap(frame, "full_gain")
    summary = {"full": full, "by_rank": {}}
    for rank in RANKS:
        column = f"rank{rank}__ci_within"
        if column not in frame.columns:
            continue
        frame["gain"] = frame[column] - frame["global__ci_within"]
        cell = paired_bootstrap(frame, "gain")
        cell["retained_fraction"] = (
            float(cell["mean"] / full["mean"]) if full["mean"] not in (0.0, None) else float("nan")
        )
        summary["by_rank"][f"rank{rank}"] = cell
    return summary


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    basis, substrate, stats = rich_basis(substrate_dir, role)
    splits = [s for s in build_splits(substrate.rows) if s.regime == "separated"]
    print(f"{len(splits)} separated splits", flush=True)

    curve_records = label_curve(basis, substrate, splits)
    rank_records, spectrum = rank_probe(basis, substrate, splits)
    curve = summarise_curve(curve_records)
    rank = summarise_rank(rank_records)

    knee = None
    for k in K_SWEEP:
        cell = curve.get(f"k{k}", {}).get("best")
        if cell and cell.get("lower95", float("nan")) > 0.005:
            knee = k
            break
    best_rank = max(
        (cell["retained_fraction"] for cell in rank["by_rank"].values()
         if np.isfinite(cell.get("retained_fraction", np.nan)) and cell["components"] > 0
         and cell is not rank["by_rank"].get("rank266")),
        default=float("nan"),
    )
    low_rank_retained = rank["by_rank"].get("rank16", {}).get("retained_fraction", float("nan"))

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED, "capacity": CAPACITY, "k_sweep": list(K_SWEEP),
            "ridges": list(K_RIDGES), "draws": DRAWS, "ranks": list(RANKS),
            "reference": "target-agnostic global head of identical capacity",
            "selection": "best ridge per k (optimistic)",
            "basis_holdout": "rank basis estimated with the evaluated homology component removed",
        },
        "label_curve": curve,
        "rank": rank,
        "head_spectrum": spectrum,
        "headline": {
            "smallest_k_with_lower95_above_mde": knee,
            "rank16_retained_fraction": low_rank_retained,
            "best_sub_full_rank_retained_fraction": best_rank,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    curve_records.to_parquet(records_path, index=False)
    rank_records.to_parquet(records_path.with_name(records_path.stem + "_rank.parquet"), index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="M2 label-budget and rank probe")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print("label-budget curve (gain over target-agnostic head):")
    for key, value in payload["label_curve"].items():
        cell = value["best"]
        print(f"  {key:6s} ridge={cell['ridge']:>7g} {cell['mean']:+.4f} "
              f"[{cell['lower95']:+.4f},{cell['upper95']:+.4f}] n={cell['components']}")
    print("rank projection (retained fraction of the full-head gain):")
    print(f"  full   {payload['rank']['full']['mean']:+.4f}")
    for key, cell in payload["rank"]["by_rank"].items():
        print(f"  {key:8s} {cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}] "
              f"retained={cell['retained_fraction']:+.3f}")
    print("head spectrum:", json.dumps(payload["head_spectrum"], indent=2))


if __name__ == "__main__":
    main()
