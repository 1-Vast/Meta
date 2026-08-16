"""M3: is the source-head subspace the *right* subspace?

M2 measured the obstruction precisely.  The per-target adaptation object is
high-dimensional -- 553 heads in 266 dimensions with a participation ratio of
115 -- and a rank-16 principal subspace of the source heads retains only 29.6 %
of the object's ranking gain.  That is almost exactly what k=5 free fitting
already achieves (+0.0055 of an available +0.0190), so the low-rank route buys
nothing *as currently constructed*.

But that construction has a defect worth taking seriously before concluding the
object is irreducibly high-dimensional.  Principal components were taken in raw
coefficient space, i.e. under the Euclidean inner product on `w`.  The quantity
that actually matters is the *function* the head computes, and two coefficient
vectors that differ along a direction the chemistry never varies in produce
identical predictions.  The correct geometry is therefore the one induced by the
design second moment

    Sigma = E[ x x^T ],        <w1, w2>_Sigma = w1^T Sigma w2,

under which distance between heads equals expected squared difference between
their predictions.  Compressing in the wrong metric spends rank on directions
that carry coefficient variance but no predictive variance.

This probe compares three subspaces at matched rank, all estimated with the
evaluated homology component held out:

  `raw`         principal components of the heads (what M2 measured)
  `whitened`    principal components of Sigma^{1/2} w -- the function-space
                geometry
  `precision`   whitened, and additionally weighted by sqrt(n_train) so that
                heads estimated from more labels count for more, which is the
                empirical-Bayes weighting for heads of unequal reliability

If `whitened` retains materially more than `raw` at low rank, the object's
apparent high dimensionality is partly an artefact of the metric, and a
meta-learned representation has somewhere to go.  If all three coincide, the
dimensionality is intrinsic and the few-shot framing is in serious trouble.

Reads the `discover` role only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits, pair_concordance, paired_bootstrap
from research.psep.psep_m0 import rich_basis
from research.psep.psep_m2 import RANKS, ridge_fit
from research.psep.psep_m4 import CAPACITY, RIDGE, fit_global_heads

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_m3_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_m3_records_2026-08-02.parquet"

GEOMETRIES = ("raw", "whitened", "precision")


def symmetric_power(matrix: np.ndarray, power: float, floor: float = 1e-8) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    values = np.clip(values, floor, None)
    return (vectors * (values ** power)) @ vectors.T


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    basis, substrate, stats = rich_basis(substrate_dir, role)
    splits = [s for s in build_splits(substrate.rows) if s.regime == "separated"]
    design = basis[:, :CAPACITY]
    documents = substrate.rows.docs.astype(str).to_numpy()
    globals_ = fit_global_heads(basis, substrate)
    folds = substrate.rows.oof_fold.to_numpy()

    heads, keep = [], []
    for split in splits:
        if len(split.train) < 20:
            continue
        heads.append(ridge_fit(design[split.train], substrate.residual[split.train], (RIDGE,))[0])
        keep.append(split)
    matrix = np.stack(heads)
    weights = np.sqrt(np.asarray([len(split.train) for split in keep], dtype=np.float64))
    components = np.asarray([split.component for split in keep])

    covariance = design.T @ design / len(design)
    root = symmetric_power(covariance, 0.5)
    inverse_root = symmetric_power(covariance, -0.5)
    print(f"design covariance condition: {np.linalg.cond(covariance):.3e}", flush=True)

    transformed = {
        "raw": matrix,
        "whitened": matrix @ root,
        "precision": (matrix @ root) * weights[:, None],
    }

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
        row: dict[str, object] = {
            "unit": split.unit, "component": split.component, "endpoint": split.endpoint,
            "base__ci_within": reference["ci_within"],
            "global__ci_within": pair_concordance(
                label, base + design[evaluation] @ global_weight, evaluation_docs
            )["ci_within"],
            "full__ci_within": pair_concordance(
                label, base + design[evaluation] @ heads[position], evaluation_docs
            )["ci_within"],
        }
        held = components != split.component
        for geometry in GEOMETRIES:
            points = transformed[geometry]
            others = points[held]
            centre = others.mean(axis=0, keepdims=True)
            _, _, directions = np.linalg.svd(others - centre, full_matrices=False)
            own = points[position]
            if geometry == "precision":
                own = own / weights[position]
                centre = centre / weights[held].mean()
            for rank in RANKS:
                if rank > directions.shape[0]:
                    continue
                projector = directions[:rank]
                projected = centre[0] + projector.T @ (projector @ (own - centre[0]))
                vector = projected if geometry == "raw" else projected @ inverse_root
                row[f"{geometry}__rank{rank}"] = pair_concordance(
                    label, base + design[evaluation] @ vector, evaluation_docs
                )["ci_within"]
        records.append(row)
        if position % 100 == 0:
            print(f"  {position}/{len(keep)}", flush=True)

    frame = pd.DataFrame.from_records(records)
    frame["full_gain"] = frame["full__ci_within"] - frame["global__ci_within"]
    full = paired_bootstrap(frame, "full_gain")

    summary: dict[str, object] = {"full": full, "geometries": {}}
    for geometry in GEOMETRIES:
        cells: dict[str, object] = {}
        for rank in RANKS:
            column = f"{geometry}__rank{rank}"
            if column not in frame.columns:
                continue
            frame["gain"] = frame[column] - frame["global__ci_within"]
            cell = paired_bootstrap(frame, "gain")
            cell["retained_fraction"] = float(cell["mean"] / full["mean"]) if full["mean"] else float("nan")
            cells[f"rank{rank}"] = cell
        summary["geometries"][geometry] = cells

    def retained(geometry: str, rank: int) -> float:
        return summary["geometries"][geometry].get(f"rank{rank}", {}).get("retained_fraction", float("nan"))

    improvement = {
        f"rank{rank}": {
            "raw": retained("raw", rank),
            "whitened": retained("whitened", rank),
            "precision": retained("precision", rank),
        }
        for rank in RANKS
    }
    metric_matters = bool(
        np.isfinite(retained("whitened", 16)) and np.isfinite(retained("raw", 16))
        and retained("whitened", 16) - retained("raw", 16) > 0.10
    )

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED, "capacity": CAPACITY, "ranks": list(RANKS),
            "geometries": list(GEOMETRIES),
            "basis_holdout": "evaluated homology component removed before the SVD",
            "reference": "target-agnostic global head of identical capacity",
        },
        "summary": summary,
        "retained_by_rank": improvement,
        "metric_materially_matters": metric_matters,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="M3 head-space geometry probe")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(f"full gain {payload['summary']['full']['mean']:+.4f}")
    print(f"{'rank':>8s} {'raw':>18s} {'whitened':>18s} {'precision':>18s}")
    for rank, cells in payload["retained_by_rank"].items():
        print(f"{rank:>8s} {cells['raw']:>18.3f} {cells['whitened']:>18.3f} {cells['precision']:>18.3f}")
    print("metric materially matters:", payload["metric_materially_matters"])


if __name__ == "__main__":
    main()
