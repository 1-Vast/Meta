"""M0: is the D0 null an absence, or a mis-parameterisation?

D0 measured one estimator -- a 26-dimensional ridge head on raw residuals,
fitted on documents disjoint from the evaluation documents -- and found no
within-document ranking gain.  That single measurement cannot distinguish three
very different worlds, and the whole mechanism programme turns on which one we
are in:

  (1) *Not expressible.*  No function of the ligand in this feature space orders
      compounds within a document at the target level.  Then there is no
      adaptation object and no mechanism can exist.

  (2) *Expressible but not transferable.*  A chemical head fitted on compounds
      from the same document ranks held-out compounds in that document, but a
      head fitted on other documents does not.  Then the object is real and
      document-local -- few-shot adaptation would have to identify a
      context-specific state, not a target-specific one.

  (3) *Transferable but obscured by nuisance.*  The object is there, but the
      estimator spends its budget on the per-context offsets that N0 showed
      carry 68 % of residual variance.  Fitting on within-document contrasts
      instead of absolute residuals would then recover it.  This is the NEA
      mechanism the A2S programme designed and never implemented.

The probe crosses three axes and grants every arm its best case:

  estimator   `cross_doc`     fit on training documents, evaluate on held-out
                              documents (this is D0's own head)
              `same_doc_cv`   fit and evaluate inside the *same* documents,
                              cross-fitted over compounds so no compound is
                              scored by a head that saw it
  target      `raw`           absolute residuals
              `centred`       residuals minus their document mean (G0-invariant)
              `studentised`   additionally divided by the document SD (affine)
  capacity    d in {10, 26, 74, 266} nested feature dimensions
  ridge       a grid, reported in full, with the *best* value taken per
              configuration -- an optimistic choice, which is the correct
              direction when the result is a null

Within-document concordance is invariant to any per-document additive constant,
so a head trained on contrasts is scored on exactly the same footing as one
trained on absolute residuals.

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

from research.psep.psep_d0 import (
    DEFAULT_SUBSTRATE,
    MIN_EVAL_ROWS,
    SEED,
    build_splits,
    load,
    pair_concordance,
    paired_bootstrap,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_m0_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_m0_records_2026-08-02.parquet"

PCA_MAX = 256
CAPACITIES = (10, 26, 74, 266)          # 10 descriptors + {0, 16, 64, 256} components
RIDGES = (1.0, 10.0, 100.0, 1000.0)
TARGETS = ("raw", "centred", "studentised")
INNER_FOLDS = 5
MIN_FIT_ROWS = 12


def rich_basis(substrate_dir: Path, role: str) -> tuple[np.ndarray, object, dict[str, object]]:
    """Label-free basis with up to PCA_MAX Morgan components, nested by capacity."""

    from scipy.sparse import load_npz

    substrate, stats = load(substrate_dir, role)
    bits = load_npz(substrate_dir / "morgan.npz").tocsr()
    descriptors = np.load(substrate_dir / "descriptors.npy")
    index = substrate.rows.structure_row.to_numpy()

    dense = np.asarray(bits[index].todense(), dtype=np.float64)
    centre = dense.mean(axis=0, keepdims=True)
    centred = dense - centre
    eigenvalues, eigenvectors = np.linalg.eigh(centred.T @ centred)
    order = np.argsort(eigenvalues)[::-1][:PCA_MAX]
    projection = eigenvectors[:, order]
    leading = np.abs(projection).argmax(axis=0)
    signs = np.sign(projection[leading, np.arange(projection.shape[1])])
    signs[signs == 0] = 1.0
    scores = centred @ (projection * signs)
    scores = (scores - scores.mean(axis=0)) / np.maximum(scores.std(axis=0), 1e-6)

    raw = descriptors[index]
    basis = np.hstack([(raw - raw.mean(axis=0)) / np.maximum(raw.std(axis=0), 1e-6), scores])
    stats["basis_dimension_max"] = int(basis.shape[1])
    stats["pca_max"] = PCA_MAX
    stats["explained_variance_ratio_max"] = float(
        eigenvalues[order].sum() / np.clip(eigenvalues, 0, None).sum()
    )
    return basis, substrate, stats


def ridge_path(design: np.ndarray, target: np.ndarray, ridges: tuple[float, ...]) -> np.ndarray:
    """Weights for every ridge from one eigendecomposition.  Returns (len(ridges), d)."""

    centre = design.mean(axis=0)
    centred = design - centre
    gram = centred.T @ centred
    values, vectors = np.linalg.eigh(gram)
    projected = vectors.T @ (centred.T @ (target - target.mean()))
    return np.stack([
        vectors @ (projected / (values + ridge)) for ridge in ridges
    ])


def contrast_targets(residual: np.ndarray, documents: np.ndarray) -> dict[str, np.ndarray]:
    frame = pd.DataFrame({"doc": documents, "r": residual})
    grouped = frame.groupby("doc", sort=False)["r"]
    means = grouped.transform("mean").to_numpy()
    deviations = grouped.transform(lambda values: values.std(ddof=1)).to_numpy()
    deviations = np.where(np.isfinite(deviations) & (deviations > 1e-6), deviations, 1.0)
    return {
        "raw": residual,
        "centred": residual - means,
        "studentised": (residual - means) / deviations,
    }


def document_folds(documents: np.ndarray, digest: int, folds: int) -> np.ndarray:
    """Fold assignment that keeps a document whole, so `same_doc_cv` never scores
    a compound with a head fitted on a document-mate of that compound... it does
    the opposite: folds split *compounds*, deliberately keeping documents shared,
    because the question is whether chemistry orders compounds inside a document."""

    rng = np.random.default_rng(digest ^ 0xA11CE)
    return rng.integers(0, folds, size=len(documents))


def evaluate(basis: np.ndarray, substrate, splits: list) -> pd.DataFrame:
    documents_all = substrate.rows.docs.astype(str).to_numpy()
    residual_all = substrate.residual
    affinity_all = substrate.affinity
    base_all = substrate.base

    records: list[dict[str, object]] = []
    separated = [split for split in splits if split.regime == "separated"]
    for number, split in enumerate(separated):
        train, evaluation = split.train, split.evaluation
        evaluation_docs = documents_all[evaluation]
        reference = pair_concordance(affinity_all[evaluation], base_all[evaluation], evaluation_docs)
        if not np.isfinite(reference["ci_within"]):
            continue

        row: dict[str, object] = {
            "unit": split.unit, "component": split.component, "endpoint": split.endpoint,
            "n_train": int(len(train)), "n_eval": int(len(evaluation)),
            "base__ci_within": reference["ci_within"], "pairs_within": reference["pairs_within"],
        }

        train_targets = contrast_targets(residual_all[train], documents_all[train])
        eval_targets = contrast_targets(residual_all[evaluation], evaluation_docs)
        digest = int(sha256(f"{SEED}:m0:{split.unit}".encode()).hexdigest()[:8], 16)
        folds = document_folds(evaluation_docs, digest, INNER_FOLDS)

        for capacity in CAPACITIES:
            design_train = basis[train][:, :capacity]
            design_eval = basis[evaluation][:, :capacity]
            for name in TARGETS:
                # --- cross-document head (D0's estimator, extended) -----------
                weights = ridge_path(design_train, train_targets[name], RIDGES)
                for ridge, weight in zip(RIDGES, weights):
                    score = pair_concordance(
                        affinity_all[evaluation],
                        base_all[evaluation] + design_eval @ weight,
                        evaluation_docs,
                    )
                    row[f"cross_doc__{name}__d{capacity}__r{ridge:g}"] = score["ci_within"]

                # --- same-document head, cross-fitted over compounds ---------
                predictions = np.zeros((len(RIDGES), len(evaluation)), dtype=np.float64)
                usable = True
                for fold in range(INNER_FOLDS):
                    held = folds == fold
                    inner = ~held
                    if inner.sum() < MIN_FIT_ROWS or held.sum() == 0:
                        usable = False
                        break
                    inner_weights = ridge_path(
                        design_eval[inner], eval_targets[name][inner], RIDGES
                    )
                    predictions[:, held] = (design_eval[held] @ inner_weights.T).T
                if usable:
                    for position, ridge in enumerate(RIDGES):
                        score = pair_concordance(
                            affinity_all[evaluation],
                            base_all[evaluation] + predictions[position],
                            evaluation_docs,
                        )
                        row[f"same_doc_cv__{name}__d{capacity}__r{ridge:g}"] = score["ci_within"]
        records.append(row)
        if number % 100 == 0:
            print(f"  {number}/{len(separated)}", flush=True)
    return pd.DataFrame.from_records(records)


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    columns = [c for c in records.columns if "__d" in c]
    for column in columns:
        estimator, target, capacity, ridge = column.split("__")
        frame = records[["component", "unit", column, "base__ci_within"]].copy()
        frame["gain"] = frame[column] - frame["base__ci_within"]
        cell = paired_bootstrap(frame, "gain")
        summary.setdefault(estimator, {}).setdefault(target, {}).setdefault(capacity, {})[ridge] = cell
    return summary


def best_of(summary: dict[str, object]) -> dict[str, object]:
    """Optimistic selection: the best ridge and capacity per (estimator, target)."""

    best: dict[str, object] = {}
    for estimator, targets in summary.items():
        for target, capacities in targets.items():
            champion = None
            for capacity, ridges in capacities.items():
                for ridge, cell in ridges.items():
                    if not np.isfinite(cell.get("mean", np.nan)):
                        continue
                    if champion is None or cell["mean"] > champion["cell"]["mean"]:
                        champion = {"capacity": capacity, "ridge": ridge, "cell": cell}
            if champion is not None:
                best[f"{estimator}__{target}"] = champion
    return best


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    basis, substrate, stats = rich_basis(substrate_dir, role)
    splits = build_splits(substrate.rows)
    print(f"built {len(splits)} splits", flush=True)
    records = evaluate(basis, substrate, splits)
    summary = summarise(records)
    best = best_of(summary)

    cross = max(
        (value["cell"]["mean"] for key, value in best.items() if key.startswith("cross_doc")),
        default=float("nan"),
    )
    same = max(
        (value["cell"]["mean"] for key, value in best.items() if key.startswith("same_doc_cv")),
        default=float("nan"),
    )
    same_lower = max(
        (value["cell"]["lower95"] for key, value in best.items() if key.startswith("same_doc_cv")),
        default=float("nan"),
    )
    if np.isfinite(same_lower) and same_lower > 0.005:
        world = "EXPRESSIBLE_LOCALLY" if cross < 0.005 else "EXPRESSIBLE_AND_TRANSFERABLE"
    else:
        world = "NOT_EXPRESSIBLE_NO_ADAPTATION_OBJECT"

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED,
            "capacities": list(CAPACITIES),
            "ridges": list(RIDGES),
            "targets": list(TARGETS),
            "inner_folds": INNER_FOLDS,
            "selection": "best ridge and capacity per configuration (optimistic)",
            "metric": "within-document pair concordance minus base, component bootstrap",
        },
        "best": best,
        "summary": summary,
        "world": world,
        "headline": {
            "best_cross_document_gain": cross,
            "best_same_document_gain": same,
            "best_same_document_lower95": same_lower,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="M0 capacity / invariance / locality probe")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({"world": payload["world"], "headline": payload["headline"]}, indent=2))
    for key, value in sorted(payload["best"].items()):
        cell = value["cell"]
        print(f"{key:28s} {value['capacity']:>5s} {value['ridge']:>6s} "
              f"{cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}] "
              f"n={cell['components']}")


if __name__ == "__main__":
    main()
