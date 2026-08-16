"""D0 at power: does a target-specific chemical adaptation object survive
simultaneous scaffold + document + assay separation?

This is the A2S programme's registered reopening test, run on the PSEP substrate
(828 homology components, 508 of them in the `discover` role) instead of the 92
components on which it originally failed.  The estimand, the split rule, the
metric and the decision thresholds are ported unchanged from
`research/a2s_nea_preconditions.py` so the two results are directly comparable.

The measurement, restated:

  base            a target-agnostic ridge on chemistry alone, cross-fitted by
                  homology component.
  own head        a 26-dimensional ridge fitted on the unit's *own* training
                  rows -- the most generous possible per-target chemical head,
                  fitted from labels, with no few-shot restriction.
  document oracle a chemistry-free predictor that knows only a compound's
                  document and that document's mean residual.

  separated       train and evaluation rows share no document, no assay and no
                  Murcko scaffold.
  scaffold_only   the conventional split, for attribution.

  headline        within-document pair concordance, where any additive
                  per-context offset cancels exactly.

`split_is_clean` is a structural self-check, not a statistical one: on a
document-disjoint split the document oracle can only predict a constant, so its
measured gain must be ~0.  If it is not, the harness is wrong and the run is void.

Only the `discover` role is read.  `validate` and `confirm` are untouched.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.sparse import load_npz

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SUBSTRATE = ROOT / "dataset" / "processed" / "psep.v1"
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_d0_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_d0_records_2026-08-02.parquet"

SEED = 20260802
BOOTSTRAP_DRAWS = 2000
PCA_COMPONENTS = 16
HEAD_RIDGE = 1.0
MIN_PAIRS = 8
MIN_EVAL_ROWS = 8
MIN_TRAIN_ROWS = 20
MAX_EVAL_ROWS = 600        # bounds pair cost on the deepest units; deterministic

D0_MIN_HEADROOM = 0.005
D0_MAX_ORACLE_LEAK = 0.005


# --------------------------------------------------------------------------- #
# Substrate
# --------------------------------------------------------------------------- #


@dataclass
class Substrate:
    rows: pd.DataFrame
    basis: np.ndarray
    residual: np.ndarray
    affinity: np.ndarray
    base: np.ndarray


def load(substrate_dir: Path, role: str) -> tuple[Substrate, dict[str, object]]:
    rows = pd.read_parquet(substrate_dir / "rows.parquet")
    manifest = json.loads((substrate_dir / "manifest.json").read_text(encoding="utf-8"))
    available = set(rows.role.unique())
    if role not in available:
        raise RuntimeError(f"role {role!r} not present in substrate")
    rows = rows.loc[rows.role == role].reset_index(drop=True)

    bits = load_npz(substrate_dir / "morgan.npz").tocsr()
    descriptors = np.load(substrate_dir / "descriptors.npy")
    index = rows.structure_row.to_numpy()

    # Label-free basis: 10 standardised descriptors + top Morgan principal
    # components, from an exact sign-fixed eigendecomposition so the basis is
    # reproducible across processes (the defect that invalidated earlier gates).
    dense = np.asarray(bits[index].todense(), dtype=np.float64)
    centre = dense.mean(axis=0, keepdims=True)
    centred = dense - centre
    eigenvalues, eigenvectors = np.linalg.eigh(centred.T @ centred)
    order = np.argsort(eigenvalues)[::-1][:PCA_COMPONENTS]
    projection = eigenvectors[:, order]
    leading = np.abs(projection).argmax(axis=0)
    signs = np.sign(projection[leading, np.arange(projection.shape[1])])
    signs[signs == 0] = 1.0
    projection = projection * signs
    scores = centred @ projection
    scores = (scores - scores.mean(axis=0)) / np.maximum(scores.std(axis=0), 1e-6)

    raw = descriptors[index]
    scale = np.maximum(raw.std(axis=0), 1e-6)
    basis = np.hstack([(raw - raw.mean(axis=0)) / scale, scores])

    affinity = rows.affinity.to_numpy(dtype=np.float64)
    base = rows.base.to_numpy(dtype=np.float64)
    stats = {
        "role": role,
        "rows": int(len(rows)),
        "units": int(rows.unit.nunique()),
        "components": int(rows.component.nunique()),
        "basis_dimension": int(basis.shape[1]),
        "basis_decomposition": "exact symmetric eigendecomposition, sign-fixed",
        "explained_variance_ratio": float(
            eigenvalues[order].sum() / np.clip(eigenvalues, 0, None).sum()
        ),
        "substrate_manifest": manifest["counts"],
    }
    return Substrate(rows, basis, affinity - base, affinity, base), stats


# --------------------------------------------------------------------------- #
# Estimators (ported unchanged)
# --------------------------------------------------------------------------- #


def fit_head(design: np.ndarray, residual: np.ndarray, ridge: float = HEAD_RIDGE) -> np.ndarray:
    centre = design.mean(axis=0)
    centred = design - centre
    gram = centred.T @ centred + ridge * np.eye(design.shape[1])
    return np.linalg.solve(gram, centred.T @ (residual - residual.mean()))


def pair_concordance(
    label: np.ndarray, prediction: np.ndarray, group: np.ndarray | None = None
) -> dict[str, float]:
    label = np.asarray(label, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if len(label) < 2:
        return {"ci": float("nan"), "ci_within": float("nan"), "pairs": 0, "pairs_within": 0}
    left, right = np.triu_indices(len(label), k=1)
    truth = np.sign(label[left] - label[right])
    guess = np.sign(prediction[left] - prediction[right])
    active = truth != 0
    correct = (guess == truth).astype(np.float64) + 0.5 * (guess == 0).astype(np.float64)

    def score(mask: np.ndarray) -> float:
        selected = mask & active
        if selected.sum() < MIN_PAIRS:
            return float("nan")
        return float(correct[selected].mean())

    if group is None:
        same = np.zeros(len(left), dtype=bool)
    else:
        codes = pd.factorize(np.asarray(group, dtype=object))[0]
        same = codes[left] == codes[right]
    return {
        "ci": score(np.ones(len(left), dtype=bool)),
        "ci_within": score(same),
        "pairs": int(active.sum()),
        "pairs_within": int((same & active).sum()),
    }


def document_offset_prediction(
    train_docs: np.ndarray, train_residual: np.ndarray, eval_docs: np.ndarray
) -> np.ndarray:
    means = pd.DataFrame({"doc": train_docs, "r": train_residual}).groupby("doc", sort=False)["r"].mean()
    return np.asarray([means.get(doc, 0.0) for doc in eval_docs], dtype=np.float64)


# --------------------------------------------------------------------------- #
# Splits
# --------------------------------------------------------------------------- #


@dataclass
class Split:
    unit: str
    component: str
    endpoint: str
    train: np.ndarray
    evaluation: np.ndarray
    regime: str


def build_splits(rows: pd.DataFrame) -> list[Split]:
    splits: list[Split] = []
    for unit, group in rows.groupby("unit", sort=True):
        index = group.index.to_numpy()
        if len(index) < MIN_TRAIN_ROWS + MIN_EVAL_ROWS:
            continue
        component = str(group.component.iloc[0])
        endpoint = str(group.endpoint.iloc[0])
        documents = group.docs.astype(str).to_numpy()
        assays = group.assays.astype(str).to_numpy()
        scaffolds = group.scaffold.astype(str).to_numpy()

        counts = pd.Series(documents).value_counts()
        eligible = [doc for doc in counts.index if counts[doc] >= MIN_EVAL_ROWS]
        if len(eligible) < 2:
            continue
        digest = int(sha256(f"{SEED}:{unit}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(digest)

        held: set[str] = set()
        for position in rng.permutation(len(eligible)):
            held.add(eligible[position])
            if np.isin(documents, list(held)).sum() >= max(MIN_EVAL_ROWS, int(0.3 * len(index))):
                break
        evaluation_mask = np.isin(documents, list(held))
        train_mask = ~evaluation_mask
        if train_mask.sum() < MIN_TRAIN_ROWS:
            continue

        clean = evaluation_mask.copy()
        clean &= ~np.isin(scaffolds, list(set(scaffolds[train_mask].tolist())))
        clean &= ~np.isin(assays, list(set(assays[train_mask].tolist())))
        if clean.sum() >= MIN_EVAL_ROWS:
            splits.append(Split(str(unit), component, endpoint,
                                index[train_mask], _cap(index[clean], digest), "separated"))

        unique = sorted(set(scaffolds))
        chosen: set[str] = set()
        for position in rng.permutation(len(unique)):
            chosen.add(unique[position])
            if np.isin(scaffolds, list(chosen)).sum() >= max(MIN_EVAL_ROWS, int(0.3 * len(index))):
                break
        scaffold_eval = np.isin(scaffolds, list(chosen))
        if scaffold_eval.sum() >= MIN_EVAL_ROWS and (~scaffold_eval).sum() >= MIN_TRAIN_ROWS:
            splits.append(Split(str(unit), component, endpoint,
                                index[~scaffold_eval], _cap(index[scaffold_eval], digest),
                                "scaffold_only"))
    return splits


def _cap(rows: np.ndarray, digest: int) -> np.ndarray:
    """Deterministically bound evaluation size so pair cost stays finite."""

    if len(rows) <= MAX_EVAL_ROWS:
        return rows
    rng = np.random.default_rng(digest ^ 0x5EED)
    return np.sort(rng.choice(rows, size=MAX_EVAL_ROWS, replace=False))


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #


def evaluate(substrate: Substrate, splits: list[Split]) -> pd.DataFrame:
    documents = substrate.rows.docs.astype(str).to_numpy()
    assays = substrate.rows.assays.astype(str).to_numpy()
    records: list[dict[str, object]] = []
    for number, split in enumerate(splits):
        train, evaluation = split.train, split.evaluation
        evaluation_docs = documents[evaluation]
        reference = pair_concordance(
            substrate.affinity[evaluation], substrate.base[evaluation], evaluation_docs
        )
        if not np.isfinite(reference["ci"]):
            continue
        weight = fit_head(substrate.basis[train], substrate.residual[train])
        own = pair_concordance(
            substrate.affinity[evaluation],
            substrate.base[evaluation] + substrate.basis[evaluation] @ weight,
            evaluation_docs,
        )
        offset = document_offset_prediction(
            documents[train], substrate.residual[train], evaluation_docs
        )
        oracle = pair_concordance(
            substrate.affinity[evaluation], substrate.base[evaluation] + offset, evaluation_docs
        )
        train_docs = set(documents[train].tolist())
        train_assays = set(assays[train].tolist())
        records.append({
            "unit": split.unit,
            "component": split.component,
            "endpoint": split.endpoint,
            "regime": split.regime,
            "n_train": int(len(train)),
            "n_eval": int(len(evaluation)),
            "eval_documents": int(len(set(evaluation_docs.tolist()))),
            "document_overlap": float(np.mean([doc in train_docs for doc in evaluation_docs])),
            "assay_overlap": float(np.mean([a in train_assays for a in assays[evaluation]])),
            "base__ci": reference["ci"],
            "base__ci_within": reference["ci_within"],
            "own__ci": own["ci"],
            "own__ci_within": own["ci_within"],
            "docoffset__ci": oracle["ci"],
            "pairs": reference["pairs"],
            "pairs_within": reference["pairs_within"],
        })
        if number % 250 == 0:
            print(f"  evaluated {number}/{len(splits)}", flush=True)
    return pd.DataFrame.from_records(records)


CONTRASTS = {
    "own_minus_base": ("own__ci", "base__ci"),
    "own_minus_base_within_document": ("own__ci_within", "base__ci_within"),
    "docoffset_minus_base": ("docoffset__ci", "base__ci"),
}


def paired_bootstrap(frame: pd.DataFrame, column: str, draws: int = BOOTSTRAP_DRAWS) -> dict[str, float]:
    usable = frame[["component", "unit", column]].dropna()
    if usable.empty:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    per_unit = usable.groupby(["component", "unit"], sort=True)[column].mean().reset_index()
    per_component = per_unit.groupby("component", sort=True)[column].mean().to_numpy(dtype=np.float64)
    per_component = per_component[np.isfinite(per_component)]
    if len(per_component) == 0:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    rng = np.random.default_rng(SEED)
    sample = rng.integers(0, len(per_component), size=(draws, len(per_component)))
    means = per_component[sample].mean(axis=1)
    return {
        "components": int(len(per_component)),
        "units": int(len(per_unit)),
        "mean": float(per_component.mean()),
        "sd_across_components": float(per_component.std(ddof=1)) if len(per_component) > 1 else float("nan"),
        "lower95": float(np.quantile(means, 0.025)),
        "upper95": float(np.quantile(means, 0.975)),
    }


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for regime in sorted(records.regime.unique()):
        frame = records.loc[records.regime == regime].copy()
        cells: dict[str, object] = {}
        for name, (left, right) in CONTRASTS.items():
            frame[name] = frame[left] - frame[right]
            cells[name] = paired_bootstrap(frame, name)
        cells["by_endpoint"] = {
            endpoint: paired_bootstrap(part, "own_minus_base_within_document")
            for endpoint, part in frame.groupby("endpoint")
        }
        cells["descriptives"] = {
            "units": int(frame.unit.nunique()),
            "components": int(frame.component.nunique()),
            "mean_document_overlap": float(frame.document_overlap.mean()),
            "mean_assay_overlap": float(frame.assay_overlap.mean()),
            "mean_eval_rows": float(frame.n_eval.mean()),
            "mean_base_ci": float(frame.base__ci.mean()),
            "mean_base_ci_within": float(frame.base__ci_within.mean(skipna=True)),
            "total_within_document_pairs": int(frame.pairs_within.sum()),
        }
        summary[regime] = cells
    return summary


def decide(summary: dict[str, object]) -> dict[str, object]:
    separated = summary.get("separated", {})

    def cell(name: str) -> dict[str, float]:
        value = separated.get(name)
        return value if isinstance(value, dict) else {"mean": float("nan"), "lower95": float("nan")}

    headroom = cell("own_minus_base_within_document")
    leak = cell("docoffset_minus_base")
    clean = bool(
        abs(leak.get("mean", 1.0)) < D0_MAX_ORACLE_LEAK
        and separated.get("descriptives", {}).get("mean_document_overlap", 1.0) < 1e-9
    )
    passed = bool(headroom.get("lower95", float("nan")) > D0_MIN_HEADROOM)
    if not clean:
        verdict = "D0_SPLIT_NOT_CLEAN_HARNESS_INVALID"
    elif passed:
        verdict = "ADAPTATION_OBJECT_SURVIVES_SEPARATION_AT_POWER"
    else:
        verdict = "NO_ADAPTATION_OBJECT_AT_POWER"
    return {
        "gate_D0": {
            "within_document_headroom": headroom,
            "all_pair_headroom": cell("own_minus_base"),
            "document_oracle_leak": leak,
            "split_is_clean": clean,
            "pass": passed,
            "threshold": D0_MIN_HEADROOM,
        },
        "verdict": verdict,
    }


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    substrate, stats = load(substrate_dir, role)
    if role != "discover":
        raise RuntimeError("D0 is registered to run on the discover role only")
    splits = build_splits(substrate.rows)
    print(f"built {len(splits)} splits", flush=True)
    records = evaluate(substrate, splits)
    summary = summarise(records)
    decision = decide(summary)

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {
            "evaluated_role": role,
            "validate_read": False,
            "confirm_read": False,
            "chembl_probe_locked_recipients_excluded_by_homology": True,
            "trains_no_gradient_model": True,
        },
        "protocol": {
            "seed": SEED,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "max_eval_rows": MAX_EVAL_ROWS,
            "estimand": "within-document pair concordance, target-macro, component bootstrap",
            "ported_from": "research/a2s_nea_preconditions.py (unchanged thresholds)",
        },
        "summary": summary,
        **decision,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="D0 at power on the PSEP substrate")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({
        "verdict": payload["verdict"],
        "gate_D0": {
            key: payload["gate_D0"][key]
            for key in ("within_document_headroom", "all_pair_headroom",
                        "document_oracle_leak", "split_is_clean", "pass")
        },
    }, indent=2))


if __name__ == "__main__":
    main()
