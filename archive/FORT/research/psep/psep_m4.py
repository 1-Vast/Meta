"""M4: is the high-capacity cross-document gain a *target-specific* object?

M0 found that the cross-document head's within-document gain rises monotonically
with capacity -- +0.002 at d=10, +0.005 at d=26 (the dimension D0 used, where it
is indistinguishable from zero), +0.008 at d=74, +0.015 at d=266.  Monotone in
capacity across every ridge and every target variant, so it is a systematic
effect rather than selection noise.

That is necessary but nowhere near sufficient.  A per-target head fitted at
d=266 can gain for a reason that has nothing to do with the target: the frozen
base is a ridge on the same features, and if it is simply *underfit*, then any
richer head -- target-specific or not -- would close the same gap.  The gain
would then be a base-capacity artefact wearing a per-target costume.

This gate separates them.  Every arm uses the identical d=266 basis, the
identical ridge, and the identical evaluation rows; only the *rows the head is
fitted on* change:

  global    all discover rows from other cross-fit folds, pooled.  Target-
            agnostic by construction: it cannot know which target it is scoring.
  own       the unit's own training rows (documents disjoint from evaluation).
  wrong     the training rows of a *different* unit in a different homology
            component -- a derangement.  Same row count, same fitting procedure,
            same capacity; only the identity of the target is wrong.
  stacked   `global`, plus an `own` head fitted on the residual that `global`
            leaves.  This is the honest form of the question "what does knowing
            the target add on top of the best target-agnostic chemistry?"

The contrast that decides the programme is **own - global**.  If it is zero,
there is no target-specific adaptation object and the whole few-shot framing is
misconceived on this data; the correct response is a better base model, not a
meta-learner.  `own - wrong` is the second, stricter form of the same question.

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
from research.psep.psep_m0 import CAPACITIES, ridge_path, rich_basis

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_m4_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_m4_records_2026-08-02.parquet"

CAPACITY = 266             # pre-specified from the M0 sweep
RIDGE = 100.0              # pre-specified from the M0 sweep
GLOBAL_RIDGE = 1000.0      # the global head sees ~10^5 rows, so it is scaled up
MIN_TRAIN = 20


def fit_global_heads(basis: np.ndarray, substrate) -> dict[int, np.ndarray]:
    """One target-agnostic head per cross-fit fold, fitted on the other folds.

    Fitted on *residual* exactly as the per-target heads are, so the two arms
    differ only in which rows they see.
    """

    folds = substrate.rows.oof_fold.to_numpy()
    design = basis[:, :CAPACITY]
    heads: dict[int, np.ndarray] = {}
    for fold in sorted(set(folds.tolist())):
        train = folds != fold
        weights = ridge_path(design[train], substrate.residual[train], (GLOBAL_RIDGE,))
        heads[int(fold)] = weights[0]
        print(f"  global head fold {fold}: {int(train.sum())} rows", flush=True)
    return heads


def derange(splits: list, seed: int) -> dict[str, int]:
    """Map each split to a donor split from a different homology component."""

    rng = np.random.default_rng(seed)
    components = np.asarray([split.component for split in splits])
    donors: dict[str, int] = {}
    for position, split in enumerate(splits):
        allowed = np.flatnonzero(components != split.component)
        if len(allowed) == 0:
            continue
        donors[split.unit] = int(rng.choice(allowed))
    return donors


def evaluate(basis: np.ndarray, substrate, splits: list) -> pd.DataFrame:
    documents = substrate.rows.docs.astype(str).to_numpy()
    design = basis[:, :CAPACITY]
    separated = [split for split in splits if split.regime == "separated"]
    globals_ = fit_global_heads(basis, substrate)
    donors = derange(separated, SEED)
    by_unit = {split.unit: split for split in separated}

    records: list[dict[str, object]] = []
    for number, split in enumerate(separated):
        train, evaluation = split.train, split.evaluation
        evaluation_docs = documents[evaluation]
        fold = int(substrate.rows.oof_fold.to_numpy()[evaluation[0]])
        base = substrate.base[evaluation]
        label = substrate.affinity[evaluation]

        reference = pair_concordance(label, base, evaluation_docs)
        if not np.isfinite(reference["ci_within"]):
            continue

        global_weight = globals_[fold]
        global_prediction = design[evaluation] @ global_weight
        own_weight = ridge_path(design[train], substrate.residual[train], (RIDGE,))[0]

        donor_index = donors.get(split.unit)
        donor = separated[donor_index] if donor_index is not None else None
        wrong_prediction = None
        if donor is not None and len(donor.train) >= MIN_TRAIN:
            wrong_weight = ridge_path(design[donor.train], substrate.residual[donor.train], (RIDGE,))[0]
            wrong_prediction = design[evaluation] @ wrong_weight

        # `stacked`: what the target adds on top of the best target-agnostic head.
        global_train = design[train] @ global_weight
        increment_weight = ridge_path(
            design[train], substrate.residual[train] - global_train, (RIDGE,)
        )[0]
        stacked_prediction = global_prediction + design[evaluation] @ increment_weight

        arms = {
            "base": np.zeros(len(evaluation)),
            "global": global_prediction,
            "own": design[evaluation] @ own_weight,
            "stacked": stacked_prediction,
        }
        if wrong_prediction is not None:
            arms["wrong"] = wrong_prediction

        row: dict[str, object] = {
            "unit": split.unit, "component": split.component, "endpoint": split.endpoint,
            "n_train": int(len(train)), "n_eval": int(len(evaluation)),
            "pairs_within": reference["pairs_within"],
            "donor_unit": donor.unit if donor is not None else None,
            "donor_n_train": int(len(donor.train)) if donor is not None else None,
        }
        for name, delta in arms.items():
            score = pair_concordance(label, base + delta, evaluation_docs)
            row[f"{name}__ci_within"] = score["ci_within"]
            row[f"{name}__ci"] = score["ci"]
        records.append(row)
        if number % 100 == 0:
            print(f"  {number}/{len(separated)}", flush=True)
    return pd.DataFrame.from_records(records)


CONTRASTS = {
    "global_minus_base": ("global__ci_within", "base__ci_within"),
    "own_minus_base": ("own__ci_within", "base__ci_within"),
    "wrong_minus_base": ("wrong__ci_within", "base__ci_within"),
    "stacked_minus_base": ("stacked__ci_within", "base__ci_within"),
    "own_minus_global": ("own__ci_within", "global__ci_within"),
    "own_minus_wrong": ("own__ci_within", "wrong__ci_within"),
    "stacked_minus_global": ("stacked__ci_within", "global__ci_within"),
}


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for name, (left, right) in CONTRASTS.items():
        if left not in records.columns or right not in records.columns:
            continue
        frame = records.copy()
        frame[name] = frame[left] - frame[right]
        summary[name] = paired_bootstrap(frame, name)
        summary[name]["by_endpoint"] = {
            endpoint: paired_bootstrap(part, name)
            for endpoint, part in frame.groupby("endpoint")
        }
    return summary


def decide(summary: dict[str, object]) -> dict[str, object]:
    def cell(name: str) -> dict[str, float]:
        value = summary.get(name, {})
        return value if isinstance(value, dict) else {}

    specific = cell("own_minus_global")
    versus_wrong = cell("own_minus_wrong")
    agnostic = cell("global_minus_base")

    target_specific = bool(specific.get("lower95", float("nan")) > 0.005)
    beats_wrong = bool(versus_wrong.get("lower95", float("nan")) > 0.005)
    if target_specific and beats_wrong:
        verdict = "TARGET_SPECIFIC_ADAPTATION_OBJECT_CONFIRMED"
    elif agnostic.get("lower95", 0.0) > 0.005:
        verdict = "GAIN_IS_BASE_CAPACITY_NOT_TARGET_ADAPTATION"
    else:
        verdict = "NO_ADAPTATION_OBJECT"
    return {
        "gate_M4": {
            "target_specific_increment": specific,
            "own_versus_wrong_target": versus_wrong,
            "target_agnostic_capacity_gain": agnostic,
            "pass_target_specific": target_specific,
            "pass_beats_wrong_target": beats_wrong,
        },
        "verdict": verdict,
    }


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    basis, substrate, stats = rich_basis(substrate_dir, role)
    splits = build_splits(substrate.rows)
    records = evaluate(basis, substrate, splits)
    summary = summarise(records)
    decision = decide(summary)

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED,
            "capacity": CAPACITY,
            "ridge": RIDGE,
            "global_ridge": GLOBAL_RIDGE,
            "capacity_and_ridge_source": "pre-specified from the M0 sweep on the same role",
            "metric": "within-document pair concordance, component bootstrap",
            "arms": ["base", "global", "own", "wrong", "stacked"],
        },
        "summary": summary,
        **decision,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="M4 counterfactual identifiability gate")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({"verdict": payload["verdict"]}, indent=2))
    for name, cell in payload["summary"].items():
        print(f"{name:26s} {cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}] "
              f"n={cell['components']}")


if __name__ == "__main__":
    main()
