"""NEA preconditions: gates D0, N0 and N1.

Nuisance-Equivariant Adaptation rests on three empirical claims that can be
measured without building the mechanism.  This runner measures them, on source
``fit`` targets only, and trains nothing beyond closed-form ridge heads.

``D0``  Does a target-specific *chemical* ranking object survive when support and
        query are separated simultaneously by scaffold, **document** and
        **assay**?  Gate T0 showed a chemistry-free document oracle outscores the
        whole chemical head under scaffold-only separation, so the size of the
        prize is unknown.  D0 measures it, and reports the same protocol under
        scaffold-only separation so the drop is attributable.

        *Built-in validity check:* under document-disjoint separation the
        document-mean oracle is structurally powerless -- every evaluation
        document is unseen, so it predicts a constant.  Its measured gain must be
        ~0.  If it is not, the split is not what it claims to be.

``N0``  Which nuisance group actually acts on measurements?  ``G0`` is a
        per-context offset (``y -> y + b_c``); the full affine group ``G`` adds a
        per-context scale (``y -> a_c y + b_c``).  Invariance to ``G`` costs one
        degree of freedom per context more than ``G0``, so the group should be
        adopted only if the scale term is real.

``N1``  Coverage.  Under the sealed *passive* support protocol, how often does a
        k <= 5 draw contain a context holding >= 2 supports (usable under ``G0``)
        or >= 3 (usable under ``G``)?  If it rarely does, NEA has no deployment
        path at k <= 5 regardless of its other merits.

Only source ``fit`` is opened.  ``probe`` is excluded because it was consumed once
by PIRS; ``locked`` and the recipient roster are never requested.
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
import torch

from research.a2s.a2s_information_gate import canonical, sha256_file
from research.a2s.a2s_trace import DEVICE, Substrate, load_substrate
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK, DEFAULT_OOF, paired_bootstrap
from research.a2s.a2s_mode_gates import build_basis, fit_head
from research.a2s.a2s_transfer_object_gate import (
    document_offset_prediction,
    pair_concordance,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_nea_preconditions_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_nea_preconditions_records_2026-08-02.parquet"

SEED = 20260802
BOOTSTRAP_DRAWS = 2000
MDE = 0.005
MIN_EVAL_ROWS = 8
MIN_TRAIN_ROWS = 20
MIN_CONTEXT_ROWS = 5
COVERAGE_DRAWS = 200
K_SWEEP = (1, 3, 5, 10)

D0_MIN_HEADROOM = 0.005          # lower 95% bound, same-document pairs, separated split
D0_MAX_ORACLE_LEAK = 0.005       # document oracle must be powerless on a clean split
N0_SCALE_MATERIAL = 0.15         # SD of log within-context dispersion to call scale real
N1_MIN_COVERAGE = 0.50           # fraction of k=5 draws with a usable context


@dataclass
class SeparatedSplit:
    target: str
    component: str
    train_rows: np.ndarray
    eval_rows: np.ndarray
    regime: str                  # "separated" | "scaffold_only"


# --------------------------------------------------------------------------- #
# D0 -- simultaneously separated within-target splits
# --------------------------------------------------------------------------- #


def separated_splits(substrate: Substrate) -> list[SeparatedSplit]:
    """Within-target splits that are scaffold-, document- and assay-disjoint.

    Documents are partitioned first because they are the coarsest provenance
    unit; scaffold and assay collisions are then *removed from the evaluation
    side*, so the guarantee is exact rather than approximate.
    """

    frame = substrate.labeled
    splits: list[SeparatedSplit] = []
    for target, group in frame.loc[frame.role == "fit"].groupby("target", sort=True):
        rows = group.index.to_numpy()
        if len(rows) < MIN_TRAIN_ROWS + MIN_EVAL_ROWS:
            continue
        component = str(group.component.iloc[0])
        documents = group.docs.astype(str).to_numpy()
        assays = group.assays.astype(str).to_numpy()
        scaffolds = group.scaffold.astype(str).to_numpy()

        counts = pd.Series(documents).value_counts()
        eligible = [d for d in counts.index if counts[d] >= MIN_EVAL_ROWS]
        if len(eligible) < 2:
            continue
        digest = int(sha256(f"{SEED}:{target}".encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(digest)

        # Grow the held-out document set until it holds enough rows.
        order = rng.permutation(len(eligible))
        held: set[str] = set()
        for index in order:
            held.add(eligible[index])
            if np.isin(documents, list(held)).sum() >= max(MIN_EVAL_ROWS, int(0.3 * len(rows))):
                break
        eval_mask = np.isin(documents, list(held))
        train_mask = ~eval_mask
        if train_mask.sum() < MIN_TRAIN_ROWS:
            continue

        # Exact scaffold and assay separation, enforced on the evaluation side.
        train_scaffolds = set(scaffolds[train_mask].tolist())
        train_assays = set(assays[train_mask].tolist())
        clean = eval_mask & ~np.isin(scaffolds, list(train_scaffolds))
        clean &= ~np.isin(assays, list(train_assays))
        if clean.sum() >= MIN_EVAL_ROWS:
            splits.append(
                SeparatedSplit(str(target), component, rows[train_mask], rows[clean], "separated")
            )

        # The same targets under scaffold-only separation, for attribution.
        unique = sorted(set(scaffolds))
        shuffled = rng.permutation(len(unique))
        chosen: set[str] = set()
        for index in shuffled:
            chosen.add(unique[index])
            if np.isin(scaffolds, list(chosen)).sum() >= max(MIN_EVAL_ROWS, int(0.3 * len(rows))):
                break
        scaffold_eval = np.isin(scaffolds, list(chosen))
        if scaffold_eval.sum() >= MIN_EVAL_ROWS and (~scaffold_eval).sum() >= MIN_TRAIN_ROWS:
            splits.append(
                SeparatedSplit(
                    str(target), component, rows[~scaffold_eval], rows[scaffold_eval], "scaffold_only"
                )
            )
    return splits


def evaluate_d0(substrate: Substrate, basis: np.ndarray, splits: list[SeparatedSplit]) -> pd.DataFrame:
    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    label = substrate.affinity.cpu().numpy().astype(np.float64)
    base = substrate.base.cpu().numpy().astype(np.float64)
    documents = substrate.labeled.docs.astype(str).to_numpy()
    assays = substrate.labeled.assays.astype(str).to_numpy()

    records: list[dict[str, object]] = []
    for split in splits:
        train_rows, eval_rows = split.train_rows, split.eval_rows
        eval_docs = documents[eval_rows]
        reference = pair_concordance(label[eval_rows], base[eval_rows], eval_docs)
        if not np.isfinite(reference["ci"]):
            continue

        weight, _ = fit_head(basis[train_rows], residual[train_rows])
        own = pair_concordance(
            label[eval_rows], base[eval_rows] + basis[eval_rows] @ weight, eval_docs
        )
        offset = document_offset_prediction(
            documents[train_rows], residual[train_rows], eval_docs
        )
        oracle = pair_concordance(label[eval_rows], base[eval_rows] + offset, eval_docs)

        records.append({
            "target": split.target,
            "component": split.component,
            "regime": split.regime,
            "n_train": int(len(train_rows)),
            "n_eval": int(len(eval_rows)),
            "eval_documents": int(len(set(eval_docs.tolist()))),
            "document_overlap": float(
                np.mean([d in set(documents[train_rows].tolist()) for d in eval_docs])
            ),
            "assay_overlap": float(
                np.mean([a in set(assays[train_rows].tolist()) for a in assays[eval_rows]])
            ),
            "base__ci": reference["ci"],
            "base__ci_within": reference["ci_within"],
            "own__ci": own["ci"],
            "own__ci_within": own["ci_within"],
            "docoffset__ci": oracle["ci"],
            "pairs_within": reference["pairs_within"],
            "pairs": reference["pairs"],
        })
    return pd.DataFrame.from_records(records)


# --------------------------------------------------------------------------- #
# N0 -- which nuisance group acts?
# --------------------------------------------------------------------------- #


def nuisance_group(substrate: Substrate) -> dict[str, object]:
    residual = (substrate.affinity - substrate.base).cpu().numpy().astype(np.float64)
    frame = substrate.labeled.loc[substrate.labeled.role == "fit"].copy()
    frame["residual"] = residual[frame.index.to_numpy()]

    report: dict[str, object] = {}
    for name in ("docs", "assays"):
        table = frame.groupby(["target", name])["residual"]
        sizes = table.transform("size")
        usable = frame.loc[sizes >= MIN_CONTEXT_ROWS]
        if usable.empty:
            continue
        grouped = usable.groupby(["target", name])["residual"]
        offsets = grouped.mean()
        scales = grouped.std(ddof=1).clip(lower=1e-6)

        total = float(np.var(usable.residual.to_numpy()))
        centred = usable.residual - grouped.transform("mean")
        after_offset = float(np.var(centred.to_numpy()))
        studentised = centred / grouped.transform(lambda v: max(float(v.std(ddof=1)), 1e-6))
        report[name] = {
            "contexts": int(len(offsets)),
            "median_rows_per_context": float(usable.groupby(["target", name]).size().median()),
            "offset_sd": float(np.std(offsets.to_numpy())),
            "scale_median": float(np.median(scales.to_numpy())),
            "log_scale_sd": float(np.std(np.log(scales.to_numpy()))),
            "variance_total": total,
            "variance_after_offset": after_offset,
            "variance_explained_by_offset": float(1.0 - after_offset / total) if total else float("nan"),
            "variance_after_studentising": float(np.var(studentised.to_numpy())),
        }
    scale_evidence = max(
        (value["log_scale_sd"] for value in report.values() if isinstance(value, dict)),
        default=float("nan"),
    )
    report["scale_is_material"] = bool(scale_evidence > N0_SCALE_MATERIAL)
    report["recommended_group"] = "G_affine" if report["scale_is_material"] else "G0_offset_only"
    return report


# --------------------------------------------------------------------------- #
# N1 -- coverage under the passive support protocol
# --------------------------------------------------------------------------- #


def coverage(substrate: Substrate) -> dict[str, object]:
    frame = substrate.labeled.loc[substrate.labeled.role == "fit"]
    rng = np.random.default_rng(SEED)
    report: dict[str, object] = {}
    for name in ("docs", "assays"):
        per_k: dict[str, dict[str, float]] = {}
        for k in K_SWEEP:
            usable_pair, usable_triple, pairs = [], [], []
            for _, group in frame.groupby("target", sort=True):
                contexts = group[name].astype(str).to_numpy()
                if len(contexts) < k:
                    continue
                for _ in range(COVERAGE_DRAWS // 10):
                    draw = contexts[rng.choice(len(contexts), size=k, replace=False)]
                    counts = pd.Series(draw).value_counts().to_numpy()
                    usable_pair.append(float((counts >= 2).any()))
                    usable_triple.append(float((counts >= 3).any()))
                    pairs.append(float(sum(c * (c - 1) / 2 for c in counts)))
            if usable_pair:
                per_k[f"k{k}"] = {
                    "fraction_with_a_context_of_2plus": float(np.mean(usable_pair)),
                    "fraction_with_a_context_of_3plus": float(np.mean(usable_triple)),
                    "mean_within_context_pairs": float(np.mean(pairs)),
                }
        report[name] = per_k
    return report


# --------------------------------------------------------------------------- #
# Summary and decision
# --------------------------------------------------------------------------- #


CONTRASTS = {
    "own_minus_base": ("own__ci", "base__ci"),
    "own_minus_base_within_document": ("own__ci_within", "base__ci_within"),
    "docoffset_minus_base": ("docoffset__ci", "base__ci"),
}


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for regime in sorted(records.regime.unique()):
        frame = records.loc[records.regime == regime].copy()
        cells: dict[str, object] = {}
        for name, (left, right) in CONTRASTS.items():
            frame[name] = frame[left] - frame[right]
            cells[name] = paired_bootstrap(frame, name, draws=BOOTSTRAP_DRAWS)
        cells["descriptives"] = {
            "targets": int(frame.target.nunique()),
            "components": int(frame.component.nunique()),
            "mean_document_overlap": float(frame.document_overlap.mean()),
            "mean_assay_overlap": float(frame.assay_overlap.mean()),
            "mean_eval_rows": float(frame.n_eval.mean()),
            "mean_eval_documents": float(frame.eval_documents.mean()),
        }
        summary[regime] = cells
    return summary


def decide(
    summary: dict[str, object], groups: dict[str, object], cover: dict[str, object]
) -> dict[str, object]:
    separated = summary.get("separated", {})

    def cell(name: str) -> dict[str, float]:
        value = separated.get(name)
        return value if isinstance(value, dict) else {"mean": float("nan"), "lower95": float("nan")}

    headroom = cell("own_minus_base_within_document")
    leak = cell("docoffset_minus_base")
    gates: dict[str, object] = {
        "D0": {
            "question": "does chemical headroom survive simultaneous separation",
            "within_document_headroom": headroom,
            "all_pair_headroom": cell("own_minus_base"),
            "document_oracle_leak": leak,
            "split_is_clean": bool(
                abs(leak.get("mean", 1.0)) < D0_MAX_ORACLE_LEAK
                and separated.get("descriptives", {}).get("mean_document_overlap", 1.0) < 1e-9
            ),
            "pass": bool(headroom.get("lower95", float("nan")) > D0_MIN_HEADROOM),
        },
        "N0": {
            "question": "which nuisance group acts",
            "report": groups,
            "pass": True,   # descriptive; it selects the group rather than gating
        },
        "N1": {
            "question": "is there a deployment path at k<=5",
            "report": cover,
            "pass": bool(
                cover.get("docs", {}).get("k5", {}).get("fraction_with_a_context_of_2plus", 0.0)
                > N1_MIN_COVERAGE
            ),
        },
    }

    if not gates["D0"]["split_is_clean"]:
        verdict = "D0_SPLIT_NOT_CLEAN_HARNESS_INVALID"
    elif not gates["D0"]["pass"]:
        verdict = "NO_CHEMICAL_ADAPTATION_OBJECT_SURVIVES_SEPARATION_STOP_PROGRAMME"
    elif not gates["N1"]["pass"]:
        verdict = "OBJECT_SURVIVES_BUT_NEA_HAS_NO_K5_DEPLOYMENT_PATH"
    else:
        verdict = "NEA_PRECONDITIONS_MET_PROCEED_TO_MECHANISM_DESIGN"
    return {"gates": gates, "verdict": verdict}


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path) -> dict[str, object]:
    started = time.time()
    substrate, provenance = load_substrate(lock_path, oof_cache)
    roles = set(substrate.labeled.role.unique())
    if "locked" in roles:
        raise RuntimeError("role firewall breached")

    basis, basis_stats = build_basis(substrate)
    splits = separated_splits(substrate)
    records = evaluate_d0(substrate, basis, splits)
    summary = summarise(records)
    groups = nuisance_group(substrate)
    cover = coverage(substrate)
    decision = decide(summary, groups, cover)

    payload: dict[str, object] = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "device": str(DEVICE),
        "firewall": {
            "evaluated_role": "fit",
            "probe_excluded": True,
            "probe_exclusion_reason": "consumed once by PIRS; not reusable for model selection",
            "locked_requested": False,
            "recipient_requested": False,
            "fits_closed_form_ridge_heads_from_labels": True,
            "trains_no_gradient_model": True,
        },
        "protocol": {
            "seed": SEED,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "mde": MDE,
            "basis": basis_stats,
            "thresholds": {
                "D0_min_headroom": D0_MIN_HEADROOM,
                "D0_max_oracle_leak": D0_MAX_ORACLE_LEAK,
                "N0_scale_material": N0_SCALE_MATERIAL,
                "N1_min_coverage": N1_MIN_COVERAGE,
            },
        },
        "summary": summary,
        "nuisance_group": groups,
        "coverage": cover,
        **decision,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    payload["artifacts"] = {
        "json": str(output),
        "records": str(records_path),
        "records_sha256": sha256_file(records_path),
        "lock": provenance["lock"].get("content_sha256") if isinstance(provenance.get("lock"), dict) else None,
    }
    payload["content_sha256"] = sha256(canonical(payload).encode()).hexdigest()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="NEA preconditions: gates D0, N0, N1")
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.lock, arguments.output, arguments.records, arguments.oof_cache)
    print(json.dumps({
        "verdict": payload["verdict"],
        "gates": {name: gate.get("pass") for name, gate in payload["gates"].items()},
    }, indent=2))


if __name__ == "__main__":
    main()
