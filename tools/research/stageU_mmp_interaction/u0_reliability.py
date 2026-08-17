"""Stage U0 — measurement-reliability audit of the governed pK supervision.

What U0 is: an estimate of how reliable the *supervision* is, at the three
provenance levels the corpus actually distinguishes, under the aggregation rule
the corpus actually used.

What U0 is **not**: a universal benchmark MSE floor. Repeated-measure variance
on the compounds that happen to be measured twice is not the irreducible error
of a prediction task, and this module refuses to present it as one. The
repeated-measure subset is small and selected, and the selection is reported
alongside every number.

Levels:

    L1  rows sharing (panel_id, assay protocol)  -> same-document repeat
    L2  rows sharing panel_id, different protocol -> within document+endpoint+target
    L3  rows of one (target, ligand) across panels -> between assay / document

Run:
    python -m tools.research.stageU_mmp_interaction.u0_reliability
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tools.research.stageU_mmp_interaction.provenance import (
    build_cache, group_by_cell, load_cache, train_allow_list,
)

HERE = Path(__file__).resolve().parent

MIN_GROUPS_FOR_ESTIMATE = 30
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20260820


def _dispersion(values: list[float], label: str) -> dict:
    array = np.asarray(values, dtype=np.float64)
    if array.size < MIN_GROUPS_FOR_ESTIMATE:
        return {
            "identifiable": False,
            "groups": int(array.size),
            "minimum_required": MIN_GROUPS_FOR_ESTIMATE,
            "reason": (f"{label}: only {array.size} independent groups; the "
                       "dispersion is not identifiable and no value is invented"),
        }
    centred = array
    return {
        "identifiable": True,
        "groups": int(array.size),
        "sd": float(np.std(centred, ddof=1)),
        "variance": float(np.var(centred, ddof=1)),
        "robust_sd": float(1.4826 * np.median(np.abs(centred - np.median(centred)))),
        "iqr": float(np.quantile(centred, 0.75) - np.quantile(centred, 0.25)),
        "median_abs": float(np.median(np.abs(centred))),
        "p95_abs": float(np.quantile(np.abs(centred), 0.95)),
        "max_abs": float(np.abs(centred).max()),
    }


def _pairwise_residuals(groups: list[list[float]]) -> list[float]:
    out: list[float] = []
    for values in groups:
        if len(values) < 2:
            continue
        centre = float(np.median(values))
        out.extend(float(value - centre) for value in values)
    return out


def _group_ranges(groups: list[list[float]]) -> list[float]:
    return [float(max(v) - min(v)) for v in groups if len(v) >= 2]


def _duplication_diagnostic(groups: list[list[float]], label: str) -> dict:
    ranges = _group_ranges(groups)
    if not ranges:
        return {"groups": 0, "identifiable": False,
                "reason": f"{label}: no multi-row groups"}
    exact = [value for value in ranges if value == 0.0]
    disagreeing = [values for values in groups
                   if len(values) >= 2 and max(values) - min(values) > 0.0]
    return {
        "groups": len(ranges),
        "exact_duplicate_groups": len(exact),
        "exact_duplicate_share": len(exact) / len(ranges),
        "disagreeing_groups": len(disagreeing),
        "residuals_on_disagreeing_only": _dispersion(
            _pairwise_residuals(disagreeing), f"{label} disagreeing residuals"),
        "range_on_disagreeing_only": _dispersion(
            _group_ranges(disagreeing), f"{label} disagreeing ranges"),
        "note": ("a zero-range group is one physical measurement curated more "
                 "than once; it is evidence of duplication, not of agreement, "
                 "and it biases the pooled dispersion downwards"),
    }


def build_levels(rows: list[dict]) -> dict:
    """L1/L2/L3 group structures and pooled dispersions."""
    by_cell = group_by_cell(rows)
    l1_groups: list[list[float]] = []
    l2_groups: list[list[float]] = []
    panel_level: dict[tuple[str, str], list[float]] = defaultdict(list)
    for cell_id, cell_rows in by_cell.items():
        by_panel: dict[str, list[dict]] = defaultdict(list)
        for row in cell_rows:
            by_panel[row["panel_id"]].append(row)
        for panel, panel_rows in by_panel.items():
            panel_level[(cell_id, panel)] = [r["pK"] for r in panel_rows]
            by_protocol: dict[tuple[str, ...], list[float]] = defaultdict(list)
            for row in panel_rows:
                key = tuple(row.get("assay_protocols") or ())
                by_protocol[key].append(row["pK"])
            for values in by_protocol.values():
                if len(values) >= 2:
                    l1_groups.append(values)
            if len(by_protocol) >= 2:
                l2_groups.append([float(np.median(v))
                                  for v in by_protocol.values()])

    l3_groups: list[list[float]] = []
    l3_detail: list[dict] = []
    for cell_id, cell_rows in by_cell.items():
        panels = sorted({row["panel_id"] for row in cell_rows})
        if len(panels) < 2:
            continue
        medians = []
        for panel in panels:
            values = [r["pK"] for r in cell_rows if r["panel_id"] == panel]
            medians.append(float(np.median(values)))
        l3_groups.append(medians)
        l3_detail.append({
            "cell_id": cell_id,
            "panels": len(panels),
            "documents": len({row["document_id"] for row in cell_rows}),
            "range": float(max(medians) - min(medians)),
        })
    return {
        "by_cell": by_cell,
        "l1_groups": l1_groups,
        "l2_groups": l2_groups,
        "l3_groups": l3_groups,
        "l3_detail": l3_detail,
    }


def _cluster_bootstrap_variance(groups: list[list[float]], draws: int,
                                seed: int) -> dict:
    """Resample repeated-measure groups (not rows) and pool residual variance."""
    usable = [values for values in groups if len(values) >= 2]
    if len(usable) < MIN_GROUPS_FOR_ESTIMATE:
        return {"identifiable": False, "groups": len(usable)}
    rng = np.random.default_rng(seed)
    samples: list[float] = []
    for _ in range(draws):
        picked = [usable[int(i)] for i in rng.integers(0, len(usable),
                                                       size=len(usable))]
        residuals: list[float] = []
        for values in picked:
            centre = float(np.median(values))
            residuals.extend(float(value - centre) for value in values)
        if len(residuals) >= 2:
            samples.append(float(np.var(residuals, ddof=1)))
    values = np.asarray(samples)
    return {
        "identifiable": True,
        "draws": draws,
        "lo": float(np.quantile(values, 0.025)),
        "median": float(np.quantile(values, 0.5)),
        "hi": float(np.quantile(values, 0.975)),
        "note": ("bootstrap unit is the repeated-measure group; residual rows "
                 "are not treated as independent"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "U0_RELIABILITY.json")
    parser.add_argument("--force-cache", action="store_true")
    args = parser.parse_args()

    cache = build_cache(force=args.force_cache)
    rows = load_cache()
    allow = train_allow_list()
    levels = build_levels(rows)

    report: dict = {
        "schema": "MetaSieve.StageU.U0Reliability.v1",
        "stage": "stageU_mmp_interaction",
        "preregistration_sha256": "fdc0a830aa92882d07b9aea50f22a4c72fc6d93f92c55a3be6bc15cd6a645c11",
        "purpose": ("measure the reliability of the governed pK supervision at "
                    "the provenance levels the corpus distinguishes"),
        "not_a_claim": (
            "these dispersions are repeated-measure variability on the subset of "
            "compounds that happen to be measured more than once. They are NOT a "
            "universal benchmark MSE floor, NOT an irreducible prediction error, "
            "and must never be quoted as either."),
        "provenance": cache,
        "aggregation_rule_actually_used": cache["authority"]["aggregation_rule"],
        "corpus_admission_filters": cache["authority"]["admission"],
        "population": {
            "meta_train_cells": len(allow.cells),
            "meta_train_source_rows": len(allow.row_ids),
            "rows_recovered": len(rows),
            "cells_with_recovered_rows": len(levels["by_cell"]),
            "recovery_fraction": len(levels["by_cell"]) / max(len(allow.cells), 1),
        },
    }

    report["levels"] = {
        "L1_same_panel_same_protocol": {
            "definition": "rows sharing (panel_id, assay protocol_sha256)",
            "groups": len(levels["l1_groups"]),
            "rows": sum(len(v) for v in levels["l1_groups"]),
            "residuals_from_group_median": _dispersion(
                _pairwise_residuals(levels["l1_groups"]), "L1 residuals"),
            "group_range": _dispersion(
                _group_ranges(levels["l1_groups"]), "L1 ranges"),
            "duplication": _duplication_diagnostic(levels["l1_groups"], "L1"),
        },
        "L2_same_panel_different_protocol": {
            "definition": "per-protocol medians inside one panel_id",
            "groups": len(levels["l2_groups"]),
            "residuals_from_group_median": _dispersion(
                _pairwise_residuals(levels["l2_groups"]), "L2 residuals"),
            "group_range": _dispersion(
                _group_ranges(levels["l2_groups"]), "L2 ranges"),
        },
        "L3_across_panels": {
            "definition": ("per-panel medians of one (target, ligand), i.e. the "
                           "inputs to the corpus's equal-panel median"),
            "groups": len(levels["l3_groups"]),
            "residuals_from_group_median": _dispersion(
                _pairwise_residuals(levels["l3_groups"]), "L3 residuals"),
            "group_range": _dispersion(
                _group_ranges(levels["l3_groups"]), "L3 ranges"),
            "duplication": _duplication_diagnostic(levels["l3_groups"], "L3"),
        },
    }

    l3 = report["levels"]["L3_across_panels"]["residuals_from_group_median"]
    l1 = report["levels"]["L1_same_panel_same_protocol"]["residuals_from_group_median"]
    between = l3.get("variance") if l3.get("identifiable") else None
    within = l1.get("variance") if l1.get("identifiable") else None
    within_boot = _cluster_bootstrap_variance(
        levels["l1_groups"], BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)
    between_boot = _cluster_bootstrap_variance(
        levels["l3_groups"], BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)

    def _twice(value, boot):
        if value is None or not boot.get("identifiable"):
            return {"identifiable": False, "point": None}
        return {
            "identifiable": True,
            "point": 2.0 * value,
            "bootstrap_ci": {
                "lo": 2.0 * boot["lo"],
                "median": 2.0 * boot["median"],
                "hi": 2.0 * boot["hi"],
                "draws": boot["draws"],
                "seed": BOOTSTRAP_SEED,
                "unit": "repeated-measure group",
            },
        }

    report["difference_label_uncertainty"] = {
        "model": ("delta_y = y_b - y_a with each y an aggregated median. A "
                  "between-assay offset shared by both cells of a same-panel "
                  "pair cancels in the difference; it does not cancel for a "
                  "cross-panel pair."),
        "between_assay_variance_per_cell": {
            "point": between,
            "bootstrap_ci": between_boot,
        },
        "within_assay_variance_per_cell": {
            "point": within,
            "bootstrap_ci": within_boot,
        },
        "same_panel_difference_variance": _twice(within, within_boot),
        "cross_panel_difference_variance": (
            {"identifiable": True,
             "point": 2.0 * (within + between),
             "bootstrap_ci": {
                 "lo": 2.0 * (within_boot["lo"] + between_boot["lo"]),
                 "median": 2.0 * (within_boot["median"] + between_boot["median"]),
                 "hi": 2.0 * (within_boot["hi"] + between_boot["hi"]),
                 "draws": BOOTSTRAP_DRAWS,
                 "seed": BOOTSTRAP_SEED,
                 "unit": "repeated-measure group",
             }}
            if (within is not None and between is not None
                and within_boot.get("identifiable")
                and between_boot.get("identifiable"))
            else {"identifiable": False, "point": None}),
        "caveat": ("both terms are estimated on small, selected subsets; the "
                   "cross-panel figure is the one that matters for the S3 "
                   "stratum and it is the weaker of the two estimates"),
    }

    report["identifiability_limits"] = {
        "technical_vs_condition_variation": (
            "NOT IDENTIFIABLE. If the projection's assay protocol_sha256 never "
            "splits a meta_train panel into two protocols (L2 groups = 0), the "
            "provenance cannot separate a technical replicate from the same "
            "paper reporting two measurements under different conditions."),
        "single_measurement_cells": (
            "96% of meta_train cells carry exactly one source row, so their "
            "measurement error is not observable at all. Every number here is "
            "conditional on the repeated-measure subset."),
        "curation_duplication": (
            "a large zero-range share means one physical measurement curated "
            "under two article DOIs; the pooled dispersion is deflated by them."),
        "not_an_mse_floor": (
            "None of this bounds achievable model error. It bounds how much of "
            "the supervision is reproducible on the compounds that happen to be "
            "measured twice."),
    }

    repeated_cells = {cell for cell, values in levels["by_cell"].items()
                      if len(values) > 1}
    cross_panel_cells = {row["cell_id"] for row in levels["l3_detail"]}
    ligand_counts = Counter(
        allow.cells[cell]["ligand_id"] for cell in repeated_cells
        if cell in allow.cells)
    report["selection_bias"] = {
        "cells_with_more_than_one_row": len(repeated_cells),
        "share_of_meta_train_cells": len(repeated_cells) / max(len(allow.cells), 1),
        "cells_measured_in_more_than_one_panel": len(cross_panel_cells),
        "distinct_ligands_in_repeated_subset": len(ligand_counts),
        "warning": (
            "a (target, ligand) pair appears twice only if two groups chose to "
            "measure it, which favours reference compounds and well-studied "
            "targets. The corpus additionally dropped rows at admission, so any "
            "disagreement estimate here is a LOWER bound on the disagreement in "
            "the raw source."),
    }

    report["confidence_strata"] = {
        "S1_same_panel_single": ("both cells share a panel and each has "
                                 "panel_count == 1"),
        "S2_same_panel_multi": ("both cells share a panel, at least one has "
                                "panel_count > 1"),
        "S3_cross_panel": ("no shared panel; weak/noise stratum, never pooled "
                           "into the primary bank"),
        "rule": ("defined by provenance only, so no label value can move a pair "
                 "between strata"),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "rows_recovered": len(rows),
        "L1_groups": len(levels["l1_groups"]),
        "L2_groups": len(levels["l2_groups"]),
        "L3_groups": len(levels["l3_groups"]),
        "same_panel_difference_variance": report[
            "difference_label_uncertainty"]["same_panel_difference_variance"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
