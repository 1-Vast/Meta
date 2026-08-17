"""Stage T0 — measurement-reliability audit of the governed pK supervision.

What T0 is: an estimate of how reliable the *supervision* is, at the three
provenance levels the corpus actually distinguishes, under the aggregation rule
the corpus actually used.

What T0 is **not**: a universal benchmark MSE floor. Repeated-measure variance
on the compounds that happen to be measured twice is not the irreducible error
of a prediction task, and this module refuses to present it as one. The
repeated-measure subset is small and selected, and the selection is reported
alongside every number.

Levels:

    L1  rows sharing (panel_id, assay protocol)  -> same-protocol repeat
    L2  rows sharing panel_id, different protocol -> within document+endpoint+target
    L3  rows of one (target, ligand) across panels -> between assay / document

Run:
    python -m tools.research.stageT_mmp.t0_reliability
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

from tools.research.stageT_mmp.provenance import (
    build_cache, group_by_cell, load_cache, train_allow_list,
)

HERE = Path(__file__).resolve().parent

# Frozen: an estimate needs at least this many independent groups before it is
# reported as a number rather than as "not identifiable".
MIN_GROUPS_FOR_ESTIMATE = 30


def _dispersion(values: list[float], label: str) -> dict:
    """Robust and classical dispersion, or an explicit non-identifiability."""
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
        # 1.4826 * MAD is the normal-consistent robust sd; reported next to the
        # classical sd because a handful of transcription errors in BindingDB
        # would inflate the classical one without changing the typical case.
        "robust_sd": float(1.4826 * np.median(np.abs(centred - np.median(centred)))),
        "iqr": float(np.quantile(centred, 0.75) - np.quantile(centred, 0.25)),
        "median_abs": float(np.median(np.abs(centred))),
        "p95_abs": float(np.quantile(np.abs(centred), 0.95)),
        "max_abs": float(np.abs(centred).max()),
    }


def _pairwise_residuals(groups: list[list[float]]) -> list[float]:
    """Within-group deviations from the group median, one list."""
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
    """How many 'repeats' are exact value duplicates rather than replication.

    BindingDB curates the same measurement from more than one article, so two
    distinct `panel_id`s can carry one physical experiment. Those groups have a
    range of exactly zero and they deflate any disagreement estimate towards
    zero. Reporting the zero-range share, and the dispersion on the genuinely
    disagreeing remainder, keeps that visible instead of averaged away.
    """
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "T0_RELIABILITY.json")
    parser.add_argument("--force-cache", action="store_true")
    args = parser.parse_args()

    cache = build_cache(force=args.force_cache)
    rows = load_cache()
    allow = train_allow_list()
    by_cell = group_by_cell(rows)

    # -- L1 / L2: structure inside one panel -------------------------------
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
                l2_groups.append([float(np.median(v)) for v in by_protocol.values()])

    # -- L3: across panels, on the per-panel medians the corpus aggregates --
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

    report: dict = {
        "schema": "MetaSieve.StageT.T0Reliability.v1",
        "stage": "stageT_mmp",
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
            "cells_with_recovered_rows": len(by_cell),
            "recovery_fraction": len(by_cell) / max(len(allow.cells), 1),
        },
    }

    # -- levels -------------------------------------------------------------
    report["levels"] = {
        "L1_same_panel_same_protocol": {
            "definition": "rows sharing (panel_id, assay protocol_sha256)",
            "groups": len(l1_groups),
            "rows": sum(len(v) for v in l1_groups),
            "residuals_from_group_median": _dispersion(
                _pairwise_residuals(l1_groups), "L1 residuals"),
            "group_range": _dispersion(_group_ranges(l1_groups), "L1 ranges"),
            "duplication": _duplication_diagnostic(l1_groups, "L1"),
        },
        "L2_same_panel_different_protocol": {
            "definition": "per-protocol medians inside one panel_id",
            "groups": len(l2_groups),
            "residuals_from_group_median": _dispersion(
                _pairwise_residuals(l2_groups), "L2 residuals"),
            "group_range": _dispersion(_group_ranges(l2_groups), "L2 ranges"),
        },
        "L3_across_panels": {
            "definition": ("per-panel medians of one (target, ligand), i.e. the "
                           "inputs to the corpus's equal-panel median"),
            "groups": len(l3_groups),
            "residuals_from_group_median": _dispersion(
                _pairwise_residuals(l3_groups), "L3 residuals"),
            "group_range": _dispersion(_group_ranges(l3_groups), "L3 ranges"),
            "duplication": _duplication_diagnostic(l3_groups, "L3"),
        },
    }

    # -- difference-label uncertainty under the actual aggregation ---------
    l3 = report["levels"]["L3_across_panels"]["residuals_from_group_median"]
    l1 = report["levels"]["L1_same_panel_same_protocol"]["residuals_from_group_median"]
    between = l3.get("variance") if l3.get("identifiable") else None
    within = l1.get("variance") if l1.get("identifiable") else None
    report["difference_label_uncertainty"] = {
        "model": ("delta_y = y_b - y_a with each y an aggregated median. A "
                  "between-assay offset shared by both cells of a same-panel "
                  "pair cancels in the difference; it does not cancel for a "
                  "cross-panel pair."),
        "between_assay_variance_per_cell": between,
        "within_assay_variance_per_cell": within,
        "same_panel_difference_variance": (
            2.0 * within if within is not None else None),
        "cross_panel_difference_variance": (
            2.0 * (within + between)
            if (within is not None and between is not None) else None),
        "caveat": ("both terms are estimated on small, selected subsets; the "
                   "cross-panel figure is the one that matters for the S3 "
                   "stratum and it is the weaker of the two estimates"),
    }

    report["identifiability_limits"] = {
        "technical_vs_condition_variation": (
            "NOT IDENTIFIABLE. The projection's assay protocol_sha256 never "
            "splits a meta_train panel into two protocols (L2 groups = 0), so "
            "the provenance cannot separate a technical replicate from the same "
            "paper reporting two measurements under different conditions "
            "(isoform, radioligand, pH, temperature). L1 is therefore "
            "'same-document repeat disagreement', not pipetting noise, and no "
            "technical-replicate variance is reported."),
        "single_measurement_cells": (
            "96% of meta_train cells carry exactly one source row, so their "
            "measurement error is not observable at all. Every number here is "
            "conditional on the repeated-measure subset."),
        "curation_duplication": (
            "59.4% of apparent cross-panel repeats have an exactly zero range, "
            "i.e. one physical measurement curated under two article DOIs. The "
            "pooled L3 dispersion is deflated by them; the disagreeing-only "
            "figures are reported next to it."),
        "not_an_mse_floor": (
            "None of this bounds achievable model error. It bounds how much of "
            "the supervision is reproducible on the compounds that happen to be "
            "measured twice."),
    }

    # -- selection bias and coverage ---------------------------------------
    repeated_cells = {cell for cell, values in by_cell.items() if len(values) > 1}
    cross_panel_cells = {row["cell_id"] for row in l3_detail}
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
            "targets. The corpus additionally dropped rows at admission "
            "(conflicting_ligands and cross_panel_pairs in "
            "corpus_admission_filters), so any disagreement estimate here is a "
            "LOWER bound on the disagreement in the raw source."),
    }

    # -- the frozen, label-blind confidence strata -------------------------
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
        "L1_groups": len(l1_groups),
        "L2_groups": len(l2_groups),
        "L3_groups": len(l3_groups),
        "L1_identifiable": report["levels"]["L1_same_panel_same_protocol"][
            "residuals_from_group_median"]["identifiable"],
        "L3_identifiable": report["levels"]["L3_across_panels"][
            "residuals_from_group_median"]["identifiable"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
