"""Stage 0 — forensic reconciliation of Stage T's transformation key.

Stage T defined

    exact_key = f"{attachment_context}|{R_a}>>{R_b}"

which **omits the shared core**, then aggregated every observation of one
`(exact_key, target)` with a median and built a descriptor that also never sees
the core. The requested estimand was

    tau = (shared core, R_a -> R_b, attachment context, stereochemistry, charge)

so Stage T's `tau` is a *coarsening* of the requested one. The consequence is
algebraic: `D(tau,t1,t2)` cancels `mu_tau` exactly only when both targets
realise the **same** tau. Under a core-blind key, two targets can realise the
same nominal edit on different scaffolds, and then

    D = [mu_{core_1} - mu_{core_2}] + [delta(t1,tau) - delta(t2,tau)] + noise

carries a generic chemical residual that has nothing to do with protein biology.

This module measures how large that residual actually is. The decisive statistic
is computed **within a single target**, where the protein is held fixed by
construction, so every bit of the measured spread is generic chemical context.

No model is trained. Stage T's frozen artifacts are read-only here.

Run:
    python -m tools.research.stageV_core_mmp.stage0_forensics
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from tools.research.stageT_mmp.observations import build_observations, load_governed

HERE = Path(__file__).resolve().parent
STAGE_T = Path(__file__).resolve().parents[1] / "stageT_mmp"


def core_inclusive_key(observation) -> str:
    """The requested estimand's key: core + attachment context + directed edit.

    Stage T's `exact_key` is already `context|R_a>>R_b`, so prefixing the core
    reconstructs the full tau without re-fragmenting anything.
    """
    payload = f"{observation.core}|{observation.exact_key}"
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _stats(values) -> dict:
    array = np.asarray(list(values), dtype=np.float64)
    if array.size == 0:
        return {"n": 0}
    return {
        "n": int(array.size),
        "mean": float(array.mean()),
        "sd": float(array.std(ddof=1)) if array.size > 1 else 0.0,
        "median": float(np.median(array)),
        "p95": float(np.quantile(array, 0.95)),
        "max": float(array.max()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "STAGE0_FORENSICS.json")
    args = parser.parse_args()

    data, seal = load_governed()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to run without the physical split view")
    fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    fit_targets = sorted(t for t, c in component_of.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component_of.items()
                              if c in set(internal))

    built = build_observations(data, fit_targets + internal_targets)
    same_panel = [o for o in built["observations"] if o.same_panel]
    fit_set = set(fit_targets)
    fit_rows = [o for o in same_panel if o.target in fit_set]
    internal_rows = [o for o in same_panel if o.target not in fit_set]

    report: dict = {
        "schema": "MetaSieve.StageV.Stage0Forensics.v1",
        "purpose": ("quantify the effect of Stage T's core-blind transformation "
                    "key on the crossed double difference"),
        "meta_test": seal,
        "stage_t_key_definition": "f'{attachment_context}|{R_a}>>{R_b}'  (no core)",
        "requested_key_definition": (
            "sha256(core_isomeric | attachment_context | R_a>>R_b)"),
        "stage_t_artifacts_read_only": True,
    }

    # -- 1. cores per Stage T key ------------------------------------------
    for label, rows in (("fit_same_panel", fit_rows),
                        ("internal_same_panel", internal_rows)):
        cores_per_key: dict[str, set[str]] = defaultdict(set)
        for row in rows:
            cores_per_key[row.exact_key].add(row.core)
        counts = [len(v) for v in cores_per_key.values()]
        observations_in_multicore = sum(
            1 for row in rows if len(cores_per_key[row.exact_key]) > 1)
        report.setdefault("cores_per_stage_t_key", {})[label] = {
            "keys": len(cores_per_key),
            "keys_with_multiple_cores": sum(1 for c in counts if c > 1),
            "keys_with_multiple_cores_fraction": (
                sum(1 for c in counts if c > 1) / len(counts) if counts else 0.0),
            "cores_per_key": _stats(counts),
            "observations": len(rows),
            "observations_under_a_multicore_key": observations_in_multicore,
            "observations_under_a_multicore_key_fraction": (
                observations_in_multicore / len(rows) if rows else 0.0),
        }

    # -- 2. target effects that aggregate multiple cores -------------------
    def effect_groups(rows):
        grouped: dict[tuple[str, str], list] = defaultdict(list)
        for row in rows:
            grouped[(row.exact_key, row.target)].append(row)
        return grouped

    for label, rows in (("fit_same_panel", fit_rows),
                        ("internal_same_panel", internal_rows)):
        grouped = effect_groups(rows)
        multi = {k: v for k, v in grouped.items()
                 if len({r.core for r in v}) > 1}
        report.setdefault("target_effects", {})[label] = {
            "effects": len(grouped),
            "effects_aggregating_multiple_cores": len(multi),
            "fraction": len(multi) / len(grouped) if grouped else 0.0,
            "cores_within_an_aggregated_effect": _stats(
                len({r.core for r in v}) for v in multi.values()),
        }

    # -- 3. the decisive statistic: within-target spread across cores -------
    # The protein is held fixed, so all of this spread is generic chemical
    # context that Stage T's key claimed had cancelled.
    within_target_core_gaps: list[float] = []
    within_target_core_ranges: list[float] = []
    for row_set in (fit_rows, internal_rows):
        grouped = effect_groups(row_set)
        for (_key, _target), rows in grouped.items():
            by_core: dict[str, list[float]] = defaultdict(list)
            for row in rows:
                by_core[row.core].append(row.delta_y)
            if len(by_core) < 2:
                continue
            medians = [float(np.median(v)) for v in by_core.values()]
            within_target_core_ranges.append(max(medians) - min(medians))
            for left, right in combinations(medians, 2):
                within_target_core_gaps.append(abs(left - right))
    report["within_target_core_effect"] = {
        "definition": ("for one target and one Stage T key realised on >= 2 "
                       "cores, the spread of delta_y across cores. The protein "
                       "is identical, so this is pure generic chemical context"),
        "pairs": _stats(within_target_core_gaps),
        "ranges": _stats(within_target_core_ranges),
        "interpretation": (
            "this is the magnitude of the nuisance term that Stage T's D "
            "carried whenever two targets realised the same nominal edit on "
            "different cores; under the requested core-inclusive key it is "
            "exactly zero by construction"),
    }

    # -- 4. D rows comparing different core distributions -------------------
    def double_difference_core_audit(rows):
        grouped = effect_groups(rows)
        by_key: dict[str, list[tuple[str, set[str], float]]] = defaultdict(list)
        for (key, target), items in grouped.items():
            by_key[key].append((target, {r.core for r in items},
                                float(np.median([r.delta_y for r in items]))))
        total = mismatched = 0
        contaminating_keys = set()
        for key, entries in by_key.items():
            for left, right in combinations(sorted(entries), 2):
                total += 1
                if not (left[1] & right[1]):
                    mismatched += 1
                    contaminating_keys.add(key)
        return {
            "double_difference_rows": total,
            "rows_with_disjoint_core_sets": mismatched,
            "fraction_with_disjoint_core_sets": (
                mismatched / total if total else 0.0),
            "keys_contributing_a_disjoint_row": len(contaminating_keys),
        }

    report["double_difference_core_mismatch"] = {
        "fit_same_panel": double_difference_core_audit(fit_rows),
        "internal_same_panel": double_difference_core_audit(internal_rows),
    }

    # -- 5. what the core-inclusive key does to the census -------------------
    def census(rows, key_function) -> dict:
        keys = defaultdict(set)
        components = defaultdict(set)
        counts: Counter = Counter()
        for row in rows:
            key = key_function(row)
            keys[key].add(row.target)
            components[key].add(row.component)
            counts[key] += 1
        rich = [k for k in counts
                if len(keys[k]) >= 3 and len(components[k]) >= 3]
        return {
            "observations": len(rows),
            "keys": len(counts),
            "targets": len({r.target for r in rows}),
            "components": len({r.component for r in rows}),
            "keys_with_ge2_targets": sum(1 for v in keys.values() if len(v) >= 2),
            "keys_with_ge3_targets_and_ge3_components": len(rich),
            "observations_in_rich_keys": int(sum(counts[k] for k in rich)),
            "singleton_key_share": (
                sum(1 for v in counts.values() if v == 1) / len(counts)
                if counts else 0.0),
            "top_key_share": (max(counts.values()) / sum(counts.values())
                              if counts else 0.0),
        }

    report["census_comparison"] = {
        "fit_same_panel": {
            "stage_t_core_blind_key": census(fit_rows, lambda r: r.exact_key),
            "requested_core_inclusive_key": census(fit_rows, core_inclusive_key),
        },
        "internal_same_panel": {
            "stage_t_core_blind_key": census(internal_rows, lambda r: r.exact_key),
            "requested_core_inclusive_key": census(internal_rows,
                                                   core_inclusive_key),
        },
    }

    # Double-difference availability under each key.
    def d_pairs(rows, key_function) -> dict:
        grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows:
            grouped[(key_function(row), row.target)].append(row.delta_y)
        by_key: dict[str, list[tuple[str, float]]] = defaultdict(list)
        for (key, target), values in grouped.items():
            by_key[key].append((target, float(np.median(values))))
        total = cross = 0
        keys = set()
        components = set()
        for key, entries in by_key.items():
            if len(entries) < 2:
                continue
            for left, right in combinations(sorted(entries), 2):
                total += 1
                keys.add(key)
                lc, rc = component_of[left[0]], component_of[right[0]]
                components.update({lc, rc})
                if lc != rc:
                    cross += 1
        return {"rows": total, "cross_component_rows": cross,
                "keys": len(keys), "components": len(components),
                "effective_independent_units": min(len(components), len(keys))}

    report["double_difference_availability"] = {
        "fit_same_panel": {
            "stage_t_core_blind_key": d_pairs(fit_rows, lambda r: r.exact_key),
            "requested_core_inclusive_key": d_pairs(fit_rows, core_inclusive_key),
        },
        "internal_same_panel": {
            "stage_t_core_blind_key": d_pairs(internal_rows, lambda r: r.exact_key),
            "requested_core_inclusive_key": d_pairs(internal_rows,
                                                    core_inclusive_key),
        },
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    fit_census = report["census_comparison"]["fit_same_panel"]
    print(json.dumps({
        "output": str(args.output),
        "within_target_core_gap_median": report["within_target_core_effect"][
            "pairs"].get("median"),
        "fit_rich_keys_core_blind": fit_census["stage_t_core_blind_key"][
            "keys_with_ge3_targets_and_ge3_components"],
        "fit_rich_keys_core_inclusive": fit_census[
            "requested_core_inclusive_key"][
                "keys_with_ge3_targets_and_ge3_components"],
        "fit_D_rows_core_inclusive": report["double_difference_availability"][
            "fit_same_panel"]["requested_core_inclusive_key"],
        "internal_D_rows_core_inclusive": report[
            "double_difference_availability"]["internal_same_panel"][
                "requested_core_inclusive_key"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
