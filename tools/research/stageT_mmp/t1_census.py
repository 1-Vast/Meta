"""Stage T1 — true-MMP census and the bipartite evidence graph.

Answers one question: **is the transformation graph identifiable at all?** A
protein x transformation interaction can only be estimated if the same
transformation is observed in several targets from several homology components.
If it is not, no model of any size can recover the interaction, and T2 must not
be trained.

Everything here is label-blind except the descriptive statistics that are
explicitly about labels (delta_y quantiles, activity-cliff counts). No label
enters the MMP definition, the transformation keys, the deduplication rule or
any split.

Run:
    python -m tools.research.stageT_mmp.t1_census
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

from scripts.internal_validation import partition_components
from tools.research.stageT_mmp.mmp import fragment, transformation
from tools.research.stageT_mmp.observations import (
    MMPObservation, build_observations, load_governed,
)

HERE = Path(__file__).resolve().parent

# Frozen in PREREGISTRATION.md section 3.
THRESHOLDS = {
    "same_panel_fit_observations": 2000,
    "fit_targets": 50,
    "keys_with_3_targets_and_3_components": 30,
    "internal_observations": 300,
    "internal_components": 10,
}
MIN_TARGETS_PER_KEY = 3
MIN_COMPONENTS_PER_KEY = 3


def _quantiles(values) -> dict:
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


def key_table(observations: list[MMPObservation], attribute: str) -> dict:
    """Per-transformation-key target and component degree."""
    targets: dict[str, set[str]] = defaultdict(set)
    components: dict[str, set[str]] = defaultdict(set)
    counts: Counter = Counter()
    for item in observations:
        key = getattr(item, attribute)
        targets[key].add(item.target)
        components[key].add(item.component)
        counts[key] += 1
    return {
        "keys": len(counts),
        "observations": int(sum(counts.values())),
        "target_degree": _quantiles(len(v) for v in targets.values()),
        "component_degree": _quantiles(len(v) for v in components.values()),
        "keys_with_ge2_targets": sum(1 for v in targets.values() if len(v) >= 2),
        "keys_with_ge3_targets": sum(1 for v in targets.values() if len(v) >= 3),
        "keys_with_ge3_targets_and_ge3_components": sum(
            1 for key in counts
            if len(targets[key]) >= MIN_TARGETS_PER_KEY
            and len(components[key]) >= MIN_COMPONENTS_PER_KEY),
        "top_key_share": (max(counts.values()) / sum(counts.values())
                          if counts else 0.0),
        "top10_key_share": (sum(v for _, v in counts.most_common(10))
                            / sum(counts.values()) if counts else 0.0),
        "singleton_key_share": (sum(1 for v in counts.values() if v == 1)
                                / len(counts) if counts else 0.0),
    }


def connected_components(observations: list[MMPObservation], attribute: str) -> dict:
    """Connected components of the bipartite (target, key) evidence graph."""
    parent: dict[tuple[str, str], tuple[str, str]] = {}

    def find(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[left] = right

    for item in observations:
        union(("target", item.target), ("key", getattr(item, attribute)))
    groups: Counter = Counter(find(node) for node in list(parent))
    sizes = sorted(groups.values(), reverse=True)
    return {
        "nodes": len(parent),
        "components": len(groups),
        "largest_component_nodes": sizes[0] if sizes else 0,
        "largest_component_share": (sizes[0] / len(parent)) if parent else 0.0,
        "size_distribution": sizes[:10],
    }


def summarize(observations: list[MMPObservation], label: str) -> dict:
    same_panel = [o for o in observations if o.same_panel]
    cross_panel = [o for o in observations if not o.same_panel]
    return {
        "label": label,
        "observations": len(observations),
        "targets": len({o.target for o in observations}),
        "components": len({o.component for o in observations}),
        "ligands": len({lig for o in observations for lig in (o.ligand_a, o.ligand_b)}),
        "same_panel": len(same_panel),
        "cross_panel": len(cross_panel),
        "strata": dict(Counter(o.stratum for o in observations)),
        "activity_cliffs": sum(1 for o in observations if o.activity_cliff),
        "stereo_edits": sum(1 for o in observations if o.stereo_edit),
        "charge_changing": sum(1 for o in observations if o.charge_change != 0),
        "abs_delta_y": _quantiles(abs(o.delta_y) for o in observations),
        "exact_key": key_table(observations, "exact_key"),
        "coarse_key": key_table(observations, "coarse_key"),
        "graph_exact": connected_components(observations, "exact_key"),
        "graph_coarse": connected_components(observations, "coarse_key"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "T1_CENSUS.json")
    args = parser.parse_args()

    data, seal = load_governed()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to census without the physical split view")
    if "meta_test" in data.tasks:
        raise SystemExit("meta_test is present; the seal is broken")

    fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    fit_targets = sorted(t for t, c in component_of.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component_of.items()
                              if c in set(internal))

    built = build_observations(data, fit_targets + internal_targets)
    observations = built["observations"]
    fit_set, internal_set = set(fit_targets), set(internal_targets)
    fit_obs = [o for o in observations if o.target in fit_set]
    internal_obs = [o for o in observations if o.target in internal_set]
    fit_same_panel = [o for o in fit_obs if o.same_panel]
    internal_same_panel = [o for o in internal_obs if o.same_panel]

    report: dict = {
        "schema": "MetaSieve.StageT.T1Census.v1",
        "stage": "stageT_mmp",
        "meta_test": seal,
        "construction": built["construction"],
        "populations": {
            "fit": summarize(fit_obs, "fit"),
            "fit_same_panel": summarize(fit_same_panel, "fit_same_panel"),
            "internal": summarize(internal_obs, "internal"),
            "internal_same_panel": summarize(internal_same_panel,
                                             "internal_same_panel"),
        },
    }

    # -- overlap between fit and internal ----------------------------------
    fit_keys = {o.exact_key for o in fit_obs}
    internal_keys = {o.exact_key for o in internal_obs}
    fit_coarse = {o.coarse_key for o in fit_obs}
    internal_coarse = {o.coarse_key for o in internal_obs}
    fit_scaffolds = {o.core for o in fit_obs}
    internal_scaffolds = {o.core for o in internal_obs}
    report["overlap"] = {
        "exact_keys_shared": len(fit_keys & internal_keys),
        "exact_keys_internal_only": len(internal_keys - fit_keys),
        "exact_key_reuse_fraction": (
            len(fit_keys & internal_keys) / len(internal_keys)
            if internal_keys else 0.0),
        "coarse_keys_shared": len(fit_coarse & internal_coarse),
        "coarse_key_reuse_fraction": (
            len(fit_coarse & internal_coarse) / len(internal_coarse)
            if internal_coarse else 0.0),
        "cores_shared": len(fit_scaffolds & internal_scaffolds),
        "core_reuse_fraction": (
            len(fit_scaffolds & internal_scaffolds) / len(internal_scaffolds)
            if internal_scaffolds else 0.0),
    }

    # -- the frozen admission decision -------------------------------------
    exact_fit_same_panel = report["populations"]["fit_same_panel"]["exact_key"]
    measured = {
        "same_panel_fit_observations": len(fit_same_panel),
        "fit_targets": len({o.target for o in fit_same_panel}),
        "keys_with_3_targets_and_3_components": exact_fit_same_panel[
            "keys_with_ge3_targets_and_ge3_components"],
        "internal_observations": len(internal_obs),
        "internal_components": len({o.component for o in internal_obs}),
    }
    checks = {name: {"measured": measured[name], "threshold": value,
                     "pass": bool(measured[name] >= value)}
              for name, value in THRESHOLDS.items()}
    report["admission"] = {
        "thresholds_frozen_in": "PREREGISTRATION.md section 3",
        "granularity_for_thresholds": "exact key",
        "checks": checks,
        "coarse_key_reference": {
            "keys_with_3_targets_and_3_components": report["populations"][
                "fit_same_panel"]["coarse_key"][
                    "keys_with_ge3_targets_and_ge3_components"],
            "note": ("reported for coverage only; a coarse-only pass does NOT "
                     "admit T2"),
        },
        "all_pass": bool(all(item["pass"] for item in checks.values())),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "admission": {name: item for name, item in checks.items()},
        "all_pass": report["admission"]["all_pass"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
