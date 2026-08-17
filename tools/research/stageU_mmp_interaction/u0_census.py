"""Stage U0 — true-MMP census, bipartite evidence graph and admission gate.

Answers one question: **is the transformation graph identifiable at all?** A
protein x transformation interaction can only be estimated if the same
transformation is observed in several targets from several homology components.

Everything here is label-blind except the descriptive statistics that are
explicitly about labels (delta_y quantiles, activity-cliff counts). No label
enters the MMP definition, the transformation keys, the deduplication rule or
any split.

Run:
    python -m tools.research.stageU_mmp_interaction.u0_census
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
from scripts.qpsmp_data import stable_seed
from tools.research.stageU_mmp_interaction.mmp import fragment, strip_stereochemistry, transformation
from tools.research.stageU_mmp_interaction.observation_cache import (
    CACHE as OBSERVATION_CACHE, cache_sha256, load_observations,
    save_observations,
)
from tools.research.stageU_mmp_interaction.observations import (
    MMPObservation, build_observations, load_governed,
)

HERE = Path(__file__).resolve().parent

# Frozen in PREREGISTRATION.md section 2.5.
THRESHOLDS = {
    "same_panel_fit_observations": 2000,
    "fit_targets": 50,
    "keys_with_3_targets_and_3_components": 30,
    "internal_observations": 300,
    "internal_components": 10,
}
DEGREE_GATES = {
    "top1_key_share": 0.05,
    "top10_key_share": 0.20,
    "top1_target_share": 0.25,
    "top5_target_share": 0.75,
    "top1_component_share": 0.25,
    "top5_component_share": 0.75,
}
MIN_TARGETS_PER_KEY = 3
MIN_COMPONENTS_PER_KEY = 3
SUPPORT_SIZES = (1, 2, 3, 5)
QUERY_SIZE = 16
DRAWS = 1
BANK_SEED = 20260820
NOVELTY_TERCILES = (0.30136987566947937, 0.5606504082679749)


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
    targets: dict[str, set[str]] = defaultdict(set)
    components: dict[str, set[str]] = defaultdict(set)
    counts: Counter = Counter()
    for item in observations:
        key = getattr(item, attribute)
        targets[key].add(item.target)
        components[key].add(item.component)
        counts[key] += 1
    total = sum(counts.values())
    return {
        "keys": len(counts),
        "observations": int(total),
        "target_degree": _quantiles(len(v) for v in targets.values()),
        "component_degree": _quantiles(len(v) for v in components.values()),
        "keys_with_ge2_targets": sum(1 for v in targets.values() if len(v) >= 2),
        "keys_with_ge3_targets": sum(1 for v in targets.values() if len(v) >= 3),
        "keys_with_ge3_targets_and_ge3_components": sum(
            1 for key in counts
            if len(targets[key]) >= MIN_TARGETS_PER_KEY
            and len(components[key]) >= MIN_COMPONENTS_PER_KEY),
        "top_key_share": (max(counts.values()) / total if counts else 0.0),
        "top10_key_share": (sum(v for _, v in counts.most_common(10))
                            / total if counts else 0.0),
        "singleton_key_share": (sum(1 for v in counts.values() if v == 1)
                                / len(counts) if counts else 0.0),
    }


def degree_concentration(observations: list[MMPObservation]) -> dict:
    counts = Counter(o.exact_key for o in observations)
    target_counts = Counter(o.target for o in observations)
    component_counts = Counter(o.component for o in observations)
    total = len(observations) or 1
    return {
        "top1_key_share": max(counts.values(), default=0) / total,
        "top10_key_share": sum(v for _, v in counts.most_common(10)) / total,
        "top1_target_share": max(target_counts.values(), default=0) / total,
        "top5_target_share": sum(v for _, v in target_counts.most_common(5)) / total,
        "top1_component_share": max(component_counts.values(), default=0) / total,
        "top5_component_share": sum(v for _, v in component_counts.most_common(5)) / total,
    }


def connected_components(observations: list[MMPObservation], attribute: str) -> dict:
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


def incidence_spectrum(observations: list[MMPObservation],
                       attribute: str = "exact_key") -> dict:
    """Empirical effective rank of the target x transformation incidence.

    Reported only as a sufficient-richness diagnostic; it is not a proof of
    persistent excitation and the artifact says so.
    """
    targets = sorted({o.target for o in observations})
    keys = sorted({getattr(o, attribute) for o in observations})
    if not targets or not keys:
        return {"identifiable": False, "targets": 0, "keys": 0}
    t_index = {name: i for i, name in enumerate(targets)}
    k_index = {name: i for i, name in enumerate(keys)}
    matrix = np.zeros((len(targets), len(keys)), dtype=np.float64)
    for item in observations:
        matrix[t_index[item.target], k_index[getattr(item, attribute)]] += 1.0
    singular = np.linalg.svd(matrix, compute_uv=False)
    total = float(singular.sum())
    top = max(float(singular[0]), 1e-300)
    tolerance = 1e-10 * top
    rank = int(np.sum(singular > tolerance))
    p = singular / max(total, 1e-300)
    entropy_rank = float(np.exp(-np.sum(p[p > 0] * np.log(p[p > 0]))))
    stable_rank = float((total * total) / max(float((singular ** 2).sum()), 1e-300))
    return {
        "identifiable": True,
        "targets": len(targets),
        "keys": len(keys),
        "numerical_rank": rank,
        "stable_rank": stable_rank,
        "entropy_effective_rank": entropy_rank,
        "condition_number_sigma1_over_sigma_rank": (
            top / float(singular[rank - 1]) if rank > 0 else float("inf")),
        "top_singular_values_normalized": [
            float(v / top) for v in singular[:10]],
        "diagnostic_only": ("empirical sufficient-richness diagnostics; not a "
                            "proof of persistent excitation"),
    }


def summarize(observations: list[MMPObservation], label: str) -> dict:
    same_panel = [o for o in observations if o.same_panel]
    cross_panel = [o for o in observations if not o.same_panel]
    return {
        "label": label,
        "observations": len(observations),
        "targets": len({o.target for o in observations}),
        "components": len({o.component for o in observations}),
        "ligands": len({lig for o in observations
                        for lig in (o.ligand_a, o.ligand_b)}),
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
        "degree_concentration": degree_concentration(observations),
        "incidence_spectrum_exact": incidence_spectrum(observations, "exact_key"),
        "incidence_spectrum_coarse": incidence_spectrum(observations, "coarse_key"),
    }


def _coarse_relation(left, right) -> bool:
    return (strip_stereochemistry(left.core) == strip_stereochemistry(right.core)
            and left.coarse_context == right.coarse_context
            and left.r_group != right.r_group)


def _pairs_form_mmp(data, support: tuple[int, ...], query_index: int,
                    relation: str) -> bool:
    query_smiles = data._ligand_smiles.get(data.cells[query_index]["ligand_id"])
    if not query_smiles:
        return False
    query_pieces = fragment(query_smiles)
    if not query_pieces:
        return False
    query_by_core = defaultdict(list)
    for piece in query_pieces:
        query_by_core[piece.core].append(piece)
    for cell in support:
        smiles = data._ligand_smiles.get(data.cells[cell]["ligand_id"])
        if not smiles:
            continue
        pieces = fragment(smiles)
        for piece in pieces:
            if relation == "exact":
                for other in query_by_core.get(piece.core, ()):
                    if transformation(piece, other) is not None:
                        return True
            else:
                for other in query_pieces:
                    if _coarse_relation(piece, other):
                        return True
    return False


def coverage_audit(data, fit_components, internal_components) -> dict:
    banks = data.fixed_nested_episode_banks(
        "meta_train", SUPPORT_SIZES, QUERY_SIZE, DRAWS, BANK_SEED)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    fit_set = set(fit_components)
    internal_set = set(internal_components)

    fit_ligands = sorted({cell["ligand_id"] for cell in data.cells
                          if cell["protein_group_40"] in fit_set})
    table = data.fingerprints
    reference = np.stack([table[key].numpy() for key in fit_ligands])
    reference_counts = reference.sum(axis=1)
    novelty_cache: dict[str, float] = {}

    def novelty(ligand_id: str) -> float:
        if ligand_id not in novelty_cache:
            row = table[ligand_id].numpy()
            intersection = reference @ row
            union = reference_counts + row.sum() - intersection
            with np.errstate(divide="ignore", invalid="ignore"):
                value = np.where(union > 0,
                                 intersection / np.maximum(union, 1e-12), 0.0)
            novelty_cache[ligand_id] = float(value.max()) if value.size else 0.0
        return novelty_cache[ligand_id]

    report: dict = {
        "definition": ("C_k = P(at least one support-query pair forms a valid "
                       "single-cut MMP), on the frozen nested episode banks"),
        "bank": {"support_sizes": list(SUPPORT_SIZES), "query_size": QUERY_SIZE,
                 "draws": DRAWS, "seed": BANK_SEED,
                 "constructor": "QPSMPData.fixed_nested_episode_banks"},
        "query_labels_read": False,
        "coverage": {},
    }
    for size in SUPPORT_SIZES:
        exact_hits: list[float] = []
        coarse_hits: list[float] = []
        by_component: dict[str, list[float]] = defaultdict(list)
        by_novelty: dict[str, list[float]] = defaultdict(list)
        by_population: dict[str, list[float]] = defaultdict(list)
        for spec in banks[size]:
            for query_index in spec.query:
                exact = _pairs_form_mmp(data, spec.support, query_index, "exact")
                coarse = exact or _pairs_form_mmp(
                    data, spec.support, query_index, "coarse")
                exact_hits.append(float(exact))
                coarse_hits.append(float(coarse))
                by_component[spec.component].append(float(exact))
                value = novelty(data.cells[query_index]["ligand_id"])
                bucket = ("novelty_low" if value < NOVELTY_TERCILES[0]
                          else "novelty_mid" if value < NOVELTY_TERCILES[1]
                          else "novelty_high")
                by_novelty[bucket].append(float(exact))
                population = ("internal"
                              if component_of[spec.target] in internal_set
                              else "fit")
                by_population[population].append(float(exact))
        component_means = [float(np.mean(v)) for v in by_component.values() if v]
        report["coverage"][str(size)] = {
            "queries_scored": len(exact_hits),
            "C_k_exact": float(np.mean(exact_hits)) if exact_hits else 0.0,
            "C_k_coarse": float(np.mean(coarse_hits)) if coarse_hits else 0.0,
            "C_k_exact_component_equal_weight": (
                float(np.mean(component_means)) if component_means else 0.0),
            "components": len(by_component),
            "by_novelty": {name: float(np.mean(v))
                           for name, v in sorted(by_novelty.items())},
            "by_population": {name: float(np.mean(v))
                              for name, v in sorted(by_population.items())},
        }
    values = [report["coverage"][str(s)]["C_k_exact"] for s in SUPPORT_SIZES]
    report["interpretation"] = {
        "max_C_k": max(values),
        "verdict": (
            "MMP can serve as a reference-based inference mechanism for at most "
            f"{max(values):.1%} of governed queries at k<=5. Below that share it "
            "is a TRAINING signal, not a universal deployment mechanism, and no "
            "artifact may present it as one."),
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "U0_CENSUS.json")
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

    if OBSERVATION_CACHE.exists():
        observations = load_observations()
        built = {
            "observations": observations,
            "construction": {
                "machinery": ("rdkit.Chem.rdMMPA.FragmentMol (Hussain-Rea), "
                              "single cut, isomeric SMILES"),
                "exact_key": ("SHA-256 of core | full attachment context | "
                              "R_a >> R_b; includes shared core"),
                "coarse_key": ("SHA-256 of stereo-stripped core | (element, "
                               "aromatic) | stereo-stripped R_a >> R_b"),
                "targets_considered": len({o.target for o in observations}),
                "ligand_slots_considered": len({
                    lig for o in observations
                    for lig in (o.ligand_a, o.ligand_b)}),
                "ligands_without_smiles": None,
                "ligands_with_no_admissible_cut": None,
                "deduplication": ("one row per (target, core, R_a, R_b); the "
                                  "lower cell index wins, never the label"),
                "canonical_direction": ("canonical SMILES sort of the two R "
                                        "groups; delta_y = y(r_b ligand) - "
                                        "y(r_a ligand)"),
                "cliff_definition": ("tanimoto >= 0.6 and |delta_y| >= 1.0"),
                "cache": str(OBSERVATION_CACHE),
                "cache_sha256": cache_sha256(),
            },
        }
    else:
        built = build_observations(data, fit_targets + internal_targets)
        observations = built["observations"]
        save_observations(observations)
    fit_set, internal_set = set(fit_targets), set(internal_targets)
    fit_obs = [o for o in observations if o.target in fit_set]
    internal_obs = [o for o in observations if o.target in internal_set]
    fit_same_panel = [o for o in fit_obs if o.same_panel]
    internal_same_panel = [o for o in internal_obs if o.same_panel]

    report: dict = {
        "schema": "MetaSieve.StageU.U0Census.v1",
        "stage": "stageU_mmp_interaction",
        "preregistration_sha256": "fdc0a830aa92882d07b9aea50f22a4c72fc6d93f92c55a3be6bc15cd6a645c11",
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

    report["deployment_coverage"] = coverage_audit(data, fit, internal)

    exact_fit_same_panel = report["populations"]["fit_same_panel"]["exact_key"]
    measured = {
        "same_panel_fit_observations": len(fit_same_panel),
        "fit_targets": len({o.target for o in fit_same_panel}),
        "keys_with_3_targets_and_3_components": exact_fit_same_panel[
            "keys_with_ge3_targets_and_ge3_components"],
        "internal_observations": len(internal_same_panel),
        "internal_components": len({o.component for o in internal_same_panel}),
    }
    concentration = report["populations"]["fit_same_panel"]["degree_concentration"]
    checks = {name: {"measured": measured[name], "threshold": value,
                     "pass": bool(measured[name] >= value)}
              for name, value in THRESHOLDS.items()}
    for name, threshold in DEGREE_GATES.items():
        checks[name] = {"measured": concentration[name],
                        "threshold": threshold,
                        "pass": bool(concentration[name] <= threshold)}
    report["admission"] = {
        "thresholds_frozen_in": "PREREGISTRATION.md section 2.5",
        "granularity_for_thresholds": "exact key",
        "checks": checks,
        "coarse_key_reference": {
            "keys_with_3_targets_and_3_components": report["populations"][
                "fit_same_panel"]["coarse_key"][
                    "keys_with_ge3_targets_and_ge3_components"],
            "note": ("reported for coverage only; a coarse-only pass does NOT "
                     "admit U1/U2"),
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
