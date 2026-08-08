"""Run the label-blind E-AFF-X0-B crossed design re-registration.

X0-B recomputes the crossed-source design under the cell-disjoint rectangle unit
and a clustered effective sample size. It reads only the label-blind panel
geometry X0 already published, selects no affinity field, and trains nothing.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-X0B_CROSSED_DESIGN_REREGISTRATION"
REQUIRED_EFFECTIVE_N = 245
RHO_GRID = (0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.00)
CAP_GRID = (1, 2, 3, 5, 8, 12, 20, 32, 50, 80, 125, 200, 320, 500, 800, 1250, 2000, None)


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def pack_cell_disjoint(target_ligands: dict[str, set[str]]) -> list[tuple]:
    """Pack cell-disjoint rectangles: greedy target matching, then ligand pairs.

    Each target joins at most one matched pair and each matched pair splits its
    common ligands into disjoint consecutive pairs, so no measured (target,
    ligand) cell is consumed twice. This is a lower bound on the maximum packing.
    """
    available = set(target_ligands)
    packed: list[tuple] = []
    while len(available) >= 2:
        targets = sorted(available)
        best = None
        for index, left in enumerate(targets):
            for right in targets[index + 1:]:
                shared = len(target_ligands[left] & target_ligands[right])
                if shared >= 2 and (best is None or shared > best[0]):
                    best = (shared, left, right)
        if best is None:
            break
        _, left, right = best
        available.discard(left)
        available.discard(right)
        common = sorted(target_ligands[left] & target_ligands[right])
        for offset in range(0, len(common) - 1, 2):
            packed.append((left, right, common[offset], common[offset + 1]))
    return packed


def pack_target_ligand_disjoint(target_ligands: dict[str, set[str]]) -> int:
    """Conservative comparator: no target and no ligand is reused."""
    remaining = {target: set(ligands) for target, ligands in target_ligands.items()}
    packed = 0
    while True:
        targets = sorted(target for target, ligands in remaining.items() if len(ligands) >= 2)
        chosen = None
        for index, left in enumerate(targets):
            for right in targets[index + 1:]:
                common = remaining[left] & remaining[right]
                if len(common) >= 2:
                    chosen = (left, right, sorted(common)[:2])
                    break
            if chosen:
                break
        if not chosen:
            return packed
        left, right, ligands = chosen
        del remaining[left]
        del remaining[right]
        for ligand in ligands:
            for available in remaining.values():
                available.discard(ligand)
        packed += 1


def verify_cell_disjoint(packed: list[tuple]) -> bool:
    used: set[tuple[str, str]] = set()
    for left, right, first, second in packed:
        cells = {(left, first), (left, second), (right, first), (right, second)}
        if used & cells:
            return False
        used |= cells
    return True


def capped_design(sizes: list[int], cap: int | None) -> tuple[int, float]:
    limited = [min(size, cap) if cap else size for size in sizes]
    total = sum(limited)
    if total == 0:
        return 0, 0.0
    return total, sum(size * size for size in limited) / total


def effective_n(total: int, m_a: float, rho: float) -> float:
    return total / (1.0 + (m_a - 1.0) * rho)


def breakeven_rho(total: int, m_a: float) -> float:
    """Largest rho at which n_eff still reaches the requirement."""
    if total < REQUIRED_EFFECTIVE_N:
        return 0.0
    if m_a <= 1.0:
        return 1.0
    return min(1.0, (total / REQUIRED_EFFECTIVE_N - 1.0) / (m_a - 1.0))


def run(args) -> dict:
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output exists: {output}")
    x0_root = Path(args.x0)

    x0_report = json.loads((x0_root / "report.json").read_text(encoding="utf-8"))
    if x0_report.get("affinity_value_fields_selected") != 0 or x0_report.get("affinity_values_read"):
        raise RuntimeError("X0-B refuses a source whose label firewall is not clean")

    panel_targets: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    panel_endpoint: dict[str, str] = {}
    panel_cluster: dict[str, str] = {}
    for cell in _read_jsonl(x0_root / "cells.jsonl"):
        panel = cell["panel_id"]
        panel_endpoint[panel] = cell["endpoint_family"]
        panel_cluster[panel] = cell["closure_component_id"]
        panel_targets[panel][cell["protein_sequence_sha256"]].add(cell["ligand_connectivity_key"])

    endpoints = {}
    verdicts = {}
    for endpoint in ("Ki", "Kd"):
        by_cluster: dict[str, int] = defaultdict(int)
        target_pairs: set[tuple[str, str]] = set()
        used_targets: set[str] = set()
        used_ligands: set[str] = set()
        conservative = 0
        packing_verified = True
        for panel, family in panel_endpoint.items():
            if family != endpoint:
                continue
            packed = pack_cell_disjoint(panel_targets[panel])
            if packed and not verify_cell_disjoint(packed):
                packing_verified = False
            if packed:
                by_cluster[panel_cluster[panel]] += len(packed)
            for left, right, first, second in packed:
                target_pairs.add((left, right))
                used_targets.update((left, right))
                used_ligands.update((first, second))
            conservative += pack_target_ligand_disjoint(panel_targets[panel])
        if not packing_verified:
            raise RuntimeError(f"cell-disjointness violated while packing {endpoint}")

        sizes = sorted(by_cluster.values(), reverse=True)
        clusters = len(sizes)
        total = sum(sizes)

        best_by_rho = {}
        for rho in RHO_GRID:
            best = max(
                (
                    {
                        "cap": cap,
                        "units": capped[0],
                        "mean_cluster_influence": round(capped[1], 4),
                        "design_effect": round(1.0 + (capped[1] - 1.0) * rho, 4),
                        "n_eff": round(effective_n(capped[0], capped[1], rho), 2),
                    }
                    for cap in CAP_GRID
                    for capped in (capped_design(sizes, cap),)
                    if capped[0] > 0
                ),
                key=lambda entry: entry["n_eff"],
                default={"cap": None, "units": 0, "mean_cluster_influence": 0.0,
                         "design_effect": 0.0, "n_eff": 0.0},
            )
            best["meets_requirement"] = best["n_eff"] >= REQUIRED_EFFECTIVE_N
            best_by_rho[f"{rho:.2f}"] = best

        rho_star, rho_star_cap = 0.0, None
        for cap in CAP_GRID:
            capped_total, capped_m_a = capped_design(sizes, cap)
            candidate = breakeven_rho(capped_total, capped_m_a)
            if candidate > rho_star:
                rho_star, rho_star_cap = candidate, cap

        feasible = total >= REQUIRED_EFFECTIVE_N and rho_star > 0.0
        verdicts[endpoint] = (
            f"X0B_CONDITIONAL_DESIGN_SUPPORTED_{endpoint.upper()}" if feasible
            else f"X0B_DESIGN_INSUFFICIENT_{endpoint.upper()}")

        endpoints[endpoint] = {
            "cell_disjoint_units": total,
            "target_and_ligand_disjoint_units": conservative,
            "distinct_target_pairs": len(target_pairs),
            "distinct_targets": len(used_targets),
            "distinct_ligands": len(used_ligands),
            "units_per_target_pair": round(total / len(target_pairs), 2) if target_pairs else 0.0,
            "x0_effective_components": x0_report["endpoints"][endpoint]["dependency_components"],
            "clusters": clusters,
            "largest_cluster": sizes[0] if sizes else 0,
            "largest_cluster_share": round(sizes[0] / total, 4) if total else 0.0,
            "cluster_sizes": sizes,
            "mean_cluster_influence_uncapped": round(capped_design(sizes, None)[1], 4),
            "required_effective_n": REQUIRED_EFFECTIVE_N,
            "units_reach_requirement_at_rho_zero": total >= REQUIRED_EFFECTIVE_N,
            "cluster_bound_rho_max": round(clusters / REQUIRED_EFFECTIVE_N, 4),
            "breakeven_rho_star": round(rho_star, 4),
            "breakeven_rho_star_cap": rho_star_cap,
            "n_eff_by_rho": best_by_rho,
            "verdict": verdicts[endpoint],
        }

    output.mkdir(parents=True)
    report = {
        "schema": "MetaSieve.EAffX0B.v1",
        "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": "|".join(verdicts[endpoint] for endpoint in ("Ki", "Kd")),
        "label_blind": True,
        "affinity_value_fields_selected": 0,
        "affinity_values_read": False,
        "davis_label_reads": 0,
        "recipient_label_reads": 0,
        "training_performed": False,
        "unit": "cell_disjoint_rectangle_within_panel",
        "cluster": "d1_homology_document_closure_component",
        "population": "e0_core_governed_panel_geometry_published_by_x0",
        "frozen_unchanged": {
            "interaction_rms_over_assay_noise": 0.5,
            "variance_ratio": 1.25,
            "alpha_one_sided": 0.05,
            "power": 0.80,
            "required_effective_n": REQUIRED_EFFECTIVE_N,
            "affinity_gate_margins": "+0.03 correct-minus-ligand and +0.03 correct-minus-deranged",
        },
        "rho_grid": list(RHO_GRID),
        "endpoints": endpoints,
        "interpretation_limits": [
            "rho is a property of measured values and is not estimable label-blind",
            "a conditional verdict is a design statement, not evidence of interaction",
            "packing counts are auditable lower bounds, not maximum packings",
            "per-cluster caps are reported as a design option and are not applied here",
            "n_eff is bounded by clusters/rho, so cluster count is the binding resource",
        ],
    }
    _write_json(output / "report.json", report)
    manifest = {
        "stage": STAGE,
        "label_blind": True,
        "inputs": {
            "x0_cells": sha256_file(x0_root / "cells.jsonl"),
            "x0_report": sha256_file(x0_root / "report.json"),
            "x0_manifest": sha256_file(x0_root / "manifest.json"),
            "feasibility_report": sha256_file(Path(args.feasibility) / "report.json"),
            "preregistration": sha256_file(
                Path(__file__).with_name("EAFF_X0B_PREREGISTRATION.md")),
        },
        "outputs": {"report.json": sha256_file(output / "report.json")},
        "label_reads": {"affinity_values": 0, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--x0", default="research/e0_identifiability/artifacts/eaff_x0_v1")
    parser.add_argument("--feasibility",
                        default="research/e0_identifiability/artifacts/eaff_x0_feas_v1")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_x0b_v1")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
