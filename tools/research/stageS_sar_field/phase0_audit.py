"""Phase 0: audit the eligible within-target pair population before any coding.

No model is built here and nothing is trained.  The audit answers, in order:

1. how many eligible within-target ligand pairs exist, per split and per
   component partition;
2. how many of them share a panel / document / assay context, and how many are
   cross-panel comparisons whose difference is partly inter-assay offset;
3. targets, components, pair counts, activity-cliff counts, ligand novelty and
   usable high-confidence pairs;
4. that the physically isolated split view is mounted and no `meta_test` label
   is present;
5. that the frozen `meta_train` fit / internal-validation component partition
   from `scripts/internal_validation.py` is used unchanged.

Run:
    python -m tools.research.stageS_sar_field.phase0_audit
"""
from __future__ import annotations

from collections import Counter, defaultdict
import argparse
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from scripts.qpsmp_data import QPSMPData
from tools.research.stageS_sar_field.pairs import (
    CLIFF_GAP, CLIFF_TANIMOTO, LOCAL_TANIMOTO, MEDIUM_TANIMOTO, PairSpec,
    build_target_pairs, component_of_target, load_data, target_cell_index,
)

HERE = Path(__file__).resolve().parent


def _quantiles(values: list[float]) -> dict[str, float]:
    if not values:
        return {"n": 0}
    array = np.asarray(values, dtype=np.float64)
    return {
        "n": int(array.size),
        "mean": float(array.mean()),
        "sd": float(array.std(ddof=1)) if array.size > 1 else 0.0,
        "p05": float(np.quantile(array, 0.05)),
        "median": float(np.quantile(array, 0.50)),
        "p95": float(np.quantile(array, 0.95)),
        "max": float(array.max()),
    }


def summarize_pairs(specs: list[PairSpec]) -> dict:
    strata = Counter(spec.stratum for spec in specs)
    same_panel = [spec for spec in specs if spec.same_panel]
    cross_panel = [spec for spec in specs if not spec.same_panel]
    same_document = [spec for spec in specs if spec.same_document]
    return {
        "pairs": len(specs),
        "targets": len({spec.target for spec in specs}),
        "components": len({spec.component for spec in specs}),
        "strata": {name: int(strata.get(name, 0))
                   for name in ("local", "medium", "cliff", "distant")},
        "same_panel": len(same_panel),
        "cross_panel": len(cross_panel),
        "same_document": len(same_document),
        "cross_document": len(specs) - len(same_document),
        "same_panel_fraction": (len(same_panel) / len(specs)) if specs else 0.0,
        "activity_cliffs": int(strata.get("cliff", 0)),
        "activity_cliffs_same_panel": sum(
            1 for spec in same_panel if spec.stratum == "cliff"),
        "abs_delta_y": _quantiles([abs(spec.delta_y) for spec in specs]),
        "abs_delta_y_same_panel": _quantiles(
            [abs(spec.delta_y) for spec in same_panel]),
        "abs_delta_y_cross_panel": _quantiles(
            [abs(spec.delta_y) for spec in cross_panel]),
        "tanimoto": _quantiles([spec.tanimoto for spec in specs]),
    }


def ligand_novelty(data: QPSMPData, reference_targets: set[str],
                   query_targets: set[str]) -> dict:
    """Chemical novelty of the query population against the reference one.

    Reports identity overlap, Murcko-scaffold overlap and the distribution of
    each query ligand's maximum Tanimoto to the reference set.  The novelty
    terciles used later as evaluation strata are cut on this quantity.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")

    def ligands_of(targets: set[str]) -> dict[str, str]:
        out: dict[str, str] = {}
        for cell in data.cells:
            if cell["target_id"] in targets:
                out.setdefault(cell["ligand_id"], data._ligand_smiles.get(
                    cell["ligand_id"]))
        return out

    reference = ligands_of(reference_targets)
    query = ligands_of(query_targets)
    table = data.fingerprints

    def scaffold(smiles: str | None) -> str | None:
        if not smiles:
            return None
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return None
        return MurckoScaffold.MurckoScaffoldSmiles(mol=molecule)

    reference_scaffolds = {scaffold(smiles) for smiles in reference.values()}
    reference_scaffolds.discard(None)
    query_scaffolds = {key: scaffold(smiles) for key, smiles in query.items()}

    reference_keys = sorted(reference)
    reference_rows = np.stack([table[key].numpy() for key in reference_keys])
    reference_counts = reference_rows.sum(axis=1)
    maxima: dict[str, float] = {}
    for key in sorted(query):
        row = table[key].numpy()
        intersection = reference_rows @ row
        union = reference_counts + row.sum() - intersection
        with np.errstate(divide="ignore", invalid="ignore"):
            value = np.where(union > 0, intersection / np.maximum(union, 1e-12), 0.0)
        maxima[key] = float(value.max()) if value.size else 0.0

    shared = set(reference).intersection(query)
    shared_scaffolds = {key for key, value in query_scaffolds.items()
                        if value is not None and value in reference_scaffolds}
    return {
        "reference_ligands": len(reference),
        "query_ligands": len(query),
        "shared_ligand_identities": len(shared),
        "shared_ligand_fraction": (len(shared) / len(query)) if query else 0.0,
        "query_ligands_with_shared_scaffold": len(shared_scaffolds),
        "shared_scaffold_fraction": (
            len(shared_scaffolds) / len(query)) if query else 0.0,
        "max_tanimoto_to_reference": _quantiles(list(maxima.values())),
        "novelty_terciles": [
            float(np.quantile(list(maxima.values()), 1 / 3)),
            float(np.quantile(list(maxima.values()), 2 / 3)),
        ] if maxima else [],
    }


def duplicate_report(data: QPSMPData, split: str) -> dict:
    seen: dict[tuple[str, str], list[float]] = defaultdict(list)
    for cell in data.cells:
        if cell["split"] == split:
            seen[(cell["target_id"], cell["ligand_id"])].append(float(cell["pK"]))
    duplicated = {key: values for key, values in seen.items() if len(values) > 1}
    spreads = [max(values) - min(values) for values in duplicated.values()]
    return {
        "target_ligand_keys": len(seen),
        "duplicated_keys": len(duplicated),
        "rows_dropped_by_first_occurrence_rule": sum(
            len(values) - 1 for values in duplicated.values()),
        "duplicate_pk_spread": _quantiles(spreads),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "PHASE0_AUDIT.json")
    args = parser.parse_args()

    data = load_data()
    seal = data.seal_record()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to audit without the physical split view")
    if "meta_test" in data.tasks or "meta_test" in data.components:
        raise SystemExit("meta_test labels are present; the seal is broken")

    fit, internal = partition_components(data)
    component_map = component_of_target(data)
    fit_targets = {target for target, component in component_map.items()
                   if component in set(fit)}
    internal_targets = {target for target, component in component_map.items()
                        if component in set(internal)}
    report: dict = {
        "schema": "MetaSieve.StageS.Phase0Audit.v1",
        "stage": "stageS_sar_field",
        "purpose": ("quantify eligible within-target ligand pairs and their "
                    "assay context before any SAR-field code is written"),
        "meta_test": seal,
        "definitions": {
            "pair": "two distinct ligands measured against ONE target",
            "same_panel": ("the two cells share at least one panel_ids entry "
                           "= same document, same endpoint, same target"),
            "same_document": "the two cells share a DOI, possibly across panels",
            "activity_cliff": (f"tanimoto >= {CLIFF_TANIMOTO} and "
                               f"|delta pK| >= {CLIFF_GAP}"),
            "strata": {
                "cliff": "activity cliff (checked first)",
                "local": f"tanimoto >= {LOCAL_TANIMOTO}",
                "medium": f"{MEDIUM_TANIMOTO} <= tanimoto < {LOCAL_TANIMOTO}",
                "distant": f"tanimoto < {MEDIUM_TANIMOTO}",
            },
        },
        "partition": {
            "module": "scripts/internal_validation.py",
            "fit_components": len(fit),
            "internal_validation_components": len(internal),
            "fit_targets": len(fit_targets),
            "internal_validation_targets": len(internal_targets),
            # The development-validation split is not opened anywhere in this
            # stage: not for training, not for selection, not for reporting,
            # and not even for a census. Only the fit / internal-validation
            # partition of meta_train is read.
            "development_validation_split_read_in_this_stage": False,
        },
        "corpus": {
            "cells_mounted": len(data.cells),
            "splits": {split: len(tasks) for split, tasks in data.tasks.items()},
            "components": {split: len(groups)
                           for split, groups in data.components.items()},
        },
        "duplicates": {"meta_train": duplicate_report(data, "meta_train")},
    }

    populations: dict[str, dict] = {}
    for name, split, targets in (
            ("meta_train_fit", "meta_train", fit_targets),
            ("meta_train_internal_validation", "meta_train", internal_targets),
            ("meta_train_all", "meta_train", None)):
        pairs_by_target = build_target_pairs(
            data, split, sorted(targets) if targets is not None else None)
        specs = [spec for group in pairs_by_target.values() for spec in group]
        summary = summarize_pairs(specs)
        index = target_cell_index(data, split)
        eligible = {target: cells for target, cells in index.items()
                    if (targets is None or target in targets) and len(cells) >= 2}
        summary["targets_with_cells"] = len(
            [t for t in index if targets is None or t in targets])
        summary["targets_eligible_for_pairs"] = len(eligible)
        summary["unique_ligands_per_target"] = _quantiles(
            [float(len(cells)) for cells in eligible.values()])
        summary["pairs_per_target"] = _quantiles(
            [float(len(group)) for group in pairs_by_target.values()])
        populations[name] = summary

    report["populations"] = populations
    report["novelty"] = {
        "internal_validation_vs_fit": ligand_novelty(
            data, fit_targets, internal_targets),
    }

    high_confidence = populations["meta_train_fit"]["same_panel"]
    report["usable_high_confidence_pairs"] = {
        "definition": "same-panel within-target pairs on the fit components",
        "count": high_confidence,
        "cross_panel_count": populations["meta_train_fit"]["cross_panel"],
        "policy": (
            "primary supervision uses same-panel pairs; cross-panel differences "
            "are trained and reported as a separate stratum and are never "
            "pooled into the high-confidence label set"),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "fit_components": len(fit),
        "internal_components": len(internal),
        "fit_pairs": populations["meta_train_fit"]["pairs"],
        "fit_same_panel": populations["meta_train_fit"]["same_panel"],
        "fit_cross_panel": populations["meta_train_fit"]["cross_panel"],
        "fit_cliffs": populations["meta_train_fit"]["activity_cliffs"],
        "internal_pairs": populations["meta_train_internal_validation"]["pairs"],
        "internal_same_panel": populations["meta_train_internal_validation"]["same_panel"],
        "internal_cliffs": populations["meta_train_internal_validation"]["activity_cliffs"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
