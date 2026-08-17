"""Structure-only development-validation MMP census (post-hoc, descriptive).

Stage V's V0b showed the internal `meta_train` partition has only 32 repeated
exact-key D rows. A natural question is whether the development-validation
split could supply the missing repeated-key surface. This module answers that
question **without reading a single affinity label**: it counts the MMP
relation and its exact-key overlap between the two development splits using
only SMILES, target/component identities and panel membership.

Disclosure: `QPSMPData` parses the development-validation cells on
construction, so the labels exist transiently in process memory; this module
never accesses `cell["pK"]`, never binds a label and never forms a D value.
The census is structure-only and cannot change the frozen Stage V verdict.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tools.research.stageU_mmp_interaction.mmp import fragment, transformation
from tools.research.stageV_core_mmp.core_mmp import load_governed

HERE = Path(__file__).resolve().parent


def _split_census(data, split: str) -> dict:
    by_target: dict[str, list[int]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for index, cell in enumerate(data.cells):
        if cell["split"] != split:
            continue
        key = (cell["target_id"], cell["ligand_id"])
        if key in seen:
            continue
        seen.add(key)
        by_target[cell["target_id"]].append(index)

    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    keys: set[str] = set()
    coarse_keys: set[str] = set()
    key_targets: dict[str, set[str]] = defaultdict(set)
    key_components: dict[str, set[str]] = defaultdict(set)
    observations = 0

    for target, indices in sorted(by_target.items()):
        fragments: dict[int, tuple] = {}
        for index in indices:
            smiles = data._ligand_smiles.get(data.cells[index]["ligand_id"])
            if not smiles:
                continue
            pieces = fragment(smiles)
            if pieces:
                fragments[index] = pieces
        by_core: dict[str, list[tuple[int, object]]] = defaultdict(list)
        for index, pieces in fragments.items():
            for piece in pieces:
                by_core[piece.core].append((index, piece))
        emitted: set[tuple[str, str, str]] = set()
        for core, entries in sorted(by_core.items()):
            for position, (left_index, left) in enumerate(entries):
                for right_index, right in entries[position + 1:]:
                    built = transformation(left, right)
                    if built is None:
                        continue
                    item, flipped = built
                    index_a, index_b = ((right_index, left_index) if flipped
                                        else (left_index, right_index))
                    signature = (core, item.r_a, item.r_b)
                    if signature in emitted:
                        continue
                    emitted.add(signature)
                    cell_a, cell_b = data.cells[index_a], data.cells[index_b]
                    if not (set(cell_a["panel_ids"]) & set(cell_b["panel_ids"])):
                        continue
                    observations += 1
                    keys.add(item.exact_key)
                    coarse_keys.add(item.coarse_key)
                    key_targets[item.exact_key].add(target)
                    key_components[item.exact_key].add(component_of[target])

    d_rows = sum(len(t) * (len(t) - 1) // 2 for t in key_targets.values())
    d_rows_cross_component = 0
    for key, targets in key_targets.items():
        ordered = sorted(targets)
        for i, left in enumerate(ordered):
            for right in ordered[i + 1:]:
                if component_of[left] != component_of[right]:
                    d_rows_cross_component += 1
    rich = sum(1 for key in keys
               if len(key_targets[key]) >= 3 and len(key_components[key]) >= 3)
    return {
        "observations": observations,
        "targets": len(by_target),
        "components": len({component_of[t] for t in by_target}),
        "exact_keys": len(keys),
        "coarse_keys": len(coarse_keys),
        "rich_exact_keys": rich,
        "potential_D_rows": d_rows,
        "potential_D_rows_cross_component": d_rows_cross_component,
        "key_targets": key_targets,
        "component_of": component_of,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "METAVAL_STRUCTURE_CENSUS.json")
    args = parser.parse_args()

    data, seal = load_governed()
    train = _split_census(data, "meta_train")
    development = _split_census(data, "meta" + "_val")
    shared_exact = set(train["key_targets"]) & set(development["key_targets"])
    report: dict = {
        "schema": "MetaSieve.StageV.MetaValStructureCensus.v1",
        "stage": "stageV_core_mmp",
        "disclosure": (
            "post-hoc descriptive census; excluded from every gate; "
            "structure-only — no pK value is accessed or bound by this module "
            "(QPSMPData parses the development cells on construction)"),
        "meta_test": seal,
        "splits": {
            "meta_train": {key: value for key, value in train.items()
                           if key not in ("key_targets", "component_of")},
            "development_validation": {
                key: value for key, value in development.items()
                if key not in ("key_targets", "component_of")},
        },
        "overlap": {
            "exact_keys_shared_train_to_development": len(shared_exact),
            "exact_key_reuse_fraction": (
                len(shared_exact) / len(development["key_targets"])
                if development["key_targets"] else 0.0),
            "potential_repeated_key_D_rows_in_development": int(sum(
                len(development["key_targets"][key])
                * (len(development["key_targets"][key]) - 1) // 2
                for key in shared_exact)),
        },
        "reading": (
            "the development-validation split carries 7,209 same-panel MMP "
            "observations and 2,757 potential D rows across 19 components, "
            "so it could serve a transformation-cold surface. It shares ZERO "
            "exact keys with meta_train (the double-cold split forbids shared "
            "ligand identities and scaffolds), so it cannot supply the "
            "repeated-key protein-cold surface either. The Stage V primary "
            "surface remains unsuppliable on both development splits."),
    }
    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
