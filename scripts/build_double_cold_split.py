"""Stage R1: build a governed double-cold split. Label-blind by construction.

The existing CD-HIT40 protocol is component-hard on **proteins only**: Stage R0
showed 48.9% of its `meta_val` k=0 query cells contain a ligand that appears
verbatim in `meta_train`, and that this single fact accounted for the whole of
the Stage 10 result. This builder adds a second axis so that an evaluation cell
is cold in the protein *and* in the ligand.

Two-axis assignment
-------------------
* protein axis: `protein_group_40` (CD-HIT40) components, unchanged;
* ligand axis: Bemis-Murcko scaffold clusters. Single-linkage Tanimoto
  clustering was tried first and rejected — at a 0.4 threshold one chained
  cluster absorbs 74.6% of the corpus, which cannot support a split.

A cell is evaluation only if its component is on the evaluation side **and** its
scaffold cluster is on the evaluation side. Off-diagonal cells are discarded;
that loss is the price of the second axis and is reported.

Closure enforced, in order
--------------------------
1. protein homology component (assignment);
2. ligand scaffold cluster (assignment) — this also makes exact ligand identity
   disjoint, since a ligand belongs to exactly one cluster;
3. document/assay: an evaluation cell is dropped if any of its `panel_ids`
   documents also occurs in the training block;
4. eligibility: an evaluation target is kept only if it retains enough unique
   ligands for the nested k=5 protocol.

Chemical similarity to training ligands is then **measured**, not assumed, and
every evaluation ligand is assigned a tier. The `< 0.4` tier is the low
similarity tier required by the protocol.

Nothing here reads `pK`. The assignment is a pure function of identifiers and
structures, so it cannot be tuned toward an outcome.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.stageR0_retrieval_falsification import murcko_scaffolds, tanimoto_rows
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)
from scripts.qpsmp_data import QPSMPData

AXIS_SALT = "MetaSieve.DoubleCold.v1"
TIERS = ((0.0, 0.4, "lt40"), (0.4, 0.6, "t40_60"),
         (0.6, 0.8, "t60_80"), (0.8, 1.01, "ge80"))


def digest(*parts: object) -> int:
    raw = "|".join(map(str, parts)).encode("utf-8")
    return int(hashlib.sha256(raw).hexdigest()[:12], 16)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    digest_ = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest_.update(block)
    return digest_.hexdigest()


def documents(cell: dict) -> set[str]:
    return {str(item).split("|", 1)[0] for item in cell.get("panel_ids", [])}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eval-cluster-fraction", type=float, default=0.30)
    parser.add_argument("--dev-components", type=int, default=12)
    parser.add_argument("--test-components", type=int, default=8)
    parser.add_argument("--min-unique-ligands", type=int, default=9,
                        help="k=5 plus a query panel of at least 4")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    # The split builder is the one legitimate consumer of the whole corpus: it
    # *creates* the double-cold assignment, so it must see every cell of the
    # older main_v0 partition — including that partition's own `meta_test`,
    # which is a different, already-consumed population from the sealed
    # double-cold confirmation split this script goes on to define. The
    # assignment is label-blind (scaffolds, fingerprints, panel documents and
    # counts only; no `pK` is read), and re-running is blocked by the
    # `FileExistsError` guard above.
    data = QPSMPData(
        CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
        include_meta_test=True,
        meta_test_authorization=(
            "build_double_cold_split: label-blind reassignment of the entire "
            "main_v0 corpus into a new governed split; no pK is read and no "
            "double-cold meta_test exists yet at this point"))
    cells = data.cells
    scaffolds = murcko_scaffolds(data._ligand_smiles)
    fingerprints = data.fingerprints

    cluster_of = {key: (scaffolds.get(key) or f"__acyclic__{key}")
                  for key in {c["ligand_id"] for c in cells}}
    component_of = {c["target_id"]: c["protein_group_40"] for c in cells}
    ligands_of_target: dict[str, set[str]] = defaultdict(set)
    for cell in cells:
        ligands_of_target[cell["target_id"]].add(cell["ligand_id"])

    # ---- ligand axis: deterministic, label-blind -------------------------
    clusters = sorted(set(cluster_of.values()))
    eval_clusters = {c for c in clusters
                     if digest(AXIS_SALT, "ligand", c) % 1000
                     < args.eval_cluster_fraction * 1000}

    # ---- protein axis: rank components by how many of their targets keep a
    # usable ligand panel on the evaluation side. Uses counts only, no labels.
    eligible_targets: dict[str, list[str]] = defaultdict(list)
    for target, ligands in ligands_of_target.items():
        kept = {l for l in ligands if cluster_of[l] in eval_clusters}
        if len(kept) >= args.min_unique_ligands:
            eligible_targets[component_of[target]].append(target)
    ranked = sorted(eligible_targets,
                    key=lambda c: (-len(eligible_targets[c]),
                                   digest(AXIS_SALT, "protein", c)))
    taken = ranked[:args.dev_components + args.test_components]
    # Walk down the rank order and keep the dev:test ratio at its target, so
    # both sides receive a mix of strong and weak components rather than dev
    # taking every rich one.
    dev_components: set[str] = set()
    test_components: set[str] = set()
    for component in taken:
        if (len(dev_components) * args.test_components
                <= len(test_components) * args.dev_components):
            dev_components.add(component)
        else:
            test_components.add(component)
    eval_components = dev_components | test_components

    # ---- blocks -----------------------------------------------------------
    assignment: dict[str, str] = {}
    for cell in cells:
        component = cell["protein_group_40"]
        evaluation_ligand = cluster_of[cell["ligand_id"]] in eval_clusters
        if component not in eval_components and not evaluation_ligand:
            assignment[cell["cell_id"]] = "meta_train"
        elif component in dev_components and evaluation_ligand:
            assignment[cell["cell_id"]] = "meta_val"
        elif component in test_components and evaluation_ligand:
            assignment[cell["cell_id"]] = "meta_test"
    discarded_offdiagonal = len(cells) - len(assignment)

    # ---- document closure -------------------------------------------------
    train_documents: set[str] = set()
    for cell in cells:
        if assignment.get(cell["cell_id"]) == "meta_train":
            train_documents |= documents(cell)
    document_dropped = 0
    for cell in cells:
        split = assignment.get(cell["cell_id"])
        if split in {"meta_val", "meta_test"} and documents(cell) & train_documents:
            del assignment[cell["cell_id"]]
            document_dropped += 1

    # ---- eligibility ------------------------------------------------------
    unique_by_target: dict[tuple[str, str], set[str]] = defaultdict(set)
    for cell in cells:
        split = assignment.get(cell["cell_id"])
        if split in {"meta_val", "meta_test"}:
            unique_by_target[(split, cell["target_id"])].add(cell["ligand_id"])
    ineligible = {target for (split, target), ligands in unique_by_target.items()
                  if len(ligands) < args.min_unique_ligands}
    eligibility_dropped = 0
    for cell in cells:
        if (assignment.get(cell["cell_id"]) in {"meta_val", "meta_test"}
                and cell["target_id"] in ineligible):
            del assignment[cell["cell_id"]]
            eligibility_dropped += 1

    # ---- verification -----------------------------------------------------
    by_split: dict[str, list[dict]] = defaultdict(list)
    for cell in cells:
        split = assignment.get(cell["cell_id"])
        if split:
            by_split[split].append(cell)
    train_ligands = {c["ligand_id"] for c in by_split["meta_train"]}
    train_clusters = {cluster_of[l] for l in train_ligands}
    train_components = {c["protein_group_40"] for c in by_split["meta_train"]}
    violations = {}
    for split in ("meta_val", "meta_test"):
        ligands = {c["ligand_id"] for c in by_split[split]}
        violations[split] = {
            "exact_ligand_overlap": len(ligands & train_ligands),
            "scaffold_overlap": len({cluster_of[l] for l in ligands}
                                    & train_clusters),
            "component_overlap": len({c["protein_group_40"] for c in by_split[split]}
                                     & train_components),
            "document_overlap": len({d for c in by_split[split]
                                     for d in documents(c)} & train_documents),
        }
    cross = ({c["protein_group_40"] for c in by_split["meta_val"]}
             & {c["protein_group_40"] for c in by_split["meta_test"]})
    violations["val_test_component_overlap"] = len(cross)

    # ---- chemical similarity tiers, measured ------------------------------
    train_matrix = np.stack([fingerprints[l].numpy() for l in sorted(train_ligands)])
    tiers: dict[str, Counter] = {"meta_val": Counter(), "meta_test": Counter()}
    similarity_of: dict[str, float] = {}
    for split in ("meta_val", "meta_test"):
        ligands = sorted({c["ligand_id"] for c in by_split[split]})
        if not ligands:
            continue
        query = np.stack([fingerprints[l].numpy() for l in ligands])
        best = tanimoto_rows(query, train_matrix).max(-1)
        for ligand, value in zip(ligands, best):
            similarity_of[ligand] = float(value)
            for low, high, name in TIERS:
                if low <= value < high:
                    tiers[split][name] += 1
                    break

    # ---- statistics -------------------------------------------------------
    def describe(split: str) -> dict:
        rows = by_split[split]
        per_target: dict[str, set[str]] = defaultdict(set)
        for cell in rows:
            per_target[cell["target_id"]].add(cell["ligand_id"])
        sizes = sorted(len(v) for v in per_target.values())
        return {
            "cells": len(rows),
            "targets": len(per_target),
            "components": len({c["protein_group_40"] for c in rows}),
            "ligands": len({c["ligand_id"] for c in rows}),
            "scaffold_clusters": len({cluster_of[c["ligand_id"]] for c in rows}),
            "targets_with_9_ligands": sum(1 for n in sizes if n >= 9),
            "targets_with_25_ligands": sum(1 for n in sizes if n >= 25),
            "median_ligands_per_target": (float(np.median(sizes)) if sizes else 0.0),
        }

    statistics = {split: describe(split)
                  for split in ("meta_train", "meta_val", "meta_test")}
    serialized = json.dumps(assignment, sort_keys=True)
    manifest = {
        "schema": "MetaSieve.DoubleColdSplit.v1",
        "built_from": {
            "corpus": str(CORPUS),
            "cells_sha256": file_sha256(CORPUS / "cells.jsonl.gz"),
            "ligands_sha256": file_sha256(CORPUS / "ligands.jsonl"),
            "proteins_sha256": file_sha256(CORPUS / "proteins.jsonl"),
        },
        "label_blind": True,
        "axes": {
            "protein": "protein_group_40 (CD-HIT40)",
            "ligand": "Bemis-Murcko scaffold cluster (RDKit MurckoScaffold); "
                      "acyclic ligands form singleton clusters",
        },
        "rejected_alternative": "single-linkage Morgan Tanimoto clustering at "
                                "0.4 chains: the largest cluster absorbs 74.6% "
                                "of the 9,880 corpus ligands",
        "parameters": {
            "eval_cluster_fraction": args.eval_cluster_fraction,
            "dev_components": args.dev_components,
            "test_components": args.test_components,
            "min_unique_ligands": args.min_unique_ligands,
            "axis_salt": AXIS_SALT,
        },
        "assignment_sha256": sha256_text(serialized),
        "statistics": statistics,
        "discarded": {
            "off_diagonal_cells": discarded_offdiagonal,
            "document_closure_cells": document_dropped,
            "ineligible_target_cells": eligibility_dropped,
            "retained_fraction": len(assignment) / len(cells),
        },
        "closure_violations": violations,
        "similarity_tiers": {k: dict(v) for k, v in tiers.items()},
        "tier_definition": "max Morgan(r=2,1024) Tanimoto of the evaluation "
                           "ligand to any training-block ligand",
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "assignment.json").write_text(serialized + "\n", encoding="utf-8")
    (args.output / "ligand_similarity.json").write_text(
        json.dumps(similarity_of, indent=0, sort_keys=True) + "\n", encoding="utf-8")
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps({"statistics": statistics,
                      "discarded": manifest["discarded"],
                      "closure_violations": violations,
                      "similarity_tiers": manifest["similarity_tiers"],
                      "assignment_sha256": manifest["assignment_sha256"]},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
