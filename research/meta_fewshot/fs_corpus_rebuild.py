"""Label-blind FS-C0/FS-C1 audit for the BindingDB Ki few-shot estimand."""
from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[2]
PROJECTION = (ROOT / "dataset" / "processed" / "crossed_interaction" /
              "bindingdb_202608" / "metadata_projection.jsonl.gz")
OUT = ROOT / "report" / "meta_fewshot"
K_VALUES = (1, 2, 3, 5)
MIN_QUERY_LIGANDS = 3
MIN_LIGANDS = {k: k + MIN_QUERY_LIGANDS for k in K_VALUES}
MIN_EVAL_TARGETS = 30
MIN_SOURCE_TARGETS = 100
MIN_K5_COMPONENTS = 5
IDENTITY_THRESHOLD = 0.40
CANDIDATE_IDENTITY_THRESHOLD = 0.30
MMSEQS = ROOT / "tools" / "mmseqs2" / "mmseqs" / "bin" / "mmseqs.exe"
Z_SUM = 1.6448536269514722 + 0.8416212335729143
MAX_MDE_D = 0.600


class Components:
    def __init__(self, nodes):
        self.parent = {node: node for node in nodes}

    def find(self, node):
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def scaffold(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    core = MurckoScaffold.GetScaffoldForMol(Chem.RemoveHs(molecule))
    value = Chem.MolToSmiles(core, canonical=True, isomericSmiles=False)
    return value or None


def admitted_rows(path: Path = PROJECTION) -> tuple[list[dict], dict]:
    """Read metadata only; the projection has no affinity-value field."""
    cells, counts = {}, Counter()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            counts["projection_rows"] += 1
            if row["chain_count"] != 1 or row["endpoint_available"] != ["Ki"]:
                continue
            counts["single_chain_exact_ki_rows"] += 1
            target, ligand = row["target_sequence_sha256"], row["ligand_inchikey"]
            row_scaffold = scaffold(row["ligand_smiles"])
            if row_scaffold is None:
                counts["invalid_or_scaffoldless_rows_excluded"] += 1
                continue
            cell = cells.setdefault((target, ligand), {
                "target": target, "ligand": ligand, "sequence": row["target_sequence"],
                "scaffold": row_scaffold, "documents": set(),
                "assays": set(), "source_rows": 0,
            })
            if cell["sequence"] != row["target_sequence"]:
                raise ValueError("target hash has inconsistent sequences")
            cell["documents"].add(row["document_id"])
            cell["assays"].update(assay["assay_id"] for assay in row["assays"])
            cell["source_rows"] += 1
    counts["canonical_cells"] = len(cells)
    return list(cells.values()), dict(counts)


def local_identity(left: str, right: str) -> float:
    import parasail

    result = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
    return result.matches / result.length if result.length else 0.0


def homology_candidates(sequences: dict[str, str]) -> set[tuple[str, str]]:
    if not MMSEQS.is_file():
        raise FileNotFoundError(f"MMseqs2 executable is required: {MMSEQS}")
    with tempfile.TemporaryDirectory(prefix="metasieve_fs_homology_") as raw:
        work = Path(raw)
        fasta, output, temporary = work / "targets.fasta", work / "hits.tsv", work / "tmp"
        with fasta.open("w", encoding="ascii", newline="\n") as handle:
            for target in sorted(sequences):
                handle.write(f">{target}\n{sequences[target]}\n")
        command = [str(MMSEQS), "easy-search", str(fasta), str(fasta), str(output), str(temporary),
                   "--min-seq-id", str(CANDIDATE_IDENTITY_THRESHOLD), "-c", "0.0",
                   "--max-seqs", "100000", "-s", "7.5", "--format-output", "query,target",
                   "--threads", "8"]
        subprocess.run(command, check=True, capture_output=True, text=True)
        candidates = set()
        for line in output.read_text(encoding="utf-8").splitlines():
            left, right = line.split("\t")[:2]
            if left != right:
                candidates.add(tuple(sorted((left, right))))
        return candidates


def closure(cells: list[dict]) -> tuple[Components, int, int]:
    sequences = {cell["target"]: cell["sequence"] for cell in cells}
    components = Components(sequences)
    by_document, by_scaffold = defaultdict(set), defaultdict(set)
    for cell in cells:
        for document in cell["documents"]:
            by_document[document].add(cell["target"])
        by_scaffold[cell["scaffold"]].add(cell["target"])
    for groups in (by_document.values(), by_scaffold.values()):
        for group in groups:
            ordered = sorted(group)
            for target in ordered[1:]:
                components.union(ordered[0], target)
    candidates = homology_candidates(sequences)
    homology_edges = 0
    for left, right in sorted(candidates):
        if local_identity(sequences[left], sequences[right]) >= IDENTITY_THRESHOLD:
            components.union(left, right)
            homology_edges += 1
    return components, len(candidates), homology_edges


def mde_d(count: int) -> float:
    return float(Z_SUM / np.sqrt(count)) if count else float("inf")


def audit_cells(cells: list[dict], components: Components, homology_candidates_count: int,
                homology_edges: int) -> dict:
    by_target = defaultdict(list)
    for cell in cells:
        by_target[cell["target"]].append(cell)
    target_stats = {}
    for target, rows in by_target.items():
        ligands = {row["ligand"] for row in rows}
        scaffolds = {row["scaffold"] for row in rows}
        documents = set().union(*(row["documents"] for row in rows))
        target_stats[target] = {
            "cells": len(rows), "ligands": len(ligands), "scaffolds": len(scaffolds),
            "documents": len(documents), "component": components.find(target),
        }
    component_targets, component_cells = defaultdict(list), Counter()
    for target, stats in target_stats.items():
        component_targets[stats["component"]].append(target)
        component_cells[stats["component"]] += stats["cells"]
    eligible = {
        k: sorted(target for target, stats in target_stats.items()
                  if stats["ligands"] >= MIN_LIGANDS[k]) for k in K_VALUES
    }
    scaffold_feasible = {
        k: sum(stats["ligands"] >= MIN_LIGANDS[k] and stats["scaffolds"] >= k + 1
               for stats in target_stats.values()) for k in K_VALUES
    }
    eligible_components = {components.find(target) for target in eligible[5]}
    component_summary = sorted(({
        "component": component, "targets": len(targets),
        "cells": component_cells[component],
        "k5_eligible_targets": sum(target in eligible[5] for target in targets),
    } for component, targets in component_targets.items()),
        key=lambda row: (-row["targets"], row["component"]))
    total_targets, total_cells = len(target_stats), len(cells)
    return {
        "canonical_cells": total_cells,
        "targets": total_targets,
        "documents": len({doc for cell in cells for doc in cell["documents"]}),
        "targets_usable_at_k": {str(k): len(eligible[k]) for k in K_VALUES},
        "targets_usable_at_k_scaffold_disjoint": {str(k): scaffold_feasible[k] for k in K_VALUES},
        "targets_ge8_ligands_ge2_documents": sum(
            stats["ligands"] >= 8 and stats["documents"] >= 2 for stats in target_stats.values()),
        "dependency_components": len(component_targets),
        "components_with_k5_eligible_target": len(eligible_components),
        "largest_component_target_share": max(map(len, component_targets.values())) / total_targets,
        "largest_component_cell_share": max(component_cells.values()) / total_cells,
        "homology_candidates_at_30pct": homology_candidates_count,
        "homology_edges_at_40pct": homology_edges,
        "component_summary": component_summary,
        "target_stats": target_stats,
        "eligible_k5_by_component": {
            component: sum(target in eligible[5] for target in targets)
            for component, targets in component_targets.items()},
    }


def freeze_split(summary: dict) -> dict:
    components = summary["component_summary"]
    source_component = min(components, key=lambda row: (-row["targets"], row["component"]))["component"]
    component_k5 = summary["eligible_k5_by_component"]
    held_components = sorted(component for component in component_k5 if component != source_component)
    held_out = sum(component_k5[component] for component in held_components)
    source = component_k5[source_component]
    return {
        "rule": "largest_dependency_component_to_source_all_other_components_to_evaluation",
        "source_component": source_component,
        "evaluation_components": held_components,
        "evaluation_component_count": len(held_components),
        "evaluation_k5_eligible_targets": held_out,
        "source_k5_eligible_targets": source,
        "mde_d_target_unit": mde_d(held_out),
    }


def run(projection: Path = PROJECTION) -> dict:
    started = time.time()
    cells, admission_counts = admitted_rows(projection)
    components, candidate_count, homology_edges = closure(cells)
    summary = audit_cells(cells, components, candidate_count, homology_edges)
    split = freeze_split(summary)
    fs_c0 = summary["components_with_k5_eligible_target"] >= MIN_K5_COMPONENTS
    fs_c1 = (split["evaluation_k5_eligible_targets"] >= MIN_EVAL_TARGETS and
             split["source_k5_eligible_targets"] >= MIN_SOURCE_TARGETS and
             split["mde_d_target_unit"] <= MAX_MDE_D)
    verdict = ("FEWSHOT_CORPUS_STRUCTURALLY_AND_POWER_FEASIBLE_FROZEN_SPLIT" if fs_c0 and fs_c1
               else "FEWSHOT_CORPUS_DEPENDENCY_NOT_IDENTIFIABLE" if not fs_c0
               else "FEWSHOT_EVALUATION_POWER_NOT_IDENTIFIABLE")
    result = {
        "schema": "MetaSieve.MetaFewshot.FSCorpusRebuild.v1",
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "audit_script_sha256": sha256(Path(__file__)),
        "label_blind": True, "affinity_label_reads": 0,
        "projection_sha256": sha256(projection),
        "admission": {"chain_count": 1, "endpoint_available_exact": ["Ki"],
                      "canonical_cell": ["target_sequence_sha256", "ligand_inchikey"],
                      "counts": admission_counts},
        "dependency_closure": {"document": True, "murcko_scaffold": True,
                               "protein_local_identity_threshold": IDENTITY_THRESHOLD,
                               "candidate_identity_threshold": CANDIDATE_IDENTITY_THRESHOLD},
        "declared_thresholds": {"min_query_ligands": MIN_QUERY_LIGANDS,
                                 "min_eval_targets": MIN_EVAL_TARGETS,
                                 "min_source_targets": MIN_SOURCE_TARGETS,
                                 "min_k5_components": MIN_K5_COMPONENTS,
                                 "max_mde_d": MAX_MDE_D},
        "summary": {key: value for key, value in summary.items()
                    if key not in {"target_stats", "eligible_k5_by_component"}},
        "frozen_split": split,
        "gates": {"FS_C0_structural": fs_c0, "FS_C1_frozen_split_power": fs_c1},
        "split_frozen": fs_c0, "training_authorized": False,
        "TERMINAL_VERDICT": verdict, "elapsed_seconds": round(time.time() - started, 2),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "FS_CORPUS_REBUILD_AUDIT.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    (OUT / "FS_CORPUS_REBUILD_SPLIT.json").write_text(json.dumps({
        "schema": "MetaSieve.MetaFewshot.FSCorpusSplit.v1", "projection_sha256": sha256(projection),
        "label_blind": True, "affinity_label_reads": 0, "split": split,
        "training_authorized": False,
    }, indent=2, sort_keys=True), encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
