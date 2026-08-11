"""Build the literature-aligned BindingDB exact-Ki protein-task corpus."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import random
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from research.crossed_interaction.bindingdb_cq_r0 import iter_projection, sha256_file
from research.crossed_interaction.bindingdb_cq_r1 import iter_labels

ROOT = Path(__file__).resolve().parents[2]
PROJECTION = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/metadata_projection.jsonl.gz"
LABELS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/exact_labels.jsonl.gz"
OUTPUT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
CDHIT = ROOT / "tools/cdhit/4.8.1/cd-hit.exe"
SPLIT_SEED = 20260811
K_VALUES = (1, 2, 3, 5)
MIN_QUERY = 3


def stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_ligand(smiles: str) -> tuple[str, str, float, int] | None:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    molecule = Chem.RemoveHs(molecule)
    if (molecule.GetNumAtoms() < 1 or molecule.GetNumAtoms() > 128
            or any(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms())):
        return None
    weight = float(Descriptors.MolWt(molecule))
    if not np.isfinite(weight) or weight > 1000.0:
        return None
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(molecule)
    scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=True)
    return canonical, scaffold, weight, molecule.GetNumAtoms()


def aggregate_observations(metadata: dict[str, dict], labels: list[dict]) -> tuple[list[dict], dict]:
    counts = Counter()
    panel_values: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    ligand_identity: dict[str, set[tuple[str, str, float, int]]] = defaultdict(set)
    target_identity: dict[str, set[str]] = defaultdict(set)
    for label in labels:
        if label.get("endpoint") != "Ki":
            continue
        counts["exact_ki_rows"] += 1
        record = metadata.get(label["source_row_id"])
        if record is None or record.get("chain_count") != 1:
            counts["missing_or_non_single_protein"] += 1
            continue
        target, ligand = label.get("target_id"), label.get("ligand_id")
        sequence = record.get("target_sequence", "")
        parsed = canonical_ligand(record.get("ligand_smiles", ""))
        if not target or not ligand or not sequence or parsed is None:
            counts["invalid_required_input"] += 1
            continue
        ligand_identity[ligand].add(parsed)
        target_identity[target].add(sequence)
        panel_values[(label["panel_id"], target, ligand)].append(label)

    bad_ligands = {key for key, values in ligand_identity.items() if len(values) != 1}
    bad_targets = {key for key, values in target_identity.items() if len(values) != 1}
    counts["conflicting_ligands"] = len(bad_ligands)
    counts["conflicting_targets"] = len(bad_targets)
    pair_panels: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for (panel, target, ligand), values in sorted(panel_values.items()):
        if target in bad_targets or ligand in bad_ligands:
            continue
        pair_panels[(target, ligand)].append({
            "panel_id": panel,
            "pK": float(np.median([float(value["pK"]) for value in values])),
            "replicates": len(values),
            "source_row_ids": sorted(value["source_row_id"] for value in values),
        })

    cells = []
    for (target, ligand), panels in sorted(pair_panels.items()):
        cells.append({
            "cell_id": stable_hash(f"main-v0|{target}|{ligand}"),
            "target_id": target,
            "ligand_id": ligand,
            "pK": float(np.median([panel["pK"] for panel in panels])),
            "panel_count": len(panels),
            "replicate_count": sum(panel["replicates"] for panel in panels),
            "panel_ids": [panel["panel_id"] for panel in panels],
            "source_row_ids": sorted({row for panel in panels for row in panel["source_row_ids"]}),
        })
    counts["pair_observations"] = len(cells)
    counts["cross_panel_pairs"] = sum(cell["panel_count"] > 1 for cell in cells)
    return cells, dict(counts)


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for target in sorted(sequences):
            handle.write(f">{target}\n{sequences[target]}\n")


def parse_cdhit_clusters(path: Path) -> dict[str, str]:
    clusters: list[list[str]] = []
    current: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">Cluster "):
            current = []
            clusters.append(current)
            continue
        match = re.search(r">([0-9a-f]{64})\.\.\.", line)
        if current is None or match is None:
            raise ValueError(f"invalid CD-HIT cluster line: {line}")
        current.append(match.group(1))
    result = {}
    for members in clusters:
        cluster_id = stable_hash("cdhit40|" + "|".join(sorted(members)))
        for target in members:
            if target in result:
                raise ValueError("target appears in multiple CD-HIT clusters")
            result[target] = cluster_id
    return result


def assign_clusters(target_to_cluster: dict[str, str], seed: int = SPLIT_SEED) -> dict[str, str]:
    members: dict[str, list[str]] = defaultdict(list)
    for target, cluster in target_to_cluster.items():
        members[cluster].append(target)
    rng = random.Random(seed)
    tie = {cluster: rng.random() for cluster in members}
    ordered = sorted(members, key=lambda cluster: (-len(members[cluster]), tie[cluster], cluster))
    names = ("meta_train", "meta_val", "meta_test")
    fractions = {"meta_train": 0.8, "meta_val": 0.1, "meta_test": 0.1}
    desired = {name: fractions[name] * len(target_to_cluster) for name in names}
    assigned = {name: 0 for name in names}
    cluster_split = {}
    for cluster in ordered:
        size = len(members[cluster])
        def cost(candidate: str) -> tuple[float, int]:
            future = dict(assigned)
            future[candidate] += size
            score = sum(((future[name] - desired[name]) / max(desired[name], 1.0)) ** 2
                        for name in names)
            return score, names.index(candidate)
        selected = min(names, key=cost)
        cluster_split[cluster] = selected
        assigned[selected] += size
    return {target: cluster_split[cluster] for target, cluster in target_to_cluster.items()}


def _gzip_rows(path: Path, rows: list[dict]) -> None:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    writer = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
    try:
        for row in rows:
            writer.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    finally:
        writer.close()
        raw.close()


def run(output: Path = OUTPUT) -> dict:
    if not CDHIT.is_file():
        raise FileNotFoundError(f"pinned CD-HIT binary missing: {CDHIT}")
    output.mkdir(parents=True, exist_ok=True)
    metadata = {row["source_row_id"]: row for row in iter_projection(PROJECTION)}
    cells, admission = aggregate_observations(metadata, list(iter_labels(LABELS)))
    sequences = {}
    ligand_data = {}
    for cell in cells:
        record = metadata[cell["source_row_ids"][0]]
        sequences[cell["target_id"]] = record["target_sequence"]
        parsed = canonical_ligand(record["ligand_smiles"])
        if parsed is None:
            raise AssertionError("admitted ligand became invalid")
        smiles, scaffold, weight, atoms = parsed
        ligand_data[cell["ligand_id"]] = {
            "drug_key": cell["ligand_id"], "smiles": smiles, "scaffold": scaffold,
            "molecular_weight": weight, "atoms": atoms,
        }

    by_target = Counter(cell["target_id"] for cell in cells)
    split_targets = {target for target, count in by_target.items() if count >= 1 + MIN_QUERY}
    fasta = output / "proteins.fasta"
    write_fasta(fasta, {target: sequences[target] for target in split_targets})
    cluster_prefix = output / "cdhit40"
    command = [str(CDHIT), "-i", str(fasta), "-o", str(cluster_prefix),
               "-c", "0.40", "-n", "2", "-d", "0", "-T", "1", "-M", "0"]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    target_cluster = parse_cdhit_clusters(cluster_prefix.with_suffix(".clstr"))
    if set(target_cluster) != split_targets:
        raise ValueError("CD-HIT output does not cover every eligible target")
    target_split = assign_clusters(target_cluster)

    retained = []
    for cell in cells:
        target = cell["target_id"]
        if target not in target_split:
            continue
        retained.append({
            **cell,
            "protein_group_40": target_cluster[target],
            "split": target_split[target],
            "legal_k": [k for k in K_VALUES if by_target[target] >= k + MIN_QUERY],
        })
    retained.sort(key=lambda row: (row["split"], row["target_id"], row["ligand_id"]))
    retained_targets = sorted({row["target_id"] for row in retained})
    retained_ligands = sorted({row["ligand_id"] for row in retained})
    proteins = [{"sequence_sha256": target, "sequence": sequences[target],
                 "protein_group_40": target_cluster[target], "split": target_split[target]}
                for target in retained_targets]
    ligands = [ligand_data[ligand] for ligand in retained_ligands]
    _gzip_rows(output / "cells.jsonl.gz", retained)
    (output / "proteins.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in proteins), encoding="utf-8")
    (output / "ligands.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in ligands), encoding="utf-8")

    split_summary = {}
    for split in ("meta_train", "meta_val", "meta_test"):
        targets = {target for target in retained_targets if target_split[target] == split}
        split_summary[split] = {
            "targets": len(targets),
            "clusters": len({target_cluster[target] for target in targets}),
            "cells": sum(row["target_id"] in targets for row in retained),
            "targets_usable_at_k": {
                str(k): sum(by_target[target] >= k + MIN_QUERY for target in targets)
                for k in K_VALUES
            },
        }
    manifest = {
        "schema": "MetaSieve.MetaFewshot.MainV0Corpus.v1",
        "projection_sha256": sha256_file(PROJECTION),
        "labels_sha256": sha256_file(LABELS),
        "cdhit_binary_sha256": sha256_file(CDHIT),
        "cdhit_command": command,
        "cdhit_stdout_tail": completed.stdout.splitlines()[-12:],
        "split_seed": SPLIT_SEED,
        "split_rule": "complete_cdhit40_clusters_target_count_balanced_8_1_1",
        "affinity_used_for_split": False,
        "cleaning": {
            "endpoint": "exact positive uncensored Ki",
            "transform": "pKi=9-log10(Ki[nM])",
            "max_molecular_weight": 1000.0,
            "max_atoms": 128,
            "within_panel_aggregation": "median",
            "cross_panel_pair_aggregation": "equal-panel median",
            "panel_id_role": "BindingDB assay proxy",
        },
        "admission": admission,
        "cells": len(retained), "targets": len(retained_targets),
        "ligands": len(retained_ligands), "clusters": len(set(target_cluster.values())),
        "splits": split_summary,
        "files": {},
        "training_authorized": True,
        "strict_confirmation_open": False,
    }
    for name in ("cells.jsonl.gz", "proteins.jsonl", "ligands.jsonl",
                 "proteins.fasta", "cdhit40", "cdhit40.clstr"):
        manifest["files"][name] = sha256_file(output / name)
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
