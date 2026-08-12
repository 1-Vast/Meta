"""Physically isolate MetaSieve v1 source optimization and meta-validation scoring."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from research.meta_fewshot.v1_source_supervision_audit import (
    LABELS,
    median_rows,
    read_gzip_jsonl,
    sha256,
    source_key,
)

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
FEATURES = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_tbasis_features.npz"
OUT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_v1_development"
SEEDS = (20260821, 20260822, 20260823, 20260824, 20260825)


def write_jsonl_gz(path: Path, rows) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    writer = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
    try:
        for row in rows:
            writer.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    finally:
        writer.close()
        raw.close()
    return sha256(path)


def stable_hash(values) -> str:
    return hashlib.sha256("|".join(sorted(values)).encode()).hexdigest()


def fingerprint_rows(cells: list[dict], ligand_smiles: dict[str, str]) -> np.ndarray:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=288)
    cache = {}
    for ligand_id in sorted({row["ligand_id"] for row in cells}):
        molecule = Chem.MolFromSmiles(ligand_smiles[ligand_id])
        if molecule is None:
            raise ValueError(f"invalid admitted ligand {ligand_id}")
        cache[ligand_id] = generator.GetFingerprintAsNumPy(molecule).astype(np.float32)
    return np.stack([cache[row["ligand_id"]] for row in cells])


def fixed_episodes(cells: list[dict], k: int = 5, draws: int = 5) -> tuple[list[dict], list[dict]]:
    by_target = defaultdict(list)
    for row in cells:
        by_target[row["target_id"]].append(row)
    episodes, truth = [], {}
    for target, rows in sorted(by_target.items()):
        if len(rows) < k + 3:
            continue
        for seed in SEEDS:
            for draw in range(draws):
                key = f"v1-dev|{seed}|{target}|{draw}"
                episode_seed = int(hashlib.sha256(key.encode()).hexdigest()[:16], 16)
                order = np.random.default_rng(episode_seed).permutation(len(rows))
                support_rows = [rows[int(index)] for index in order[:k]]
                query_rows = [rows[int(index)] for index in order[k:]]
                support_ids = [row["cell_id"] for row in support_rows]
                query_ids = [row["cell_id"] for row in query_rows]
                episodes.append({
                    "seed": seed,
                    "draw": draw,
                    "target_id": target,
                    "protein_group_40": rows[0]["protein_group_40"],
                    "support_cell_ids": support_ids,
                    "support_pK": [float(row["pK"]) for row in support_rows],
                    "query_cell_ids": query_ids,
                    "support_hash": stable_hash(support_ids),
                    "query_hash": stable_hash(query_ids),
                })
                for row in query_rows:
                    truth[row["cell_id"]] = {
                        "cell_id": row["cell_id"], "target_id": target,
                        "pK": float(row["pK"]),
                    }
    return episodes, [truth[key] for key in sorted(truth)]


def contrast_groups(source_cells: list[dict], exact_labels: Path = LABELS) -> list[dict]:
    admitted = {}
    cell_by_pair = {}
    target_group = {}
    for cell in source_cells:
        cell_by_pair[(cell["target_id"], cell["ligand_id"])] = cell["cell_id"]
        target_group[cell["target_id"]] = cell["protein_group_40"]
        for row_id in cell["source_row_ids"]:
            admitted[source_key(row_id, cell["target_id"], cell["ligand_id"])] = "meta_train"
    rows = []
    for row in read_gzip_jsonl(exact_labels):
        if row["endpoint"] != "Ki":
            continue
        if source_key(row["source_row_id"], row["target_id"], row["ligand_id"]) in admitted:
            rows.append({**row, "split": "meta_train"})
    medians = median_rows(rows)
    within, partner = defaultdict(list), defaultdict(list)
    for row in medians:
        member = {
            "cell_id": cell_by_pair[(row["target_id"], row["ligand_id"])],
            "target_id": row["target_id"], "ligand_id": row["ligand_id"],
            "protein_group_40": target_group[row["target_id"]], "pK": row["pK"],
        }
        within[(row["panel_id"], row["target_id"])].append(member)
        partner[(row["panel_id"], row["ligand_id"])].append(member)
    result = []
    for (panel_id, target_id), members in sorted(within.items()):
        if len({row["ligand_id"] for row in members}) >= 2:
            result.append({"kind": "within_panel_ligand", "panel_id": panel_id,
                           "target_id": target_id, "members": members})
    for (panel_id, ligand_id), members in sorted(partner.items()):
        if len({row["protein_group_40"] for row in members}) >= 2:
            result.append({"kind": "measured_partner", "panel_id": panel_id,
                           "ligand_id": ligand_id, "members": members})
    return result


def seal(corpus: Path = CORPUS, features: Path = FEATURES, output: Path = OUT) -> dict:
    cells = list(read_gzip_jsonl(corpus / "cells.jsonl.gz"))
    source = [row for row in cells if row["split"] == "meta_train"]
    validation = [row for row in cells if row["split"] == "meta_val"]
    forbidden = [row for row in cells if row["split"] == "meta_test"]
    forbidden_targets = {row["target_id"] for row in forbidden}
    forbidden_cells = {row["cell_id"] for row in forbidden}
    if ({row["target_id"] for row in source + validation} & forbidden_targets or
            {row["cell_id"] for row in source + validation} & forbidden_cells):
        raise ValueError("v1 development overlaps main-v0 meta-test")

    ligand_smiles = {
        row["drug_key"]: row["smiles"]
        for row in map(json.loads, (corpus / "ligands.jsonl").read_text().splitlines())
    }
    with np.load(features, allow_pickle=False) as stored:
        cell_ids = stored["cell_id"].tolist()
        if cell_ids != [row["cell_id"] for row in cells]:
            raise ValueError("feature and corpus order differ")
        position = {cell_id: index for index, cell_id in enumerate(cell_ids)}
        source_index = np.asarray([position[row["cell_id"]] for row in source])
        val_index = np.asarray([position[row["cell_id"]] for row in validation])
        source_correct = stored["correct"][source_index].astype(np.float32)
        val_correct = stored["correct"][val_index].astype(np.float32)

    output.mkdir(parents=True, exist_ok=True)
    source_cells_hash = write_jsonl_gz(output / "source_cells.jsonl.gz", source)
    val_public = [{key: value for key, value in row.items() if key != "pK"} for row in validation]
    val_cells_hash = write_jsonl_gz(output / "metaval_cells_without_labels.jsonl.gz", val_public)
    episodes, truth = fixed_episodes(validation)
    episodes_hash = write_jsonl_gz(output / "metaval_episodes.jsonl.gz", episodes)
    truth_hash = write_jsonl_gz(output / "metaval_query_truth.jsonl.gz", truth)
    groups = contrast_groups(source)
    groups_hash = write_jsonl_gz(output / "source_contrast_groups.jsonl.gz", groups)

    np.savez_compressed(output / "source_features.npz",
                        cell_id=np.asarray([row["cell_id"] for row in source]),
                        correct=source_correct,
                        ligand=fingerprint_rows(source, ligand_smiles))
    np.savez_compressed(output / "metaval_features.npz",
                        cell_id=np.asarray([row["cell_id"] for row in validation]),
                        correct=val_correct,
                        ligand=fingerprint_rows(validation, ligand_smiles))
    blacklist = {
        "targets": sorted(forbidden_targets), "cells": sorted(forbidden_cells),
        "target_hash": stable_hash(forbidden_targets), "cell_hash": stable_hash(forbidden_cells),
    }
    (output / "main_v0_test_blacklist.json").write_text(
        json.dumps(blacklist, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "schema": "MetaSieve.V1DevelopmentSeal.v1",
        "corpus_manifest_sha256": sha256(corpus / "manifest.json"),
        "input_features_sha256": sha256(features),
        "source_cells": len(source), "source_targets": len({row["target_id"] for row in source}),
        "metaval_cells": len(validation),
        "metaval_targets": len({row["target_id"] for row in validation}),
        "main_v0_test_cells_excluded": len(forbidden_cells),
        "main_v0_test_targets_excluded": len(forbidden_targets),
        "main_v0_test_values_used": 0,
        "metaval_episode_rows": len(episodes), "metaval_query_truth_rows": len(truth),
        "source_contrast_groups": len(groups),
        "files": {
            "source_cells.jsonl.gz": source_cells_hash,
            "metaval_cells_without_labels.jsonl.gz": val_cells_hash,
            "metaval_episodes.jsonl.gz": episodes_hash,
            "metaval_query_truth.jsonl.gz": truth_hash,
            "source_contrast_groups.jsonl.gz": groups_hash,
            "source_features.npz": sha256(output / "source_features.npz"),
            "metaval_features.npz": sha256(output / "metaval_features.npz"),
            "main_v0_test_blacklist.json": sha256(output / "main_v0_test_blacklist.json"),
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    print(json.dumps(seal(), indent=2, sort_keys=True))

