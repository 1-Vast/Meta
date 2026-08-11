"""Build immutable exact residue-by-ligand-atom geometry for R0-B.

This stage reads structure coordinates but no affinity labels and fits no
parameters. It preserves the frozen CCD canonical atom order and full sequence
residue identity required by the R0-B preregistration.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from contracts.mechanism import DISTANCE_BINS_ANGSTROM, MECHANISM_RESIDUE_SLOTS
from research.e0_identifiability.run_tdir_pilot import _coordinate_bundle
from scripts.build_structure_supervision import _sha256_json
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


FROZEN_PANEL_SHA256 = (
    "a1f3d29a3b5d876f81a23819f82ac1fa07d681ea500ddb1dbc0f81cc95d89a65")
FROZEN_LIGAND_BANK_SHA256 = (
    "823815c2e436d28403b14bc79fa5e86b99b511fdfda50b104995ad64abbe7012")
FROZEN_SPLIT_COUNTS = {"train": 2516, "val": 185, "heldout_a": 144}
R0C_PANEL_SHA256 = (
    "fb7232ae70c974fc5c4e6ff8ef5517cd5a982cebf10cbcb80e651e8e043b48b9")
R0C_LIGAND_BANK_SHA256 = (
    "56883ded880903e984d620445d093c19e0892f5ddb22263a91cab070add6ec3e")
R0C_SPLIT_COUNTS = {"heldout_b": 219}


def _validate_frozen_panel(records: list[dict], panel_path: Path,
                           ligand_bank_path: Path, contract: str = "r0b") -> None:
    if contract == "r0b":
        panel_sha = FROZEN_PANEL_SHA256
        ligand_sha = FROZEN_LIGAND_BANK_SHA256
        split_counts = FROZEN_SPLIT_COUNTS
        component_field = "homology_group_id"
    elif contract == "r0c":
        panel_sha = R0C_PANEL_SHA256
        ligand_sha = R0C_LIGAND_BANK_SHA256
        split_counts = R0C_SPLIT_COUNTS
        component_field = "r0c_final_component_id"
    else:
        raise ValueError(f"unknown exact-geometry contract: {contract}")
    if sha256_file(panel_path) != panel_sha:
        raise ValueError(f"panel is not the frozen {contract.upper()} input")
    if sha256_file(ligand_bank_path) != ligand_sha:
        raise ValueError(f"ligand bank is not the frozen {contract.upper()} graph input")
    if len(records) != sum(split_counts.values()):
        raise ValueError(f"frozen {contract.upper()} panel record count differs")
    if len({record["source_entry_id"] for record in records}) != len(records):
        raise ValueError("frozen R0-B panel has duplicate source entries")
    if Counter(record["r0_split"] for record in records) != split_counts:
        raise ValueError(f"frozen {contract.upper()} split counts differ")
    if {record["schema"] for record in records} != {"MetaSieve.HoloComplex.v1"}:
        raise ValueError("frozen R0-B panel schema differs")

    dependency_fields = (
        component_field, "ccd_sha256", "connectivity_sha256",
        "murcko_scaffold")
    for field in dependency_fields:
        split_of: dict[str, str] = {}
        for record in records:
            value = record[field]
            if field == "murcko_scaffold" and not value:
                continue
            prior = split_of.setdefault(value, record["r0_split"])
            if prior != record["r0_split"]:
                raise ValueError(f"{field} dependency straddles {contract.upper()} splits")
    if contract == "r0c" and len({
            record[component_field] for record in records}) != len(records):
        raise ValueError("R0-C final components must be one-to-one with records")


def _validate_geometry_arrays(metadata: dict,
                              arrays: dict[str, np.ndarray]) -> None:
    atoms, residues = int(metadata["atoms"]), int(metadata["residues"])
    distance = np.asarray(arrays["distance_angstrom"])
    distance_bin = np.asarray(arrays["distance_bin"])
    slots = np.asarray(arrays["slot_of_residue"])
    if distance.shape != (atoms, residues) or distance_bin.shape != distance.shape:
        raise ValueError("exact distance arrays disagree with atom/residue dimensions")
    if slots.shape != (residues,):
        raise ValueError("exact residue slot array has the wrong shape")
    if not np.isfinite(distance).all() or (distance < 0).any() or \
            float(distance.max(initial=0.0)) >= float(DISTANCE_BINS_ANGSTROM[-1]):
        raise ValueError("exact distances violate the frozen bin contract")
    expected_bin = np.digitize(distance, DISTANCE_BINS_ANGSTROM[1:-1])
    if distance_bin.dtype.kind not in "iu" or not np.array_equal(
            distance_bin, expected_bin):
        raise ValueError("stored distance bins disagree with exact distances")
    if slots.dtype.kind not in "iu" or (slots >= MECHANISM_RESIDUE_SLOTS).any():
        raise ValueError("exact residue slots violate the frozen slot contract")


def _exact_geometry(record: dict) -> tuple[dict, dict[str, np.ndarray]]:
    structure = Path(record["structure_path"])
    ccd = Path(record["ccd_path"])
    if sha256_file(structure) != record["structure_sha256"]:
        raise ValueError("mmCIF hash differs from governed record")
    if sha256_file(ccd) != record["ccd_sha256"]:
        raise ValueError("CCD hash differs from governed record")
    if float(record["protein_mapping_coverage"]) < 0.999999:
        raise ValueError("R0-B exact geometry requires complete protein mapping")

    bundle = _coordinate_bundle(record)
    ligand_names = list(bundle["chemistry"]["heavy_atom_names"])
    coordinate_names = [row["label_atom_id"] for row in bundle["ligand_rows"]]
    if len(coordinate_names) != len(set(coordinate_names)) or \
            len(ligand_names) != len(set(ligand_names)) or \
            coordinate_names != ligand_names:
        raise ValueError("ligand coordinate-to-CCD atom mapping is not one-to-one")
    ligand = np.asarray([
        [float(row[axis]) for axis in ("Cartn_x", "Cartn_y", "Cartn_z")]
        for row in bundle["ligand_rows"]
    ], dtype=np.float32)
    if len(ligand) != int(record["ligand_heavy_atoms"]):
        raise ValueError("exact ligand atom count differs from governed record")

    sequence_length = len(record["sequence"])
    residue_atoms: list[list[list[float]]] = [
        [] for _ in range(sequence_length)]
    for row in bundle["protein_rows"]:
        index = int(bundle["label_to_sequence"][row["label_seq_id"]])
        residue_atoms[index].append([
            float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])
        ])
    if any(not values for values in residue_atoms):
        raise ValueError("one or more exact sequence residues lack heavy-atom coordinates")

    distances = np.empty((len(ligand), sequence_length), dtype=np.float32)
    for residue_index, values in enumerate(residue_atoms):
        protein = np.asarray(values, dtype=np.float32)
        difference = ligand[:, None, :] - protein[None, :, :]
        distances[:, residue_index] = np.sqrt(
            np.square(difference).sum(axis=-1)).min(axis=1)
    if (not np.isfinite(distances).all() or (distances < 0).any()
            or float(distances.max()) >= float(DISTANCE_BINS_ANGSTROM[-1])):
        raise ValueError("exact distances violate the frozen bin contract")
    distance_bin = np.digitize(
        distances, DISTANCE_BINS_ANGSTROM[1:-1]).astype(np.uint8)
    slot_of_residue = np.floor_divide(
        np.arange(sequence_length, dtype=np.int64) * MECHANISM_RESIDUE_SLOTS,
        sequence_length).astype(np.uint8)
    metadata = {
        "source_entry_id": record["source_entry_id"],
        "r0_split": record["r0_split"],
        "homology_group_id": record.get(
            "homology_group_id", record.get("r0c_final_component_id")),
        "sequence_sha256": record["sequence_sha256"],
        "ccd_sha256": record["ccd_sha256"],
        "ligand_comp_id": record["ligand_comp_id"],
        "structure_sha256": record["structure_sha256"],
        "atom_mapping_hash": _sha256_json(ligand_names),
        "atoms": len(ligand),
        "residues": sequence_length,
    }
    arrays = {
        "distance_angstrom": distances,
        "distance_bin": distance_bin,
        "slot_of_residue": slot_of_residue,
    }
    _validate_geometry_arrays(metadata, arrays)
    return metadata, arrays


def _geometry_worker(record: dict) -> tuple[str, object, object]:
    try:
        metadata, arrays = _exact_geometry(record)
        return "ok", metadata, arrays
    except Exception as error:  # quarantined and surfaced in the fail-closed audit
        return "error", {
            "source_entry_id": record.get("source_entry_id"),
            "r0_split": record.get("r0_split"),
            "error_type": type(error).__name__,
            "message": str(error),
        }, None


def pack_exact_geometry_rows(
        rows: list[tuple[dict, dict[str, np.ndarray]]]) -> dict[str, np.ndarray]:
    if not rows:
        raise ValueError("cannot pack an empty exact-geometry shard")
    for metadata, arrays in rows:
        _validate_geometry_arrays(metadata, arrays)
    pair_lengths = np.asarray([
        value[1]["distance_bin"].size for value in rows], dtype=np.int64)
    residue_lengths = np.asarray([
        len(value[1]["slot_of_residue"]) for value in rows], dtype=np.int32)
    pair_offsets = np.concatenate(([0], np.cumsum(pair_lengths, dtype=np.int64)))
    residue_offsets = np.concatenate((
        [0], np.cumsum(residue_lengths, dtype=np.int64)))
    return {
        "entry_ids": np.asarray([value[0]["source_entry_id"] for value in rows]),
        "atoms": np.asarray([value[0]["atoms"] for value in rows], dtype=np.int16),
        "residues": residue_lengths,
        "pair_offsets": pair_offsets,
        "residue_offsets": residue_offsets,
        "distance_angstrom": np.concatenate([
            value[1]["distance_angstrom"].reshape(-1) for value in rows
        ]).astype(np.float32, copy=False),
        "distance_bin": np.concatenate([
            value[1]["distance_bin"].reshape(-1) for value in rows
        ]).astype(np.uint8, copy=False),
        "slot_of_residue": np.concatenate([
            value[1]["slot_of_residue"] for value in rows
        ]).astype(np.uint8, copy=False),
    }


def unpack_exact_geometry_row(
        payload: dict[str, np.ndarray], index: int) -> dict[str, np.ndarray]:
    if index < 0 or index >= len(payload["entry_ids"]):
        raise IndexError("exact-geometry row is outside the shard")
    atom_count = int(payload["atoms"][index])
    residue_count = int(payload["residues"][index])
    pair_left, pair_right = (
        int(payload["pair_offsets"][index]),
        int(payload["pair_offsets"][index + 1]))
    residue_left, residue_right = (
        int(payload["residue_offsets"][index]),
        int(payload["residue_offsets"][index + 1]))
    if pair_right - pair_left != atom_count * residue_count:
        raise ValueError("exact-geometry pair offsets disagree with dimensions")
    if residue_right - residue_left != residue_count:
        raise ValueError("exact-geometry residue offsets disagree with dimensions")
    return {
        "distance_angstrom": payload["distance_angstrom"][
            pair_left:pair_right].reshape(atom_count, residue_count).copy(),
        "distance_bin": payload["distance_bin"][
            pair_left:pair_right].reshape(atom_count, residue_count).copy(),
        "slot_of_residue": payload["slot_of_residue"][
            residue_left:residue_right].copy(),
    }


def build_exact_geometry_bank(
        panel_path: str | Path, ligand_bank_path: str | Path,
        output_dir: str | Path, *, shard_size: int = 64,
        workers: int = 1, contract: str = "r0b") -> dict:
    if min(shard_size, workers) < 1:
        raise ValueError("shard size and workers must be positive")
    panel_path = Path(panel_path)
    ligand_bank_path = Path(ligand_bank_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"exact geometry output already exists: {output}")
    records = read_jsonl(panel_path)
    if not records:
        raise ValueError("R0-B exact geometry received an empty panel")
    _validate_frozen_panel(records, panel_path, ligand_bank_path, contract)
    ligand_bank = torch.load(
        ligand_bank_path, map_location="cpu", weights_only=False)
    graph_contract = {
        key: {
            "atom_mapping_hash": value["atom_mapping_hash"],
            "ligand_comp_id": value["ligand_comp_id"],
            "atoms": int(value["mask"].sum()),
        }
        for key, value in ligand_bank.items()
    }
    output.mkdir(parents=True)
    if workers == 1:
        built = [_geometry_worker(record) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            built = list(pool.map(_geometry_worker, records, chunksize=8))

    exclusions = []
    accepted: list[tuple[dict, dict[str, np.ndarray]]] = []
    for status, metadata, arrays in built:
        if status == "error":
            exclusions.append(metadata)
            continue
        graph = graph_contract.get(metadata["ccd_sha256"])
        reasons = []
        if graph is None:
            reasons.append("missing_frozen_ligand_graph")
        else:
            if graph["atom_mapping_hash"] != metadata["atom_mapping_hash"]:
                reasons.append("atom_mapping_hash_mismatch")
            if graph["ligand_comp_id"] != metadata["ligand_comp_id"]:
                reasons.append("ligand_component_mismatch")
            if graph["atoms"] != metadata["atoms"]:
                reasons.append("ligand_atom_count_mismatch")
        if reasons:
            exclusions.append({
                "source_entry_id": metadata["source_entry_id"],
                "r0_split": metadata["r0_split"],
                "error_type": "FrozenGraphContractError",
                "message": "+".join(reasons),
            })
        else:
            accepted.append((metadata, arrays))

    index_rows = []
    shards = []
    bin_counts = np.zeros(len(DISTANCE_BINS_ANGSTROM) - 1, dtype=np.int64)
    split_counts = Counter()
    split_components: dict[str, set[str]] = defaultdict(set)
    movable, residues_total, pair_total = 0, 0, 0
    for start in range(0, len(accepted), shard_size):
        shard_rows = accepted[start:start + shard_size]
        filename = f"shard_{start:06d}.npz"
        path = output / filename
        payload = pack_exact_geometry_rows(shard_rows)
        np.savez_compressed(path, **payload)
        for row_index, (metadata, arrays) in enumerate(shard_rows):
            counts = np.bincount(
                arrays["slot_of_residue"], minlength=MECHANISM_RESIDUE_SLOTS)
            movable += int((counts[arrays["slot_of_residue"]] >= 2).sum())
            residues_total += metadata["residues"]
            pair_total += metadata["atoms"] * metadata["residues"]
            bin_counts += np.bincount(
                arrays["distance_bin"].reshape(-1),
                minlength=len(bin_counts))[:len(bin_counts)]
            split_counts[metadata["r0_split"]] += 1
            split_components[metadata["r0_split"]].add(
                metadata["homology_group_id"])
            index_rows.append({
                **metadata,
                "shard": filename,
                "row": row_index,
            })
        shards.append({
            "path": filename,
            "records": len(shard_rows),
            "sha256": sha256_file(path),
        })
    write_jsonl(output / "index.jsonl", index_rows)
    write_jsonl(output / "exclusions.jsonl", exclusions)
    heldout_split = "heldout_a" if contract == "r0b" else "heldout_b"
    heldout_sizes = Counter(
        row[0]["homology_group_id"] for row in accepted
        if row[0]["r0_split"] == heldout_split)
    largest = max(heldout_sizes.values(), default=0)
    heldout_records = int(split_counts[heldout_split])
    manifest = {
        "schema": ("MetaSieve.R0BExactGeometryBank.v1" if contract == "r0b"
                   else "MetaSieve.R0CExactGeometryBank.v1"),
        "contract": contract,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "panel_path": str(panel_path.resolve()),
        "panel_sha256": sha256_file(panel_path),
        "ligand_bank_path": str(ligand_bank_path.resolve()),
        "ligand_bank_sha256": sha256_file(ligand_bank_path),
        "selected_records": len(records),
        "built_records": len(accepted),
        "excluded_records": len(exclusions),
        "split_records": dict(split_counts),
        "split_components": {
            key: len(value) for key, value in split_components.items()},
        "heldout_largest_component_share": (
            largest / heldout_records if heldout_records else 1.0),
        "pair_cells": pair_total,
        "distance_bins_angstrom": list(DISTANCE_BINS_ANGSTROM),
        "distance_bin_counts": bin_counts.tolist(),
        "movable_residues": movable,
        "exact_residues": residues_total,
        "movable_residue_fraction": (
            movable / residues_total if residues_total else 0.0),
        "index_sha256": sha256_file(output / "index.jsonl"),
        "exclusions_sha256": sha256_file(output / "exclusions.jsonl"),
        "shards": shards,
        "affinity_labels_used": False,
        "trainable_parameters": 0,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    return manifest


def verify_exact_geometry_bank(output_dir: str | Path,
                               contract: str = "r0b") -> dict[str, int]:
    """Verify an already written bank without reading affinity labels."""
    output = Path(output_dir)
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    expected_schema = ("MetaSieve.R0BExactGeometryBank.v1" if contract == "r0b"
                       else "MetaSieve.R0CExactGeometryBank.v1")
    expected_panel = FROZEN_PANEL_SHA256 if contract == "r0b" else R0C_PANEL_SHA256
    expected_ligand = (FROZEN_LIGAND_BANK_SHA256 if contract == "r0b"
                       else R0C_LIGAND_BANK_SHA256)
    if manifest.get("schema") != expected_schema:
        raise ValueError("exact geometry manifest schema differs")
    if manifest.get("panel_sha256") != expected_panel or \
            manifest.get("ligand_bank_sha256") != expected_ligand:
        raise ValueError("exact geometry manifest does not reference frozen inputs")
    index_rows = read_jsonl(output / "index.jsonl")
    if sha256_file(output / "index.jsonl") != manifest["index_sha256"]:
        raise ValueError("exact geometry index hash differs")
    if sha256_file(output / "exclusions.jsonl") != manifest["exclusions_sha256"]:
        raise ValueError("exact geometry exclusions hash differs")
    if len(index_rows) != manifest["built_records"] or \
            len({row["source_entry_id"] for row in index_rows}) != len(index_rows):
        raise ValueError("exact geometry index count or uniqueness differs")

    by_shard: dict[str, list[dict]] = defaultdict(list)
    for row in index_rows:
        by_shard[row["shard"]].append(row)
    verified_pairs = 0
    verified_records = 0
    for shard in manifest["shards"]:
        path = output / shard["path"]
        if sha256_file(path) != shard["sha256"]:
            raise ValueError(f"exact geometry shard hash differs: {path}")
        rows = sorted(by_shard[shard["path"]], key=lambda row: row["row"])
        if len(rows) != shard["records"] or \
                [row["row"] for row in rows] != list(range(len(rows))):
            raise ValueError("exact geometry shard index rows are incomplete")
        with np.load(path, allow_pickle=False) as payload:
            if payload["entry_ids"].tolist() != [
                    row["source_entry_id"] for row in rows]:
                raise ValueError("exact geometry shard entry order differs")
            for index, metadata in enumerate(rows):
                arrays = unpack_exact_geometry_row(payload, index)
                _validate_geometry_arrays(metadata, arrays)
                verified_pairs += int(metadata["atoms"]) * int(metadata["residues"])
                verified_records += 1
    if verified_records != manifest["built_records"] or \
            verified_pairs != manifest["pair_cells"]:
        raise ValueError("exact geometry verified totals differ from manifest")
    return {"records": verified_records, "pair_cells": verified_pairs,
            "shards": len(manifest["shards"])}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--shard-size", type=int, default=64)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--contract", choices=("r0b", "r0c"), default="r0b")
    args = parser.parse_args()
    if args.verify_only:
        if len(args.paths) != 1:
            parser.error("--verify-only requires exactly one output path")
        result = verify_exact_geometry_bank(args.paths[0], contract=args.contract)
    else:
        if len(args.paths) != 3:
            parser.error("build requires panel, ligand_bank, and output paths")
        result = build_exact_geometry_bank(
            args.paths[0], args.paths[1], args.paths[2],
            shard_size=args.shard_size, workers=args.workers,
            contract=args.contract)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
