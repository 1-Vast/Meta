"""Compile canonical holo-complex records into immutable geometry sidecars."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from contracts.ligand_graph import MAX_ATOMS
from contracts.mechanism import (CONTACT_THRESHOLD_ANGSTROM,
    DISTANCE_BINS_ANGSTROM, HOLO_COMPLEX_SCHEMA, MECHANISM_RESIDUE_SLOTS,
    STRUCTURE_SUPERVISION_SCHEMA)
from scripts.build_holo_complex_index import (_atom_rows, _canonical_altloc_rows,
    _ccd_molecule, _protein_sequence_mapping)
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


def _sha256_json(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _ligand_coordinates(record: dict, rows: list[dict[str, str]]) -> tuple[np.ndarray, list[str]]:
    chemistry = _ccd_molecule(Path(record["ccd_path"]))
    canonical_names = chemistry["heavy_atom_names"]
    selected = [row for row in rows if row["group_PDB"] == "HETATM" and
                row["label_asym_id"] == record["ligand_asym_id"] and
                row["label_comp_id"].upper() == record["ligand_comp_id"].upper() and
                row["label_seq_id"] == record["ligand_label_seq_id"] and
                row["type_symbol"].upper() != "H"]
    selected = _canonical_altloc_rows(
        selected, ("label_asym_id", "label_seq_id", "label_atom_id"))
    by_name = {row["label_atom_id"]: row for row in selected}
    if len(by_name) != len(selected) or set(by_name) != set(canonical_names):
        raise ValueError("canonical ligand mapping is not one-to-one and complete")
    coordinates = np.asarray([[float(by_name[name][axis]) for axis in
                               ("Cartn_x", "Cartn_y", "Cartn_z")]
                              for name in canonical_names], dtype=np.float32)
    return coordinates, canonical_names


def _protein_residue_atoms(record: dict, rows: list[dict[str, str]]) -> tuple[list, list[str], list[int]]:
    selected = [row for row in rows if row["group_PDB"] == "ATOM" and
                row["label_asym_id"] == record["protein_asym_id"] and
                row["label_seq_id"] not in {"", ".", "?"}]
    selected = _canonical_altloc_rows(
        selected, ("label_asym_id", "label_seq_id", "label_atom_id"))
    residues: dict[str, list[list[float]]] = {}
    for row in selected:
        residue_id = row["label_seq_id"]
        residues.setdefault(residue_id, []).append([
            float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])])
    if not residues:
        raise ValueError("no resolved protein coordinates for canonical chain")
    mapped_labels, sequence_indices, coverage = _protein_sequence_mapping(
        selected, record["sequence"])
    if coverage < 0.9:
        raise ValueError("protein mapping coverage fell below canonical QC")
    arrays = [np.asarray(residues[value], dtype=np.float32) for value in mapped_labels]
    return arrays, mapped_labels, sequence_indices


def _geometry(record: dict, threshold: float,
              distance_bins: tuple[float, ...]) -> tuple[dict, dict]:
    required = {"schema", "source_entry_id", "pdb_id", "protein_asym_id",
                "ligand_asym_id", "ligand_label_seq_id", "ligand_comp_id",
                "sequence", "structure_path", "ccd_path", "structure_sha256",
                "ccd_sha256"}
    missing = required - set(record)
    if missing:
        raise ValueError(f"holo record missing fields: {sorted(missing)}")
    if record["schema"] != HOLO_COMPLEX_SCHEMA:
        raise ValueError(f"unsupported holo record schema: {record['schema']!r}")
    structure_path, ccd_path = Path(record["structure_path"]), Path(record["ccd_path"])
    if sha256_file(structure_path) != record["structure_sha256"]:
        raise ValueError("mmCIF hash differs from canonical holo record")
    if sha256_file(ccd_path) != record["ccd_sha256"]:
        raise ValueError("CCD hash differs from canonical holo record")

    import gemmi
    block = gemmi.cif.read(str(structure_path)).sole_block()
    rows = _atom_rows(block)
    ligand, ligand_names = _ligand_coordinates(record, rows)
    residues, label_seq_ids, sequence_indices = _protein_residue_atoms(record, rows)
    sequence_length = len(record["sequence"])
    if len(ligand) > MAX_ATOMS or len(residues) > sequence_length * 1.1:
        raise ValueError("canonical structure exceeds declared dimensions")

    slot_atoms: list[list[np.ndarray]] = [[] for _ in range(MECHANISM_RESIDUE_SLOTS)]
    for sequence_index, atoms in zip(sequence_indices, residues):
        slot = min(MECHANISM_RESIDUE_SLOTS - 1,
                   sequence_index * MECHANISM_RESIDUE_SLOTS // sequence_length)
        slot_atoms[slot].append(atoms)
    atom_mask = np.zeros(MAX_ATOMS, dtype=np.uint8)
    residue_mask = np.asarray([bool(values) for values in slot_atoms], dtype=np.uint8)
    distances = np.full((MAX_ATOMS, MECHANISM_RESIDUE_SLOTS), np.inf, dtype=np.float32)
    atom_mask[:len(ligand)] = 1
    for slot, values in enumerate(slot_atoms):
        if values:
            protein_atoms = np.concatenate(values, axis=0)
            delta = ligand[:, None, :] - protein_atoms[None, :, :]
            distances[:len(ligand), slot] = np.sqrt(np.square(delta).sum(-1)).min(-1)
    pair_mask = atom_mask[:, None] * residue_mask[None, :]
    contact = ((distances <= threshold) * pair_mask).astype(np.uint8)
    distance_bin = np.digitize(distances, distance_bins[1:-1]).astype(np.uint8)
    distance_bin[pair_mask == 0] = np.iinfo(np.uint8).max
    mapping = {"ligand_atom_names": ligand_names, "protein_label_seq_ids": label_seq_ids,
               "sequence_length": sequence_length,
               "residue_slots": MECHANISM_RESIDUE_SLOTS,
               "slot_policy": "compact_resolved_label_seq_id_linear_bins"}
    metadata = {
        "source_entry_id": record["source_entry_id"], "pdb_id": record["pdb_id"],
        "protein_asym_id": record["protein_asym_id"],
        "ligand_asym_id": record["ligand_asym_id"],
        "ligand_comp_id": record["ligand_comp_id"],
        "atom_mapping_hash": _sha256_json(ligand_names),
        "residue_mapping_hash": _sha256_json(mapping),
        "contact_definition": f"distance<={threshold:g}A",
        "distance_bins": list(distance_bins),
        "structure_sha256": record["structure_sha256"],
        "ccd_sha256": record["ccd_sha256"],
        "source_split": record.get("source_split", "unassigned_homology_group"),
    }
    arrays = {"contact": contact, "distance_bin": distance_bin,
              "atom_mask": atom_mask, "residue_mask": residue_mask}
    return metadata, arrays


def _geometry_worker(payload: tuple[dict, float, tuple[float, ...]]) -> tuple[dict, dict]:
    return _geometry(*payload)


def build_structure_supervision(records_path: str | Path, output_dir: str | Path, *,
                                shard_size: int = 128,
                                workers: int = 1,
                                contact_threshold: float = CONTACT_THRESHOLD_ANGSTROM,
                                distance_bins=DISTANCE_BINS_ANGSTROM) -> dict:
    records = read_jsonl(records_path)
    if not records:
        raise ValueError("structure supervision requires canonical holo records")
    bins = tuple(float(value) for value in distance_bins)
    if shard_size < 1 or workers < 1 or contact_threshold <= 0 or sorted(bins) != list(bins):
        raise ValueError("invalid structure supervision configuration")
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"structure supervision output already exists: {output}")
    if workers == 1:
        built = [_geometry(record, contact_threshold, bins) for record in records]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            built = list(pool.map(
                _geometry_worker,
                ((record, contact_threshold, bins) for record in records),
                chunksize=16))
    index_manifest_path = Path(records_path).parent / "manifest.json"
    index_manifest = (json.loads(index_manifest_path.read_text(encoding="utf-8"))
                      if index_manifest_path.is_file() else {})
    output.mkdir(parents=True)
    metadata, shards = [], []
    for start in range(0, len(built), shard_size):
        batch = built[start:start + shard_size]
        filename = f"shard_{start:06d}.npz"
        np.savez_compressed(output / filename,
            contact=np.stack([value[1]["contact"] for value in batch]),
            distance_bin=np.stack([value[1]["distance_bin"] for value in batch]),
            atom_mask=np.stack([value[1]["atom_mask"] for value in batch]),
            residue_mask=np.stack([value[1]["residue_mask"] for value in batch]))
        for offset, (item, _) in enumerate(batch):
            metadata.append(item | {"shard": filename, "shard_index": offset})
        shards.append({"path": filename, "sha256": sha256_file(output / filename),
                       "records": len(batch)})
    write_jsonl(output / "pairs.jsonl", metadata)
    manifest = {
        "schema": STRUCTURE_SUPERVISION_SCHEMA, "pairs": len(metadata),
        "record_schema": HOLO_COMPLEX_SCHEMA,
        "residue_slots": MECHANISM_RESIDUE_SLOTS, "max_atoms": MAX_ATOMS,
        "slot_policy": "compact_resolved_label_seq_id_linear_bins",
        "contact_threshold_angstrom": contact_threshold,
        "distance_bins_angstrom": list(bins),
        "contact_dtype": "uint8", "distance_bin_dtype": "uint8",
        "invalid_distance_bin": 255,
        "coordinate_source": index_manifest.get("coordinate_source"),
        "coordinate_license": index_manifest.get("coordinate_license"),
        "biological_relevance_source": index_manifest.get("biological_relevance_source"),
        "pdb_release_cutoff": index_manifest.get("pdb_release_cutoff"),
        "pdb_snapshot_utc": index_manifest.get("pdb_snapshot_utc"),
        "quality_filter_sha256": index_manifest.get("quality_filter_sha256"),
        "homology_governance_manifest_sha256": index_manifest.get(
            "homology_governance_manifest_sha256"),
        "ligand_overlap_audit_sha256": index_manifest.get("ligand_overlap_audit_sha256"),
        "ccd_snapshot_sha256": index_manifest.get("ccd_snapshot_sha256"),
        "records_sha256": sha256_file(records_path), "shards": shards,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def audit_plinder_mlsb(directory: str | Path) -> dict:
    """Verify whether local MLSB inputs include holo ligand coordinates."""
    root = Path(directory)
    pdbs = list(root.glob("*/input.pdb"))
    ligand_coordinate_files = list(root.glob("*/ligand.sdf")) + list(root.glob("*/ligand.mol2"))
    return {"schema": "MetaSieve.StructureSourceAudit.v1", "source": "PLINDER-MLSB",
            "systems": len(pdbs), "ligand_coordinate_files": len(ligand_coordinate_files),
            "admissible_for_geometry_supervision": bool(pdbs and ligand_coordinate_files),
            "reason": None if ligand_coordinate_files else
                "MLSB contains receptor input PDBs and ligand SMILES, not holo ligand coordinates"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records", nargs="?")
    parser.add_argument("output", nargs="?")
    parser.add_argument("--shard-size", type=int, default=128)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--audit-plinder-mlsb")
    args = parser.parse_args()
    if args.audit_plinder_mlsb:
        result = audit_plinder_mlsb(args.audit_plinder_mlsb)
        print(json.dumps(result, indent=2))
        return 0 if result["admissible_for_geometry_supervision"] else 2
    if not args.records or not args.output:
        parser.error("records and output are required when not auditing PLINDER MLSB")
    print(json.dumps(build_structure_supervision(
        args.records, args.output, shard_size=args.shard_size,
        workers=args.workers), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
