"""Compile BioLiP2 annotations and RCSB PDB/CCD files into canonical holo records."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from contracts.mechanism import HOLO_COMPLEX_SCHEMA
from scripts.data_contract import write_jsonl
from scripts.structure_sources.biolip import (BioLiPEntry, pilot_candidates,
    regular_ligand_ids)
from scripts.structure_sources.rcsb import sha256_file


ATOM_SITE_TAGS = (
    "_atom_site.group_PDB", "_atom_site.label_atom_id", "_atom_site.label_comp_id",
    "_atom_site.label_asym_id", "_atom_site.label_seq_id",
    "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
    "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z",
    "_atom_site.type_symbol", "_atom_site.label_alt_id", "_atom_site.occupancy",
    "_atom_site.pdbx_PDB_model_num",
)

QUALITY_FILTER = {
    "experimental_method": "X-RAY DIFFRACTION",
    "resolution_angstrom_max": 3.0,
    "protein_length": [50, 1022],
    "ligand_heavy_atoms": [6, 96],
    "regular_ligand": "HETATM instance with label_seq_id='.' and non-peptide CCD type",
    "ligand_mapping_coverage": 1.0,
    "protein_mapping_coverage_min": 0.9,
    "covalent_complexes": "exclude",
}


def record_set_sha256(records: list[dict]) -> str:
    digest = hashlib.sha256()
    for record in records:
        line = json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n"
        digest.update(line.encode("utf-8"))
    return digest.hexdigest()


def _values(block, tag: str) -> list[str]:
    return [str(value) for value in block.find_values(tag)]


def _ccd_molecule(path: Path) -> dict:
    import gemmi
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    block = gemmi.cif.read(str(path)).sole_block()
    component_type = str(block.find_value("_chem_comp.type")).strip("'\"").upper()
    atom_names = _values(block, "_chem_comp_atom.atom_id")
    elements = _values(block, "_chem_comp_atom.type_symbol")
    charges = _values(block, "_chem_comp_atom.charge")
    if not atom_names or not (len(atom_names) == len(elements) == len(charges)):
        raise ValueError("incomplete CCD atom table")
    if len(set(atom_names)) != len(atom_names):
        raise ValueError("duplicate CCD atom names")

    molecule = Chem.RWMol()
    name_to_index: dict[str, int] = {}
    for name, element, charge in zip(atom_names, elements, charges):
        atom = Chem.Atom(element.title())
        if charge not in {"", ".", "?"}:
            atom.SetFormalCharge(int(charge))
        atom.SetProp("_CCDAtomName", name)
        name_to_index[name] = molecule.AddAtom(atom)
    bond_types = {
        "SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
        "TRIP": Chem.BondType.TRIPLE, "QUAD": Chem.BondType.QUADRUPLE,
        "AROM": Chem.BondType.AROMATIC, "DELO": Chem.BondType.AROMATIC,
    }
    left = _values(block, "_chem_comp_bond.atom_id_1")
    right = _values(block, "_chem_comp_bond.atom_id_2")
    orders = _values(block, "_chem_comp_bond.value_order")
    if not (len(left) == len(right) == len(orders)):
        raise ValueError("incomplete CCD bond table")
    for atom_1, atom_2, order in zip(left, right, orders):
        bond_type = bond_types.get(order.upper())
        if bond_type is None:
            raise ValueError(f"unsupported CCD bond order {order!r}")
        molecule.AddBond(name_to_index[atom_1], name_to_index[atom_2], bond_type)
        if bond_type == Chem.BondType.AROMATIC:
            bond = molecule.GetBondBetweenAtoms(name_to_index[atom_1], name_to_index[atom_2])
            bond.SetIsAromatic(True)
            molecule.GetAtomWithIdx(name_to_index[atom_1]).SetIsAromatic(True)
            molecule.GetAtomWithIdx(name_to_index[atom_2]).SetIsAromatic(True)
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    heavy = [atom for atom in result.GetAtoms() if atom.GetAtomicNum() > 1]
    ranks = list(Chem.CanonicalRankAtoms(result, breakTies=True))
    heavy_indices = sorted((atom.GetIdx() for atom in heavy), key=lambda index: ranks[index])
    heavy_names = [result.GetAtomWithIdx(index).GetProp("_CCDAtomName")
                   for index in heavy_indices]
    smiles = Chem.MolToSmiles(Chem.RemoveHs(result), canonical=True, isomericSmiles=True)
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(Chem.RemoveHs(result))
    scaffold = Chem.MolToSmiles(scaffold_mol, canonical=True, isomericSmiles=False)
    return {"component_type": component_type, "molecule": result,
            "heavy_atom_names": heavy_names, "heavy_atoms": len(heavy_names),
            "canonical_smiles": smiles, "murcko_scaffold": scaffold}


def _atom_rows(block) -> list[dict[str, str]]:
    table = block.find(list(ATOM_SITE_TAGS))
    return [{tag.rsplit(".", 1)[1]: str(value) for tag, value in zip(ATOM_SITE_TAGS, row)}
            for row in table]


def _altloc_ok(value: str) -> bool:
    return value in {"", ".", "?", "A"}


def _canonical_altloc_rows(rows: list[dict[str, str]], key_fields: tuple[str, ...]) -> list[dict[str, str]]:
    """Choose one coordinate per atom identity, preferring blank then occupancy."""
    groups: dict[tuple[str, ...], list[dict[str, str]]] = {}
    for row in rows:
        if row["pdbx_PDB_model_num"] not in {"", ".", "?", "1"}:
            continue
        groups.setdefault(tuple(row[field] for field in key_fields), []).append(row)
    selected = []
    for values in groups.values():
        values.sort(key=lambda row: (
            row["label_alt_id"] not in {"", ".", "?"},
            -float(row["occupancy"] if row["occupancy"] not in {"", ".", "?"} else 0),
            row["label_alt_id"],
        ))
        selected.append(values[0])
    return selected


def _protein_sequence_mapping(protein_rows: list[dict[str, str]],
                              sequence: str) -> tuple[list[str], list[int], float]:
    import gemmi
    import parasail

    residue_names: dict[str, str] = {}
    for row in protein_rows:
        residue_names.setdefault(row["label_seq_id"], row["label_comp_id"])
    try:
        label_seq_ids = sorted(residue_names, key=lambda value: int(value))
    except ValueError as error:
        raise ValueError("non-integer_protein_label_seq_id") from error
    structure_sequence = ""
    for label_seq_id in label_seq_ids:
        code = gemmi.find_tabulated_residue(residue_names[label_seq_id]).one_letter_code
        structure_sequence += code if code and code != " " else "X"
    alignment = parasail.nw_trace_striped_16(
        structure_sequence, sequence, 10, 1, parasail.blosum62).traceback
    structure_index = sequence_index = 0
    mapped_labels, mapped_indices = [], []
    for structure_code, sequence_code in zip(alignment.query, alignment.ref):
        if structure_code != "-" and sequence_code != "-":
            mapped_labels.append(label_seq_ids[structure_index])
            mapped_indices.append(sequence_index)
        if structure_code != "-":
            structure_index += 1
        if sequence_code != "-":
            sequence_index += 1
    coverage = len(mapped_labels) / max(len(label_seq_ids), len(sequence))
    return mapped_labels, mapped_indices, coverage


def _has_covalent_link(block, ligand_asym_id: str, ligand_comp_id: str,
                       ligand_seq_id: str) -> bool:
    tags = ["_struct_conn.conn_type_id"]
    for side in ("ptnr1", "ptnr2"):
        tags.extend((f"_struct_conn.{side}_label_asym_id",
                     f"_struct_conn.{side}_label_comp_id",
                     f"_struct_conn.{side}_label_seq_id"))
    table = block.find(tags)
    for row in table:
        values = [str(value) for value in row]
        if not values[0].lower().startswith("covale"):
            continue
        for offset in (1, 4):
            if (values[offset] == ligand_asym_id and
                    values[offset + 1].upper() == ligand_comp_id.upper() and
                    values[offset + 2] == ligand_seq_id):
                return True
    return False


def _compile_entry(entry: BioLiPEntry, mmcif_path: Path, ccd_path: Path) -> dict:
    import gemmi

    chemistry = _ccd_molecule(ccd_path)
    component_type = chemistry["component_type"]
    if "PEPTIDE" in component_type:
        raise ValueError("ligand_not_regular_small_molecule")
    if not 6 <= chemistry["heavy_atoms"] <= 96:
        raise ValueError("ligand_heavy_atom_count")

    block = gemmi.cif.read(str(mmcif_path)).sole_block()
    method = str(block.find_value("_exptl.method")).strip("'\"").upper()
    if method != "X-RAY DIFFRACTION":
        raise ValueError("not_xray")
    resolution_text = str(block.find_value("_refine.ls_d_res_high"))
    resolution = float(resolution_text)
    if resolution > 3.0:
        raise ValueError("resolution_above_3A")
    rows = _atom_rows(block)

    protein_rows = [row for row in rows if row["group_PDB"] == "ATOM" and
                    row["auth_asym_id"] == entry.receptor_auth_asym_id and
                    row["label_seq_id"] not in {"", ".", "?"}]
    protein_rows = _canonical_altloc_rows(
        protein_rows, ("label_asym_id", "label_seq_id", "label_atom_id"))
    if not protein_rows:
        raise ValueError("protein_chain_not_found")
    protein_asym_ids = {row["label_asym_id"] for row in protein_rows}
    if len(protein_asym_ids) != 1:
        raise ValueError("ambiguous_protein_label_asym_id")
    _, _, mapping_coverage = _protein_sequence_mapping(protein_rows, entry.sequence)
    if mapping_coverage < 0.9:
        raise ValueError("protein_mapping_coverage")

    ligand_rows = [row for row in rows if row["group_PDB"] == "HETATM" and
                   row["label_comp_id"].upper() == entry.ligand_comp_id and
                   row["auth_asym_id"] == entry.ligand_auth_asym_id and
                   row["auth_seq_id"] == entry.ligand_auth_seq_id and
                   row["type_symbol"].upper() != "H"]
    ligand_rows = _canonical_altloc_rows(
        ligand_rows, ("label_asym_id", "label_seq_id", "label_atom_id"))
    if not ligand_rows:
        raise ValueError("ligand_instance_not_found")
    ligand_asym_ids = {row["label_asym_id"] for row in ligand_rows}
    ligand_seq_ids = {row["label_seq_id"] for row in ligand_rows}
    if len(ligand_asym_ids) != 1 or len(ligand_seq_ids) != 1:
        raise ValueError("ambiguous_ligand_instance")
    if ligand_seq_ids != {"."}:
        raise ValueError("ligand_is_polymer_not_regular")
    coordinate_names = [row["label_atom_id"] for row in ligand_rows]
    if len(coordinate_names) != len(set(coordinate_names)):
        raise ValueError("duplicate_ligand_atom_names_or_altloc")
    required_names = chemistry["heavy_atom_names"]
    if set(coordinate_names) != set(required_names):
        raise ValueError("ligand_mapping_not_100_percent")
    ligand_asym_id = next(iter(ligand_asym_ids))
    ligand_seq_id = next(iter(ligand_seq_ids))
    if _has_covalent_link(block, ligand_asym_id, entry.ligand_comp_id, ligand_seq_id):
        raise ValueError("covalent_complex")

    return {
        "schema": HOLO_COMPLEX_SCHEMA,
        "source": {"annotation": "BioLiP2", "coordinates": "RCSB_PDB",
                   "chemistry": "RCSB_CCD", "coordinate_license": "CC0-1.0"},
        "source_entry_id": entry.source_entry_id,
        "pdb_id": entry.pdb_id,
        "protein_asym_id": next(iter(protein_asym_ids)),
        "protein_auth_asym_id": entry.receptor_auth_asym_id,
        "ligand_asym_id": ligand_asym_id,
        "ligand_auth_asym_id": entry.ligand_auth_asym_id,
        "ligand_label_seq_id": ligand_seq_id,
        "ligand_auth_seq_id": entry.ligand_auth_seq_id,
        "ligand_comp_id": entry.ligand_comp_id,
        "sequence": entry.sequence,
        "sequence_sha256": hashlib.sha256(entry.sequence.encode("ascii")).hexdigest(),
        "canonical_smiles": chemistry["canonical_smiles"],
        "connectivity_sha256": hashlib.sha256(
            chemistry["canonical_smiles"].encode("utf-8")).hexdigest(),
        "murcko_scaffold": chemistry["murcko_scaffold"],
        "structure_path": str(mmcif_path.resolve()),
        "ccd_path": str(ccd_path.resolve()),
        "experimental_method": method,
        "resolution_angstrom": resolution,
        "protein_mapping_coverage": mapping_coverage,
        "ligand_mapping_coverage": 1.0,
        "ligand_heavy_atoms": chemistry["heavy_atoms"],
        "structure_sha256": sha256_file(mmcif_path),
        "ccd_sha256": sha256_file(ccd_path),
        "biologically_relevant": True,
        "noncovalent": True,
    }


def _compile_candidate_worker(payload: tuple[BioLiPEntry, str]) -> tuple[BioLiPEntry, dict | None, str | None]:
    entry, root_value = payload
    root = Path(root_value)
    mmcif = root / "mmcif" / f"{entry.pdb_id}.cif.gz"
    ccd = root / "ccd" / f"{entry.ligand_comp_id}.cif"
    if not mmcif.is_file():
        return entry, None, "missing_mmcif"
    if not ccd.is_file():
        return entry, None, "missing_ccd"
    try:
        return entry, _compile_entry(entry, mmcif, ccd), None
    except Exception as error:
        return entry, None, str(error)


def build_holo_complex_index(annotation_path: str | Path, raw_root: str | Path,
                             output_dir: str | Path, *, candidate_limit: int,
                             governance_manifest: str | Path | None = None,
                             workers: int = 1) -> dict:
    root, output = Path(raw_root), Path(output_dir)
    if output.exists():
        raise FileExistsError(f"holo index output already exists: {output}")
    ligand_summary = root / "biolip2" / "ligand.tsv.gz"
    if not ligand_summary.is_file():
        raise FileNotFoundError(f"BioLiP ligand summary not found: {ligand_summary}")
    candidates = pilot_candidates(annotation_path, limit=candidate_limit,
                                  allowed_ligands=regular_ligand_ids(ligand_summary))
    acquisition_path = root / "acquisition_manifest.json"
    acquisition = (json.loads(acquisition_path.read_text(encoding="utf-8"))
                   if acquisition_path.is_file() else None)
    provenance_pass = bool(
        acquisition and acquisition.get("source") == "BioLiP2+RCSB_PDB_CCD" and
        acquisition.get("coordinate_license") == "CC0-1.0" and
        acquisition.get("files_failed") == 0)
    accepted: list[dict] = []
    rejected: Counter[str] = Counter()
    rejection_examples: dict[str, list[str]] = {}
    if workers < 1:
        raise ValueError("workers must be positive")

    if workers == 1:
        compiled = map(_compile_candidate_worker,
                       ((entry, str(root)) for entry in candidates))
        pool = None
    else:
        pool = ProcessPoolExecutor(max_workers=workers)
        compiled = pool.map(_compile_candidate_worker,
                            ((entry, str(root)) for entry in candidates), chunksize=16)
    try:
        for entry, record, reason in compiled:
            if record is not None:
                accepted.append(record)
                continue
            assert reason is not None
            rejected[reason] += 1
            examples = rejection_examples.setdefault(reason, [])
            if len(examples) < 10:
                examples.append(entry.source_entry_id)
    finally:
        if pool is not None:
            pool.shutdown()

    ungoverned_record_sha256 = record_set_sha256(accepted)
    ungoverned_entry_ids = {record["source_entry_id"] for record in accepted}
    governance = None
    if governance_manifest:
        governance = json.loads(Path(governance_manifest).read_text(encoding="utf-8"))
    governance_pass = bool(
        governance and governance.get("gate_status") == "PASS" and
        governance.get("identity_threshold") == 0.4 and
        governance.get("structure_records") == len(accepted) and
        governance.get("structure_records_sha256") == ungoverned_record_sha256 and
        set(governance.get("excluded_source_entry_ids", [])) <= ungoverned_entry_ids)
    excluded_entries = set(governance.get("excluded_source_entry_ids", [])) if governance else set()
    if governance_pass:
        accepted = [record for record in accepted
                    if record["source_entry_id"] not in excluded_entries]
    unique_sequences = len({record["sequence"] for record in accepted})
    unique_chemotypes = len({record["canonical_smiles"] for record in accepted})
    quantitative_pass = (len(accepted) >= 10000 and unique_sequences >= 2000 and
                         unique_chemotypes >= 2000)
    gate_status = ("PASS" if quantitative_pass and governance_pass and provenance_pass
                   else "NOT_RUN_FAIL_CLOSED")
    ccd_snapshot = sorted({(record["ligand_comp_id"], record["ccd_sha256"])
                           for record in accepted})
    output.mkdir(parents=True)
    write_jsonl(output / "complexes.jsonl", accepted)
    manifest = {
        "schema": "MetaSieve.HoloComplexIndexManifest.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "record_schema": HOLO_COMPLEX_SCHEMA,
        "annotation_path": str(Path(annotation_path).resolve()),
        "annotation_sha256": sha256_file(annotation_path),
        "acquisition_manifest_sha256": (sha256_file(acquisition_path)
                                         if acquisition_path.is_file() else None),
        "coordinate_source": "RCSB_PDB",
        "coordinate_license": "CC0-1.0",
        "biological_relevance_source": "BioLiP2",
        "pdb_release_cutoff": None,
        "pdb_snapshot_utc": acquisition.get("created_utc") if acquisition else None,
        "quality_filter": QUALITY_FILTER,
        "quality_filter_sha256": hashlib.sha256(json.dumps(
            QUALITY_FILTER, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest(),
        "ccd_snapshot_sha256": hashlib.sha256(json.dumps(
            ccd_snapshot, separators=(",", ":")).encode("utf-8")).hexdigest(),
        "candidate_limit": candidate_limit,
        "candidates": len(candidates),
        "ungoverned_record_sha256": ungoverned_record_sha256,
        "valid_holo_complexes": len(accepted),
        "unique_receptor_sequences": unique_sequences,
        "unique_ligand_chemotypes": unique_chemotypes,
        "rejections": dict(sorted(rejected.items())),
        "rejection_examples": dict(sorted(rejection_examples.items())),
        "homology_governance": governance,
        "homology_governance_manifest_sha256": (
            sha256_file(governance_manifest) if governance_manifest else None),
        "ligand_overlap_audit_sha256": (hashlib.sha256(json.dumps(
            governance.get("ligand_overlap_report"), sort_keys=True,
            separators=(",", ":")).encode("utf-8")).hexdigest()
            if governance else None),
        "gate_requirements": {"valid_holo_complexes": 10000,
                              "unique_receptor_sequences": 2000,
                              "unique_ligand_chemotypes": 2000,
                              "benchmark_identity_exclusion": 0.4},
        "quantitative_pass": quantitative_pass,
        "governance_pass": governance_pass,
        "provenance_pass": provenance_pass,
        "gate_status": gate_status,
        "gate_reason": None if gate_status == "PASS" else
            "P1A requires corpus thresholds, CC0 acquisition provenance, and hash-bound 40% benchmark homology exclusion",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation")
    parser.add_argument("raw_root")
    parser.add_argument("output")
    parser.add_argument("--candidate-limit", type=int, default=15000)
    parser.add_argument("--governance-manifest")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = build_holo_complex_index(
        args.annotation, args.raw_root, args.output,
        candidate_limit=args.candidate_limit,
        governance_manifest=args.governance_manifest, workers=args.workers)
    print(json.dumps(result, indent=2))
    return 0 if result["gate_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
