"""Run the preregistered research-only T-DIR-P0 feasibility pilot."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import platform
import sys
import warnings

import numpy as np
import torch

from scripts.build_holo_complex_index import (
    _atom_rows,
    _canonical_altloc_rows,
    _ccd_molecule,
    _protein_sequence_mapping,
)
from scripts.pretrain_mechanistic_bridge import MechanismPretrainer, TrainConfig
from scripts.structure_sources.rcsb import sha256_file


STAGE = "P1R2B-TDIR-P0_ANNOTATION_AND_LEARNABILITY_PILOT"
SEED = 1701
SPLIT_COUNTS = {"train": 24, "val": 8, "test": 8}
PAIR_CUTOFF_ANGSTROM = 8.0
CHANNELS = (
    "hbond_protein_donor",
    "hbond_ligand_donor",
    "hydrophobic",
    "salt_ligand_negative",
    "salt_protein_negative",
    "pi_stacking_group",
    "cation_pi_group",
    "halogen",
)
PRIMARY_CHANNEL = "hydrophobic"
RESIDUE_CLASSES = {
    "A": 0, "V": 0, "I": 0, "L": 0, "M": 0, "C": 0,
    "F": 1, "W": 1, "Y": 1,
    "S": 2, "T": 2, "N": 2, "Q": 2,
    "K": 3, "R": 3, "H": 3,
    "D": 4, "E": 4,
    "G": 5, "P": 5, "X": 5,
}


def _canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict]) -> None:
    path.write_text("".join(_canonical_json(row) + "\n" for row in rows), encoding="utf-8")


def _selection_key(record: dict) -> str:
    return hashlib.sha256(
        f"TDIR-P0|{record['source_entry_id']}".encode("utf-8")
    ).hexdigest()


def select_records(records: list[dict]) -> tuple[list[dict], dict]:
    from rdkit import Chem

    selected: list[dict] = []
    used_groups: set[str] = set()
    used_pdb: set[str] = set()
    used_sequences: set[str] = set()
    earlier_scaffolds: set[str] = set()
    exclusions = Counter()
    for split in ("train", "val", "test"):
        candidates = sorted(
            (row for row in records if row["source_split"] == split),
            key=_selection_key,
        )
        split_rows = []
        for row in candidates:
            scaffold = row.get("murcko_scaffold", "")
            reasons = []
            if not scaffold:
                reasons.append("empty_scaffold")
            if row["homology_group_id"] in used_groups:
                reasons.append("homology_group_overlap")
            if row["pdb_id"] in used_pdb:
                reasons.append("pdb_overlap")
            if row["sequence_sha256"] in used_sequences:
                reasons.append("exact_sequence_overlap")
            if split != "train" and scaffold in earlier_scaffolds:
                reasons.append("earlier_split_scaffold_overlap")
            if len(row["ligand_comp_id"]) > 3:
                reasons.append("pdb_resname_too_long")
            molecule = Chem.MolFromSmiles(row["canonical_smiles"])
            if molecule is None:
                reasons.append("canonical_smiles_parse_failure")
            elif any(atom.GetAtomicNum() not in {
                    5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53
                    } for atom in molecule.GetAtoms()):
                reasons.append("metal_containing_ligand")
            if reasons:
                exclusions.update(reasons)
                continue
            split_rows.append(row)
            used_groups.add(row["homology_group_id"])
            used_pdb.add(row["pdb_id"])
            used_sequences.add(row["sequence_sha256"])
            earlier_scaffolds.add(scaffold)
            if len(split_rows) == SPLIT_COUNTS[split]:
                break
        if len(split_rows) != SPLIT_COUNTS[split]:
            raise RuntimeError(f"insufficient label-blind candidates for {split}")
        selected.extend(split_rows)
    audit = {
        "selection_rule": "sha256(TDIR-P0|source_entry_id)",
        "counts": dict(Counter(row["source_split"] for row in selected)),
        "unique_homology_groups": len({row["homology_group_id"] for row in selected}),
        "unique_pdb_ids": len({row["pdb_id"] for row in selected}),
        "unique_sequences": len({row["sequence_sha256"] for row in selected}),
        "exact_connectivity_overlap": _cross_split_overlap(selected, "connectivity_sha256"),
        "scaffold_overlap": _cross_split_overlap(selected, "murcko_scaffold"),
        "exclusions_before_quota": dict(exclusions),
        "p1b_exposure": {
            "train": "checkpoint training split",
            "val": "checkpoint held-out validation split",
            "test": "checkpoint held-out test split",
        },
    }
    if audit["unique_homology_groups"] != sum(SPLIT_COUNTS.values()):
        raise RuntimeError("cross-split homology group overlap")
    if audit["scaffold_overlap"]:
        raise RuntimeError("cross-split scaffold overlap")
    return selected, audit


def _cross_split_overlap(records: list[dict], field: str) -> list[dict]:
    owners: dict[str, set[str]] = defaultdict(set)
    for row in records:
        owners[str(row[field])].add(row["source_split"])
    return [
        {field: key, "splits": sorted(splits)}
        for key, splits in sorted(owners.items()) if len(splits) > 1
    ]


def _coordinate_bundle(record: dict) -> dict:
    import gemmi

    block = gemmi.cif.read(record["structure_path"]).sole_block()
    rows = _atom_rows(block)
    protein_rows = [
        row for row in rows
        if row["group_PDB"] == "ATOM"
        and row["auth_asym_id"] == record["protein_auth_asym_id"]
        and row["label_seq_id"] not in {"", ".", "?"}
        and row["type_symbol"].upper() != "H"
    ]
    protein_rows = _canonical_altloc_rows(
        protein_rows, ("label_asym_id", "label_seq_id", "label_atom_id")
    )
    labels, indices, coverage = _protein_sequence_mapping(protein_rows, record["sequence"])
    if coverage < 0.9:
        raise ValueError("protein mapping coverage below governed threshold")
    label_to_sequence = dict(zip(labels, indices))
    mapped_protein = [row for row in protein_rows if row["label_seq_id"] in label_to_sequence]

    ligand_rows = [
        row for row in rows
        if row["group_PDB"] == "HETATM"
        and row["label_comp_id"].upper() == record["ligand_comp_id"].upper()
        and row["auth_asym_id"] == record["ligand_auth_asym_id"]
        and row["auth_seq_id"] == record["ligand_auth_seq_id"]
        and row["type_symbol"].upper() != "H"
    ]
    ligand_rows = _canonical_altloc_rows(
        ligand_rows, ("label_asym_id", "label_seq_id", "label_atom_id")
    )
    chemistry = _ccd_molecule(Path(record["ccd_path"]))
    by_name = {row["label_atom_id"]: row for row in ligand_rows}
    if set(by_name) != set(chemistry["heavy_atom_names"]):
        raise ValueError("canonical ligand atom mapping is not exact")
    ordered_ligand = [by_name[name] for name in chemistry["heavy_atom_names"]]
    return {
        "protein_rows": mapped_protein,
        "label_to_sequence": label_to_sequence,
        "ligand_rows": ordered_ligand,
        "chemistry": chemistry,
    }


def _pdb_atom_line(group: str, serial: int, atom_name: str, residue: str,
                   chain: str, residue_number: int, row: dict) -> str:
    name = atom_name[:4]
    return (
        f"{group:<6}{serial:5d} {name:^4s} {residue[:3]:>3s} {chain}"
        f"{residue_number:4d}    {float(row['Cartn_x']):8.3f}"
        f"{float(row['Cartn_y']):8.3f}{float(row['Cartn_z']):8.3f}"
        f"{1.0:6.2f}{0.0:6.2f}          {row['type_symbol'][:2]:>2s}\n"
    )


def write_plip_input(record: dict, destination: Path) -> dict:
    bundle = _coordinate_bundle(record)
    serial = 1
    lines: list[str] = []
    residue_coordinates: dict[int, list[list[float]]] = defaultdict(list)
    for row in sorted(bundle["protein_rows"], key=lambda value: (
            bundle["label_to_sequence"][value["label_seq_id"]], value["label_atom_id"])):
        sequence_index = bundle["label_to_sequence"][row["label_seq_id"]]
        lines.append(_pdb_atom_line(
            "ATOM", serial, row["label_atom_id"], row["label_comp_id"],
            "A", sequence_index + 1, row,
        ))
        residue_coordinates[sequence_index].append([
            float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])
        ])
        serial += 1
    lines.append("TER\n")
    ligand_serial_to_atom: dict[int, int] = {}
    atom_to_serial: dict[int, int] = {}
    ligand_coordinates = []
    for atom_index, row in enumerate(bundle["ligand_rows"]):
        lines.append(_pdb_atom_line(
            "HETATM", serial, row["label_atom_id"], record["ligand_comp_id"],
            "Z", 1, row,
        ))
        ligand_serial_to_atom[serial] = atom_index
        atom_to_serial[atom_index] = serial
        ligand_coordinates.append([
            float(row["Cartn_x"]), float(row["Cartn_y"]), float(row["Cartn_z"])
        ])
        serial += 1
    heavy_names = bundle["chemistry"]["heavy_atom_names"]
    name_to_heavy = {name: index for index, name in enumerate(heavy_names)}
    molecule = bundle["chemistry"]["molecule"]
    for bond in molecule.GetBonds():
        left = molecule.GetAtomWithIdx(bond.GetBeginAtomIdx())
        right = molecule.GetAtomWithIdx(bond.GetEndAtomIdx())
        if left.GetAtomicNum() <= 1 or right.GetAtomicNum() <= 1:
            continue
        left_name = left.GetProp("_CCDAtomName")
        right_name = right.GetProp("_CCDAtomName")
        if left_name in name_to_heavy and right_name in name_to_heavy:
            lines.append(
                f"CONECT{atom_to_serial[name_to_heavy[left_name]]:5d}"
                f"{atom_to_serial[name_to_heavy[right_name]]:5d}\n"
            )
    lines.append("END\n")
    destination.write_text("".join(lines), encoding="ascii")
    return {
        "ligand_serial_to_atom": ligand_serial_to_atom,
        "ligand_coordinates": np.asarray(ligand_coordinates, dtype=np.float32),
        "residue_coordinates": {
            index: np.asarray(coords, dtype=np.float32)
            for index, coords in residue_coordinates.items()
        },
    }


def _collect_ligand_serials(value, ligand_serials: set[int], key: str = "") -> set[int]:
    result: set[int] = set()
    if isinstance(value, (int, np.integer)) and "orig_idx" in key and int(value) in ligand_serials:
        result.add(int(value))
    elif isinstance(value, dict):
        for child_key, child in value.items():
            result.update(_collect_ligand_serials(child, ligand_serials, str(child_key)))
    elif hasattr(value, "_asdict"):
        result.update(_collect_ligand_serials(value._asdict(), ligand_serials, key))
    elif isinstance(value, (list, tuple)):
        for child in value:
            result.update(_collect_ligand_serials(child, ligand_serials, key))
    return result


def extract_plip_labels(record: dict, pdb_path: Path, mapping: dict) -> tuple[dict, list[dict]]:
    from openbabel import pybel
    from plip.structure.preparation import PDBComplex

    original_write = pybel.Molecule.write

    def metadata_safe_write(self, format="smi", filename=None, overwrite=False, opt=None):
        if format == "inchikey":
            return "INCHIKEY_UNAVAILABLE_IN_LOCAL_OPENBABEL_BUILD"
        return original_write(self, format, filename, overwrite, opt)

    pybel.Molecule.write = metadata_safe_write
    try:
        complex_value = PDBComplex()
        complex_value.load_pdb(str(pdb_path))
        matches = [
            ligand for ligand in complex_value.ligands
            if ligand.hetid.upper() == record["ligand_comp_id"].upper()
            and ligand.chain == "Z" and int(ligand.position) == 1
        ]
        if len(matches) != 1:
            raise ValueError(f"PLIP target ligand match count is {len(matches)}")
        complex_value.characterize_complex(matches[0])
        site = complex_value.interaction_sets[
            f"{record['ligand_comp_id']}:Z:1"
        ]
    finally:
        pybel.Molecule.write = original_write

    source_lists = {
        "hbond_protein_donor": site.hbonds_pdon,
        "hbond_ligand_donor": site.hbonds_ldon,
        "hydrophobic": site.hydrophobic_contacts,
        "salt_ligand_negative": site.saltbridge_lneg,
        "salt_protein_negative": site.saltbridge_pneg,
        "pi_stacking_group": site.pistacking,
        "cation_pi_group": list(site.pication_paro) + list(site.pication_laro),
        "halogen": site.halogen_bonds,
    }
    ligand_serials = set(mapping["ligand_serial_to_atom"])
    labels = {channel: set() for channel in CHANNELS}
    raw_rows = []
    unmapped = Counter()
    for channel, interactions in source_lists.items():
        for interaction in interactions:
            values = interaction._asdict()
            residue_number = int(values["resnr"])
            sequence_index = residue_number - 1
            serial_values = _collect_ligand_serials(values, ligand_serials)
            if sequence_index not in mapping["residue_coordinates"]:
                unmapped[f"{channel}:residue"] += 1
                continue
            if not serial_values:
                unmapped[f"{channel}:ligand_atom"] += 1
                continue
            atom_indices = sorted(mapping["ligand_serial_to_atom"][value] for value in serial_values)
            for atom_index in atom_indices:
                labels[channel].add((atom_index, sequence_index))
            raw_rows.append({
                "channel": channel,
                "atom_indices": atom_indices,
                "sequence_index": sequence_index,
                "residue_type": str(values.get("restype", "")),
            })
    audit = {
        "counts": {channel: len(labels[channel]) for channel in CHANNELS},
        "raw_interaction_counts": {
            channel: len(source_lists[channel]) for channel in CHANNELS
        },
        "unmapped": dict(unmapped),
    }
    return {channel: sorted(values) for channel, values in labels.items()}, raw_rows + [
        {"_mapping_audit": audit}
    ]


def _candidate_pairs(mapping: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    atoms, residues, distances = [], [], []
    ligand = mapping["ligand_coordinates"]
    for sequence_index, residue_coords in sorted(mapping["residue_coordinates"].items()):
        difference = ligand[:, None, :] - residue_coords[None, :, :]
        minimum = np.sqrt(np.square(difference).sum(axis=2)).min(axis=1)
        for atom_index in np.flatnonzero(minimum <= PAIR_CUTOFF_ANGSTROM):
            atoms.append(int(atom_index))
            residues.append(int(sequence_index))
            distances.append(float(minimum[atom_index]))
    return (
        np.asarray(atoms, dtype=np.int16),
        np.asarray(residues, dtype=np.int16),
        np.asarray(distances, dtype=np.float32),
    )


def _load_protein_rows(root: Path, keys: set[str]) -> dict[str, dict[str, np.ndarray]]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    unresolved = set(keys)
    result = {}
    for item in manifest["shards"]:
        if not unresolved:
            break
        path = root / item["path"]
        with np.load(path, allow_pickle=False) as shard:
            shard_keys = [str(value) for value in shard["keys"]]
            wanted = [(key, shard_keys.index(key)) for key in unresolved if key in shard_keys]
            if not wanted:
                continue
            if sha256_file(path) != item["sha256"]:
                raise ValueError(f"protein cache hash mismatch: {path}")
            for key, index in wanted:
                result[key] = {
                    name: shard[name][index].copy()
                    for name in ("pooled", "residues", "mask")
                }
                unresolved.remove(key)
    if unresolved:
        raise KeyError(f"missing protein cache rows: {sorted(unresolved)}")
    return result


def _load_frozen_model(checkpoint_path: Path, protein_dim: int, device: str):
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    config = TrainConfig(**checkpoint["config"])
    model = MechanismPretrainer(protein_dim, config)
    model.load_state_dict(checkpoint["model_state"])
    model.eval().to(device)
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model, checkpoint


def _frozen_features(model, graph: dict, protein: dict, device: str) -> dict[str, np.ndarray]:
    from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM, MAX_ATOMS

    X = torch.zeros(1, MAX_ATOMS, ATOM_FEAT_DIM, device=device)
    A = torch.zeros(1, MAX_ATOMS, MAX_ATOMS, BOND_FEAT_DIM, device=device)
    mask = graph["mask"].float().unsqueeze(0).to(device)
    X[0] = graph["X"].float().to(device)
    edge = graph["edge_index"].long().to(device)
    A[0, edge[0], edge[1]] = graph["edge_attr"].float().to(device)
    pooled = torch.from_numpy(protein["pooled"]).float().unsqueeze(0).to(device)
    residues = torch.from_numpy(protein["residues"]).float().unsqueeze(0).to(device)
    residue_mask = torch.from_numpy(protein["mask"]).float().unsqueeze(0).to(device)
    with torch.inference_mode():
        _, atom_states = model.ligand(X, A, mask)
        _, residue_states = model.protein(pooled, residues)
        prediction = model.bridge(atom_states, mask, residue_states, residue_mask)
    return {
        "raw_atom": X[0].cpu().numpy(),
        "atom_state": atom_states[0].cpu().numpy(),
        "residue_state": residue_states[0].cpu().numpy(),
        "contact": torch.sigmoid(prediction.contact_logits[0]).cpu().numpy(),
        "distance": torch.softmax(prediction.distance_logits[0], dim=-1).cpu().numpy(),
        "atom_mask": mask[0].cpu().numpy(),
        "residue_mask": residue_mask[0].cpu().numpy(),
    }


def _residue_one_hot(sequence_code: str) -> np.ndarray:
    result = np.zeros(6, dtype=np.float32)
    result[RESIDUE_CLASSES.get(sequence_code, 5)] = 1.0
    return result


def build_pair_dataset(record: dict, mapping: dict, labels: dict,
                       frozen: dict) -> dict[str, np.ndarray]:
    atom_indices, residue_indices, actual_distance = _candidate_pairs(mapping)
    slots = np.floor(residue_indices.astype(np.float64) * 128 / len(record["sequence"])).astype(np.int16)
    if len(atom_indices) == 0:
        raise ValueError("no oracle near pairs")
    if atom_indices.max() >= int(frozen["atom_mask"].sum()):
        raise ValueError("ligand canonical atom index exceeds frozen graph mask")
    if np.any(frozen["residue_mask"][slots] == 0):
        raise ValueError("mapped explicit residue enters a masked P1B slot")
    d0, d1, d2, target = [], [], [], []
    label_sets = {channel: set(map(tuple, values)) for channel, values in labels.items()}
    for atom_index, residue_index, slot in zip(atom_indices, residue_indices, slots):
        geometry = np.concatenate((
            [frozen["contact"][atom_index, slot]],
            frozen["distance"][atom_index, slot],
        )).astype(np.float32)
        chemistry = np.concatenate((
            geometry,
            frozen["raw_atom"][atom_index],
            _residue_one_hot(record["sequence"][residue_index]),
        )).astype(np.float32)
        local = np.concatenate((
            chemistry,
            frozen["atom_state"][atom_index],
            frozen["residue_state"][slot],
        )).astype(np.float32)
        d0.append(geometry)
        d1.append(chemistry)
        d2.append(local)
        target.append([
            int((int(atom_index), int(residue_index)) in label_sets[channel])
            for channel in CHANNELS
        ])
    multiplicity = Counter(slots.tolist())
    return {
        "d0": np.asarray(d0, dtype=np.float32),
        "d1": np.asarray(d1, dtype=np.float32),
        "d2": np.asarray(d2, dtype=np.float32),
        "target": np.asarray(target, dtype=np.uint8),
        "atom_index": atom_indices,
        "residue_index": residue_indices,
        "slot": slots,
        "actual_min_distance": actual_distance,
        "slot_multiplicity_mean": np.asarray([np.mean(list(multiplicity.values()))]),
        "slot_multiplicity_max": np.asarray([max(multiplicity.values())]),
    }


def _complex_weights(complex_ids: np.ndarray) -> np.ndarray:
    counts = Counter(complex_ids.tolist())
    scale = len(complex_ids) / len(counts)
    return np.asarray([scale / counts[value] for value in complex_ids], dtype=np.float64)


def _metrics(y: np.ndarray, prediction: np.ndarray, complex_ids: np.ndarray) -> dict:
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

    weights = _complex_weights(complex_ids)
    prevalence = float(np.average(y, weights=weights))
    ap = float(average_precision_score(y, prediction, sample_weight=weights))
    per_complex = []
    directions = 0
    for value in sorted(set(complex_ids.tolist())):
        mask = complex_ids == value
        if len(np.unique(y[mask])) < 2:
            continue
        local_prevalence = float(y[mask].mean())
        local_ap = float(average_precision_score(y[mask], prediction[mask]))
        directions += int(local_ap > local_prevalence)
        per_complex.append(local_ap)
    return {
        "prevalence": prevalence,
        "average_precision": ap,
        "normalized_ap_lift": (ap - prevalence) / max(1.0 - prevalence, 1e-12),
        "auroc": float(roc_auc_score(y, prediction, sample_weight=weights)),
        "brier": float(brier_score_loss(y, prediction, sample_weight=weights)),
        "prediction_sd": float(np.std(prediction)),
        "eligible_complexes": len(per_complex),
        "per_complex_ap_median": float(np.median(per_complex)) if per_complex else None,
        "per_complex_ap_iqr": (
            [float(np.quantile(per_complex, 0.25)), float(np.quantile(per_complex, 0.75))]
            if per_complex else None
        ),
        "complexes_ap_above_prevalence": directions,
    }


def _channel_counts(target: np.ndarray, complex_ids: np.ndarray, channel_index: int) -> dict:
    y = target[:, channel_index]
    positives = {value for value in set(complex_ids.tolist()) if y[complex_ids == value].sum() > 0}
    return {"pairs": len(y), "positives": int(y.sum()), "positive_complexes": len(positives),
            "both_classes": len(np.unique(y)) == 2}


def _evaluable(counts: dict[str, dict]) -> bool:
    return (
        counts["train"]["both_classes"] and counts["val"]["both_classes"]
        and counts["test"]["both_classes"]
        and counts["train"]["positives"] >= 50
        and counts["train"]["positive_complexes"] >= 8
        and counts["val"]["positives"] >= 10
        and counts["val"]["positive_complexes"] >= 3
        and counts["test"]["positives"] >= 10
        and counts["test"]["positive_complexes"] >= 3
    )


def fit_and_evaluate(arrays: dict[str, dict[str, np.ndarray]]) -> dict:
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    report = {"channels": {}, "primary_channel": PRIMARY_CHANNEL}
    for channel_index, channel in enumerate(CHANNELS):
        counts = {
            split: _channel_counts(value["target"], value["complex_id"], channel_index)
            for split, value in arrays.items()
        }
        channel_report = {"counts": counts, "evaluable": _evaluable(counts), "arms": {}}
        if not channel_report["evaluable"]:
            channel_report["verdict"] = "NOT_EVALUABLE_PILOT_CHANNEL"
            report["channels"][channel] = channel_report
            continue
        for arm in ("d0", "d1", "d2"):
            model = make_pipeline(
                StandardScaler(),
                LogisticRegression(
                    C=1.0, class_weight="balanced", max_iter=500,
                    penalty="l2", random_state=SEED, solver="lbfgs",
                ),
            )
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", ConvergenceWarning)
                model.fit(
                    arrays["train"][arm], arrays["train"]["target"][:, channel_index],
                    logisticregression__sample_weight=_complex_weights(arrays["train"]["complex_id"]),
                )
            if any(issubclass(item.category, ConvergenceWarning) for item in caught):
                raise RuntimeError(f"logistic solver did not converge: {channel}/{arm}")
            arm_report = {
                "feature_dimension": int(arrays["train"][arm].shape[1]),
                "iterations": int(model[-1].n_iter_[0]),
            }
            for split in ("train", "val", "test"):
                prediction = model.predict_proba(arrays[split][arm])[:, 1]
                if not np.isfinite(prediction).all() or np.std(prediction) == 0:
                    raise RuntimeError(f"invalid predictions: {channel}/{arm}/{split}")
                arm_report[split] = _metrics(
                    arrays[split]["target"][:, channel_index], prediction,
                    arrays[split]["complex_id"],
                )
            channel_report["arms"][arm] = arm_report
        channel_report["verdict"] = "DESCRIPTIVE_MATRIX_COMPLETE"
        report["channels"][channel] = channel_report

    primary = report["channels"][PRIMARY_CHANNEL]
    if not primary["evaluable"]:
        report["verdict"] = "PILOT_ANNOTATION_INFEASIBLE"
        return report
    d0, d2 = primary["arms"]["d0"], primary["arms"]["d2"]
    criteria = {
        "val_d2_above_prevalence": d2["val"]["average_precision"] > d2["val"]["prevalence"],
        "test_d2_above_prevalence": d2["test"]["average_precision"] > d2["test"]["prevalence"],
        "val_d2_above_d0": d2["val"]["average_precision"] > d0["val"]["average_precision"],
        "test_d2_above_d0": d2["test"]["average_precision"] > d0["test"]["average_precision"],
        "test_d2_ap_minus_prevalence_ge_0_10": (
            d2["test"]["average_precision"] - d2["test"]["prevalence"] >= 0.10
        ),
        "test_complex_direction_count_ge_5": d2["test"]["complexes_ap_above_prevalence"] >= 5,
    }
    report["feasibility_criteria"] = criteria
    report["verdict"] = (
        "PILOT_FEASIBILITY_SIGNAL_OBSERVED"
        if all(criteria.values()) else "PILOT_LEARNABILITY_SIGNAL_NOT_OBSERVED"
    )
    return report


def _stack_datasets(per_record: list[dict]) -> dict[str, dict[str, np.ndarray]]:
    result = {}
    for split in ("train", "val", "test"):
        values = [item for item in per_record if item["split"] == split]
        result[split] = {
            key: np.concatenate([item[key] for item in values], axis=0)
            for key in ("d0", "d1", "d2", "target")
        }
        result[split]["complex_id"] = np.concatenate([
            np.repeat(item["source_entry_id"], len(item["target"])) for item in values
        ])
    return result


def run(args) -> dict:
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("T-DIR-P0 frozen P1B inference requires CUDA")
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.mkdir(parents=True)
    input_dir = output / "plip_inputs"
    input_dir.mkdir()

    records = [json.loads(line) for line in Path(args.records).read_text(encoding="utf-8").splitlines()]
    selected, selection_audit = select_records(records)
    selection_rows = [{
        "source_entry_id": row["source_entry_id"], "source_split": row["source_split"],
        "pdb_id": row["pdb_id"], "homology_group_id": row["homology_group_id"],
        "sequence_sha256": row["sequence_sha256"],
        "connectivity_sha256": row["connectivity_sha256"],
        "murcko_scaffold": row["murcko_scaffold"], "selection_key": _selection_key(row),
    } for row in selected]
    _write_jsonl(output / "selection.jsonl", selection_rows)
    selection_audit["selection_sha256"] = sha256_file(output / "selection.jsonl")
    _write_json(output / "selection_audit.json", selection_audit)

    protein_rows = _load_protein_rows(
        Path(args.protein_bank), {row["sequence_sha256"] for row in selected}
    )
    protein_dim = len(next(iter(protein_rows.values()))["pooled"])
    model, checkpoint = _load_frozen_model(Path(args.checkpoint), protein_dim, args.device)
    ligands = torch.load(args.ligand_bank, map_location="cpu", weights_only=False)

    per_record, annotation_rows, failures = [], [], []
    for number, record in enumerate(selected, start=1):
        try:
            pdb_path = input_dir / f"{number:03d}_{record['pdb_id']}.pdb"
            mapping = write_plip_input(record, pdb_path)
            labels, raw_interactions = extract_plip_labels(record, pdb_path, mapping)
            frozen = _frozen_features(
                model, ligands[record["ccd_sha256"]],
                protein_rows[record["sequence_sha256"]], args.device,
            )
            dataset = build_pair_dataset(record, mapping, labels, frozen)
            dataset.update({"split": record["source_split"],
                            "source_entry_id": record["source_entry_id"]})
            per_record.append(dataset)
            annotation_rows.append({
                "source_entry_id": record["source_entry_id"],
                "source_split": record["source_split"],
                "pdb_input": str(pdb_path.resolve()),
                "pdb_sha256": sha256_file(pdb_path),
                "candidate_pairs": len(dataset["target"]),
                "slot_multiplicity_mean": float(dataset["slot_multiplicity_mean"][0]),
                "slot_multiplicity_max": int(dataset["slot_multiplicity_max"][0]),
                "channel_positive_pairs": {
                    channel: int(dataset["target"][:, index].sum())
                    for index, channel in enumerate(CHANNELS)
                },
                "plip_interactions": raw_interactions,
            })
        except Exception as error:
            failures.append({"source_entry_id": record["source_entry_id"],
                             "source_split": record["source_split"],
                             "error": f"{type(error).__name__}: {error}"})
    _write_jsonl(output / "annotations.jsonl", annotation_rows)
    _write_jsonl(output / "failures.jsonl", failures)

    annotation_success = len(per_record)
    if annotation_success < 36 or any(
            sum(item["split"] == split for item in per_record) == 0
            for split in SPLIT_COUNTS):
        report = {
            "stage": STAGE, "verdict": "PILOT_DATA_OR_MAPPING_CONTRACT_FAIL_CLOSED",
            "selected": len(selected), "annotation_success": annotation_success,
            "failures": failures,
        }
    else:
        arrays = _stack_datasets(per_record)
        report = fit_and_evaluate(arrays)
        report.update({"stage": STAGE, "selected": len(selected),
                       "annotation_success": annotation_success,
                       "annotation_failures": len(failures)})

    try:
        import plip
        plip_version = getattr(plip, "__version__", "3.0.1-package-metadata")
    except Exception:
        plip_version = "unknown"
    manifest = {
        "stage": STAGE,
        "research_only": True,
        "oracle_near_pair_conditional": True,
        "forbidden_label_reads": {"chembl_affinity": 0, "davis": 0, "recipient": 0},
        "versions": {
            "python": sys.version, "platform": platform.platform(),
            "torch": torch.__version__, "numpy": np.__version__, "plip": plip_version,
        },
        "inputs": {
            "records": {"path": str(Path(args.records).resolve()),
                        "sha256": sha256_file(args.records)},
            "checkpoint": {"path": str(Path(args.checkpoint).resolve()),
                           "sha256": sha256_file(args.checkpoint),
                           "epoch": checkpoint["epoch"]},
            "ligand_bank": {"path": str(Path(args.ligand_bank).resolve()),
                            "sha256": sha256_file(args.ligand_bank)},
            "protein_manifest": {"path": str((Path(args.protein_bank) / 'manifest.json').resolve()),
                                 "sha256": sha256_file(Path(args.protein_bank) / "manifest.json")},
        },
        "plip_openbabel_compatibility": (
            "local Open Babel lacks the inchikey writer; only PLIP ligand metadata "
            "serialization is replaced with a fixed unavailable marker"
        ),
        "selection_sha256": selection_audit["selection_sha256"],
        "annotations_sha256": sha256_file(output / "annotations.jsonl"),
        "failures_sha256": sha256_file(output / "failures.jsonl"),
    }
    _write_json(output / "report.json", report)
    manifest["report_sha256"] = sha256_file(output / "report.json")
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", default="dataset/processed/open_structures/pilot20k_homology_split_v2/complexes.jsonl")
    parser.add_argument("--protein-bank", default="dataset/processed/open_structures/pilot20k_esm2_t30_slots128_v1")
    parser.add_argument("--ligand-bank", default="dataset/processed/open_structures/pilot20k_mechanism_ligands_v1.pt")
    parser.add_argument("--checkpoint", default="report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/tdir_p0_v1")
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    result = run(parse_args())
    print(json.dumps(result, indent=2, sort_keys=True))
