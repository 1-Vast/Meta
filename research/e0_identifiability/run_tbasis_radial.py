"""Run T-BASIS-R0 fixed radial basis recoverability on a fresh structure panel."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from contracts.mechanism import DISTANCE_BINS_ANGSTROM, MECHANISM_RESIDUE_SLOTS
from research.e0_identifiability.run_tdir_pilot import (
    RESIDUE_CLASSES,
    _coordinate_bundle,
    _frozen_features,
    _load_frozen_model,
    _load_protein_rows,
)
from scripts.structure_sources.rcsb import sha256_file


STAGE = "P1R2B-T-BASIS-R0_FIXED_RADIAL_BASIS_RECOVERABILITY"
SEED = 2718
SPLIT_COUNTS = {"train": 192, "val": 64, "test": 64}
CENTERS = np.asarray([2.0, 3.5, 5.0, 6.5, 8.0, 9.5], dtype=np.float64)
RBF_SIGMA = 1.0
CUTOFF = 10.0
ATOM_CHANNELS = (
    "hydrophobe", "aromatic", "donor", "acceptor",
    "positive", "negative", "halogen", "other",
)


def _canonical_json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict]) -> None:
    path.write_text("".join(_canonical_json(row) + "\n" for row in rows), encoding="utf-8")


def _selection_key(record: dict) -> str:
    return hashlib.sha256(
        f"T-BASIS-R0|{record['source_entry_id']}".encode("utf-8")
    ).hexdigest()


def select_panel(records: list[dict], excluded_entries: set[str]) -> tuple[list[dict], dict]:
    selected = []
    used_groups, used_pdb, used_sequences, earlier_scaffolds = set(), set(), set(), set()
    exclusions = Counter()
    for split in ("train", "val", "test"):
        split_rows = []
        for row in sorted(
                (item for item in records if item["source_split"] == split),
                key=_selection_key):
            scaffold = row.get("murcko_scaffold", "")
            reasons = []
            if row["source_entry_id"] in excluded_entries:
                reasons.append("tdir_p0_panel")
            if float(row.get("protein_mapping_coverage", 0.0)) < 0.999:
                reasons.append("incomplete_protein_mapping")
            if not scaffold:
                reasons.append("empty_scaffold")
            if row["homology_group_id"] in used_groups:
                reasons.append("homology_group_overlap")
            if row["pdb_id"] in used_pdb:
                reasons.append("pdb_overlap")
            if row["sequence_sha256"] in used_sequences:
                reasons.append("sequence_overlap")
            if split != "train" and scaffold in earlier_scaffolds:
                reasons.append("earlier_split_scaffold_overlap")
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
            raise RuntimeError(f"insufficient eligible records for {split}")
        selected.extend(split_rows)
    audit = {
        "counts": dict(Counter(row["source_split"] for row in selected)),
        "unique_homology_groups": len({row["homology_group_id"] for row in selected}),
        "unique_pdb_ids": len({row["pdb_id"] for row in selected}),
        "unique_sequences": len({row["sequence_sha256"] for row in selected}),
        "unique_scaffolds": len({row["murcko_scaffold"] for row in selected}),
        "excluded_tdir_p0": sum(row["source_entry_id"] in excluded_entries for row in records),
        "exclusions_before_quota": dict(exclusions),
        "p1b_exposure": {
            "train": "P1B training split",
            "val": "P1B held-out validation split",
            "test": "P1B held-out test split",
        },
    }
    expected = sum(SPLIT_COUNTS.values())
    if any(audit[key] != expected for key in (
            "unique_homology_groups", "unique_pdb_ids", "unique_sequences")):
        raise RuntimeError("fresh panel global closure failed")
    return selected, audit


def sequence_identity(left: str, right: str) -> float:
    import parasail

    trace = parasail.nw_trace_striped_16(left, right, 10, 1, parasail.blosum62).traceback
    matches = sum(a == b for a, b in zip(trace.query, trace.ref) if a != "-" and b != "-")
    return matches / min(len(left), len(right))


def build_derangement(selected: list[dict]) -> tuple[dict[str, dict], list[dict]]:
    from scipy.optimize import linear_sum_assignment

    mapping, audit = {}, []
    for split in ("train", "val", "test"):
        rows = sorted(
            (row for row in selected if row["source_split"] == split),
            key=lambda row: row["source_entry_id"],
        )
        size = len(rows)
        ineligible_cost = 1e20
        cost = np.full((size, size), ineligible_cost, dtype=np.float64)
        identities = np.full((size, size), np.nan, dtype=np.float64)
        ratios = np.full((size, size), np.nan, dtype=np.float64)
        for left, correct in enumerate(rows):
            for right, wrong in enumerate(rows):
                if wrong["homology_group_id"] == correct["homology_group_id"]:
                    continue
                ratio = len(wrong["sequence"]) / len(correct["sequence"])
                if not 0.5 <= ratio <= 2.0:
                    continue
                identity = sequence_identity(correct["sequence"], wrong["sequence"])
                if identity >= 0.40:
                    continue
                key = hashlib.sha256(
                    f"T-BASIS-R0-WRONG|{correct['source_entry_id']}|{wrong['source_entry_id']}".encode()
                ).hexdigest()
                cost[left, right] = int(key[:13], 16)
                identities[left, right] = identity
                ratios[left, right] = ratio
        left_indices, right_indices = linear_sum_assignment(cost)
        if np.any(cost[left_indices, right_indices] >= ineligible_cost):
            raise RuntimeError(f"no complete eligible derangement matching for {split}")
        for left, right in zip(left_indices, right_indices):
            correct, wrong = rows[left], rows[right]
            mapping[correct["source_entry_id"]] = wrong
            audit.append({
                "correct_entry": correct["source_entry_id"],
                "wrong_entry": wrong["source_entry_id"],
                "split": split,
                "identity": float(identities[left, right]),
                "length_ratio": float(ratios[left, right]),
                "correct_group": correct["homology_group_id"],
                "wrong_group": wrong["homology_group_id"],
            })
    if len(mapping) != len(selected) or len({row["wrong_entry"] for row in audit}) != len(selected):
        raise RuntimeError("derangement is not one-to-one")
    return mapping, audit


def radial_basis(distance: np.ndarray) -> np.ndarray:
    value = np.asarray(distance, dtype=np.float64)[..., None]
    gaussian = np.exp(-0.5 * np.square((value - CENTERS) / RBF_SIGMA))
    clipped = np.clip(value / CUTOFF, 0.0, 1.0)
    cutoff = 0.5 * (np.cos(np.pi * clipped) + 1.0)
    cutoff[value >= CUTOFF] = 0.0
    return gaussian * cutoff


def bin_rbf_expectation(points: int = 4096) -> np.ndarray:
    edges = np.asarray(DISTANCE_BINS_ANGSTROM, dtype=np.float64)
    result = np.zeros((len(edges) - 1, len(CENTERS)), dtype=np.float64)
    for index, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        if left >= CUTOFF:
            continue
        upper = min(right, CUTOFF)
        grid = np.linspace(left, upper, points, endpoint=False) + (upper - left) / (2 * points)
        result[index] = radial_basis(grid).mean(axis=0)
    return result


def slot_composition(sequence: str) -> np.ndarray:
    counts = np.zeros((MECHANISM_RESIDUE_SLOTS, 6), dtype=np.float64)
    for index, code in enumerate(sequence):
        slot = min(MECHANISM_RESIDUE_SLOTS - 1,
                   index * MECHANISM_RESIDUE_SLOTS // len(sequence))
        counts[slot, RESIDUE_CLASSES.get(code, 5)] += 1.0
    totals = counts.sum(axis=1, keepdims=True)
    return np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0)


def ligand_channels(record: dict) -> np.ndarray:
    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures
    from scripts.build_holo_complex_index import _ccd_molecule

    chemistry = _ccd_molecule(Path(record["ccd_path"]))
    molecule = chemistry["molecule"]
    name_to_heavy = {name: index for index, name in enumerate(chemistry["heavy_atom_names"])}
    atom_to_heavy = {
        atom.GetIdx(): name_to_heavy[atom.GetProp("_CCDAtomName")]
        for atom in molecule.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetProp("_CCDAtomName") in name_to_heavy
    }
    result = np.zeros((len(name_to_heavy), len(ATOM_CHANNELS)), dtype=np.float64)
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    family_to_index = {
        "Hydrophobe": 0, "LumpedHydrophobe": 0, "Donor": 2, "Acceptor": 3,
        "PosIonizable": 4, "NegIonizable": 5,
    }
    for feature in factory.GetFeaturesForMol(molecule):
        output_index = family_to_index.get(feature.GetFamily())
        if output_index is None:
            continue
        for atom_index in feature.GetAtomIds():
            if atom_index in atom_to_heavy:
                result[atom_to_heavy[atom_index], output_index] = 1.0
    for atom_index, heavy_index in atom_to_heavy.items():
        atom = molecule.GetAtomWithIdx(atom_index)
        result[heavy_index, 1] = float(atom.GetIsAromatic())
        result[heavy_index, 6] = float(atom.GetAtomicNum() in {9, 17, 35, 53})
    result[:, 7] = (result[:, :7].sum(axis=1) == 0).astype(np.float64)
    if np.any(result.sum(axis=1) == 0):
        raise ValueError("ligand channel mapping left an untyped heavy atom")
    return result


def exact_slot_distances(record: dict) -> tuple[np.ndarray, np.ndarray]:
    bundle = _coordinate_bundle(record)
    ligand = np.asarray([
        [float(row[axis]) for axis in ("Cartn_x", "Cartn_y", "Cartn_z")]
        for row in bundle["ligand_rows"]
    ], dtype=np.float64)
    slot_atoms: list[list[list[float]]] = [[] for _ in range(MECHANISM_RESIDUE_SLOTS)]
    for row in bundle["protein_rows"]:
        sequence_index = bundle["label_to_sequence"][row["label_seq_id"]]
        slot = min(MECHANISM_RESIDUE_SLOTS - 1,
                   sequence_index * MECHANISM_RESIDUE_SLOTS // len(record["sequence"]))
        slot_atoms[slot].append([float(row[axis]) for axis in ("Cartn_x", "Cartn_y", "Cartn_z")])
    distance = np.full((len(ligand), MECHANISM_RESIDUE_SLOTS), np.inf, dtype=np.float64)
    mask = np.zeros(MECHANISM_RESIDUE_SLOTS, dtype=bool)
    for slot, atoms in enumerate(slot_atoms):
        if not atoms:
            continue
        mask[slot] = True
        protein = np.asarray(atoms, dtype=np.float64)
        delta = ligand[:, None, :] - protein[None, :, :]
        distance[:, slot] = np.sqrt(np.square(delta).sum(axis=2)).min(axis=1)
    return distance, mask


def aggregate_basis(atom_channels: np.ndarray, residue_composition: np.ndarray,
                    radial: np.ndarray, slot_mask: np.ndarray) -> np.ndarray:
    if radial.shape[:2] != (len(atom_channels), len(residue_composition)):
        raise ValueError("radial tensor does not match atom and slot axes")
    masked_radial = radial * slot_mask[None, :, None]
    return np.einsum(
        "ia,sr,isk->ark", atom_channels, residue_composition, masked_radial,
        optimize=True,
    ) / len(atom_channels)


def compute_record_basis(record: dict, wrong: dict, model, proteins: dict,
                         ligands: dict, bin_moments: np.ndarray, device: str) -> dict:
    graph = ligands[record["ccd_sha256"]]
    atom_channels = ligand_channels(record)
    if len(atom_channels) != int(graph["mask"].sum()):
        raise ValueError("CCD channel order differs from ligand graph")
    distance, slot_mask = exact_slot_distances(record)
    teacher = aggregate_basis(
        atom_channels, slot_composition(record["sequence"]), radial_basis(distance), slot_mask,
    )
    correct = _frozen_features(model, graph, proteins[record["sequence_sha256"]], device)
    wrong_value = _frozen_features(model, graph, proteins[wrong["sequence_sha256"]], device)
    correct_radial = np.einsum("isb,bk->isk", correct["distance"][:len(atom_channels)], bin_moments)
    wrong_radial = np.einsum("isb,bk->isk", wrong_value["distance"][:len(atom_channels)], bin_moments)
    raw_correct = aggregate_basis(
        atom_channels, slot_composition(record["sequence"]), correct_radial,
        correct["residue_mask"].astype(bool),
    )
    raw_deranged = aggregate_basis(
        atom_channels, slot_composition(wrong["sequence"]), wrong_radial,
        wrong_value["residue_mask"].astype(bool),
    )
    return {"teacher": teacher, "raw_correct": raw_correct, "raw_deranged": raw_deranged}


def fit_radial_calibrator(train_values: list[dict]):
    from sklearn.linear_model import Ridge

    source = np.concatenate([value["raw_correct"].reshape(-1, len(CENTERS))
                             for value in train_values], axis=0)
    target = np.concatenate([value["teacher"].reshape(-1, len(CENTERS))
                             for value in train_values], axis=0)
    model = Ridge(alpha=1e-3, fit_intercept=True, solver="svd")
    model.fit(source, target)
    return model


def apply_calibrator(model, value: np.ndarray) -> np.ndarray:
    shape = value.shape
    return model.predict(value.reshape(-1, shape[-1])).reshape(shape)


def _complex_errors(values: list[dict], mean: np.ndarray, scale: np.ndarray,
                    active: np.ndarray, calibrator) -> dict[str, np.ndarray]:
    result = defaultdict(list)
    for value in values:
        teacher = ((value["teacher"].reshape(-1) - mean) / scale)[active]
        raw = ((value["raw_correct"].reshape(-1) - mean) / scale)[active]
        calibrated = ((apply_calibrator(calibrator, value["raw_correct"]).reshape(-1) - mean) / scale)[active]
        deranged = ((apply_calibrator(calibrator, value["raw_deranged"]).reshape(-1) - mean) / scale)[active]
        result["mean"].append(float(np.mean(np.square(teacher))))
        result["raw_correct"].append(float(np.mean(np.square(teacher - raw))))
        result["correct"].append(float(np.mean(np.square(teacher - calibrated))))
        result["deranged"].append(float(np.mean(np.square(teacher - deranged))))
    return {key: np.asarray(items, dtype=np.float64) for key, items in result.items()}


def _gain_summary(errors: dict[str, np.ndarray], bootstrap: int = 2000) -> dict:
    baseline = errors["mean"]
    reconstruction_delta = baseline - errors["correct"]
    partner_delta = errors["deranged"] - errors["correct"]
    reconstruction = float(reconstruction_delta.mean() / baseline.mean())
    partner = float(partner_delta.mean() / baseline.mean())
    generator = np.random.default_rng(SEED)
    reconstruction_samples, partner_samples = [], []
    for _ in range(bootstrap):
        indices = generator.integers(0, len(baseline), len(baseline))
        denominator = baseline[indices].mean()
        reconstruction_samples.append(reconstruction_delta[indices].mean() / denominator)
        partner_samples.append(partner_delta[indices].mean() / denominator)
    return {
        "complexes": len(baseline),
        "mse": {key: float(value.mean()) for key, value in errors.items()},
        "reconstruction_gain": reconstruction,
        "reconstruction_gain_ci95": [float(np.quantile(reconstruction_samples, 0.025)),
                                      float(np.quantile(reconstruction_samples, 0.975))],
        "partner_gain": partner,
        "partner_gain_ci95": [float(np.quantile(partner_samples, 0.025)),
                               float(np.quantile(partner_samples, 0.975))],
    }


def evaluate(values: list[dict], calibrator) -> dict:
    train = [value for value in values if value["split"] == "train"]
    matrix = np.stack([value["teacher"].reshape(-1) for value in train])
    mean, scale = matrix.mean(axis=0), matrix.std(axis=0)
    active = scale > 1e-8
    safe_scale = np.where(active, scale, 1.0)
    report = {"basis_dimensions": len(mean), "active_dimensions": int(active.sum()), "splits": {}}
    for split in ("train", "val", "test"):
        split_values = [value for value in values if value["split"] == split]
        report["splits"][split] = _gain_summary(
            _complex_errors(split_values, mean, safe_scale, active, calibrator)
        )
    val, test = report["splits"]["val"], report["splits"]["test"]
    conditions = {
        "val_reconstruction_gain_positive": val["reconstruction_gain"] > 0,
        "val_partner_gain_positive": val["partner_gain"] > 0,
        "test_reconstruction_gain_ge_0_10": test["reconstruction_gain"] >= 0.10,
        "test_reconstruction_lcb_positive": test["reconstruction_gain_ci95"][0] > 0,
        "test_partner_gain_ge_0_10": test["partner_gain"] >= 0.10,
        "test_partner_lcb_positive": test["partner_gain_ci95"][0] > 0,
    }
    report["gate_conditions"] = conditions
    if not conditions["test_reconstruction_gain_ge_0_10"] or not conditions["test_reconstruction_lcb_positive"]:
        report["verdict"] = "RADIAL_BASIS_RECOVERY_NOT_IDENTIFIED"
    elif not all((conditions["val_partner_gain_positive"],
                  conditions["test_partner_gain_ge_0_10"],
                  conditions["test_partner_lcb_positive"])):
        report["verdict"] = "RADIAL_BASIS_PARTNER_DEPENDENCE_NOT_IDENTIFIED"
    elif all(conditions.values()):
        report["verdict"] = "RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED"
    else:
        report["verdict"] = "RADIAL_BASIS_RECOVERY_NOT_IDENTIFIED"
    return report, mean, safe_scale, active


def run(args) -> dict:
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("T-BASIS-R0 requires CUDA for frozen P1B inference")
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.mkdir(parents=True)
    records = [json.loads(line) for line in Path(args.records).read_text(encoding="utf-8").splitlines()]
    excluded = {
        row["source_entry_id"] for row in
        (json.loads(line) for line in Path(args.exclude_panel).read_text(encoding="utf-8").splitlines())
    }
    selected, selection_audit = select_panel(records, excluded)
    derangement, derangement_rows = build_derangement(selected)
    selection_rows = [{
        "source_entry_id": row["source_entry_id"], "source_split": row["source_split"],
        "pdb_id": row["pdb_id"], "homology_group_id": row["homology_group_id"],
        "sequence_sha256": row["sequence_sha256"], "murcko_scaffold": row["murcko_scaffold"],
        "selection_key": _selection_key(row),
    } for row in selected]
    _write_jsonl(output / "selection.jsonl", selection_rows)
    _write_jsonl(output / "derangement.jsonl", derangement_rows)
    selection_audit.update({
        "selection_sha256": sha256_file(output / "selection.jsonl"),
        "derangement_sha256": sha256_file(output / "derangement.jsonl"),
        "derangement_max_identity": max(row["identity"] for row in derangement_rows),
        "derangement_wrong_reuse": len(derangement_rows) - len({row["wrong_entry"] for row in derangement_rows}),
    })
    _write_json(output / "selection_audit.json", selection_audit)

    proteins = _load_protein_rows(
        Path(args.protein_bank), {row["sequence_sha256"] for row in selected}
    )
    protein_dim = len(next(iter(proteins.values()))["pooled"])
    model, checkpoint = _load_frozen_model(Path(args.checkpoint), protein_dim, args.device)
    ligands = torch.load(args.ligand_bank, map_location="cpu", weights_only=False)
    bin_moments = bin_rbf_expectation()
    values, failures = [], []
    for record in selected:
        try:
            basis = compute_record_basis(
                record, derangement[record["source_entry_id"]], model, proteins,
                ligands, bin_moments, args.device,
            )
            basis.update({"source_entry_id": record["source_entry_id"],
                          "split": record["source_split"]})
            values.append(basis)
        except Exception as error:
            failures.append({"source_entry_id": record["source_entry_id"],
                             "split": record["source_split"],
                             "error": f"{type(error).__name__}: {error}"})
    _write_jsonl(output / "failures.jsonl", failures)
    if failures:
        report = {"stage": STAGE, "verdict": "T_BASIS_DATA_OR_MAPPING_FAIL_CLOSED",
                  "selected": len(selected), "completed": len(values), "failures": failures}
        _write_json(output / "report.json", report)
        return report

    calibrator = fit_radial_calibrator([value for value in values if value["split"] == "train"])
    report, mean, scale, active = evaluate(values, calibrator)
    report.update({
        "stage": STAGE, "selected": len(selected), "completed": len(values),
        "affinity_label_reads": 0, "davis_label_reads": 0, "recipient_label_reads": 0,
        "teacher": {"atom_channels": list(ATOM_CHANNELS), "residue_channels": 6,
                    "radial_centers": CENTERS.tolist(), "sigma": RBF_SIGMA,
                    "cutoff": CUTOFF},
        "calibrator": {"type": "shared_radial_ridge", "alpha": 1e-3,
                       "coefficient_shape": list(calibrator.coef_.shape)},
    })
    _write_json(output / "report.json", report)
    np.savez_compressed(
        output / "basis_values.npz",
        entry=np.asarray([value["source_entry_id"] for value in values]),
        split=np.asarray([value["split"] for value in values]),
        teacher=np.stack([value["teacher"] for value in values]).astype(np.float32),
        raw_correct=np.stack([value["raw_correct"] for value in values]).astype(np.float32),
        raw_deranged=np.stack([value["raw_deranged"] for value in values]).astype(np.float32),
        calibration_coef=calibrator.coef_.astype(np.float64),
        calibration_intercept=calibrator.intercept_.astype(np.float64),
        train_mean=mean, train_scale=scale, active=active,
        bin_rbf_expectation=bin_moments,
    )
    manifest = {
        "stage": STAGE, "research_only": True,
        "inputs": {
            "records": {"sha256": sha256_file(args.records)},
            "excluded_panel": {"sha256": sha256_file(args.exclude_panel)},
            "checkpoint": {"sha256": sha256_file(args.checkpoint), "epoch": checkpoint["epoch"]},
            "protein_manifest": {"sha256": sha256_file(Path(args.protein_bank) / "manifest.json")},
            "ligand_bank": {"sha256": sha256_file(args.ligand_bank)},
        },
        "outputs": {
            name: sha256_file(output / name) for name in (
                "selection.jsonl", "selection_audit.json", "derangement.jsonl",
                "failures.jsonl", "basis_values.npz", "report.json",
            )
        },
        "forbidden_reads": {"affinity": 0, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", default="dataset/processed/open_structures/pilot20k_homology_split_v2/complexes.jsonl")
    parser.add_argument("--exclude-panel", default="research/e0_identifiability/artifacts/tdir_p0_v1/selection.jsonl")
    parser.add_argument("--protein-bank", default="dataset/processed/open_structures/pilot20k_esm2_t30_slots128_v1")
    parser.add_argument("--ligand-bank", default="dataset/processed/open_structures/pilot20k_mechanism_ligands_v1.pt")
    parser.add_argument("--checkpoint", default="report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/tbasis_r0_v1")
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
