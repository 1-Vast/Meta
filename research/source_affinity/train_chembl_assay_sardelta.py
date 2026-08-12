"""Train a ChEMBL assay-level SAR-delta matched-pair observable.

F-151 tests one local chemistry axis only: within a single assay, same-scaffold
ligand pairs define delta supervision for a ridge model on ligand-delta
features. The evaluation is leave-one-assay-out, so no assay is used for both
pair construction and test scoring.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_cq_observable import bootstrap_contrast
from research.crossed_interaction.train_seqchem_cq_observable import ligand_descriptor, protein_descriptor


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "dataset/processed/source_affinity/chembl37_f0_rehydrated_v1/accepted_assays"
OUT = ROOT / "report/source_affinity/chembl_assay_sardelta_gate1"
PAIR_SIMILARITY = 0.50


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def stable_id(*parts: str) -> str:
    import hashlib
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()


def canonical_ligand_key(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    molecule = Chem.RemoveHs(molecule)
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def murcko_scaffold(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    scaffold = MurckoScaffold.GetScaffoldForMol(Chem.RemoveHs(molecule))
    return Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)


def assay_rows(path: Path) -> list[dict]:
    rows = read_jsonl(path)
    by_identity: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        key = str(row["standard_inchi_key"])
        by_identity[key].append(row)
    merged = []
    for key, values in sorted(by_identity.items()):
        first = values[0]
        merged.append({
            "assay_chembl_id": first["assay_chembl_id"],
            "document_chembl_id": first["document_chembl_id"],
            "target_chembl_id": first["target_chembl_id"],
            "target_accession": first["target_accession"],
            "protein_sequence": first.get("protein_sequence", ""),
            "ligand_connectivity_key": str(first["ligand_connectivity_key"]),
            "standard_inchi_key": key,
            "canonical_smiles": first["canonical_smiles"],
            "p_value": float(np.median([float(v["p_value"]) for v in values])),
            "endpoint_family": first["endpoint_family"],
            "source_row_ids": [int(v["activity_id"]) for v in values],
            "scaffold": murcko_scaffold(first["canonical_smiles"]),
        })
    return merged


def ligand_features(smiles: str) -> np.ndarray:
    return ligand_descriptor(smiles).astype(np.float64)


def pair_feature(left: np.ndarray, right: np.ndarray, *, target: np.ndarray | None,
                 mode: str) -> np.ndarray:
    if mode == "delta":
        return left - right
    if mode == "delta_target":
        if target is None:
            raise ValueError("target descriptor is required for delta_target")
        return np.concatenate([target, left - right])
    if mode == "concat":
        return np.concatenate([left, right])
    raise ValueError(f"unknown pair feature mode: {mode}")


def build_matched_pairs(rows: list[dict], *, feature_mode: str,
                        max_pairs_per_scaffold: int | None = None) -> list[dict]:
    rows_by_scaffold: dict[str, list[dict]] = defaultdict(list)
    feature_cache = {}
    target_cache = {}
    for row in rows:
        feature_cache[row["ligand_connectivity_key"]] = ligand_features(row["canonical_smiles"])
        if row.get("protein_sequence"):
            target_cache[row["assay_chembl_id"]] = protein_descriptor(row["protein_sequence"])
        row.setdefault("standard_inchi_key", row.get("ligand_connectivity_key", ""))
        rows_by_scaffold[row["scaffold"]].append(row)
    pairs = []
    for scaffold, scaffold_rows in sorted(rows_by_scaffold.items()):
        if len(scaffold_rows) < 2:
            continue
        ordered = sorted(scaffold_rows, key=lambda row: (row["p_value"], row["ligand_connectivity_key"]))
        seen = set()
        candidate_indices = list(range(len(ordered)))
        scaffold_pairs = []
        for i in range(len(candidate_indices)):
            for j in range(i + 1, len(candidate_indices)):
                left = ordered[i]
                right = ordered[j]
                pair_key = tuple(sorted((left.get("standard_inchi_key", left["ligand_connectivity_key"]),
                                         right.get("standard_inchi_key", right["ligand_connectivity_key"]))))
                if pair_key in seen:
                    continue
                left_fp = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(left["canonical_smiles"]), 2, 1024)
                right_fp = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(right["canonical_smiles"]), 2, 1024)
                similarity = DataStructs.TanimotoSimilarity(left_fp, right_fp)
                if similarity < PAIR_SIMILARITY:
                    continue
                seen.add(pair_key)
                scaffold_pairs.append({
                    "assay_chembl_id": left["assay_chembl_id"],
                    "document_chembl_id": left["document_chembl_id"],
                    "target_chembl_id": left["target_chembl_id"],
                    "target_accession": left["target_accession"],
                    "scaffold": scaffold,
                    "left_key": left["ligand_connectivity_key"],
                    "right_key": right["ligand_connectivity_key"],
                    "left_inchi": left["standard_inchi_key"],
                    "right_inchi": right["standard_inchi_key"],
                    "left_smiles": left["canonical_smiles"],
                    "right_smiles": right["canonical_smiles"],
                    "left_p_value": float(left["p_value"]),
                    "right_p_value": float(right["p_value"]),
                    "delta_p_value": float(left["p_value"] - right["p_value"]),
                    "feature": pair_feature(
                        feature_cache[left["ligand_connectivity_key"]],
                        feature_cache[right["ligand_connectivity_key"]],
                        target=target_cache.get(left["assay_chembl_id"]),
                        mode=feature_mode),
                })
        if max_pairs_per_scaffold is not None and len(scaffold_pairs) > max_pairs_per_scaffold:
            scaffold_pairs = scaffold_pairs[:max_pairs_per_scaffold]
        pairs.extend(scaffold_pairs)
    return pairs


def fit_ridge(x: np.ndarray, y: np.ndarray, ridge: float) -> dict:
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    mean = x.mean(axis=0)
    scale = x.std(axis=0)
    scale[scale < 1e-6] = 1.0
    x_scaled = (x - mean) / scale
    y_mean = float(y.mean())
    y_centered = y - y_mean
    identity = np.eye(x_scaled.shape[1], dtype=np.float64)
    weights = np.linalg.solve(x_scaled.T @ x_scaled + ridge * identity, x_scaled.T @ y_centered)
    prediction = x_scaled @ weights + y_mean
    return {
        "mean": mean,
        "scale": scale,
        "y_mean": y_mean,
        "weights": weights,
        "ridge": ridge,
        "train_mse": float(np.square(y - prediction).mean()),
        "feature_dim": int(x.shape[1]),
    }


def predict(model: dict, x: np.ndarray) -> np.ndarray:
    return ((x - model["mean"]) / model["scale"]) @ model["weights"] + model["y_mean"]


def leave_one_assay_out(assays: dict[str, list[dict]], *, feature_mode: str,
                        ridge: float, max_pairs_per_scaffold: int | None = None) -> tuple[list[dict], dict]:
    assay_ids = sorted(assays)
    rows = []
    source_metadata = {
        "assays": len(assays),
        "assay_ids": assay_ids,
        "pair_similarity_threshold": PAIR_SIMILARITY,
        "feature_mode": feature_mode,
        "ridge": ridge,
        "source_manifest_sha256": sha256_file(SOURCE / "CHEMBL1000360.jsonl"),
    }
    for held_out in assay_ids:
        train_pairs = []
        for assay_id, assay_rows_ in assays.items():
            if assay_id == held_out:
                continue
            train_pairs.extend(build_matched_pairs(
                assay_rows_, feature_mode=feature_mode,
                max_pairs_per_scaffold=max_pairs_per_scaffold))
        test_pairs = build_matched_pairs(
            assays[held_out], feature_mode=feature_mode,
            max_pairs_per_scaffold=max_pairs_per_scaffold)
        if len(train_pairs) < 2 or len(test_pairs) < 2:
            continue
        x_train = np.stack([pair["feature"] for pair in train_pairs]).astype(np.float64)
        y_train = np.asarray([pair["delta_p_value"] for pair in train_pairs], dtype=np.float64)
        model = fit_ridge(x_train, y_train, ridge)
        x_test = np.stack([pair["feature"] for pair in test_pairs]).astype(np.float64)
        y_test = np.asarray([pair["delta_p_value"] for pair in test_pairs], dtype=np.float64)
        pred = predict(model, x_test)
        zero_pred = np.zeros_like(pred)
        for pair, true, estimate, control in zip(test_pairs, y_test, pred, zero_pred):
            rows.append({
                "held_out_assay": held_out,
                "scaffold": pair["scaffold"],
                "pair_id": stable_id(held_out, pair["left_key"], pair["right_key"]),
                "delta_p_value": float(true),
                "prediction": float(estimate),
                "zero_delta": float(control),
                "squared_error": float((true - estimate) ** 2),
                "zero_squared_error": float((true - control) ** 2),
            })
    return rows, source_metadata


def summarize(rows: list[dict], arm: str) -> dict[str, float]:
    return {
        "pairs": len(rows),
        "mse": float(np.mean([row["squared_error"] if arm == "correct"
                              else row["zero_squared_error"] for row in rows])),
    }


def component_metric(rows: list[dict], arm: str) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        grouped[row["held_out_assay"]].append(
            row["squared_error"] if arm == "correct" else row["zero_squared_error"])
    return {key: float(np.mean(values)) for key, values in sorted(grouped.items())}


def run(
        source_dir: Path = SOURCE, output: Path = OUT, ridge: float = 100.0,
        feature_mode: str = "delta_target", max_pairs_per_scaffold: int | None = None,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    assays = {
        path.stem: assay_rows(path)
        for path in sorted(source_dir.glob("*.jsonl"))
    }
    rows, metadata = leave_one_assay_out(
        assays, feature_mode=feature_mode, ridge=ridge,
        max_pairs_per_scaffold=max_pairs_per_scaffold)
    if len(rows) < 2:
        raise ValueError("insufficient matched-pair rows for evaluation")
    correct = component_metric(rows, "correct")
    zero = component_metric(rows, "zero")
    assays_shared = sorted(set(correct) & set(zero))
    delta = np.asarray([zero[key] - correct[key] for key in assays_shared], dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = delta[rng.integers(0, len(delta), size=(bootstrap_draws, len(delta)))].mean(axis=1)
    lcb = float(np.quantile(samples, 0.05))
    gate = {
        "correct_beats_zero_delta": bool(lcb > 0.0),
        "assays_ge_3": len(assays_shared) >= 3,
    }
    verdict = "CHEMBL_ASSAY_SARDELTA_GATE1_PASS" if all(gate.values()) else "CHEMBL_ASSAY_SARDELTA_GATE1_FAIL_CLOSED"
    result = {
        "schema": "MetaSieve.ChEMBLAssaySARDeltaGate1.v1",
        "hypothesis": (
            "Within-assay same-scaffold matched-pair delta supervision can "
            "learn a transferable SAR residual better than zero delta."),
        "literature_mechanism": {
            "mmpa": "matched molecular pairs expose local ligand transformations",
            "activity_cliff": "small structural changes can induce large potency deltas",
            "fsmol_boundary": "each assay is a target-as-task split unit",
        },
        "source": metadata,
        "config": {
            "ridge": ridge,
            "feature_mode": feature_mode,
            "max_pairs_per_scaffold": max_pairs_per_scaffold,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
        },
        "development_summary": {
            "correct": summarize(rows, "correct"),
            "zero": summarize(rows, "zero"),
        },
        "development_contrast": {
            "control": "zero_delta",
            "component_macro_reduction": float(delta.mean()),
            "one_sided_95_lcb": lcb,
            "pass": bool(lcb > 0.0),
        },
        "gates": gate,
        "development_training_authorized": verdict.endswith("PASS"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    (output / "development_rows.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=100.0)
    parser.add_argument(
        "--feature-mode", choices=("delta", "delta_target", "concat"),
        default="delta_target")
    parser.add_argument("--max-pairs-per-scaffold", type=int, default=None)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        source_dir=args.source_dir, output=args.output, ridge=args.ridge,
        feature_mode=args.feature_mode,
        max_pairs_per_scaffold=args.max_pairs_per_scaffold,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
