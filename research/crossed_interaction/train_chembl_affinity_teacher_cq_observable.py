"""Train a ChEMBL affinity-teacher observable on crossed Ki quotients.

F-149 tests one axis only: an external, affinity-aligned proteochemometric
ridge teacher is trained on ChEMBL37 EnergyPilot rows after exact BindingDB CQ
target and ligand-connectivity exclusion. The frozen teacher then gates the
same quotient-space BindingDB admission test used by prior CQ observables.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_cq_observable import (
    ARMS,
    OUT as CQ_OUT,
    QuotientBlock,
    additive_residual,
    bootstrap_contrast,
    fit_ridge,
    score_blocks,
)
from research.crossed_interaction.train_physchem_cq_observable import (
    ligand_pharmacophore_descriptor,
)
from research.crossed_interaction.train_seqchem_cq_observable import (
    donor_maps,
    product_feature,
    protein_descriptor,
    read_jsonl,
    read_jsonl_gz,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
SOURCE = ROOT / "dataset/processed/source_affinity/energy_pilot_v1"
GOVERNANCE = ROOT / "dataset/processed/source_affinity/energy_pilot_v1_governance"
OUT = CQ_OUT.parent / "chembl_affinity_teacher_cq_observable_gate1"
FEATURE_MODES = ("weighted_product", "prediction")
LABEL_MODES = ("affinity", "task_ligand_residual")


def read_jsonl_stream(path: Path):
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def bindingdb_ligand_connectivity_keys(corpus: Path) -> set[str]:
    keys = set()
    for row in read_jsonl(corpus / "ligands.jsonl"):
        drug_key = str(row["drug_key"])
        connectivity = drug_key.split("-", 1)[0]
        if len(connectivity) == 14:
            keys.add(connectivity)
    return keys


def bindingdb_target_sequence_keys(corpus: Path) -> set[str]:
    return {row["sequence_sha256"] for row in read_jsonl(corpus / "proteins.jsonl")}


def load_governed_task_ids(path: Path) -> set[str]:
    return {row["task_id"] for row in read_jsonl_stream(path)}


def fit_affinity_teacher_from_arrays(
        x: np.ndarray, y: np.ndarray, *, ridge: float) -> dict:
    if ridge <= 0:
        raise ValueError("teacher ridge must be strictly positive")
    if x.ndim != 2 or y.ndim != 1 or x.shape[0] != y.shape[0]:
        raise ValueError("teacher arrays have incompatible shapes")
    if x.shape[0] < 2:
        raise ValueError("at least two teacher rows are required")
    x_mean = x.mean(axis=0)
    x_scale = x.std(axis=0)
    x_scale[x_scale < 1e-6] = 1.0
    y_mean = float(y.mean())
    x_scaled = (x - x_mean) / x_scale
    y_centered = y - y_mean
    identity = np.eye(x_scaled.shape[1], dtype=np.float64)
    weights = np.linalg.solve(
        x_scaled.T @ x_scaled + ridge * identity,
        x_scaled.T @ y_centered)
    prediction = x_scaled @ weights + y_mean
    return {
        "ridge": ridge,
        "x_mean": x_mean,
        "x_scale": x_scale,
        "y_mean": y_mean,
        "weights": weights,
        "train_mse": float(np.square(y - prediction).mean()),
        "train_rows": int(x.shape[0]),
        "feature_dim": int(x.shape[1]),
        "p_affinity_mean": y_mean,
        "p_affinity_std": float(y.std()),
    }


def two_way_residual_labels(
        task_ids: list[str], ligand_ids: list[str], y: np.ndarray, *,
        iterations: int = 50, tolerance: float = 1e-10) -> np.ndarray:
    if len(task_ids) != len(ligand_ids) or len(task_ids) != y.shape[0]:
        raise ValueError("residual label inputs have incompatible lengths")
    intercept = float(y.mean())
    task_effect = {key: 0.0 for key in set(task_ids)}
    ligand_effect = {key: 0.0 for key in set(ligand_ids)}
    for _ in range(iterations):
        max_delta = 0.0
        task_sums = {key: 0.0 for key in task_effect}
        task_counts = {key: 0 for key in task_effect}
        for task, ligand, label in zip(task_ids, ligand_ids, y):
            task_sums[task] += float(label - intercept - ligand_effect[ligand])
            task_counts[task] += 1
        for task, value in task_sums.items():
            updated = value / task_counts[task]
            max_delta = max(max_delta, abs(updated - task_effect[task]))
            task_effect[task] = updated

        ligand_sums = {key: 0.0 for key in ligand_effect}
        ligand_counts = {key: 0 for key in ligand_effect}
        for task, ligand, label in zip(task_ids, ligand_ids, y):
            ligand_sums[ligand] += float(label - intercept - task_effect[task])
            ligand_counts[ligand] += 1
        for ligand, value in ligand_sums.items():
            updated = value / ligand_counts[ligand]
            max_delta = max(max_delta, abs(updated - ligand_effect[ligand]))
            ligand_effect[ligand] = updated
        if max_delta < tolerance:
            break
    residual = np.asarray([
        label - intercept - task_effect[task] - ligand_effect[ligand]
        for task, ligand, label in zip(task_ids, ligand_ids, y)
    ], dtype=np.float64)
    return residual


def predict_affinity_teacher(teacher: dict, x: np.ndarray) -> float:
    scaled = (x - teacher["x_mean"]) / teacher["x_scale"]
    return float(scaled @ teacher["weights"] + teacher["y_mean"])


def teacher_feature(product: np.ndarray, teacher: dict, *, mode: str) -> np.ndarray:
    if mode == "weighted_product":
        return product * teacher["weights"]
    if mode == "prediction":
        return np.asarray([predict_affinity_teacher(teacher, product)], dtype=np.float64)
    raise ValueError(f"unknown teacher feature mode: {mode}")


def load_external_teacher_training_arrays(
        *, source_dir: Path, governed_task_ids: set[str],
        blocked_targets: set[str], blocked_ligands: set[str],
        max_source_rows: int, label_mode: str = "affinity") -> tuple[np.ndarray, np.ndarray, dict]:
    if max_source_rows <= 0:
        raise ValueError("max_source_rows must be positive")
    if label_mode not in LABEL_MODES:
        raise ValueError(f"label_mode must be one of {LABEL_MODES}")
    protein_cache: dict[str, np.ndarray] = {}
    ligand_cache: dict[str, np.ndarray] = {}
    features = []
    labels = []
    task_labels = []
    ligand_labels = []
    counters = {
        "source_rows_scanned": 0,
        "source_rows_used": 0,
        "skipped_ungoverned_task": 0,
        "skipped_blocked_target": 0,
        "skipped_blocked_ligand": 0,
        "skipped_invalid_descriptor": 0,
        "blocked_bindingdb_targets": len(blocked_targets),
        "blocked_bindingdb_ligand_connectivities": len(blocked_ligands),
        "governed_tasks": len(governed_task_ids),
    }
    for row in read_jsonl_stream(source_dir / "canonical_rows.jsonl"):
        counters["source_rows_scanned"] += 1
        if row["task_id"] not in governed_task_ids:
            counters["skipped_ungoverned_task"] += 1
            continue
        protein_key = row["protein_sequence_sha256"]
        ligand_key = row["ligand_connectivity_key"]
        if protein_key in blocked_targets:
            counters["skipped_blocked_target"] += 1
            continue
        if ligand_key in blocked_ligands:
            counters["skipped_blocked_ligand"] += 1
            continue
        try:
            protein = protein_cache.get(protein_key)
            if protein is None:
                protein = protein_descriptor(row["protein_sequence"])
                protein_cache[protein_key] = protein
            ligand = ligand_cache.get(ligand_key)
            if ligand is None:
                ligand = ligand_pharmacophore_descriptor(row["canonical_smiles"])
                ligand_cache[ligand_key] = ligand
            label = float(row["p_affinity"])
            if not np.isfinite(label):
                raise ValueError("non-finite p_affinity")
            features.append(product_feature(protein, ligand))
            labels.append(label)
            task_labels.append(row["task_id"])
            ligand_labels.append(ligand_key)
        except (KeyError, TypeError, ValueError):
            counters["skipped_invalid_descriptor"] += 1
            continue
        counters["source_rows_used"] += 1
        if counters["source_rows_used"] >= max_source_rows:
            break
    if counters["source_rows_used"] < 10:
        raise ValueError("too few leakage-excluded ChEMBL rows for teacher training")
    x = np.stack(features).astype(np.float64)
    y_raw = np.asarray(labels, dtype=np.float64)
    y = (
        two_way_residual_labels(task_labels, ligand_labels, y_raw)
        if label_mode == "task_ligand_residual"
        else y_raw
    )
    counters["teacher_label_mode"] = label_mode
    counters["teacher_raw_p_affinity_std"] = float(y_raw.std())
    counters["teacher_training_label_std"] = float(y.std())
    counters["protein_descriptor_dim"] = int(protein_cache[next(iter(protein_cache))].shape[0])
    counters["ligand_descriptor_dim"] = int(ligand_cache[next(iter(ligand_cache))].shape[0])
    counters["teacher_raw_product_dim"] = int(x.shape[1])
    return x, y, counters


def train_external_teacher(
        *, corpus: Path, source_dir: Path, governance_dir: Path,
        max_source_rows: int, teacher_ridge: float,
        label_mode: str = "affinity") -> tuple[dict, dict]:
    blocked_targets = bindingdb_target_sequence_keys(corpus)
    blocked_ligands = bindingdb_ligand_connectivity_keys(corpus)
    governed_task_ids = load_governed_task_ids(governance_dir / "split_assignments.jsonl")
    x, y, counters = load_external_teacher_training_arrays(
        source_dir=source_dir,
        governed_task_ids=governed_task_ids,
        blocked_targets=blocked_targets,
        blocked_ligands=blocked_ligands,
        max_source_rows=max_source_rows,
        label_mode=label_mode)
    teacher = fit_affinity_teacher_from_arrays(x, y, ridge=teacher_ridge)
    metadata = {
        **counters,
        "source_schema": "MetaSieve.AffinityEnergyCorpus.v1",
        "source_manifest_sha256": sha256_file(source_dir / "corpus_manifest.json"),
        "governance_manifest_sha256": sha256_file(
            governance_dir / "governance_manifest.json"),
        "governance_split_sha256": sha256_file(governance_dir / "split_assignments.jsonl"),
        "teacher_train_mse": teacher["train_mse"],
        "teacher_p_affinity_mean": teacher["p_affinity_mean"],
        "teacher_p_affinity_std": teacher["p_affinity_std"],
    }
    return teacher, metadata


def materialize_features(
        corpus: Path, teacher: dict, *, feature_mode: str) -> tuple[
            dict[str, dict[str, np.ndarray]], dict]:
    if feature_mode not in FEATURE_MODES:
        raise ValueError(f"feature_mode must be one of {FEATURE_MODES}")
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins = {
        row["sequence_sha256"]: protein_descriptor(row["sequence"])
        for row in read_jsonl(corpus / "proteins.jsonl")
    }
    ligands = {
        row["drug_key"]: ligand_pharmacophore_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    protein_donor, ligand_donor = donor_maps(cells)
    features = {}
    for cell in cells:
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        features[cell["cell_id"]] = {
            "correct": teacher_feature(
                product_feature(proteins[target], ligands[ligand]),
                teacher, mode=feature_mode),
            "deranged_protein": teacher_feature(
                product_feature(proteins[protein_donor[target]], ligands[ligand]),
                teacher, mode=feature_mode),
            "foreign_ligand": teacher_feature(
                product_feature(proteins[target], ligands[ligand_donor[ligand]]),
                teacher, mode=feature_mode),
        }
    first = next(iter(features.values()))["correct"]
    metadata = {
        "feature_source": "chembl37_affinity_teacher_weighted_pcm_product",
        "feature_mode": feature_mode,
        "feature_dim": int(first.shape[0]),
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return features, metadata


def load_blocks(
        corpus: Path, teacher: dict, *, feature_mode: str) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    features, metadata = materialize_features(corpus, teacher, feature_mode=feature_mode)
    blocks = []
    max_orthogonality = 0.0
    for panel in panels:
        ordered = [cells[cell_id] for cell_id in panel["cell_ids"]]
        target_ids = [cell["target_id"] for cell in ordered]
        ligand_ids = [cell["ligand_id"] for cell in ordered]
        y_raw = np.asarray([cell["pK"] for cell in ordered], dtype=np.float64)
        y, retained_rank, y_orthogonality = additive_residual(target_ids, ligand_ids, y_raw)
        if retained_rank <= 0:
            continue
        arm_features = {}
        for arm in ARMS:
            raw = np.stack([features[cell["cell_id"]][arm] for cell in ordered])
            residual, feature_rank, feature_orthogonality = additive_residual(
                target_ids, ligand_ids, raw)
            if raw.ndim == 2 and residual.ndim == 1:
                residual = residual[:, None]
            if feature_rank != retained_rank:
                raise ValueError(f"feature rank mismatch in panel {panel['panel_id']}")
            arm_features[arm] = residual
            max_orthogonality = max(max_orthogonality, feature_orthogonality)
        max_orthogonality = max(max_orthogonality, y_orthogonality)
        blocks.append(QuotientBlock(
            panel_id=panel["panel_id"],
            split=panel["split"],
            dependency_component=panel["dependency_component"],
            retained_rank=retained_rank,
            y=np.asarray(y, dtype=np.float64),
            features=arm_features,
        ))
    metadata["blocks"] = len(blocks)
    metadata["max_projection_orthogonality"] = max_orthogonality
    return blocks, metadata


def run(
        corpus: Path = CORPUS, source_dir: Path = SOURCE,
        governance_dir: Path = GOVERNANCE, output: Path = OUT,
        ridge: float = 10000.0, teacher_ridge: float = 1000.0,
        max_source_rows: int = 20000, bootstrap_draws: int = 9999,
        seed: int = 20260812, feature_mode: str = "weighted_product",
        label_mode: str = "affinity") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    teacher, teacher_metadata = train_external_teacher(
        corpus=corpus, source_dir=source_dir, governance_dir=governance_dir,
        max_source_rows=max_source_rows, teacher_ridge=teacher_ridge,
        label_mode=label_mode)
    blocks, corpus_metadata = load_blocks(corpus, teacher, feature_mode=feature_mode)
    train_blocks = [block for block in blocks if block.split == "train"]
    development_blocks = [block for block in blocks if block.split == "development"]
    models = {arm: fit_ridge(train_blocks, arm, ridge) for arm in ARMS}
    train_rows, train_summary = score_blocks(train_blocks, models)
    development_rows, development_summary = score_blocks(development_blocks, models)
    controls = ("zero_additive", "deranged_protein", "foreign_ligand")
    contrasts = [
        bootstrap_contrast(
            development_rows, "correct", control,
            draws=bootstrap_draws, seed=seed + index)
        for index, control in enumerate(controls)
    ]
    gates = {
        "projection_orthogonality": corpus_metadata["max_projection_orthogonality"] <= 1e-7,
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "CHEMBL_AFFINITY_TEACHER_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "CHEMBL_AFFINITY_TEACHER_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.ChEMBLAffinityTeacherCQObservableGate1.v1",
        "hypothesis": (
            "A leakage-excluded ChEMBL37 affinity-aligned PCM ridge teacher can "
            "select protein-ligand product coordinates that beat additive, "
            "wrong-target, and wrong-ligand quotient controls on BindingDB CQ."),
        "literature_mechanism": {
            "proteochemometrics": (
                "joint target-ligand descriptors can transfer bioactivity "
                "structure across assays when leakage is controlled"),
            "deepdta_family": (
                "affinity-supervised sequence/ligand representations are a "
                "standard DTA transfer mechanism"),
            "cycle_space_gate": (
                "BindingDB CQ residuals remove target-only and ligand-only effects "
                "before positive-ridge admission"),
            "adambind_boundary": (
                "target-as-task meta-learning remains downstream and unchanged"),
        },
        "source_teacher": teacher_metadata,
        "corpus": corpus_metadata,
        "config": {
            "ridge": ridge,
            "teacher_ridge": teacher_ridge,
            "max_source_rows": max_source_rows,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "feature_mode": feature_mode,
            "label_mode": label_mode,
            "arms": list(ARMS),
            "train_split": "train",
            "evaluation_split": "development",
        },
        "train_summary": train_summary,
        "development_summary": development_summary,
        "development_contrasts": contrasts,
        "gates": gates,
        "development_training_authorized": verdict.endswith("PASS_DEVELOPMENT"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    (output / "development_panel_metrics.json").write_text(
        json.dumps(development_rows, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--source-dir", type=Path, default=SOURCE)
    parser.add_argument("--governance-dir", type=Path, default=GOVERNANCE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--teacher-ridge", type=float, default=1000.0)
    parser.add_argument("--max-source-rows", type=int, default=20000)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--feature-mode", choices=FEATURE_MODES, default="weighted_product")
    parser.add_argument("--label-mode", choices=LABEL_MODES, default="affinity")
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, source_dir=args.source_dir,
        governance_dir=args.governance_dir, output=args.output,
        ridge=args.ridge, teacher_ridge=args.teacher_ridge,
        max_source_rows=args.max_source_rows,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        feature_mode=args.feature_mode, label_mode=args.label_mode)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
