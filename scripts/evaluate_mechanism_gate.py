"""Evaluate the frozen P1B bridge on structural test data and deranged controls."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random

import numpy as np
import torch
from sklearn.metrics import (average_precision_score, f1_score, precision_score,
                             recall_score, roc_auc_score)

from contracts.mechanism import DISTANCE_BINS_ANGSTROM
from scripts.data_contract import write_jsonl
from scripts.pretrain_mechanistic_bridge import (MechanismCorpus, MechanismPretrainer,
                                                 TrainConfig)
from scripts.structure_sources.rcsb import sha256_file


DISTANCE_BIN_CENTERS_ANGSTROM = (2.0, 5.0, 7.0, 9.0, 12.0)
CONTROL_SEED = 17
MIN_CONTROL_COVERAGE = 0.99
BOOTSTRAP_REPLICATES = 5000


def _stable_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _control_indices(corpus: MechanismCorpus, indices: list[int], seed: int,
                     pool_indices: list[int] | None = None
                     ) -> tuple[list[int], list[int], list[int], list[dict], list[dict]]:
    rng = random.Random(seed)
    pool_indices = list(range(len(corpus.records))) if pool_indices is None else pool_indices
    proteins_by_slots: dict[int, list[int]] = defaultdict(list)
    ligands_by_atoms: dict[int, list[int]] = defaultdict(list)
    for index in pool_indices:
        geometry = corpus.geometry(index)
        proteins_by_slots[int(geometry["residue_mask"].sum())].append(index)
        ligands_by_atoms[int(geometry["atom_mask"].sum())].append(index)

    eligible, protein_controls, ligand_controls, audit, excluded = [], [], [], [], []
    for index in indices:
        record = corpus.records[index]
        geometry = corpus.geometry(index)
        protein_candidates = [candidate for candidate in proteins_by_slots[
            int(geometry["residue_mask"].sum())]
            if corpus.records[candidate]["homology_group_id"] != record["homology_group_id"]]
        ligand_candidates = [candidate for candidate in ligands_by_atoms[
            int(geometry["atom_mask"].sum())]
            if corpus.records[candidate]["ccd_sha256"] != record["ccd_sha256"]]
        if not protein_candidates or not ligand_candidates:
            reasons = []
            if not protein_candidates:
                reasons.append("no exact-slot homology-disjoint protein")
            if not ligand_candidates:
                reasons.append("no exact-heavy-atom different-CCD ligand")
            excluded.append({"source_entry_id": record["source_entry_id"],
                             "reasons": reasons})
            continue
        protein_candidates.sort(key=lambda candidate: (
            abs(len(corpus.records[candidate]["sequence"]) - len(record["sequence"])),
            corpus.records[candidate]["source_entry_id"]))
        nearest_distance = abs(len(corpus.records[protein_candidates[0]]["sequence"])
                               - len(record["sequence"]))
        nearest = [candidate for candidate in protein_candidates
                   if abs(len(corpus.records[candidate]["sequence"])
                          - len(record["sequence"])) == nearest_distance]
        protein_control = rng.choice(nearest)
        ligand_control = rng.choice(sorted(ligand_candidates,
                                           key=lambda candidate: corpus.records[candidate][
                                               "source_entry_id"]))
        eligible.append(index)
        protein_controls.append(protein_control)
        ligand_controls.append(ligand_control)
        audit.append({
            "source_entry_id": record["source_entry_id"],
            "deranged_protein_entry_id": corpus.records[protein_control]["source_entry_id"],
            "deranged_protein_group_id": corpus.records[protein_control]["homology_group_id"],
            "deranged_protein_split": corpus.records[protein_control]["source_split"],
            "deranged_ligand_entry_id": corpus.records[ligand_control]["source_entry_id"],
            "deranged_ligand_ccd_sha256": corpus.records[ligand_control]["ccd_sha256"],
            "deranged_ligand_split": corpus.records[ligand_control]["source_split"],
        })
    return eligible, protein_controls, ligand_controls, audit, excluded


def _metrics(labels: np.ndarray, contact_prob: np.ndarray,
             distance_labels: np.ndarray, distance_prob: np.ndarray,
             top_l: list[tuple[np.ndarray, np.ndarray, int]]) -> dict:
    contact_pred = contact_prob >= 0.5
    centers = np.asarray(DISTANCE_BIN_CENTERS_ANGSTROM, dtype=np.float64)
    expected = distance_prob @ centers
    true_centers = centers[distance_labels]
    top_l_precision = []
    for sample_labels, sample_prob, count in top_l:
        selected = np.argpartition(sample_prob, -count)[-count:]
        top_l_precision.append(float(sample_labels[selected].mean()))
    clipped = np.clip(distance_prob, 1e-12, 1.0)
    return {
        "contact_auprc": float(average_precision_score(labels, contact_prob)),
        "contact_auroc": float(roc_auc_score(labels, contact_prob)),
        "contact_precision_at_top_l": float(np.mean(top_l_precision)),
        "contact_precision_0_5": float(precision_score(labels, contact_pred,
                                                        zero_division=0)),
        "contact_recall_0_5": float(recall_score(labels, contact_pred, zero_division=0)),
        "contact_f1_0_5": float(f1_score(labels, contact_pred, zero_division=0)),
        "distance_cross_entropy": float(-np.log(clipped[
            np.arange(len(distance_labels)), distance_labels]).mean()),
        "distance_bin_accuracy": float((distance_prob.argmax(axis=1)
                                         == distance_labels).mean()),
        "expected_distance_mae_angstrom": float(np.abs(expected - true_centers).mean()),
        "valid_pairs": int(len(labels)),
    }


def _evaluate_arm(model, corpus, label_indices, protein_indices, ligand_indices,
                  batch_size, device) -> dict:
    labels, contact_prob, distance_labels, distance_prob, top_l, per_complex = [], [], [], [], [], []
    centers = np.asarray(DISTANCE_BIN_CENTERS_ANGSTROM, dtype=np.float64)
    model.eval()
    with torch.inference_mode():
        for start in range(0, len(label_indices), batch_size):
            labels_batch = label_indices[start:start + batch_size]
            batch = corpus.batch(
                labels_batch, device,
                protein_indices=protein_indices[start:start + batch_size],
                ligand_indices=ligand_indices[start:start + batch_size])
            with torch.autocast(device_type="cuda", dtype=torch.float16,
                                enabled=device.startswith("cuda")):
                output = model(batch["X"], batch["A"], batch["atom_mask"],
                               batch["protein_pooled"], batch["protein_residues"],
                               batch["residue_mask"])
            probabilities = torch.sigmoid(output.contact_logits).float().cpu().numpy()
            distance_probabilities = torch.softmax(
                output.distance_logits, dim=-1).float().cpu().numpy()
            contact = batch["contact"].cpu().numpy()
            distance = batch["distance"].cpu().numpy()
            valid = output.pair_mask.bool().cpu().numpy()
            for row in range(len(labels_batch)):
                sample_labels = contact[row][valid[row]].astype(np.uint8)
                sample_prob = probabilities[row][valid[row]]
                labels.append(sample_labels)
                contact_prob.append(sample_prob)
                distance_labels.append(distance[row][valid[row]].astype(np.uint8))
                distance_prob.append(distance_probabilities[row][valid[row]])
                top_l.append((sample_labels, sample_prob,
                              int(batch["residue_mask"][row].sum().item())))
                sample_distance = distance_labels[-1]
                sample_distance_prob = distance_prob[-1]
                per_complex.append({
                    "source_entry_id": corpus.records[labels_batch[row]]["source_entry_id"],
                    "contact_positive_pairs": int(sample_labels.sum()),
                    "contact_auprc": float(average_precision_score(sample_labels, sample_prob)
                                           if sample_labels.any() else 0.0),
                    "expected_distance_mae_angstrom": float(np.abs(
                        sample_distance_prob @ centers - centers[sample_distance]).mean()),
                })
    metrics = _metrics(np.concatenate(labels), np.concatenate(contact_prob),
                       np.concatenate(distance_labels), np.concatenate(distance_prob), top_l)
    return metrics, per_complex


def _paired_bootstrap(correct: list[dict], control: list[dict], key: str, *,
                      higher_is_better: bool, seed: int) -> dict:
    correct_map = {row["source_entry_id"]: row[key] for row in correct}
    control_map = {row["source_entry_id"]: row[key] for row in control}
    if set(correct_map) != set(control_map):
        raise ValueError("P1B paired complex sets differ")
    ids = sorted(correct_map)
    differences = np.asarray([
        correct_map[value] - control_map[value] if higher_is_better
        else control_map[value] - correct_map[value] for value in ids])
    rng = np.random.default_rng(seed)
    bootstrap = np.empty(BOOTSTRAP_REPLICATES, dtype=np.float64)
    for index in range(BOOTSTRAP_REPLICATES):
        bootstrap[index] = differences[rng.integers(0, len(differences),
                                                    len(differences))].mean()
    return {"point": float(differences.mean()),
            "lower_95": float(np.quantile(bootstrap, 0.025)),
            "upper_95": float(np.quantile(bootstrap, 0.975)),
            "favorable_direction": "positive", "bootstrap_unit": "structure complex",
            "replicates": BOOTSTRAP_REPLICATES}


def evaluate_mechanism_gate(records_path: str | Path, supervision_dir: str | Path,
                            protein_bank_dir: str | Path, ligand_bank_path: str | Path,
                            checkpoint_path: str | Path, output_dir: str | Path, *,
                            batch_size: int = 8, device: str = "cuda",
                            max_records: int | None = None, seed: int = CONTROL_SEED) -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("P1B gate is registered for CUDA and fails closed without it")
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"P1B evaluation output already exists: {output}")
    output.mkdir(parents=True)
    corpus = MechanismCorpus(records_path, supervision_dir, protein_bank_dir, ligand_bank_path)
    indices = list(corpus.split_indices["test"])
    if max_records is not None:
        indices = indices[:max_records]
    if len(indices) < 2:
        raise ValueError("P1B evaluation requires at least two test records")
    total_test_records = len(indices)
    indices, protein_controls, ligand_controls, control_audit, excluded = _control_indices(
        corpus, indices, seed)
    coverage = len(indices) / total_test_records
    if len(indices) < 2:
        raise ValueError("fewer than two test records have exact-shape controls")
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    config = TrainConfig(**checkpoint["config"])
    model = MechanismPretrainer(int(checkpoint["protein_dim"]), config).to(device)
    model.load_state_dict(checkpoint["model_state"])

    evaluated = {
        "correct": _evaluate_arm(model, corpus, indices, indices, indices, batch_size, device),
        "deranged_protein": _evaluate_arm(model, corpus, indices, protein_controls, indices,
                                          batch_size, device),
        "deranged_ligand": _evaluate_arm(model, corpus, indices, indices, ligand_controls,
                                         batch_size, device),
    }
    arms = {name: value[0] for name, value in evaluated.items()}
    per_complex = {name: value[1] for name, value in evaluated.items()}
    zero_contact_complexes = sum(row["contact_positive_pairs"] == 0
                                 for row in per_complex["correct"])
    contrasts = {}
    for offset, control in enumerate(("deranged_protein", "deranged_ligand")):
        contrasts[f"correct_minus_{control}_contact_auprc"] = _paired_bootstrap(
            per_complex["correct"], per_complex[control], "contact_auprc",
            higher_is_better=True, seed=seed + offset)
        contrasts[f"correct_minus_{control}_distance_mae_improvement"] = _paired_bootstrap(
            per_complex["correct"], per_complex[control],
            "expected_distance_mae_angstrom", higher_is_better=False,
            seed=seed + 10 + offset)
    significant = all(value["lower_95"] > 0.0 for value in contrasts.values())
    passed = (coverage >= MIN_CONTROL_COVERAGE
              and significant
              and arms["correct"]["contact_auprc"]
              > arms["deranged_protein"]["contact_auprc"]
              and arms["correct"]["contact_auprc"]
              > arms["deranged_ligand"]["contact_auprc"]
              and arms["correct"]["expected_distance_mae_angstrom"]
              < arms["deranged_protein"]["expected_distance_mae_angstrom"]
              and arms["correct"]["expected_distance_mae_angstrom"]
              < arms["deranged_ligand"]["expected_distance_mae_angstrom"])
    write_jsonl(output / "control_mapping.jsonl", control_audit)
    write_jsonl(output / "control_exclusions.jsonl", excluded)
    write_jsonl(output / "per_complex_metrics.jsonl", [
        {"arm": arm, **row} for arm, values in per_complex.items() for row in values])
    report = {
        "schema": "MetaSieve.MechanismGeometryGate.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "gate": "P1B", "gate_status": "PASS" if passed else "FAIL",
        "primary_rule": "coverage >=0.99; correct global AUPRC/MAE beat both controls; all paired complex-bootstrap lower bounds >0",
        "records": len(indices), "test_records_total": total_test_records,
        "control_coverage": coverage, "minimum_control_coverage": MIN_CONTROL_COVERAGE,
        "control_exclusions": len(excluded), "split": "test", "seed": seed,
        "zero_positive_contact_complexes": zero_contact_complexes,
        "distance_bin_edges_angstrom": list(DISTANCE_BINS_ANGSTROM),
        "distance_bin_centers_angstrom": list(DISTANCE_BIN_CENTERS_ANGSTROM),
        "precision_at_top_l_policy": "L=number of valid mechanism residue slots per complex",
        "control_policy": {
            "candidate_pool": "all records in the governed structural corpus; labels and metrics remain test-only",
            "protein": "different registered homology group, exact residue-slot count, nearest sequence length",
            "ligand": "different CCD hash, exact heavy-atom count",
            "mapping_sha256": _stable_hash(control_audit),
        },
        "arms": arms, "paired_contrasts": contrasts,
        "inputs": {
            "records_sha256": sha256_file(records_path),
            "supervision_manifest_sha256": sha256_file(Path(supervision_dir) / "manifest.json"),
            "protein_manifest_sha256": sha256_file(Path(protein_bank_dir) / "manifest.json"),
            "ligand_bank_sha256": sha256_file(ligand_bank_path),
            "checkpoint_sha256": sha256_file(checkpoint_path),
        },
        "affinity_labels_used": False, "csmo_used": False,
    }
    (output / "gate_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records")
    parser.add_argument("supervision_dir")
    parser.add_argument("protein_bank_dir")
    parser.add_argument("ligand_bank")
    parser.add_argument("checkpoint")
    parser.add_argument("output")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-records", type=int)
    parser.add_argument("--seed", type=int, default=CONTROL_SEED)
    args = parser.parse_args()
    result = evaluate_mechanism_gate(
        args.records, args.supervision_dir, args.protein_bank_dir, args.ligand_bank,
        args.checkpoint, args.output, batch_size=args.batch_size, device=args.device,
        max_records=args.max_records, seed=args.seed)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
