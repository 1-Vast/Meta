"""Train a target-function annotation observable on crossed Ki quotients.

This source-only Gate uses BindingDB curator target names as external
biological metadata. Target names are mapped to a fixed, small vocabulary of
functional classes; target ids and labels are never used as features.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
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
    read_jsonl,
    read_jsonl_gz,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
PROJECTION = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/metadata_projection.jsonl.gz"
OUT = CQ_OUT.parent / "function_context_cq_observable_gate1"

FUNCTION_KEYWORDS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("kinase", ("kinase", "tyrosine-protein", "serine/threonine")),
    ("gpcr_receptor", ("g-protein coupled", "gpcr", "adrenergic", "dopamine", "serotonin", "muscarinic")),
    ("nuclear_receptor", ("nuclear receptor", "estrogen receptor", "androgen receptor", "retinoic acid receptor")),
    ("ion_channel", ("channel", "transporter channel", "voltage-gated")),
    ("transporter", ("transporter", "uptake", "pump", "exchanger")),
    ("protease", ("protease", "proteinase", "peptidase", "caspase", "thrombin", "trypsin")),
    ("phosphatase", ("phosphatase",)),
    ("phosphodiesterase", ("phosphodiesterase", "pde")),
    ("polymerase", ("polymerase", "reverse transcriptase")),
    ("transferase", ("transferase", "acetyltransferase", "methyltransferase", "prenyltransferase")),
    ("hydrolase", ("hydrolase", "esterase", "lipase", "amylase")),
    ("oxidoreductase", ("oxidoreductase", "dehydrogenase", "reductase", "oxidase", "oxygenase")),
    ("epigenetic", ("histone", "bromodomain", "hdac", "deacetylase", "sirtuin")),
    ("immune_cytokine", ("interleukin", "cytokine", "chemokine", "integrin", "cd")),
    ("viral", ("virus", "viral", "hiv", "hepatitis", "influenza")),
    ("bacterial", ("bacterial", "mycobacterium", "staphylococcus", "streptococcus", "escherichia")),
    ("metabolic_enzyme", ("synthase", "synthetase", "lyase", "isomerase", "carboxylase")),
    ("unknown_or_other", ()),
)


def normalize_name(value: str) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", value.lower()))


def target_function_descriptor(names: list[str]) -> np.ndarray:
    text = " ; ".join(normalize_name(name) for name in names)
    descriptor = np.zeros(len(FUNCTION_KEYWORDS), dtype=np.float64)
    for index, (_, keywords) in enumerate(FUNCTION_KEYWORDS[:-1]):
        if any(keyword in text for keyword in keywords):
            descriptor[index] = 1.0
    if descriptor[:-1].sum() == 0:
        descriptor[-1] = 1.0
    return descriptor


def read_projection_names(path: Path, target_ids: set[str]) -> tuple[dict[str, list[str]], dict]:
    names: dict[str, set[str]] = {target: set() for target in target_ids}
    rows = 0
    matched_rows = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            rows += 1
            record = json.loads(line)
            target = record.get("target_sequence_sha256", "")
            if target in names:
                matched_rows += 1
                name = record.get("target_name", "").strip()
                if name:
                    names[target].add(name)
    result = {target: sorted(values) for target, values in names.items()}
    metadata = {
        "projection_rows_scanned": rows,
        "projection_rows_matched": matched_rows,
        "targets_with_curator_name": sum(1 for values in result.values() if values),
        "targets_without_curator_name": sum(1 for values in result.values() if not values),
        "metadata_projection_sha256": sha256_file(path),
    }
    return result, metadata


def materialize_features(
        corpus: Path, projection: Path,
        function_mode: str = "function_class") -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    if function_mode not in {"function_class", "unknown_only"}:
        raise ValueError(f"unknown function mode: {function_mode}")
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    target_ids = {cell["target_id"] for cell in cells}
    names_by_target, projection_metadata = read_projection_names(projection, target_ids)
    target_features = {}
    for target in sorted(target_ids):
        descriptor = target_function_descriptor(names_by_target[target])
        if function_mode == "unknown_only":
            descriptor = np.zeros_like(descriptor)
            descriptor[-1] = 1.0
        target_features[target] = descriptor
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
            "correct": product_feature(target_features[target], ligands[ligand]),
            "deranged_protein": product_feature(
                target_features[protein_donor[target]], ligands[ligand]),
            "foreign_ligand": product_feature(
                target_features[target], ligands[ligand_donor[ligand]]),
        }
    counts = np.sum(np.stack(list(target_features.values())), axis=0).astype(int).tolist()
    metadata = {
        **projection_metadata,
        "feature_source": "bindingdb_curator_target_function_class_x_ligand_estate",
        "function_mode": function_mode,
        "function_classes": [name for name, _ in FUNCTION_KEYWORDS],
        "function_class_target_counts": counts,
        "feature_dim": int(next(iter(features.values()))["correct"].shape[0]),
        "target_descriptor_dim": len(FUNCTION_KEYWORDS),
        "ligand_descriptor_dim": int(next(iter(ligands.values())).shape[0]),
        "cells": len(cells),
        "targets": len(target_ids),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return features, metadata


def load_blocks(
        corpus: Path, features: dict[str, dict[str, np.ndarray]]) -> tuple[list[QuotientBlock], float]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
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
    return blocks, max_orthogonality


def run(
        corpus: Path = CORPUS, projection: Path = PROJECTION,
        output: Path = OUT, ridge: float = 1000.0,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        function_mode: str = "function_class") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    features, metadata = materialize_features(corpus, projection, function_mode=function_mode)
    blocks, max_projection_orthogonality = load_blocks(corpus, features)
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
        "projection_orthogonality": max_projection_orthogonality <= 1e-7,
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "FUNCTION_CONTEXT_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "FUNCTION_CONTEXT_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.FunctionContextCQObservableGate1.v1",
        "hypothesis": (
            "Curator target-function annotations crossed with ligand E-state "
            "chemistry provide a transferable biological context for quotient "
            "interaction residuals."),
        "literature_mechanism": {
            "proteochemometrics": (
                "target annotations and ligand descriptors can form cross-target "
                "interaction features"),
            "domain_function_dti": (
                "protein function or family annotations are external biological "
                "context beyond sequence aggregation"),
            "hodge_cycle_space": (
                "final scoring removes target and ligand main effects with the "
                "same additive quotient Gate"),
        },
        "corpus": {
            **metadata,
            "blocks": len(blocks),
            "max_projection_orthogonality": max_projection_orthogonality,
        },
        "config": {
            "ridge": ridge,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
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
    parser.add_argument("--projection", type=Path, default=PROJECTION)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=1000.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument(
        "--function-mode",
        choices=("function_class", "unknown_only"),
        default="function_class")
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, projection=args.projection, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        function_mode=args.function_mode)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
