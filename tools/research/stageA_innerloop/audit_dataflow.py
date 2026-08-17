"""Stage 0 audit: trace the ten-step data flow and measure the adaptable scope.

No training, no meta_test. Every claim the preregistration makes about the
current pipeline is produced here as a number rather than asserted from a
reading of the source.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    resolve_architecture, training_label_scale,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
OUT = Path(__file__).resolve().parent / "AUDIT_DATAFLOW.json"

# Candidate adaptable subsets, smallest first. The screening experiment adapts
# exactly one of these; the rest are recorded so the choice is visible.
CANDIDATE_SCOPES = {
    "interaction_head_last": ("interaction_head.2.weight", "interaction_head.2.bias"),
    "interaction_head_last_bias_only": ("interaction_head.2.bias",),
    "interaction_head_full": (
        "interaction_head.0.weight", "interaction_head.0.bias",
        "interaction_head.2.weight", "interaction_head.2.bias"),
    "contact_weight": ("contact_weight.weight",),
    "interaction_head_last_plus_contact": (
        "interaction_head.2.weight", "interaction_head.2.bias",
        "contact_weight.weight"),
}


def episode_contract(data: QPSMPData, draws: int = 400) -> dict:
    """Step 2-3: one target per episode, unique and disjoint support/query."""
    rng = np.random.default_rng(20260816)
    violations = {"multi_target": 0, "cell_overlap": 0, "ligand_overlap": 0,
                  "support_duplicate_ligand": 0, "query_duplicate_ligand": 0}
    sizes = []
    for _ in range(draws):
        k = int(rng.integers(0, 6))
        spec = data.draw_episode("meta_train", k, 16, rng)
        rows = [data.cells[i] for i in (*spec.support, *spec.query)]
        if len({row["target_id"] for row in rows}) != 1:
            violations["multi_target"] += 1
        if set(spec.support) & set(spec.query):
            violations["cell_overlap"] += 1
        support = [data.cells[i]["ligand_id"] for i in spec.support]
        query = [data.cells[i]["ligand_id"] for i in spec.query]
        if set(support) & set(query):
            violations["ligand_overlap"] += 1
        if len(set(support)) != len(support):
            violations["support_duplicate_ligand"] += 1
        if len(set(query)) != len(query):
            violations["query_duplicate_ligand"] += 1
        sizes.append((len(spec.support), len(spec.query)))
    return {"draws": draws, "violations": violations,
            "support_sizes_seen": sorted({s for s, _ in sizes}),
            "query_size_range": [min(q for _, q in sizes),
                                 max(q for _, q in sizes)]}


def nested_bank_contract(data: QPSMPData) -> dict:
    """Step 10: are the evaluation banks nested and query-stable across k?"""
    banks = data.fixed_nested_episode_banks(
        "meta_val", (0, 1, 2, 3, 5), 16, 1, TrainConfig.evaluation_seed, 1)
    sizes = sorted(banks)
    index = {k: {(e.target, e.query): e for e in banks[k]} for k in sizes}
    keys = set(index[sizes[0]])
    nested, stable_query, disjoint = True, True, True
    for k in sizes:
        if set(index[k]) != keys:
            stable_query = False
    for key in keys:
        prefixes = [index[k][key].support for k in sizes]
        for smaller, larger in zip(prefixes, prefixes[1:]):
            if larger[:len(smaller)] != smaller:
                nested = False
        for k in sizes:
            episode = index[k][key]
            support_ligands = {data.cells[i]["ligand_id"] for i in episode.support}
            query_ligands = {data.cells[i]["ligand_id"] for i in episode.query}
            if support_ligands & query_ligands:
                disjoint = False
    return {"support_sizes": sizes,
            "episodes_per_k": {str(k): len(banks[k]) for k in sizes},
            "support_is_nested_prefix": nested,
            "query_panel_identical_across_k": stable_query,
            "support_query_ligand_disjoint": disjoint,
            "targets": len({e.target for e in banks[sizes[0]]}),
            "components": len({e.component for e in banks[sizes[0]]})}


def parameter_scopes(data: QPSMPData) -> dict:
    """Step 5: which parameters produce the zero-shot endpoint, and how many."""
    config = TrainConfig(arch="similarity_only")
    model = resolve_architecture(config.arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks,
        adapter_scale=config.adapter_scale, use_cartesian=config.use_cartesian)
    named = dict(model.named_parameters())
    total = sum(p.numel() for p in named.values() if p.requires_grad)
    branches: dict[str, int] = {}
    for name, parameter in named.items():
        if not parameter.requires_grad:
            continue
        branches[name.split(".")[0]] = (
            branches.get(name.split(".")[0], 0) + parameter.numel())
    scopes = {}
    for label, names in CANDIDATE_SCOPES.items():
        missing = [n for n in names if n not in named]
        scopes[label] = {
            "parameters": names, "missing": missing,
            "count": int(sum(named[n].numel() for n in names if n in named)),
        }
        scopes[label]["fraction_of_trainable"] = (
            scopes[label]["count"] / total if total else float("nan"))
    frozen = sorted(name for name, p in named.items() if not p.requires_grad)
    return {"trainable_total": int(total), "by_branch": branches,
            "frozen_parameters": frozen, "candidate_scopes": scopes}


def selection_and_leakage(data: QPSMPData) -> dict:
    """Steps 4, 9, 10: what the trainer reads, and from which split."""
    scale = training_label_scale(data)
    train_pk = [c["pK"] for c in data.cells if c["split"] == "meta_train"]
    val_pk = [c["pK"] for c in data.cells if c["split"] == "meta_val"]
    return {
        "label_scale_fitted_on": "meta_train",
        "label_scale": {"mean": scale.mean, "scale": scale.scale},
        "label_scale_excludes_meta_val": bool(
            abs(scale.mean - float(np.mean(train_pk))) < 1e-9),
        "cells": {"meta_train": len(train_pk), "meta_val": len(val_pk),
                  "meta_test_withheld": data.sealed_cell_count},
        "training_episode_split": "meta_train (scripts/train_qpsmp.py::train)",
        "checkpoint_selection_split": "meta_val",
        "checkpoint_selection_uses_meta_val_labels": True,
        "checkpoint_selection_note": (
            "train() evaluates val_banks built from meta_val every "
            "val_interval steps and keeps the state with the best admission "
            "score. meta_val labels therefore inform model selection, though "
            "never the training gradient. This is matched across all arms of "
            "this experiment, so it cannot manufacture a between-arm "
            "difference, but it does make every reported meta_val number an "
            "optimistic development estimate rather than a held-out one."),
        "meta_test": data.seal_record(),
    }


def main() -> int:
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    payload = {
        "schema": "MetaSieve.StageA.DataflowAudit.v1",
        "date": "2026-08-16",
        "split_directory": str(SPLIT.relative_to(ROOT)),
        "targets": {split: len(tasks) for split, tasks in data.tasks.items()},
        "components": {split: len(c) for split, c in data.components.items()},
        "episode_contract": episode_contract(data),
        "nested_bank_contract": nested_bank_contract(data),
        "parameter_scopes": parameter_scopes(data),
        "selection_and_leakage": selection_and_leakage(data),
    }
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(json.dumps(payload["episode_contract"], indent=1))
    print(json.dumps(payload["nested_bank_contract"], indent=1))
    print(json.dumps({k: v["count"] for k, v in
                      payload["parameter_scopes"]["candidate_scopes"].items()},
                     indent=1))
    print("trainable total:", payload["parameter_scopes"]["trainable_total"])
    print("by branch:", json.dumps(payload["parameter_scopes"]["by_branch"]))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
