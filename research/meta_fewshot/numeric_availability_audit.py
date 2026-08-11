"""Identity-only audit of numeric Ki and legacy-feature coverage for frozen O1."""
from __future__ import annotations

import gzip
import json
import subprocess
import time
from collections import Counter, defaultdict
from pathlib import Path

from research.meta_fewshot.fs_corpus_rebuild import (
    MIN_LIGANDS,
    PROJECTION,
    ROOT,
    admitted_rows,
    closure,
    mde_d,
    sha256,
)

EXACT_LABELS = (ROOT / "dataset" / "processed" / "crossed_interaction" /
                "bindingdb_202608" / "exact_labels.jsonl.gz")
LEGACY_CELLS = (ROOT / "dataset" / "processed" / "crossed_interaction" /
                "bindingdb_202608" / "cq_ki_corpus" / "cells.jsonl.gz")
FROZEN_SPLIT = ROOT / "report" / "meta_fewshot" / "FS_CORPUS_REBUILD_SPLIT.json"
OUTPUT = ROOT / "report" / "meta_fewshot" / "FS_NUMERIC_AVAILABILITY_AUDIT.json"


def identity_pairs(path: Path, endpoint: str | None = None) -> tuple[set[tuple[str, str]], Counter]:
    """Return pair identities and multiplicities; affinity values are never selected."""
    counts: Counter = Counter()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if endpoint is not None and row.get("endpoint") != endpoint:
                continue
            counts[(row["target_id"], row["ligand_id"])] += 1
    return set(counts), counts


def summarize(
    cells: list[dict],
    target_to_component: dict[str, str],
    source_component: str,
    evaluation_components: set[str],
    available_pairs: set[tuple[str, str]],
) -> dict:
    by_side_target: dict[str, dict[str, set[str]]] = {
        "source": defaultdict(set),
        "evaluation": defaultdict(set),
    }
    admitted_pairs = {(cell["target"], cell["ligand"]) for cell in cells}
    for target, ligand in admitted_pairs & available_pairs:
        component = target_to_component[target]
        if component == source_component:
            by_side_target["source"][target].add(ligand)
        elif component in evaluation_components:
            by_side_target["evaluation"][target].add(ligand)
    return {
        "admitted_pair_overlap": len(admitted_pairs & available_pairs),
        "source_k5_targets": sum(
            len(ligands) >= MIN_LIGANDS[5]
            for ligands in by_side_target["source"].values()
        ),
        "evaluation_k5_targets": sum(
            len(ligands) >= MIN_LIGANDS[5]
            for ligands in by_side_target["evaluation"].values()
        ),
    }


def run() -> dict:
    started = time.time()
    frozen = json.loads(FROZEN_SPLIT.read_text(encoding="utf-8"))
    if sha256(PROJECTION) != frozen["projection_sha256"]:
        raise ValueError("projection hash does not match the frozen split")

    cells, _ = admitted_rows(PROJECTION)
    components, _, _ = closure(cells)
    target_to_component = {
        cell["target"]: components.find(cell["target"])
        for cell in cells
    }
    observed_roots = set(target_to_component.values())
    split = frozen["split"]
    frozen_roots = {split["source_component"], *split["evaluation_components"]}
    if observed_roots != frozen_roots:
        raise ValueError("reconstructed dependency roots do not match the frozen split")

    exact_pairs, exact_counts = identity_pairs(EXACT_LABELS, endpoint="Ki")
    legacy_pairs, _ = identity_pairs(LEGACY_CELLS)
    numeric = summarize(
        cells, target_to_component, split["source_component"],
        set(split["evaluation_components"]), exact_pairs,
    )
    numeric_legacy = summarize(
        cells, target_to_component, split["source_component"],
        set(split["evaluation_components"]), exact_pairs & legacy_pairs,
    )
    admitted_pairs = {(cell["target"], cell["ligand"]) for cell in cells}
    result = {
        "schema": "MetaSieve.MetaFewshot.NumericAvailability.v1",
        "execution_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "audit_script_sha256": sha256(Path(__file__)),
        "projection_sha256": sha256(PROJECTION),
        "exact_labels_sha256": sha256(EXACT_LABELS),
        "legacy_cells_sha256": sha256(LEGACY_CELLS),
        "identity_fields_used": ["endpoint", "target_id", "ligand_id"],
        "affinity_values_used": 0,
        "aggregation_performed": False,
        "structural_cells": len(cells),
        "exact_numeric": numeric,
        "exact_numeric_and_legacy_feature_pair": numeric_legacy,
        "structural_pairs_without_exact_numeric": len(cells) - numeric["admitted_pair_overlap"],
        "exact_pairs_with_replicates": sum(
            exact_counts[pair] > 1 for pair in admitted_pairs & exact_pairs
        ),
        "numeric_evaluation_mde_d": mde_d(numeric["evaluation_k5_targets"]),
        "frozen_minimum_evaluation_targets": 30,
        "numeric_gate_pass": numeric["evaluation_k5_targets"] >= 30,
        "training_authorized": False,
        "terminal_verdict": "NUMERIC_FEWSHOT_CORPUS_GATE_REQUIRED",
        "elapsed_seconds": round(time.time() - started, 2),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
