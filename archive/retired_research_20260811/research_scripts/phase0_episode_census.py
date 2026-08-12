"""Phase 0 label-blind few-shot episode feasibility census.

Reads the governed BindingDB Ki corpus metadata ONLY. The `pK` field is dropped
on load and never enters any statistic here; the census is a design audit, not a
measurement. Ki is the sole endpoint; no Kd, Kdapp, IC50, inhibition or
displacement value is touched.

Feasibility thresholds are declared below, before the census is run, so the
outcome cannot be produced by choosing a threshold afterwards.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = Path(__file__).resolve().parent
CORPUS = (ROOT / "dataset" / "processed" / "crossed_interaction" /
          "bindingdb_202608" / "cq_ki_corpus")
FEATURES = (ROOT / "dataset" / "processed" / "crossed_interaction" /
            "bindingdb_202608" / "cq_ki_tbasis_features.npz")
FEATURE_MANIFEST = (ROOT / "dataset" / "processed" / "crossed_interaction" /
                    "bindingdb_202608" / "cq_ki_tbasis_features.manifest.json")
OUT = ROOT / "report" / "meta_fewshot"

# ---- declared BEFORE running (Phase 0 feasibility contract) ----------------
K_VALUES = (1, 2, 3, 5)
MIN_QUERY_PER_EPISODE = 3          # a query set smaller than this is not a task
MIN_LIGANDS_FOR_K = {k: k + MIN_QUERY_PER_EPISODE for k in K_VALUES}
MIN_EVAL_TARGETS = 30              # held-out targets usable at the largest k
MIN_SOURCE_TARGETS = 100           # source targets usable at the largest k
MIN_EVAL_COMPONENTS = 5            # dependency components in the eval split
# Standardized minimum detectable effect at 80% power, one-sided 95%:
#   MDE_d = (1.645 + 0.842) / sqrt(N)
Z_SUM = 1.6448536269514722 + 0.8416212335729143
MAX_ACCEPTABLE_MDE_D = 0.60        # larger than this is an underpowered panel


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str),
                    encoding="utf-8")


LABEL_FIELDS = {"pK"}


def load_cells() -> list[dict]:
    """Load cell metadata with every label field stripped on read."""
    rows = []
    with gzip.open(CORPUS / "cells.jsonl.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            rows.append({k: v for k, v in row.items() if k not in LABEL_FIELDS})
    return rows


def mde_d(n_units: int) -> float:
    return float(Z_SUM / np.sqrt(n_units)) if n_units > 0 else float("inf")


def census_split(cells: list[dict], split: str) -> dict:
    rows = [c for c in cells if c["split"] == split]
    by_target = defaultdict(list)
    for row in rows:
        by_target[row["target_id"]].append(row)

    per_target = {}
    for target, items in by_target.items():
        ligands = {i["ligand_id"] for i in items}
        scaffolds = {i["scaffold"] for i in items}
        documents = {i["document_id"] for i in items}
        panels = {i["panel_id"] for i in items}
        groups = {i["protein_group_40"] for i in items}
        per_target[target] = {
            "cells": len(items), "ligands": len(ligands),
            "scaffolds": len(scaffolds), "documents": len(documents),
            "panels": len(panels), "protein_group_40": sorted(groups),
        }

    counts = np.array([v["ligands"] for v in per_target.values()])
    usable = {k: int((counts >= MIN_LIGANDS_FOR_K[k]).sum()) for k in K_VALUES}
    # scaffold-disjoint feasibility: need at least k+1 distinct scaffolds so a
    # support set can be drawn without sharing a scaffold with every query
    scaffold_counts = np.array([v["scaffolds"] for v in per_target.values()])
    usable_scaffold_disjoint = {
        k: int(((counts >= MIN_LIGANDS_FOR_K[k]) &
                (scaffold_counts >= k + 1)).sum()) for k in K_VALUES}

    return {
        "cells": len(rows),
        "targets": len(per_target),
        "distinct_ligands": len({r["ligand_id"] for r in rows}),
        "distinct_scaffolds": len({r["scaffold"] for r in rows}),
        "distinct_documents": len({r["document_id"] for r in rows}),
        "distinct_panels": len({r["panel_id"] for r in rows}),
        "distinct_protein_group_40": len({r["protein_group_40"] for r in rows}),
        "ligands_per_target": {
            "min": int(counts.min()) if counts.size else 0,
            "median": float(np.median(counts)) if counts.size else 0.0,
            "mean": float(counts.mean()) if counts.size else 0.0,
            "max": int(counts.max()) if counts.size else 0,
        },
        "targets_usable_at_k": usable,
        "targets_usable_at_k_scaffold_disjoint": usable_scaffold_disjoint,
        "per_target": per_target,
    }


def run() -> dict:
    started = time.time()
    manifest = json.loads((CORPUS / "manifest.json").read_text(encoding="utf-8"))
    feature_manifest = json.loads(FEATURE_MANIFEST.read_text(encoding="utf-8"))
    cells = load_cells()
    if len(cells) != manifest["cells"]:
        raise RuntimeError("corpus cell count does not match its manifest")

    train = census_split(cells, "train")
    development = census_split(cells, "development")

    # ---- leakage audits between the source and evaluation splits
    def keyset(split, field):
        return {c[field] for c in cells if c["split"] == split}

    leakage = {}
    for field, name in (("target_id", "target"), ("scaffold", "scaffold"),
                        ("document_id", "document"),
                        ("protein_group_40", "protein_homology_40"),
                        ("ligand_id", "ligand")):
        a, b = keyset("train", field), keyset("development", field)
        shared = a & b
        leakage[name] = {
            "train": len(a), "development": len(b), "shared": len(shared),
            "development_share_leaked": (len(shared) / len(b)) if b else 0.0,
        }

    # ---- dependency components per split, using the panel closure already frozen
    comps = defaultdict(set)
    for cell in cells:
        comps[cell["split"]].add(cell["panel_id"])
    eval_components = manifest["splits"]["development"]["components"]
    source_components = manifest["splits"]["train"]["components"]

    eval_targets_k5 = development["targets_usable_at_k"][5]
    power = {
        "primary_inference_unit": "held-out target",
        "secondary_inference_unit": "dependency component",
        "formula": "MDE_d = (z_0.95 + z_0.80) / sqrt(N)",
        "eval_targets_at_k5": eval_targets_k5,
        "mde_d_targets_k5": mde_d(eval_targets_k5),
        "eval_components": eval_components,
        "mde_d_components": mde_d(eval_components),
        "max_acceptable_mde_d": MAX_ACCEPTABLE_MDE_D,
    }

    checks = {
        "eval_targets_at_k5_sufficient": {
            "observed": eval_targets_k5, "required_at_least": MIN_EVAL_TARGETS,
            "pass": eval_targets_k5 >= MIN_EVAL_TARGETS},
        "source_targets_at_k5_sufficient": {
            "observed": train["targets_usable_at_k"][5],
            "required_at_least": MIN_SOURCE_TARGETS,
            "pass": train["targets_usable_at_k"][5] >= MIN_SOURCE_TARGETS},
        "eval_components_sufficient": {
            "observed": eval_components, "required_at_least": MIN_EVAL_COMPONENTS,
            "pass": eval_components >= MIN_EVAL_COMPONENTS},
        "target_leakage_zero": {
            "observed": leakage["target"]["shared"], "required": 0,
            "pass": leakage["target"]["shared"] == 0},
        "powered_at_primary_unit": {
            "observed": power["mde_d_targets_k5"],
            "required_at_most": MAX_ACCEPTABLE_MDE_D,
            "pass": power["mde_d_targets_k5"] <= MAX_ACCEPTABLE_MDE_D},
    }
    feasible = all(c["pass"] for c in checks.values())

    verdict = ("EPISODE_DESIGN_FEASIBLE_DEVELOPMENT_ONLY" if feasible
               else "FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE")

    result = {
        "schema": "MetaSieve.MetaFewshot.Phase0.Census.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "label_blind": True,
        "label_fields_stripped_on_read": sorted(LABEL_FIELDS),
        "endpoint": "Ki only; no Kd, Kdapp, IC50, inhibition or displacement",
        "corpus_manifest_sha256": sha_file(CORPUS / "manifest.json"),
        "corpus_files_sha256": manifest["files"],
        "features_sha256": feature_manifest["features_sha256"],
        "feature_dimensions": feature_manifest["dimensions"],
        "feature_arms": feature_manifest["arms"],
        "declared_thresholds": {
            "k_values": list(K_VALUES),
            "min_query_per_episode": MIN_QUERY_PER_EPISODE,
            "min_ligands_for_k": {str(k): v for k, v in MIN_LIGANDS_FOR_K.items()},
            "min_eval_targets": MIN_EVAL_TARGETS,
            "min_source_targets": MIN_SOURCE_TARGETS,
            "min_eval_components": MIN_EVAL_COMPONENTS,
            "max_acceptable_mde_d": MAX_ACCEPTABLE_MDE_D,
        },
        "source_split": {k: v for k, v in train.items() if k != "per_target"},
        "evaluation_split": {k: v for k, v in development.items()
                             if k != "per_target"},
        "source_components": source_components,
        "evaluation_components": eval_components,
        "largest_component_share": manifest["largest_component_share"],
        "leakage_audit": leakage,
        "power": power,
        "feasibility_checks": checks,
        "corpus_status": "DEVELOPMENT_ONLY_CLOSED_COMPONENT_CORPUS",
        "TERMINAL_VERDICT": verdict,
        "elapsed_seconds": round(time.time() - started, 2),
    }
    write_json(OUT / "PHASE0_EPISODE_CENSUS.json", result)
    write_json(OUT / "PHASE0_PER_TARGET.json",
               {"train": train["per_target"], "development": development["per_target"]})
    return result


def main(argv=None) -> int:
    argparse.ArgumentParser().parse_args(argv)
    result = run()
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {"corpus_files_sha256", "leakage_audit"}},
                     indent=2, default=str), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
