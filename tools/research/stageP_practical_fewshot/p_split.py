"""Stage P-line target-cold split builder (BindingDB-Ki governed corpus).

Read-only over dataset/processed/meta_fewshot/bindingdb_ki_main_v0.
Builds a NEW P-line split: cdhit40 clusters assigned p_train/p_val/p_test
60/20/20 by target count (cluster-balanced, greedy largest-first), whole
clusters move together => protein-component cold. Ligands are never split
by scaffold (P1/P2 allow same-series support/query).

Seed: SHA-256 stable. No Python hash(). Output:
  tools/research/stageP_practical_fewshot/artifacts/P_SPLIT.json
plus per-cell split and k-eligibility census.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset" / "processed" / "meta_fewshot" / "bindingdb_ki_main_v0"
HERE = Path(__file__).resolve().parent
OUT = HERE / "artifacts"
SCHEMA = "MetaSieve.StageP.PSplit.v1"
TRAIN_FRAC, VAL_FRAC, TEST_FRAC = 0.60, 0.20, 0.20


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_rng(*parts) -> "np.random.Generator":
    import numpy as np
    raw = "|".join(str(p) for p in parts).encode("utf-8")
    return np.random.default_rng(int(hashlib.sha256(raw).hexdigest()[:16], 16))


def load_cells():
    rows = []
    with gzip.open(CORPUS / "cells.jsonl.gz", "rt", encoding="utf-8") as f:
        for line in f:
            rows.append(json.loads(line))
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    corpus_manifest_sha = sha256_file(CORPUS / "manifest.json")
    cells = load_cells()
    by_cluster: dict[str, list[dict]] = {}
    for c in cells:
        by_cluster.setdefault(c["protein_group_40"], []).append(c)
    clusters = sorted(by_cluster)
    # target count per cluster
    cluster_targets = {
        cl: sorted({c["target_id"] for c in by_cluster[cl]}) for cl in clusters
    }
    rng = stable_rng("stageP", "psplit", "seed", 20260861)
    total_targets = sum(len(v) for v in cluster_targets.values())
    quotas = {"p_train": round(total_targets * TRAIN_FRAC),
              "p_val": round(total_targets * VAL_FRAC),
              "p_test": round(total_targets * TEST_FRAC)}
    buckets: dict[str, list[str]] = {"p_train": [], "p_val": [], "p_test": []}
    counts: dict[str, int] = {"p_train": 0, "p_val": 0, "p_test": 0}
    # seeded cluster permutation; assign each cluster to the most
    # under-filled bucket relative to its target-count quota
    for cl in rng.permutation(clusters):
        dst = min(buckets, key=lambda b: (counts[b] / quotas[b], b))
        buckets[dst].append(cl)
        counts[dst] += len(cluster_targets[cl])
    for name in buckets:
        frac = counts[name] / total_targets
        print(f"{name}: {len(buckets[name])} clusters, {counts[name]} targets "
              f"({frac:.3f})")
    split_of = {cl: name for name, cls in buckets.items() for cl in cls}
    # per-cell records
    cell_split = {}
    for c in cells:
        cell_split[c["cell_id"]] = {
            "split": split_of[c["protein_group_40"]],
            "target_id": c["target_id"],
            "ligand_id": c["ligand_id"],
            "cluster": c["protein_group_40"],
        }
    # k-eligibility census: unique ligands per target per split
    from collections import defaultdict
    lig = defaultdict(set)
    for c in cells:
        lig[(split_of[c["protein_group_40"]], c["target_id"])].add(c["ligand_id"])
    census = {}
    for split_name in ("p_train", "p_val", "p_test"):
        vals = sorted(
            (len(v) for (s, _), v in lig.items() if s == split_name),
            reverse=True)
        census[split_name] = {
            "targets": len(vals),
            "ligand_count_min": min(vals) if vals else 0,
            "ligand_count_median": vals[len(vals) // 2] if vals else 0,
            "ligand_count_max": max(vals) if vals else 0,
            "targets_usable_at_k": {
                str(k): sum(1 for v in vals if v >= k) for k in (1, 2, 3, 5, 10, 20, 40, 41)
            },
        }
    artifact = {
        "schema": SCHEMA,
        "corpus": str(CORPUS.relative_to(ROOT)),
        "corpus_manifest_sha256": corpus_manifest_sha,
        "seed_key": "stageP|psplit|seed|20260861",
        "fractions": {"p_train": TRAIN_FRAC, "p_val": VAL_FRAC, "p_test": TEST_FRAC},
        "rule": "cdhit40 clusters whole-cluster assignment, greedy largest-first "
                "by target count into the smallest bucket; 60/20/20 targets",
        "clusters": {name: buckets[name] for name in buckets},
        "cell_split": cell_split,
        "census": census,
    }
    text = json.dumps(artifact, indent=1, sort_keys=True)
    path = OUT / "P_SPLIT.json"
    path.write_text(text, encoding="utf-8")
    art_sha = sha256_file(path)
    (OUT / "P_SPLIT.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest",
        "file": "P_SPLIT.json",
        "sha256": art_sha,
        "corpus_manifest_sha256": corpus_manifest_sha,
        "cells": len(cells),
        "clusters_total": len(clusters),
    }, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
