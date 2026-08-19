"""Stage P1/P2 episode bank builder (frozen bank, shared by all arms).

Inputs (read-only): P_SPLIT.json + corpus cells.jsonl.gz.
Rules: P1_BAKEOFF_PREREGISTRATION.md + ADDENDUM AD1
- one record per (split, target, draw, k) for every eligible k in
  {0,1,2,3,5,10,20,40}; eligibility: n_unique_ligands >= k + Q;
- ligand-unique rng ordering per (split, target, draw), keyed
  stage|pbank|SEED|split|target|draw (never by arm or k);
- support(k) = first k cells of that ordering; query = next Q cells after
  k (k=0 -> first Q cells);
- donor target = p_train target in a DIFFERENT cdhit40 cluster, drawn per
  (split, target, draw) with the same stream.
Output: artifacts/P_BANK.json + manifest. SHA-256 pinned.
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
SCHEMA = "MetaSieve.StageP.PBank.v1"
K_LIST = (0, 1, 2, 3, 5, 10, 20, 40)
Q = 8
DRAWS = 8
SEED = 20260861
AD1_SHA = "a675b0eca9d69d2e96d869da1ce61ce26cfd5b9bbb3ed1d30986edfa81f57a3b"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_rng(*parts):
    import numpy as np
    raw = "|".join(str(p) for p in parts).encode("utf-8")
    return np.random.default_rng(int(hashlib.sha256(raw).hexdigest()[:16], 16))


def main() -> int:
    split_art = json.loads((OUT / "P_SPLIT.json").read_text(encoding="utf-8"))
    cells_by_target: dict[tuple[str, str], list[dict]] = {}
    lig_by_target: dict[tuple[str, str], dict[str, list[str]]] = {}
    with gzip.open(CORPUS / "cells.jsonl.gz", "rt", encoding="utf-8") as f:
        for line in f:
            c = json.loads(line)
            split = split_art["cell_split"][c["cell_id"]]["split"]
            key = (split, c["target_id"])
            cells_by_target.setdefault(key, []).append(c)
            lig_by_target.setdefault(key, {}).setdefault(c["ligand_id"], []).append(
                c["cell_id"])
    train_targets = {}
    for (split, t), cells in cells_by_target.items():
        if split == "p_train":
            train_targets[t] = split_art["cell_split"][cells[0]["cell_id"]][
                "cluster"]
    records = []
    for split_name in ("p_val", "p_test"):
        targets = sorted({k[1] for k in cells_by_target if k[0] == split_name})
        for target in targets:
            key = (split_name, target)
            ligs = sorted(lig_by_target[key])
            cluster = split_art["cell_split"][
                cells_by_target[key][0]["cell_id"]]["cluster"]
            foreign_pool = [t for t, cl in train_targets.items() if cl != cluster]
            if not foreign_pool:
                raise ValueError(f"no foreign donor pool for {target}")
            for draw in range(DRAWS):
                rng = stable_rng("stageP", "pbank", SEED, split_name, target, draw)
                order = rng.permutation(len(ligs))
                ordered = [ligs[i] for i in order]
                donor = foreign_pool[int(rng.integers(len(foreign_pool)))]
                for k in K_LIST:
                    if len(ligs) < k + Q:
                        continue
                    support = [lig_by_target[key][l][0] for l in ordered[:k]]
                    query = [lig_by_target[key][l][0] for l in ordered[k:k + Q]]
                    records.append({
                        "split": split_name,
                        "cluster": cluster,
                        "target_id": target,
                        "draw": draw,
                        "k": k,
                        "support_cell_ids": support,
                        "query_cell_ids": query,
                        "donor_target_id": donor,
                        "donor_cluster": train_targets[donor],
                        "n_ligands": len(ligs),
                    })
    artifact = {
        "schema": SCHEMA,
        "split_artifact_sha256": sha256_file(OUT / "P_SPLIT.json"),
        "corpus_manifest_sha256": sha256_file(CORPUS / "manifest.json"),
        "addendum_ad1_sha256": AD1_SHA,
        "seed_key": f"stageP|pbank|{SEED}",
        "k_list": list(K_LIST),
        "Q": Q,
        "draws": DRAWS,
        "rule": "per (split,target,draw,k) records; eligibility n_ligands>=k+Q; "
                "support=first k, query=next Q; donor=p_train different cluster",
        "records": records,
    }
    text = json.dumps(artifact, indent=1, sort_keys=True)
    path = OUT / "P_BANK.json"
    path.write_text(text, encoding="utf-8")
    art_sha = sha256_file(path)
    (OUT / "P_BANK.manifest.json").write_text(json.dumps({
        "schema": SCHEMA + ".Manifest",
        "file": "P_BANK.json",
        "sha256": art_sha,
        "records": len(records),
        "splits": {"p_val": sum(r["split"] == "p_val" for r in records),
                   "p_test": sum(r["split"] == "p_test" for r in records)},
    }, indent=1), encoding="utf-8")
    print("wrote", path, "sha", art_sha, "records", len(records))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
