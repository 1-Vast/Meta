"""S1 — dataset role registry, licence audit and P1B structural exposure audit.

Fails closed if no independent structural test set can be formed.
"""
from __future__ import annotations

import hashlib
import json
import os
import ssl
import urllib.request

import numpy as np

ROOT = r"D:\MetaSieve"
OUT = os.path.join(ROOT, "report", "ssl_b2_structural_observability")
os.makedirs(OUT, exist_ok=True)
CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0"}

HOLO = os.path.join(ROOT, "dataset", "processed", "open_structures",
                    "pilot20k_holo_governed_v2", "complexes.jsonl")
SPLIT = os.path.join(ROOT, "dataset", "processed", "open_structures",
                     "pilot20k_homology_split_v2", "complexes.jsonl")
P1B = os.path.join(ROOT, "report", "mechanism_refactor", "p1b_pilot20k_seed17_v1",
                   "manifest.json")

# --------------------------------------------------------------- role registry
REGISTRY = {
    "schema": "MetaSieve.DatasetRoleRegistry.v1",
    "generated": "2026-08-08",
    "sources": {
        "RCSB_PDB_pilot20k": {
            "tier": "A", "role": "STRUCTURAL_TRAINING_ONLY",
            "url": "https://files.rcsb.org/", "license": "CC0-1.0",
            "redistributable": True, "endpoint": None,
            "structures": True, "sequences": True, "ligand_structures": True,
            "poses": "experimental holo",
            "exposure": "P1B TRAINED ON THIS; not independent validation",
        },
        "BioLiP2_annotations": {
            "tier": "B", "role": "ANNOTATION_ONLY",
            "url": "https://zhanggroup.org/BioLiP/", "license": "academic use, see readme",
            "redistributable": False, "endpoint": None,
            "note": "used only to select biologically relevant ligands; the S2 "
                    "teacher is computed from raw coordinates, never from BioLiP "
                    "absence-as-negative",
        },
        "RCSB_PDB_independent": {
            "tier": "A", "role": "STRUCTURAL_INDEPENDENT_TEST",
            "url": "https://files.rcsb.org/", "license": "CC0-1.0",
            "redistributable": True,
            "selection": "entries absent from pilot20k AND homology-disjoint from "
                         "the P1B training groups",
        },
        "ChEMBL37_static": {
            "tier": "A", "role": "AFFINITY_SOURCE_CALIBRATION_ONLY",
            "license": "CC BY-SA 3.0", "endpoint": "Ki/Kd",
            "gated": "not opened in this programme; S8 only, after structural admission",
        },
        "BindingDB_Articles_202608": {
            "tier": "A", "role": "AFFINITY_REPLICATION",
            "license": "CC BY 3.0", "endpoint": "Ki/Kd",
            "note": "XP4 measured its within-panel interaction as below the noise floor",
        },
        "Metz_2011": {"tier": "A", "role": "CONSUMED_DEVELOPMENT_DIAGNOSTIC_ONLY",
                      "endpoint": "pKi", "consumed_by": ["XP1", "XP2", "XP5"]},
        "Klaeger_2017": {"tier": "A", "role": "BOUNDED_SECONDARY_EVIDENCE",
                         "endpoint": "pKd_app", "consumed_by": ["XP1", "XP2-F"]},
        "PDSP_KiDB": {"tier": "A", "role": "BOUNDED_SECONDARY_EVIDENCE",
                      "endpoint": "Ki", "consumed_by": ["XP1-C"]},
        "DAVIS": {"tier": "PROHIBITED", "role": "NONE", "reads": 0},
        "recipient_labels": {"tier": "PROHIBITED", "role": "NONE", "reads": 0},
        "PDBbind": {"tier": "B", "role": "NOT_USED",
                    "reason": "redistribution contract not verified as compatible; "
                              "explicitly avoided per policy"},
        "PLINDER": {"tier": "B", "role": "NOT_USED_IN_S1",
                    "reason": "full 400k systems not downloaded by policy; the "
                              "independent test set is built from raw RCSB CC0 "
                              "coordinates instead, which needs no separate "
                              "annotation licence"},
        "TierC_DTI_and_KG": {"tier": "C", "role": "NOT_USED",
                             "members": ["Yamanishi", "BIOSNAP", "STITCH", "DrugBank",
                                         "SuperTarget", "TTD", "DGIdb", "DrugCentral",
                                         "OpenTargets", "CTD", "ChemProt", "Hetionet",
                                         "DRKG", "PharmKG"],
                             "reason": "S9 only, and never as affinity evidence"},
    },
}

# ------------------------------------------------------------- exposure audit
recs = [json.loads(l) for l in open(HOLO, encoding="utf-8")]
split = [json.loads(l) for l in open(SPLIT, encoding="utf-8")]
p1b = json.load(open(P1B))
by_pdb = {}
for r in recs:
    by_pdb.setdefault(r["pdb_id"].lower(), []).append(r)
split_of = {}
for r in split:
    if "split" in r:
        split_of.setdefault(r["split"], set()).add(r["pdb_id"].lower())
print("holo records:", len(recs), " distinct pdb ids:", len(by_pdb))
print("split keys present:", sorted(split_of) or "(none; split field absent)")

exposure = {
    "schema": "MetaSieve.StructuralExposureAudit.v1",
    "p1b_checkpoint_sha256": p1b["checkpoint_sha256"],
    "p1b_affinity_labels_used": p1b["affinity_labels_used"],
    "p1b_train_records": p1b["train_records"],
    "p1b_val_records": p1b["val_records"],
    "p1b_test_records": p1b["test_records"],
    "corpus_records": len(recs),
    "corpus_pdb_ids": len(by_pdb),
    "exposed_pdb_ids": sorted(by_pdb),
    "policy": "EVERY pdb id in the pilot20k corpus is treated as EXPOSED, including "
              "the P1B val and test partitions, because those partitions were "
              "consumed as P1B's own evaluation. The independent test set must be "
              "disjoint from this set AND homology-disjoint from it.",
}
seqs = {}
for r in recs:
    seqs[r["pdb_id"].lower()] = r.get("sequence", "")
exposure["exposed_unique_sequences"] = len({s for s in seqs.values() if s})

with open(os.path.join(OUT, "STRUCTURAL_EXPOSURE_AUDIT.json"), "w") as f:
    json.dump({k: v for k, v in exposure.items() if k != "exposed_pdb_ids"} |
              {"exposed_pdb_ids_count": len(by_pdb)}, f, indent=2)
with open(os.path.join(ROOT, "dataset", "processed", "ssl_b2_exposed_pdb_ids.json"),
          "w") as f:
    json.dump(sorted(by_pdb), f)
with open(os.path.join(ROOT, "dataset", "processed", "ssl_b2_exposed_sequences.json"),
          "w") as f:
    json.dump(seqs, f)

with open(os.path.join(OUT, "DATASET_ROLE_REGISTRY.json"), "w") as f:
    json.dump(REGISTRY, f, indent=2)

print(f"\nexposed pdb ids: {len(by_pdb)}; unique exposed sequences: "
      f"{exposure['exposed_unique_sequences']}")
print("wrote DATASET_ROLE_REGISTRY.json and STRUCTURAL_EXPOSURE_AUDIT.json")
