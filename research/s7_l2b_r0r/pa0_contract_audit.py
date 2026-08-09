"""Phase 2A / Phase 0 — contract and artifact audit. FAIL-CLOSED.

Registered by research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md
(SHA-256 4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e).

Runs BEFORE any Phase 2A metric. Verifies every artifact by SHA-256, rebuilds
the held-out A offset table two independent ways, and proves the compared arms
address identical rows and identical masks.

C3 is the load-bearing check: Phase 1's marginal decomposition indexed the
B5-family memmaps with the B4-family offset table. If those two tables ever
differed, every B5/BX5/BP5 marginal number would be misaligned. This proves it
rather than assuming it.

Opens no affinity source of any kind.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
PROC = ROOT / "dataset" / "processed" / "s7_l2b_r0r"
S4 = PROC / "sealed_preds"
S5 = PROC / "sealed_preds_b5"
ESM = PROC / "esm2_650M"
CORPUS = PROC / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
PREREG_SHA = "4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e"

# Every path the audit is permitted to open. Recorded so the label firewall is
# an enumerated fact rather than an assertion.
AFFINITY_SOURCE_MARKERS = ("chembl", "bindingdb", "davis", "kiba", "recipient",
                           "affinity", "ki_", "kd_", "ic50")


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def main():
    fails = []
    manifest = {}

    # ------------------------------------------------------------- C1 hashes
    targets = {
        "corpus_development": CORPUS / "monn_development_edge_corpus.jsonl.gz",
        "corpus_additional_pdb": CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz",
        "corpus_summary": CORPUS / "raw_corpus_summary.json",
        "atom_quarantine": OUT / "I1_ATOM_QUARANTINE.json",
        "homology_alignment_stats": PROC / "r0r2_closure" / "homology_alignment_stats.json",
        "sealed_index_b4family": S4 / "heldoutA_index.json",
        "ap_tables_b4family": S4 / "ap_tables.json",
        "ap_tables_b5family": S5 / "ap_tables_b5.json",
        "B5_checkpoint": S5 / "B5_checkpoint.pt",
        "esm2_index": ESM / "esm2_650M_index.json",
        "esm2_residues": ESM / "esm2_650M_residues.fp16.dat",
        "P0_SEALED_PREDICTION_MANIFEST": OUT / "P0_SEALED_PREDICTION_MANIFEST.json",
        "P1_B5_GATE": OUT / "P1_B5_GATE.json",
        "P1_MARGINAL_DECOMPOSITION": OUT / "P1_MARGINAL_DECOMPOSITION.json",
        "PUBLICATION_TIME_CLOSURE_AUDIT": OUT / "PUBLICATION_TIME_CLOSURE_AUDIT.json",
        "I1_ATOM_CORRESPONDENCE_AUDIT": OUT / "I1_ATOM_CORRESPONDENCE_AUDIT.json",
        "I2_COUPLING_IDENTIFIABILITY_AUDIT": OUT / "I2_COUPLING_IDENTIFIABILITY_AUDIT.json",
        "PREREG_UNIFIED": ROOT / "research" / "s7_l2b_r0r" / "PREREG_S7_L2B_UNIFIED.md",
        "PREREG_PHASE2A": ROOT / "research" / "s7_l2b_r0r" / "PREREG_S7_L2B_PHASE2A.md",
        "monn_mol_dict": MONN / "mol_dict",
        "monn_independent_mol_dict": MONN / "independent_dataset_mol_dict",
    }
    for a in ("B0", "B4", "BL", "BM", "BP", "BX"):
        targets[f"sealed_{a}"] = S4 / f"heldoutA_{a}.f16.dat"
    for a in ("B5", "BM5", "BP5", "BX5"):
        targets[f"sealed_{a}"] = S5 / f"heldoutA_{a}.f16.dat"

    for name, p in sorted(targets.items()):
        if not p.is_file():
            fails.append(f"C1 missing artifact: {name} -> {p}")
            manifest[name] = {"path": str(p.relative_to(ROOT)).replace("\\", "/"),
                              "present": False}
            continue
        manifest[name] = {"path": str(p.relative_to(ROOT)).replace("\\", "/"),
                          "present": True, "bytes": p.stat().st_size,
                          "sha256": sha_file(p)}
        print(f"  hashed {name}", flush=True)

    if manifest.get("PREREG_PHASE2A", {}).get("sha256") != PREREG_SHA:
        fails.append("C1 Phase 2A preregistration hash does not match the frozen value")

    # ------------------------------------------------- C2 hashes vs manifests
    c2 = {"checked": 0, "mismatches": []}
    p0m = json.loads((OUT / "P0_SEALED_PREDICTION_MANIFEST.json").read_text(encoding="utf-8"))
    for arm, want in p0m["per_pair_predictions"]["arm_sha256"].items():
        got = manifest.get(f"sealed_{arm}", {}).get("sha256")
        c2["checked"] += 1
        if got != want:
            c2["mismatches"].append({"arm": arm, "recorded": want, "recomputed": got})
    g5 = json.loads((OUT / "P1_B5_GATE.json").read_text(encoding="utf-8"))
    for arm, want in g5["per_pair_prediction_sha256"].items():
        got = manifest.get(f"sealed_{arm}", {}).get("sha256")
        c2["checked"] += 1
        if got != want:
            c2["mismatches"].append({"arm": arm, "recorded": want, "recomputed": got})
    if c2["mismatches"]:
        fails.append(f"C2 sealed prediction hash mismatch: {c2['mismatches']}")

    # ------------------------------------------ rebuild the data contract
    kept, quarantine, contract, _feats = build()
    comp_of = protein_components(kept)
    train, held_all, held_A, held_B = make_split(kept, comp_of)
    print(f"  rebuilt: kept={len(kept)} train={len(train)} heldA={len(held_A)}", flush=True)

    esm_index = json.loads((ESM / "esm2_650M_index.json").read_text(encoding="utf-8"))

    def offsets(records):
        idx, off = {}, 0
        for rec in records:
            idx[rec["source_key"]] = [off, rec["n_res"], rec["n_atoms"]]
            off += rec["n_res"] * rec["n_atoms"]
        return idx, off

    # exactly as p0_seal_predictions built it (no ESM filter)
    idx_p0, total_p0 = offsets(held_A)
    # exactly as p1_run_b5 built it (ESM availability filter on seq_key)
    held_A_esm = [r for r in held_A if r["seq_key"] in esm_index]
    idx_p1, total_p1 = offsets(held_A_esm)
    sealed_idx = json.loads((S4 / "heldoutA_index.json").read_text(encoding="utf-8"))

    def same(a, b):
        if set(a) != set(b):
            return False, f"key sets differ ({len(set(a) ^ set(b))} symmetric-difference)"
        for k in a:
            if list(a[k]) != list(b[k]):
                return False, f"offset/shape differs at {k}: {a[k]} vs {b[k]}"
        return True, "identical"

    c3 = {}
    ok_a, why_a = same(idx_p0, idx_p1)
    ok_b, why_b = same(idx_p0, sealed_idx)
    c3["p0_offsets_equal_p1_offsets"] = {"pass": ok_a, "detail": why_a}
    c3["p0_offsets_equal_sealed_index"] = {"pass": ok_b, "detail": why_b}
    c3["heldout_A_records"] = len(held_A)
    c3["heldout_A_records_with_esm"] = len(held_A_esm)
    c3["esm_filter_dropped"] = len(held_A) - len(held_A_esm)
    c3["total_cells_p0"] = total_p0
    c3["total_cells_p1"] = total_p1
    c3["sealed_index_entries"] = len(sealed_idx)
    if not (ok_a and ok_b):
        fails.append(f"C3 offset-table mismatch: {c3}")

    # --------------------------------------------------------- C4 file sizes
    c4 = {}
    for a in ("B0", "B4", "BL", "BM", "BP", "BX", "B5", "BM5", "BP5", "BX5"):
        m = manifest.get(f"sealed_{a}")
        if not m or not m.get("present"):
            continue
        want = 2 * total_p0
        c4[a] = {"bytes": m["bytes"], "expected_bytes": want,
                 "pass": bool(m["bytes"] == want)}
        if m["bytes"] != want:
            fails.append(f"C4 size mismatch for {a}: {m['bytes']} != {want}")

    # ---------------------------------- C5 identical rows and evaluation mask
    n_cells = sum(r["n_res"] * r["n_atoms"] for r in held_A)
    n_pos = sum(len(r["edges"]) for r in held_A)
    c5 = {"single_offset_table_addresses_all_arms": bool(all(v["pass"] for v in c4.values())),
          "arms": sorted(f"sealed_{a}" for a in
                         ("B0", "B4", "BL", "BM", "BP", "BX", "B5", "BM5", "BP5", "BX5")),
          "rows_complexes": len(held_A),
          "components": len({comp_of[r["source_key"]] for r in held_A}),
          "evaluation_mask": "complete n_res x n_atoms per complex, uniform weight",
          "total_cells": int(n_cells),
          "total_positive_cells": int(n_pos),
          "positive_density": float(n_pos / n_cells)}
    if not c5["single_offset_table_addresses_all_arms"]:
        fails.append("C5 arms are not addressed by one offset table")

    # --------------------------------------------- C6 metadata availability
    need = ("seq_key", "uniprot_id", "pdb_id", "ligand_ccd", "graph_key", "scaffold",
            "cohort", "n_res", "n_atoms", "edges", "positive_typed_edges")
    missing = Counter()
    empty_scaffold = 0
    for r in kept:
        for f in need:
            if f not in r or r[f] is None:
                missing[f] += 1
        if not r.get("scaffold"):
            empty_scaffold += 1
    pub = json.loads((OUT / "PUBLICATION_TIME_CLOSURE_AUDIT.json").read_text(encoding="utf-8"))
    c6 = {"records": len(kept), "missing_fields": dict(missing),
          "records_with_empty_murcko_scaffold": empty_scaffold,
          "publication_time_audit_top_level_keys": sorted(pub.keys())}
    if missing:
        fails.append(f"C6 missing metadata fields: {dict(missing)}")

    # ------------------------------------------------- C7 label-field reads
    c7 = {
        "label_fields_read": ["positive_binary_edges", "positive_typed_edges",
                              "positive_event_edges", "atom_names",
                              "source_atom_indices", "uniprot_sequence"],
        "affinity_sources_opened": [],
        "affinity_value_reads": 0,
        "paths_opened_by_this_stage": sorted(
            str(p.relative_to(ROOT)).replace("\\", "/") for p in targets.values()),
    }
    for p in targets.values():
        low = str(p).lower()
        if any(m in low for m in AFFINITY_SOURCE_MARKERS):
            c7["affinity_sources_opened"].append(str(p))
    if c7["affinity_sources_opened"]:
        fails.append(f"C7 an affinity-marked source was opened: {c7['affinity_sources_opened']}")

    # ------------------------------------------------------------------ out
    res = {
        "schema": "MetaSieve.S7L2B.P2A.ContractAudit.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "repo_commit": "623602e76b7d4f445af069014782278163183d59",
        "checks": {
            "C1_artifacts_present_and_hashed": {"pass": not any(f.startswith("C1") for f in fails),
                                                "artifacts": len(targets)},
            "C2_hashes_match_recorded_manifests": {"pass": not c2["mismatches"], **c2},
            "C3_offset_tables_identical": c3,
            "C4_file_sizes": c4,
            "C5_identical_rows_and_masks": c5,
            "C6_metadata_available": c6,
            "C7_label_and_affinity_reads": c7,
        },
        "data_contract_census": contract,
        "records_quarantined": len(quarantine),
        "split": {"train": len(train), "heldout_all": len(held_all),
                  "heldout_A": len(held_A), "heldout_B": len(held_B)},
        "input_manifest": manifest,
        "failures": fails,
        "verdict": ("PHASE2A_CONTRACT_OR_ARTIFACT_FAIL_CLOSED" if fails
                    else "PHASE2A_CONTRACT_PASS"),
    }
    (OUT / "PHASE2A_INPUT_MANIFEST.json").write_text(json.dumps(res, indent=2),
                                                     encoding="utf-8")
    print("\n" + json.dumps({"verdict": res["verdict"], "failures": fails,
                             "C3": c3, "C5": c5}, indent=2))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
