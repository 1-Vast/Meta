"""Label-free exact-overlap audit for local R2 crossed-panel supplies."""
from __future__ import annotations

import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from research.meta_fewshot.train_main_v0 import sha256

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
REPORT = ROOT / "report/meta_fewshot/r2_crossed_overlap_audit.json"
SEQUENCE_FILES = {
    "BLK-BDB-PANELS": ROOT / "dataset/processed/multipanel/bdb_sequences.json",
    "PDSP-CORE": ROOT / "dataset/processed/crossed_panels/pdsp_sequences.json",
    "BLK-METZ-XP2": ROOT / "dataset/processed/crossed_panels_xp2/uniprot_seq_xp2.json",
}


def read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def exact_sequence_overlap(panel_path: Path, main_by_sequence: dict) -> dict:
    panel = json.loads(panel_path.read_text(encoding="utf-8"))
    sequences = set(panel.values())
    overlapping = sequences & set(main_by_sequence)
    by_split = Counter()
    for sequence in overlapping:
        by_split.update(main_by_sequence[sequence])
    return {
        "panel_unique_sequences": len(sequences),
        "exact_sequence_overlap": len(overlapping),
        "main_split_memberships": dict(sorted(by_split.items())),
    }


def normalize_pmid(value) -> str:
    text = str(value)
    return text[:-2] if text.endswith(".0") else text


def run() -> dict:
    proteins = read_jsonl(CORPUS / "proteins.jsonl")
    main_by_sequence = defaultdict(set)
    for row in proteins:
        main_by_sequence[row["sequence"]].add(row["split"])
    main_ligands = {row["smiles"] for row in read_jsonl(CORPUS / "ligands.jsonl")}

    main_pmids = set()
    with gzip.open(CORPUS / "r2_structural_index.jsonl.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            for panel_id in row["panel_ids"]:
                match = re.search(r"(?:pmid:)?(\d{7,9})", panel_id, flags=re.IGNORECASE)
                if match:
                    main_pmids.add(match.group(1))

    panels = {
        name: exact_sequence_overlap(path, main_by_sequence)
        for name, path in SEQUENCE_FILES.items()
    }
    import numpy as np
    bdb_path = ROOT / "dataset/processed/multipanel/blk_bdb_panels.npz"
    metz_path = ROOT / "dataset/processed/crossed_panels_xp2/blk_metz_xp2.npz"
    with np.load(bdb_path, allow_pickle=True) as stored:
        bdb_ligands = set(stored["smiles"].tolist())
        bdb_pmids = {normalize_pmid(value) for value in stored["pmid"].tolist()}
    with np.load(metz_path, allow_pickle=True) as stored:
        metz_ligands = set(stored["smiles"].tolist())
        pocket_count = len(stored["pocket"])
        pocket_lengths = sorted(set(map(len, stored["pocket"].tolist())))
    panels["BLK-BDB-PANELS"].update({
        "panel_unique_ligands": len(bdb_ligands),
        "exact_smiles_overlap": len(bdb_ligands & main_ligands),
        "panel_unique_pmids": len(bdb_pmids),
        "direct_pmid_overlap": len(bdb_pmids & main_pmids),
    })
    panels["BLK-METZ-XP2"].update({
        "panel_unique_ligands": len(metz_ligands),
        "exact_smiles_overlap": len(metz_ligands & main_ligands),
        "pocket_strings": pocket_count,
        "pocket_string_lengths": pocket_lengths,
        "existing_pocket_embedding_rows": 82,
    })
    payload = {
        "schema": "MetaSieve.R2CrossedOverlapAudit.v1",
        "declared_role": "LABEL_FREE_EXACT_OVERLAP_AND_ASSET_COVERAGE_AUDIT",
        "outcome_values_read": 0,
        "panels": panels,
        "fresh_confirmation_supply_identified": False,
        "verdict": "LOCAL_CROSSED_PANELS_OVERLAP_CONSUMED_MAIN_V0",
        "inputs": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                CORPUS / "proteins.jsonl",
                CORPUS / "ligands.jsonl",
                CORPUS / "r2_structural_index.jsonl.gz",
                *SEQUENCE_FILES.values(),
                bdb_path,
                metz_path,
            ]
        },
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    return payload


def main() -> int:
    print(json.dumps(run(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
