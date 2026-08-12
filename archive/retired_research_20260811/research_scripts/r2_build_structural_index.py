"""Build a physically label-redacted structural index for R2-E2.

This one-time derivation reads the already-consumed main-v0 cell file.  The
downstream E2 audit reads only the redacted artifact and therefore never opens
or deserializes an affinity value itself.
"""
from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

from research.meta_fewshot.train_main_v0 import sha256

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
SOURCE = CORPUS / "cells.jsonl.gz"
OUTPUT = CORPUS / "r2_structural_index.jsonl.gz"
MANIFEST = ROOT / "report/meta_fewshot/r2_structural_index_manifest.json"
FIELDS = ("cell_id", "target_id", "ligand_id", "protein_group_40", "split", "panel_ids")


def build() -> dict:
    rows = []
    with gzip.open(SOURCE, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            source_row = json.loads(line)
            rows.append({field: source_row[field] for field in FIELDS})
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as text:
                for row in rows:
                    text.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    payload = {
        "schema": "MetaSieve.R2StructuralIndex.v1",
        "declared_role": "LABEL_REDACTION_FROM_ALREADY_CONSUMED_MAIN_V0_CORPUS",
        "source_affinity_values_deserialized_during_build": True,
        "downstream_affinity_fields_present": False,
        "fields": list(FIELDS),
        "rows": len(rows),
        "source_sha256": sha256(SOURCE),
        "structural_index_sha256": sha256(OUTPUT),
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")
    return payload


def main() -> int:
    print(json.dumps(build(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
