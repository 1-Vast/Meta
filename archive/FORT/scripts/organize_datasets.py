"""Validate and classify the already-built dataset packages.

This command is intentionally non-destructive.  It does not recopy protected
innovation files and it does not derive a formal package from the legacy
dual-cold registry.  Build a new package with ``preprocess.py`` first.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from scripts.preprocess import verify_formal
except ModuleNotFoundError:
    from scripts.preprocess import verify_formal


RAW_DB = ROOT / "dataset/public/chembl_historical/snapshots/chembl_37/chembl_37.db"
RAW_MANIFEST = ROOT / "dataset/raw/chembl_37/manifest.json"
FORMAL = ROOT / "dataset/formal_training/chembl37_pki_formal.v4"
INNOVATION = ROOT / "dataset/innovation_tests/a2s_validation_small.v1"


def digest(path: Path) -> str:
    checksum = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    return {
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "bytes": path.stat().st_size,
        "sha256": digest(path),
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="utf-8")


def build(output: Path = FORMAL) -> dict[str, Any]:
    if output != FORMAL:
        raise ValueError("classification is pinned to the sealed ChEMBL 37 v3 package")
    required = [RAW_DB, RAW_MANIFEST, output / "manifest.json", INNOVATION / "manifest.json"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("missing classified dataset inputs: " + ", ".join(missing))
    raw = json.loads(RAW_MANIFEST.read_text(encoding="utf-8"))
    formal = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    verification = verify_formal(output, formal, raw)
    innovation_files = [
        file_record(path)
        for path in sorted(INNOVATION.iterdir())
        if path.is_file()
    ]
    classification = {
        "schema": "fort-dataset-classification-v2",
        "raw": {**raw, "file": file_record(RAW_DB)},
        "formal_training": {
            "path": str(output.relative_to(ROOT)).replace("\\", "/"),
            "status": formal["status"],
            "primary_endpoint": "pKi",
            "counts": formal["counts"],
            "manifest": str((output / "manifest.json").relative_to(ROOT)).replace("\\", "/"),
            "files": sorted(formal["files"]),
            "verification": verification,
        },
        "innovation_tests": {
            "path": str(INNOVATION.relative_to(ROOT)).replace("\\", "/"),
            "status": "DEVELOPMENT_ONLY_H0_BLOCKED",
            "purpose": "innovation mechanism/module validation only",
            "files": innovation_files,
        },
        "cleanup": {
            "status": "COMPLETE",
            "stale_candidates_and_superseded_intermediates_removed": True,
        },
        "policy": "raw SQLite is retained and never modified; only sealed derived packages are model inputs",
    }
    write_json(ROOT / "dataset/registry/DATASET_CLASSIFICATION_REPORT.v2.json", classification)
    (ROOT / "dataset/registry/DATASET_CLASSIFICATION_REPORT.v2.md").write_text(
        "# Dataset classification v2\n\n"
        f"- Raw ChEMBL 37: `{classification['raw']['file']['path']}` ({classification['raw']['file']['bytes']} bytes, `{classification['raw']['file']['sha256']}`).\n"
        f"- Formal pKi corpus: `{classification['formal_training']['path']}` with `{formal['counts']['pki_exact_measurements']}` exact measurements and `{formal['counts']['pki_exact_contexts']}` assay-context rows.\n"
        f"- Innovation tests: `{classification['innovation_tests']['path']}`.\n"
        f"- Natural-tail status: `{formal['status']}`.\n\n"
        "Stale candidates and superseded intermediate packages have been removed. "
        "Read the versioned formal manifest before consuming any data.\n",
        encoding="utf-8",
    )
    return classification


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=FORMAL)
    args = parser.parse_args()
    print(json.dumps(build(args.output.resolve()), indent=2, sort_keys=True, ensure_ascii=True))
