"""Pre-aggregation measurement provenance, restricted to `meta_train`.

The governed corpus stores one aggregated `pK` per (target, ligand) cell. The
pre-aggregation rows survive in two artifacts whose SHA-256 digests are recorded
in the corpus manifest itself, so their authority is verifiable rather than
assumed:

    exact_labels.jsonl.gz        -> manifest["labels_sha256"]
    metadata_projection.jsonl.gz -> manifest["projection_sha256"]

**Disclosed isolation exception.** Those two files are single all-label
artifacts. The allow-list of `source_row_id`s is derived from the physically
isolated split view, `meta_train` only; the artifacts are then streamed and a
row is retained *in the same pass* only if its id is allow-listed. A
non-allow-listed value is never bound to a retained structure. This is **logical
exclusion after parsing**, strictly weaker than the model path's physical
isolation, and it is reported that way. The sealed and development-validation
rows are both excluded by the allow-list.

The filtered result is cached once under this stage; every downstream module
reads only the cache.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
SOURCE = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608"
EXACT_LABELS = SOURCE / "exact_labels.jsonl.gz"
METADATA_PROJECTION = SOURCE / "metadata_projection.jsonl.gz"
SPLIT_VIEW = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views"
TRAIN_CELLS = SPLIT_VIEW / "meta_train/cells.jsonl.gz"

HERE = Path(__file__).resolve().parent
CACHE = HERE / "U0_PROVENANCE_meta_train.jsonl.gz"

ENDPOINT = "Ki"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_source_authority() -> dict:
    manifest = json.loads((CORPUS / "manifest.json").read_text(encoding="utf-8"))
    actual = {
        "exact_labels": file_sha256(EXACT_LABELS),
        "metadata_projection": file_sha256(METADATA_PROJECTION),
    }
    expected = {
        "exact_labels": manifest["labels_sha256"],
        "metadata_projection": manifest["projection_sha256"],
    }
    for key, value in expected.items():
        if actual[key] != value:
            raise ValueError(f"{key} does not match the corpus manifest")
    return {
        "verified": True,
        "expected": expected,
        "actual": actual,
        "corpus_manifest_sha256": file_sha256(CORPUS / "manifest.json"),
        "aggregation_rule": manifest["cleaning"],
        "admission": manifest["admission"],
    }


@dataclass(frozen=True)
class TrainAllowList:
    row_ids: frozenset[str]
    cell_of_row: dict[str, str]
    cells: dict[str, dict]


def train_allow_list() -> TrainAllowList:
    cells: dict[str, dict] = {}
    cell_of_row: dict[str, str] = {}
    with gzip.open(TRAIN_CELLS, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            cell = json.loads(line)
            cells[cell["cell_id"]] = cell
            for row_id in cell["source_row_ids"]:
                cell_of_row[str(row_id)] = cell["cell_id"]
    return TrainAllowList(frozenset(cell_of_row), cell_of_row, cells)


def build_cache(force: bool = False) -> dict:
    authority = verify_source_authority()
    allow = train_allow_list()
    if CACHE.exists() and not force:
        return {"cache": str(CACHE), "rebuilt": False, "authority": authority}

    kept: dict[str, dict] = {}
    scanned = 0
    dropped_endpoint = 0
    with gzip.open(EXACT_LABELS, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            scanned += 1
            row_id = str(row["source_row_id"])
            if row_id not in allow.row_ids:
                continue
            if row["endpoint"] != ENDPOINT:
                dropped_endpoint += 1
                continue
            kept[row_id] = {
                "source_row_id": row_id,
                "cell_id": allow.cell_of_row[row_id],
                "target_id": row["target_id"],
                "ligand_id": row["ligand_id"],
                "panel_id": row["panel_id"],
                "document_id": row["document_id"],
                "endpoint": row["endpoint"],
                "pK": float(row["pK"]),
            }

    assay_scanned = 0
    with gzip.open(METADATA_PROJECTION, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            assay_scanned += 1
            row_id = str(row["source_row_id"])
            if row_id not in kept:
                continue
            protocols = sorted({str(assay.get("protocol_sha256"))
                                for assay in row.get("assays", [])})
            kept[row_id]["assay_protocols"] = protocols
            kept[row_id]["assay_names"] = sorted({
                str(assay.get("assay_name_norm"))
                for assay in row.get("assays", [])})
            kept[row_id]["publication_date"] = row.get("publication_date")

    with gzip.open(CACHE, "wt", encoding="utf-8") as handle:
        for row_id in sorted(kept, key=lambda value: (len(value), value)):
            handle.write(json.dumps(kept[row_id], sort_keys=True) + "\n")

    return {
        "cache": str(CACHE),
        "cache_sha256": file_sha256(CACHE),
        "rebuilt": True,
        "authority": authority,
        "isolation": {
            "level": "logical_exclusion_after_parsing",
            "physically_isolated": False,
            "allow_list_source": "physically isolated split view, meta_train only",
            "allow_listed_row_ids": len(allow.row_ids),
            "label_rows_scanned": scanned,
            "metadata_rows_scanned": assay_scanned,
            "rows_retained": len(kept),
            "rows_dropped_wrong_endpoint": dropped_endpoint,
            "development_validation_and_sealed_rows_retained": 0,
            "why": ("exact_labels.jsonl.gz is a single all-label artifact; the "
                    "allow-list is applied in the same pass that parses it, so "
                    "no non-meta_train value is ever bound to a retained "
                    "structure. This is weaker than the model path's physical "
                    "isolation and is reported as such."),
        },
    }


def load_cache() -> list[dict]:
    if not CACHE.exists():
        raise FileNotFoundError(f"{CACHE} is missing; run u0_reliability first")
    with gzip.open(CACHE, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def group_by_cell(rows: list[dict]) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[row["cell_id"]].append(row)
    return dict(grouped)
