"""Outcome-blind H0-S inventory for historical A2S-DTA sources."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import tempfile
import time
from typing import Any, Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zipfile import ZipFile

import pyarrow.parquet as pq

from research.a2s.a2s_historical_projection import SNAPSHOTS


CHEMBL_ROOT = Path("dataset/public/chembl_37/processed/dualcold")
CHEMBL_REGISTRY = CHEMBL_ROOT / "registry.parquet"
CHEMBL_MANIFEST = CHEMBL_ROOT / "manifest.json"
CHEMBL_DOCUMENTS = CHEMBL_ROOT / "pcic_o0_document_metadata.json"
CHEMBL_DATED_DOCUMENTS = CHEMBL_ROOT / "historical_document_metadata.v1.json"
CHEMBL_HISTORICAL_ROOT = Path("dataset/public/chembl_historical")
BINDINGDB_ARCHIVE = Path(
    "dataset/public/open_s/BindingDB_BindingDB_Articles_202607_tsv.zip"
)
OPEN_SOURCE_MANIFEST = Path("manifests/open_sources.json")
GTOPDB_ROOT = Path("dataset/public/gtopdb_2026_2/metadata")
DEFAULT_OUTPUT = Path("dataset/processed/a2s_h0_metadata_inventory.v3.json")
FROZEN_H0_COMPANION_EVIDENCE = Path(
    "dataset/processed/a2s_h0_metadata_inventory.v2.json"
)

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_RELEASES = {
    "2018-12-31": "chembl_24_1",
    "2020-12-31": "chembl_27",
    "2022-12-31": "chembl_31",
    "2026-05-29": "chembl_37",
}
CHEMBL_PUBLIC_FILE_DATES = {
    "chembl_24_1": "2018-06-18",
    "chembl_27": "2020-05-21",
    "chembl_31": "2022-08-15",
    "chembl_37": "2026-05-29",
}
CHEMBL_METADATA_FILES = (
    "LICENSE",
    "REQUIRED.ATTRIBUTION",
    "checksums.txt",
)
GTOPDB_METADATA_URLS = {
    "gtp2026.2.ttl": (
        "https://www.guidetopharmacology.org/DATA/rdf/2026.2/gtp2026.2.ttl"
    ),
    "file_descriptions.txt": (
        "https://www.guidetopharmacology.org/DATA/file_descriptions.txt"
    ),
}
BINDINGDB_COMPANIONS = {
    "reaction_set_to_assay": (
        "https://www.bindingdb.org/rwd/bind/downloads/"
        "BindingDB_rsid_eaids_202607_tsv.zip"
    ),
    "assay_descriptions": (
        "https://www.bindingdb.org/rwd/bind/downloads/"
        "BindingDB_Assays_202607_tsv.zip"
    ),
}


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _request(url: str, *, method: str = "GET", timeout: int = 60) -> Any:
    request = Request(url, method=method, headers={"User-Agent": "FORT-H0-S/1.0"})
    last_error: Exception | None = None
    for attempt in range(3):
        try:
            return urlopen(request, timeout=timeout)
        except Exception as error:  # pragma: no cover - exercised only by network faults
            last_error = error
            if attempt < 2:
                time.sleep(attempt + 1)
    assert last_error is not None
    raise last_error


def fetch_bytes(url: str) -> bytes:
    with _request(url) as response:
        return response.read()


def immutable_write(path: Path, data: bytes) -> None:
    if path.exists():
        if path.read_bytes() != data:
            raise FileExistsError(f"refusing to overwrite different frozen file: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(data)
    temporary.replace(path)


def acquire_remote_metadata() -> None:
    for name, url in GTOPDB_METADATA_URLS.items():
        immutable_write(GTOPDB_ROOT / name, fetch_bytes(url))

    for release in CHEMBL_RELEASES.values():
        base = (
            "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
            f"{release}/"
        )
        for name in CHEMBL_METADATA_FILES:
            immutable_write(
                CHEMBL_HISTORICAL_ROOT / "metadata" / release / name,
                fetch_bytes(base + name),
            )


def _chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for offset in range(0, len(values), size):
        yield values[offset : offset + size]


def acquire_dated_chembl_documents(
    source: Path = CHEMBL_DOCUMENTS,
    output: Path = CHEMBL_DATED_DOCUMENTS,
) -> None:
    if output.exists():
        return
    payload = json.loads(source.read_text(encoding="utf-8"))
    local_records = payload["documents"]
    local_by_id = {str(record["document_chembl_id"]): record for record in local_records}
    if len(local_by_id) != len(local_records):
        raise ValueError("local ChEMBL document IDs are not unique")

    fields = (
        "document_chembl_id,chembl_release,doc_type,doi,patent_id,"
        "pubmed_id,src_id,year"
    )
    remote_by_id: dict[str, dict[str, Any]] = {}
    for batch in _chunks(sorted(local_by_id), 100):
        query = urlencode(
            {
                "document_chembl_id__in": ",".join(batch),
                "limit": 1000,
                "only": fields,
            }
        )
        response = json.loads(fetch_bytes(f"{CHEMBL_API}/document.json?{query}"))
        if response["page_meta"]["next"] is not None:
            raise ValueError("unexpected ChEMBL document pagination")
        for record in response["documents"]:
            document = str(record["document_chembl_id"])
            if document in remote_by_id:
                raise ValueError(f"duplicate API document: {document}")
            remote_by_id[document] = record

    missing = sorted(set(local_by_id).difference(remote_by_id))
    unexpected = sorted(set(remote_by_id).difference(local_by_id))
    if missing or unexpected:
        raise ValueError(
            f"ChEMBL API document mismatch: missing={missing[:3]}, "
            f"unexpected={unexpected[:3]}"
        )

    dated: list[dict[str, Any]] = []
    for document in sorted(local_by_id):
        local = local_by_id[document]
        remote = remote_by_id[document]
        release_value = remote.get("chembl_release")
        if not isinstance(release_value, dict):
            raise ValueError(f"document lacks expanded release metadata: {document}")
        release = release_value.get("chembl_release")
        release_date = release_value.get("creation_date")
        if release != local.get("chembl_release"):
            raise ValueError(f"release mismatch for document {document}")
        for field in ("doc_type", "doi", "patent_id", "pubmed_id", "src_id"):
            if remote.get(field) != local.get(field):
                raise ValueError(f"{field} mismatch for document {document}")
        dated.append(
            {
                "document_chembl_id": document,
                "doc_type": remote.get("doc_type"),
                "doi": remote.get("doi"),
                "patent_id": remote.get("patent_id"),
                "pubmed_id": remote.get("pubmed_id"),
                "src_id": remote.get("src_id"),
                "publication_year": remote.get("year"),
                "first_document_release": release,
                "first_document_release_date": release_date,
            }
        )

    frozen = {
        "schema": "a2s-chembl-document-time-v1",
        "source": f"{CHEMBL_API}/document.json",
        "source_document_metadata_sha256": file_sha256(source),
        "documents": dated,
    }
    immutable_write(
        output,
        (json.dumps(frozen, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )


def read_bindingdb_header(path: Path) -> tuple[list[str], dict[str, Any]]:
    with ZipFile(path) as archive:
        entries = archive.infolist()
        if len(entries) != 1:
            raise ValueError("BindingDB archive must contain exactly one TSV")
        entry = entries[0]
        header = bytearray()
        with archive.open(entry) as stream:
            while True:
                value = stream.read(1)
                if not value:
                    raise ValueError("BindingDB TSV has no header newline")
                if value == b"\n":
                    break
                if value != b"\r":
                    header.extend(value)
        fields = header.decode("utf-8-sig").split("\t")
        metadata = {
            "entry": entry.filename,
            "uncompressed_bytes": entry.file_size,
            "compressed_bytes": entry.compress_size,
            "header_bytes": len(header),
        }
    return fields, metadata


def remote_head(url: str) -> dict[str, Any]:
    try:
        with _request(url, method="HEAD") as response:
            return {
                "available": response.status == 200,
                "url": response.geturl(),
                "content_length": int(response.headers.get("Content-Length", 0)),
                "etag": response.headers.get("ETag"),
                "last_modified": response.headers.get("Last-Modified"),
            }
    except Exception as error:  # pragma: no cover - network status is reported
        return {"available": False, "url": url, "error": type(error).__name__}


def frozen_bindingdb_companion_metadata(path: Path) -> dict[str, Any]:
    """Reuse the already frozen HEAD evidence instead of issuing live requests."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    metadata = payload.get("sources", {}).get("bindingdb", {}).get(
        "companion_metadata"
    )
    if not isinstance(metadata, dict) or set(metadata) != set(BINDINGDB_COMPANIONS):
        raise ValueError("frozen BindingDB companion evidence is incomplete")
    for key, expected_url in BINDINGDB_COMPANIONS.items():
        record = metadata[key]
        if not isinstance(record, dict) or record.get("url") != expected_url:
            raise ValueError(f"frozen BindingDB companion URL mismatch: {key}")
    return {
        "evidence_path": str(path),
        "evidence_sha256": file_sha256(path),
        "network_requests_in_current_run": 0,
        "responses": metadata,
    }


def inspect_projection_certificate(
    release: str, database_files: list[str]
) -> dict[str, Any]:
    output_root = Path("dataset/processed/a2s_historical_projection") / release
    projection = output_root / "activity_identity.parquet"
    documents = output_root / "activity_identity.documents.parquet"
    manifest_path = output_root / "activity_identity.manifest.json"
    replay_path = output_root / "activity_identity.replay.json"
    blockers: list[str] = []
    if len(database_files) != 1:
        blockers.append("exactly one extracted SQLite database is required")
    for path in (projection, documents, manifest_path, replay_path):
        if not path.exists():
            blockers.append(f"missing certificate artifact: {path}")
    if blockers:
        return {"certified": False, "blockers": blockers}

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    replay = json.loads(replay_path.read_text(encoding="utf-8"))
    spec = SNAPSHOTS[release]
    database_sha256 = file_sha256(Path(database_files[0]))
    archive_path = CHEMBL_HISTORICAL_ROOT / "archives" / spec.archive_name
    checks = {
        "release": manifest.get("release") == release,
        "archive_sha256": (
            manifest.get("archive", {}).get("sha256") == spec.archive_sha256
        ),
        "database_integrity": (
            manifest.get("database", {})
            .get("certificate", {})
            .get("integrity_check")
            == "ok"
        ),
        "database_sha256": (
            database_sha256
            == manifest.get("database", {}).get("sha256_before_open")
            == replay.get("database_sha256")
        ),
        "archive_current": (
            archive_path.exists()
            and archive_path.stat().st_size == spec.archive_bytes
            and file_sha256(archive_path) == spec.archive_sha256
        ),
        "activity_projection_sha256": (
            file_sha256(projection)
            == manifest.get("projection", {}).get("sha256")
        ),
        "document_projection_sha256": (
            file_sha256(documents)
            == manifest.get("document_projection", {}).get("sha256")
        ),
        "outcome_firewall": (
            manifest.get("firewall", {}).get("numeric_affinity_columns_read") == []
            and manifest.get("firewall", {}).get("affinity_values_materialized") == 0
            and manifest.get("firewall", {}).get("denied_actions_during_projection")
            == []
        ),
        "publisher_metadata": (
            set(
                manifest.get("archive", {})
                .get("publisher_metadata", {})
                .keys()
            )
            == {"checksums.txt", "LICENSE", "REQUIRED.ATTRIBUTION"}
        ),
        "replay": (
            replay.get("status") == "PASS"
            and all(replay.get("comparisons", {}).values())
            and replay.get("primary_manifest_sha256") == file_sha256(manifest_path)
        ),
    }
    blockers.extend(name for name, passed in checks.items() if not passed)
    return {
        "certified": not blockers,
        "checks": checks,
        "blockers": blockers,
        "manifest": str(manifest_path),
        "manifest_sha256": file_sha256(manifest_path),
        "replay": str(replay_path),
        "replay_sha256": file_sha256(replay_path),
    }


def inspect_presence_certificate() -> dict[str, Any]:
    output = Path(
        "dataset/processed/a2s_historical_presence/activity_presence.v1.parquet"
    )
    manifest_path = output.with_suffix(".manifest.json")
    replay_path = output.with_suffix(".replay.json")
    if not output.exists() or not manifest_path.exists() or not replay_path.exists():
        return {
            "certified": False,
            "blockers": ["cross-release presence audit and replay are missing"],
        }
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    replay = json.loads(replay_path.read_text(encoding="utf-8"))
    checks = {
        "identity_status": manifest.get("status") == "H0_IDENTITY_PASS",
        "parquet_sha256": (
            file_sha256(output) == manifest.get("output", {}).get("parquet_sha256")
        ),
        "outcome_firewall": (
            manifest.get("firewall", {}).get("affinity_values_read") == 0
            and manifest.get("firewall", {}).get("affinity_values_materialized") == 0
        ),
        "replay": (
            replay.get("status") == "PASS"
            and all(replay.get("comparisons", {}).values())
            and replay.get("primary_manifest_sha256") == file_sha256(manifest_path)
        ),
    }
    return {
        "certified": all(checks.values()),
        "checks": checks,
        "blockers": [name for name, passed in checks.items() if not passed],
        "manifest": str(manifest_path),
        "manifest_sha256": file_sha256(manifest_path),
    }


def inspect_chembl() -> dict[str, Any]:
    base_payload = json.loads(CHEMBL_DOCUMENTS.read_text(encoding="utf-8"))
    base_documents = base_payload["documents"]
    dated_documents: list[dict[str, Any]] = []
    if CHEMBL_DATED_DOCUMENTS.exists():
        dated_payload = json.loads(CHEMBL_DATED_DOCUMENTS.read_text(encoding="utf-8"))
        dated_documents = dated_payload["documents"]

    schema = pq.ParquetFile(CHEMBL_REGISTRY).schema_arrow.names
    identity_fields = {
        "doi": sum(record.get("doi") is not None for record in base_documents),
        "pubmed_id": sum(
            record.get("pubmed_id") is not None for record in base_documents
        ),
        "patent_id": sum(
            record.get("patent_id") is not None for record in base_documents
        ),
        "src_id": sum(record.get("src_id") is not None for record in base_documents),
    }
    snapshots: dict[str, Any] = {}
    for index_date, release in CHEMBL_RELEASES.items():
        snapshot_root = CHEMBL_HISTORICAL_ROOT / "snapshots" / release
        database_files = (
            sorted(str(path) for path in snapshot_root.rglob("*.db"))
            + sorted(str(path) for path in snapshot_root.rglob("*.sqlite"))
        ) if snapshot_root.exists() else []
        checksum_path = CHEMBL_HISTORICAL_ROOT / "metadata" / release / "checksums.txt"
        projection_certificate = inspect_projection_certificate(release, database_files)
        snapshots[index_date] = {
            "release": release,
            "public_file_date": CHEMBL_PUBLIC_FILE_DATES[release],
            "selection_rule": "latest official database file available by index date",
            "database_present": bool(database_files),
            "database_files": database_files,
            "projection_certificate": projection_certificate,
            "snapshot_certified": projection_certificate["certified"],
            "checksum_manifest_present": checksum_path.exists(),
            "checksum_manifest_sha256": (
                file_sha256(checksum_path) if checksum_path.exists() else None
            ),
        }

    return {
        "version": "ChEMBL 37",
        "registry": {
            "path": str(CHEMBL_REGISTRY),
            "sha256": file_sha256(CHEMBL_REGISTRY),
            "schema_fields": schema,
            "row_values_read": 0,
            "activity_first_seen_release_field": (
                "activity_first_seen_release" in schema
            ),
        },
        "manifest": {
            "path": str(CHEMBL_MANIFEST),
            "sha256": file_sha256(CHEMBL_MANIFEST),
            "license": json.loads(CHEMBL_MANIFEST.read_text(encoding="utf-8"))[
                "license"
            ],
        },
        "documents": {
            "path": str(CHEMBL_DOCUMENTS),
            "sha256": file_sha256(CHEMBL_DOCUMENTS),
            "count": len(base_documents),
            "identity_non_null": identity_fields,
            "provider_counts": dict(
                sorted(Counter(str(record["src_id"]) for record in base_documents).items())
            ),
            "dated_path": str(CHEMBL_DATED_DOCUMENTS),
            "dated_sha256": (
                file_sha256(CHEMBL_DATED_DOCUMENTS)
                if CHEMBL_DATED_DOCUMENTS.exists()
                else None
            ),
            "publication_year_non_null": sum(
                record.get("publication_year") is not None
                for record in dated_documents
            ),
            "document_release_date_non_null": sum(
                record.get("first_document_release_date") is not None
                for record in dated_documents
            ),
            "activity_time_warning": (
                "document release is not activity-row first-seen evidence"
            ),
        },
        "historical_snapshots": snapshots,
    }


def inspect_bindingdb() -> dict[str, Any]:
    fields, entry = read_bindingdb_header(BINDINGDB_ARCHIVE)
    field_set = set(fields)
    source_manifest = json.loads(OPEN_SOURCE_MANIFEST.read_text(encoding="utf-8"))
    recorded = source_manifest["sources"]["bindingdb_native_articles"]
    required = {
        "record_id": "BindingDB Reactant_set_id",
        "parent_structure": "Ligand InChI Key",
        "upstream_source": "Curation/DataSource",
        "doi": "Article DOI",
        "pmid": "PMID",
        "patent": "Patent Number",
        "publication_date": "Date of publication",
        "curation_date": "Date in BindingDB",
        "institution": "Institution",
        "pubchem_assay": "PubChem AID",
    }
    return {
        "version": recorded["release"],
        "path": str(BINDINGDB_ARCHIVE),
        "bytes": BINDINGDB_ARCHIVE.stat().st_size,
        "sha256": file_sha256(BINDINGDB_ARCHIVE),
        "recorded_sha256_matches": (
            file_sha256(BINDINGDB_ARCHIVE).lower() == recorded["sha256"].lower()
        ),
        "license": recorded["licence"],
        "zip_entry": entry,
        "header_field_count": len(fields),
        "required_field_presence": {
            key: field in field_set for key, field in required.items()
        },
        "data_rows_read": 0,
        "generic_assay_lineage_in_main_tsv": False,
        "companion_metadata": frozen_bindingdb_companion_metadata(
            FROZEN_H0_COMPANION_EVIDENCE
        ),
    }


def inspect_gtopdb() -> dict[str, Any]:
    ttl_path = GTOPDB_ROOT / "gtp2026.2.ttl"
    description_path = GTOPDB_ROOT / "file_descriptions.txt"
    ttl = ttl_path.read_text(encoding="utf-8") if ttl_path.exists() else ""
    descriptions = (
        description_path.read_text(encoding="utf-8")
        if description_path.exists()
        else ""
    )
    issued = re.search(r'dct:issued\s+"([0-9-]+)"', ttl)
    return {
        "version": "2026.2" if "Version 2026.2" in ttl else None,
        "issued": issued.group(1) if issued else None,
        "license": (
            "ODbL data; CC BY-SA 4.0 contents"
            if "opendatacommons.org/licenses/odbl" in ttl
            and "Creative Commons Attribution-ShareAlike 4.0" in ttl
            else None
        ),
        "metadata_files": {
            str(path): file_sha256(path) if path.exists() else None
            for path in (ttl_path, description_path)
        },
        "schema_evidence": {
            "pubmed": "pubmed_id" in descriptions.lower(),
            "patent": "patent numbers" in descriptions.lower(),
            "assay_description": "assay_description" in descriptions.lower(),
            "inchikey": "inchikey" in descriptions.lower(),
        },
        "affinity_rows_read": 0,
        "per_interaction_first_seen_release": False,
    }


def historical_admission(
    *,
    dated_documents: bool,
    activity_first_seen: bool,
    snapshots_present: bool,
    presence_certified: bool = False,
) -> tuple[str, list[str]]:
    blockers: list[str] = []
    if not dated_documents:
        blockers.append("publication year and document-release date are not frozen")
    if not activity_first_seen and not snapshots_present:
        blockers.append(
            "ChEMBL-37 has neither activity-row first-seen metadata nor the "
            "24.1/27/31/37 certified databases needed to detect backfill"
        )
    if snapshots_present and not activity_first_seen and not presence_certified:
        blockers.append(
            "certified snapshots exist but cross-release stable identity/presence "
            "audit has not passed"
        )
    return ("H0_S_PASS" if not blockers else "DATA_NOT_READY", blockers)


def build_report() -> dict[str, Any]:
    chembl = inspect_chembl()
    bindingdb = inspect_bindingdb()
    gtopdb = inspect_gtopdb()
    snapshots_present = all(
        record["snapshot_certified"]
        for record in chembl["historical_snapshots"].values()
    )
    presence_certificate = inspect_presence_certificate()
    dated_documents = (
        chembl["documents"]["publication_year_non_null"] > 0
        and chembl["documents"]["document_release_date_non_null"]
        == chembl["documents"]["count"]
    )
    status, blockers = historical_admission(
        dated_documents=dated_documents,
        activity_first_seen=chembl["registry"]["activity_first_seen_release_field"],
        snapshots_present=snapshots_present,
        presence_certified=presence_certificate["certified"],
    )
    return {
        "schema": "a2s-h0-metadata-inventory-v2",
        "stage": "H0-S",
        "firewall": {
            "numeric_affinity_values_read": 0,
            "affinity_rows_materialized": 0,
            "parquet_access": "footer schema only",
            "bindingdb_access": "ZIP central directory and first TSV line only",
            "gtopdb_access": "RDF dataset metadata and file descriptions only",
        },
        "sources": {
            "chembl": chembl,
            "bindingdb": bindingdb,
            "gtopdb": gtopdb,
        },
        "historical_presence": presence_certificate,
        "cross_database_deduplication": {
            "document_keys": ["normalized DOI", "PMID", "patent identifier"],
            "compound_keys": ["parent InChIKey", "connectivity key"],
            "target_keys": ["UniProt accession", "sequence hash"],
            "feasible": True,
            "limitation": (
                "matching identifiers support deduplication but do not prove "
                "independent assay or campaign lineage"
            ),
        },
        "decision": {
            "status": status,
            "hist_s_authorized_for_h1": status == "H0_S_PASS",
            "hist_l_training_authorized": False,
            "blockers": blockers,
            "next_action": (
                "construct the H1-S affinity-blind historical roster, dependency "
                "graph, cutflow, and MDE envelope"
                if status == "H0_S_PASS"
                else "acquire and checksum ChEMBL 24.1, 27, 31, and 37; freeze "
                "integrity-checked, independently replayed activity/document "
                "identity projections without selecting or materializing affinity values"
            ),
        },
        "code_sha256": file_sha256(Path(__file__)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--acquire-remote-metadata",
        action="store_true",
        help="freeze small official metadata files and ChEMBL document dates",
    )
    return parser.parse_args()


def main() -> None:
    from research.shared.dataset_processing import require_run_context

    require_run_context(
        [str(Path.cwd() / "main.py"), "h0-inventory", *sys.argv[1:]]
    )
    args = parse_args()
    if args.acquire_remote_metadata:
        acquire_remote_metadata()
        acquire_dated_chembl_documents()
    report = build_report()
    immutable_write(
        args.output,
        (json.dumps(report, indent=2, sort_keys=True) + "\n").encode("utf-8"),
    )
    print(json.dumps(report["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
