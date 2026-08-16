"""Verify, extract, and outcome-blindly project historical ChEMBL snapshots."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import platform
import re
import shutil
import sqlite3
import sys
import tarfile
import tempfile
from typing import Any, Iterable

import pyarrow as pa
import pyarrow.parquet as pq


HISTORICAL_ROOT = Path("dataset/public/chembl_historical")
DEFAULT_OUTPUT_ROOT = Path("dataset/processed/a2s_historical_projection")


@dataclass(frozen=True)
class SnapshotSpec:
    release: str
    index_date: str
    public_file_date: str
    archive_bytes: int
    archive_sha256: str
    database_version_names: tuple[str, ...] = ()

    @property
    def archive_name(self) -> str:
        return f"{self.release}_sqlite.tar.gz"

    @property
    def archive_url(self) -> str:
        return (
            "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
            f"{self.release}/{self.archive_name}"
        )


SNAPSHOTS = {
    spec.release: spec
    for spec in (
        SnapshotSpec(
            "chembl_24_1",
            "2018-12-31",
            "2018-06-18",
            3_659_492_620,
            "6bb1030408c68b26ad8e9e9ae34ca7226e55d4214a6a057d93be987bbba5ea8c",
            ("ChEMBL_24",),
        ),
        SnapshotSpec(
            "chembl_27",
            "2020-12-31",
            "2020-05-21",
            3_742_560_551,
            "5abce60db823266834312a712b95fe1f81003b6d15be09474667b9586f1854b4",
            ("ChEMBL_27",),
        ),
        SnapshotSpec(
            "chembl_31",
            "2022-12-31",
            "2022-08-15",
            4_505_413_744,
            "cbd28e4b0e955e562d975d972d262298286d36e3992292f783dce8f498b04e26",
            ("ChEMBL_31",),
        ),
        SnapshotSpec(
            "chembl_37",
            "2026-05-29",
            "2026-05-29",
            5_764_252_857,
            "33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281",
            ("ChEMBL_37",),
        ),
    )
}


# The SQLite authorizer enforces this allowlist for every executed projection query.
READ_ALLOWLIST = {
    "activities": {
        "activity_id",
        "assay_id",
        "doc_id",
        "record_id",
        "molregno",
        "standard_relation",
        "standard_units",
        "standard_type",
        "src_id",
    },
    "assays": {
        "assay_id",
        "doc_id",
        "tid",
        "src_id",
        "src_assay_id",
        "chembl_id",
        "aidx",
    },
    "docs": {
        "doc_id",
        "chembl_id",
        "year",
        "pubmed_id",
        "doi",
        "patent_id",
        "doc_type",
        "src_id",
        "ridx",
    },
    "compound_records": {
        "record_id",
        "molregno",
        "doc_id",
        "src_id",
        "src_compound_id",
        "cidx",
    },
    "molecule_hierarchy": {"molregno", "parent_molregno"},
    "molecule_dictionary": {"molregno", "chembl_id"},
    "compound_structures": {"molregno", "standard_inchi_key"},
    "target_dictionary": {
        "tid",
        "target_type",
        "tax_id",
        "organism",
        "chembl_id",
    },
    "target_components": {"tid", "component_id", "homologue"},
    "component_sequences": {
        "component_id",
        "component_type",
        "accession",
        "sequence_md5sum",
        "tax_id",
        "db_source",
        "db_version",
    },
    "source": {"src_id", "src_short_name"},
    "version": {"name", "creation_date"},
}

AFFINITY_VALUE_COLUMNS = {
    "value",
    "published_value",
    "standard_value",
    "pchembl_value",
    "upper_value",
    "standard_upper_value",
    "text_value",
    "standard_text_value",
    "activity_comment",
}

if any(
    AFFINITY_VALUE_COLUMNS.intersection(columns)
    for columns in READ_ALLOWLIST.values()
):  # pragma: no cover - import-time contract
    raise RuntimeError("affinity value column entered the projection allowlist")


PROJECTION_SQL = """
SELECT
    a.activity_id,
    a.standard_type AS endpoint_type,
    a.standard_relation,
    a.standard_units,
    a.src_id AS activity_src_id,
    src.src_short_name AS activity_source,
    ass.assay_id,
    ass.chembl_id AS assay_chembl_id,
    ass.src_assay_id AS assay_source_id,
    ass.aidx AS assay_index,
    ass.src_id AS assay_src_id,
    ass.tid AS target_id,
    td.chembl_id AS target_chembl_id,
    td.target_type,
    td.tax_id AS target_tax_id,
    td.organism AS target_organism,
    a.doc_id AS activity_doc_id,
    ass.doc_id AS assay_doc_id,
    cr.doc_id AS record_doc_id,
    d.chembl_id AS document_chembl_id,
    d.year AS publication_year,
    d.doi,
    d.pubmed_id,
    d.patent_id,
    d.doc_type,
    d.src_id AS document_source_id,
    d.ridx AS document_index,
    a.record_id,
    cr.cidx AS compound_index,
    cr.src_compound_id AS source_compound_id,
    cr.src_id AS record_source_id,
    a.molregno,
    cr.molregno AS record_molregno,
    COALESCE(mh.parent_molregno, a.molregno, cr.molregno) AS parent_molregno,
    pmd.chembl_id AS parent_molecule_chembl_id,
    pcs.standard_inchi_key AS parent_inchikey
FROM activities AS a
JOIN assays AS ass ON ass.assay_id = a.assay_id
JOIN target_dictionary AS td ON td.tid = ass.tid
LEFT JOIN compound_records AS cr ON cr.record_id = a.record_id
LEFT JOIN molecule_hierarchy AS mh
    ON mh.molregno = COALESCE(a.molregno, cr.molregno)
LEFT JOIN molecule_dictionary AS pmd
    ON pmd.molregno = COALESCE(mh.parent_molregno, a.molregno, cr.molregno)
LEFT JOIN compound_structures AS pcs ON pcs.molregno = pmd.molregno
LEFT JOIN docs AS d ON d.doc_id = COALESCE(a.doc_id, ass.doc_id, cr.doc_id)
LEFT JOIN source AS src ON src.src_id = a.src_id
WHERE UPPER(TRIM(a.standard_type)) = 'KI'
  AND TRIM(a.standard_relation) = '='
  AND UPPER(TRIM(td.target_type)) = 'SINGLE PROTEIN'
  AND UPPER(TRIM(td.organism)) = 'HOMO SAPIENS'
ORDER BY a.activity_id
"""


TARGET_COMPONENT_SQL = """
SELECT tc.tid, tc.component_id, tc.homologue, cs.component_type,
       cs.accession, cs.sequence_md5sum, cs.tax_id, cs.db_source, cs.db_version
FROM target_components AS tc
JOIN component_sequences AS cs ON cs.component_id = tc.component_id
ORDER BY tc.tid, tc.component_id
"""


DOCUMENT_SQL = """
SELECT d.doc_id, d.chembl_id AS document_chembl_id,
       d.year AS publication_year, d.doi, d.pubmed_id, d.patent_id,
       d.doc_type, d.src_id AS document_source_id,
       src.src_short_name AS document_source, d.ridx AS document_index
FROM docs AS d
LEFT JOIN source AS src ON src.src_id = d.src_id
ORDER BY d.doc_id
"""


PROJECTION_SCHEMA = pa.schema(
    [
        ("release", pa.string()),
        ("index_date", pa.string()),
        ("public_file_date", pa.string()),
        ("activity_id", pa.int64()),
        ("native_activity_key", pa.string()),
        ("measurement_identity_sha256", pa.string()),
        ("endpoint_type", pa.string()),
        ("standard_relation", pa.string()),
        ("standard_units", pa.string()),
        ("activity_src_id", pa.int64()),
        ("activity_source", pa.string()),
        ("assay_id", pa.int64()),
        ("assay_chembl_id", pa.string()),
        ("assay_source_id", pa.string()),
        ("assay_index", pa.string()),
        ("assay_src_id", pa.int64()),
        ("target_id", pa.int64()),
        ("target_chembl_id", pa.string()),
        ("target_type", pa.string()),
        ("target_tax_id", pa.int64()),
        ("target_organism", pa.string()),
        ("target_component_count", pa.int64()),
        ("target_uniprot", pa.string()),
        ("target_sequence_hash", pa.string()),
        ("target_sequence_hash_algorithm", pa.string()),
        ("target_homologue", pa.int64()),
        ("target_component_type", pa.string()),
        ("target_component_tax_id", pa.int64()),
        ("target_db_source", pa.string()),
        ("target_db_version", pa.string()),
        ("activity_doc_id", pa.int64()),
        ("assay_doc_id", pa.int64()),
        ("record_doc_id", pa.int64()),
        ("document_id_consistent", pa.bool_()),
        ("document_chembl_id", pa.string()),
        ("publication_year", pa.int64()),
        ("doi", pa.string()),
        ("normalized_doi", pa.string()),
        ("pubmed_id", pa.int64()),
        ("patent_id", pa.string()),
        ("doc_type", pa.string()),
        ("document_source_id", pa.int64()),
        ("document_index", pa.string()),
        ("record_id", pa.int64()),
        ("compound_index", pa.string()),
        ("source_compound_id", pa.string()),
        ("record_source_id", pa.int64()),
        ("molregno", pa.int64()),
        ("record_molregno", pa.int64()),
        ("molecule_id_consistent", pa.bool_()),
        ("parent_molregno", pa.int64()),
        ("parent_molecule_chembl_id", pa.string()),
        ("parent_inchikey", pa.string()),
        ("connectivity_key", pa.string()),
    ]
)


DOCUMENT_SCHEMA = pa.schema(
    [
        ("release", pa.string()),
        ("public_file_date", pa.string()),
        ("doc_id", pa.int64()),
        ("document_chembl_id", pa.string()),
        ("publication_year", pa.int64()),
        ("doi", pa.string()),
        ("normalized_doi", pa.string()),
        ("pubmed_id", pa.int64()),
        ("patent_id", pa.string()),
        ("doc_type", pa.string()),
        ("document_source_id", pa.int64()),
        ("document_source", pa.string()),
        ("document_index", pa.string()),
    ]
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_archive(path: Path, spec: SnapshotSpec) -> dict[str, Any]:
    size = path.stat().st_size
    if size != spec.archive_bytes:
        raise ValueError(
            f"{spec.release} archive size mismatch: {size} != {spec.archive_bytes}"
        )
    digest = file_sha256(path)
    if digest.lower() != spec.archive_sha256.lower():
        raise ValueError(
            f"{spec.release} archive SHA-256 mismatch: {digest} != "
            f"{spec.archive_sha256}"
        )
    return {"path": str(path), "bytes": size, "sha256": digest}


def verify_publisher_metadata(spec: SnapshotSpec) -> dict[str, Any]:
    metadata_root = HISTORICAL_ROOT / "metadata" / spec.release
    paths = {
        name: metadata_root / name
        for name in ("checksums.txt", "LICENSE", "REQUIRED.ATTRIBUTION")
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"publisher metadata is missing: {missing}")
    checksum_text = paths["checksums.txt"].read_text(encoding="utf-8-sig")
    entries = {
        name: digest.lower()
        for digest, name in re.findall(r"(?im)^([0-9a-f]{64})\s+([^\s]+)\s*$", checksum_text)
    }
    published = entries.get(spec.archive_name)
    if published != spec.archive_sha256.lower():
        raise ValueError(
            f"publisher checksum mismatch for {spec.archive_name}: {published}"
        )
    return {
        name: {"path": str(path), "sha256": file_sha256(path)}
        for name, path in paths.items()
    }


def extract_single_sqlite(archive: Path, destination: Path) -> Path:
    """Stream one SQLite member from a verified tar archive without path traversal."""
    destination.mkdir(parents=True, exist_ok=True)
    existing = sorted((*destination.glob("*.db"), *destination.glob("*.sqlite")))
    if len(existing) == 1:
        return existing[0]
    if existing:
        raise ValueError(f"multiple SQLite files already present in {destination}")

    temporary: Path | None = None
    member_name: str | None = None
    try:
        with tarfile.open(archive, mode="r|gz") as source:
            for member in source:
                suffix = Path(member.name).suffix.lower()
                if not member.isfile() or suffix not in {".db", ".sqlite"}:
                    continue
                if member_name is not None:
                    raise ValueError("archive contains more than one SQLite database")
                member_name = member.name
                extracted = source.extractfile(member)
                if extracted is None:
                    raise ValueError(f"cannot read SQLite member {member.name}")
                with tempfile.NamedTemporaryFile(
                    dir=destination, prefix=".extract-", delete=False
                ) as output:
                    temporary = Path(output.name)
                    shutil.copyfileobj(extracted, output, length=1024 * 1024)
        if member_name is None or temporary is None:
            raise ValueError("archive contains no SQLite database")
        final = destination / Path(member_name).name
        temporary.replace(final)
        return final
    except Exception:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def _schema_columns(connection: sqlite3.Connection) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {}
    for table in sorted(READ_ALLOWLIST):
        rows = connection.execute(f'PRAGMA table_info("{table}")').fetchall()
        if not rows:
            raise ValueError(f"required ChEMBL table is missing: {table}")
        result[table] = [str(row[1]).lower() for row in rows]
        missing = READ_ALLOWLIST[table].difference(result[table])
        if missing:
            raise ValueError(f"{table} lacks required columns: {sorted(missing)}")
    return result


def _database_certificate(
    connection: sqlite3.Connection,
    spec: SnapshotSpec,
    *,
    full_integrity: bool,
) -> dict[str, Any]:
    if full_integrity:
        integrity = [str(row[0]) for row in connection.execute("PRAGMA integrity_check")]
        if integrity != ["ok"]:
            raise ValueError(f"SQLite integrity check failed: {integrity[:3]}")
        integrity_status = "ok"
    else:
        integrity_status = "inherited_from_hash-matched_primary_projection"
    master = connection.execute(
        "SELECT type, name, tbl_name, sql FROM sqlite_master "
        "WHERE name NOT LIKE 'sqlite_%' ORDER BY type, name"
    ).fetchall()
    canonical_master = json.dumps(master, ensure_ascii=True, separators=(",", ":"))
    return {
        "integrity_check": integrity_status,
        "sqlite_master_sha256": sha256(canonical_master.encode("utf-8")).hexdigest(),
        "sqlite_master_entries": len(master),
        "expected_release": spec.release,
    }


def _validate_version(
    connection: sqlite3.Connection, spec: SnapshotSpec
) -> list[dict[str, Any]]:
    rows = connection.execute(
        "SELECT name, creation_date FROM version ORDER BY name"
    ).fetchall()
    expected_names = spec.database_version_names or (
        spec.release.replace("chembl_", "ChEMBL_", 1),
    )
    expected_keys = {
        re.sub(r"[^a-z0-9]+", "", name.lower()) for name in expected_names
    }
    matches = [
        row
        for row in rows
        if re.sub(r"[^a-z0-9]+", "", str(row[0]).lower()) in expected_keys
    ]
    if len(matches) != 1 or matches[0][1] is None:
        raise ValueError(
            f"version table does not certify {spec.release}: "
            f"expected {list(expected_names)}, observed "
            f"{[(row[0], row[1]) for row in rows]}"
        )
    return [{"name": str(row[0]), "creation_date": str(row[1])} for row in rows]


def _install_read_firewall(
    connection: sqlite3.Connection,
    reads: set[tuple[str, str]],
    denied: set[tuple[int, str, str]],
) -> None:
    allowed_functions = {"coalesce", "trim", "upper"}

    def authorize(
        action: int,
        table: str | None,
        column: str | None,
        database: str | None,
        _source: str | None,
    ) -> int:
        if action == sqlite3.SQLITE_SELECT:
            return sqlite3.SQLITE_OK
        if action == sqlite3.SQLITE_FUNCTION:
            function = (column or table or "").lower()
            if function in allowed_functions:
                return sqlite3.SQLITE_OK
            denied.add((action, "function", function))
            return sqlite3.SQLITE_DENY
        if action != sqlite3.SQLITE_READ:
            denied.add((action, table or "", column or ""))
            return sqlite3.SQLITE_DENY
        table_name = (table or "").lower()
        column_name = (column or "").lower()
        if database == "main" and column_name in READ_ALLOWLIST.get(table_name, set()):
            reads.add((table_name, column_name))
            return sqlite3.SQLITE_OK
        denied.add((action, table_name, column_name))
        return sqlite3.SQLITE_DENY

    connection.set_authorizer(authorize)


def _target_components(connection: sqlite3.Connection) -> dict[int, dict[str, Any]]:
    grouped: dict[int, list[tuple[Any, ...]]] = {}
    for row in connection.execute(TARGET_COMPONENT_SQL):
        grouped.setdefault(int(row[0]), []).append(row)

    result: dict[int, dict[str, Any]] = {}
    for target_id, rows in grouped.items():
        single = rows[0] if len(rows) == 1 else None
        sequence_hash = single[5] if single is not None else None
        result[target_id] = {
            "target_component_count": len(rows),
            "target_uniprot": single[4] if single is not None else None,
            "target_sequence_hash": sequence_hash,
            "target_sequence_hash_algorithm": "MD5" if sequence_hash else None,
            "target_homologue": single[2] if single is not None else None,
            "target_component_type": single[3] if single is not None else None,
            "target_component_tax_id": single[6] if single is not None else None,
            "target_db_source": single[7] if single is not None else None,
            "target_db_version": single[8] if single is not None else None,
        }
    return result


def _normalized_doi(value: Any) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :].strip()
    return normalized or None


def _document_ids_consistent(row: dict[str, Any]) -> bool:
    values = {
        int(value)
        for value in (
            row["activity_doc_id"],
            row["assay_doc_id"],
            row["record_doc_id"],
        )
        if value is not None
    }
    return bool(values) and len(values) == 1


def _molecule_ids_consistent(row: dict[str, Any]) -> bool:
    values = {
        int(value)
        for value in (row["molregno"], row["record_molregno"])
        if value is not None
    }
    return bool(values) and len(values) == 1


def _measurement_identity(row: dict[str, Any]) -> str:
    identity = [
        row["activity_src_id"],
        row["activity_source"],
        row["document_source_id"],
        row["document_chembl_id"],
        row["document_index"],
        row["assay_chembl_id"],
        row["assay_src_id"],
        row["assay_index"],
        row["assay_source_id"],
        row["compound_index"],
        row["record_source_id"],
        row["source_compound_id"],
        row["parent_molecule_chembl_id"],
        row["target_chembl_id"],
        row["endpoint_type"],
        row["standard_relation"],
        row["standard_units"],
    ]
    payload = json.dumps(identity, ensure_ascii=True, separators=(",", ":"))
    return sha256(payload.encode("utf-8")).hexdigest()


def _projected_rows(
    cursor: sqlite3.Cursor,
    spec: SnapshotSpec,
    components: dict[int, dict[str, Any]],
    batch_size: int,
) -> Iterable[list[dict[str, Any]]]:
    columns = [description[0] for description in cursor.description]
    previous_activity_id: int | None = None
    while rows := cursor.fetchmany(batch_size):
        projected: list[dict[str, Any]] = []
        for values in rows:
            source = dict(zip(columns, values))
            activity_id = int(source["activity_id"])
            if previous_activity_id is not None and activity_id <= previous_activity_id:
                raise ValueError(
                    "activity projection is not one-to-one and strictly ordered; "
                    f"encountered {activity_id} after {previous_activity_id}"
                )
            previous_activity_id = activity_id
            target = components.get(int(source["target_id"]), {})
            parent_inchikey = source["parent_inchikey"]
            source.update(target)
            source.update(
                {
                    "release": spec.release,
                    "index_date": spec.index_date,
                    "public_file_date": spec.public_file_date,
                    "native_activity_key": f"chembl:activity:{source['activity_id']}",
                    "document_id_consistent": _document_ids_consistent(source),
                    "molecule_id_consistent": _molecule_ids_consistent(source),
                    "normalized_doi": _normalized_doi(source["doi"]),
                    "connectivity_key": (
                        str(parent_inchikey)[:14]
                        if parent_inchikey is not None
                        else None
                    ),
                }
            )
            source["measurement_identity_sha256"] = _measurement_identity(source)
            for name in (
                "target_component_count",
                "target_uniprot",
                "target_sequence_hash",
                "target_sequence_hash_algorithm",
                "target_homologue",
                "target_component_type",
                "target_component_tax_id",
                "target_db_source",
                "target_db_version",
            ):
                source.setdefault(name, None)
            projected.append({name: source.get(name) for name in PROJECTION_SCHEMA.names})
        yield projected


def _document_rows(
    cursor: sqlite3.Cursor, spec: SnapshotSpec, batch_size: int
) -> Iterable[list[dict[str, Any]]]:
    columns = [description[0] for description in cursor.description]
    previous_doc_id: int | None = None
    while rows := cursor.fetchmany(batch_size):
        projected: list[dict[str, Any]] = []
        for values in rows:
            source = dict(zip(columns, values))
            doc_id = int(source["doc_id"])
            if previous_doc_id is not None and doc_id <= previous_doc_id:
                raise ValueError(
                    "document projection is not one-to-one and strictly ordered; "
                    f"encountered {doc_id} after {previous_doc_id}"
                )
            previous_doc_id = doc_id
            source.update(
                {
                    "release": spec.release,
                    "public_file_date": spec.public_file_date,
                    "normalized_doi": _normalized_doi(source["doi"]),
                }
            )
            projected.append({name: source.get(name) for name in DOCUMENT_SCHEMA.names})
        yield projected


def project_snapshot(
    database: Path,
    output: Path,
    spec: SnapshotSpec,
    archive_record: dict[str, Any],
    *,
    batch_size: int = 50_000,
    full_integrity: bool = True,
) -> dict[str, Any]:
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = output.with_suffix(".manifest.json")
    document_output = output.with_name(f"{output.stem}.documents.parquet")
    if output.exists() or document_output.exists() or manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite frozen projection: {output}")

    database_digest = file_sha256(database)
    uri = f"file:{database.resolve().as_posix()}?mode=ro&immutable=1"
    connection = sqlite3.connect(uri, uri=True)
    reads: set[tuple[str, str]] = set()
    denied: set[tuple[int, str, str]] = set()
    temporary = output.with_name(f".{output.name}.partial")
    temporary_documents = document_output.with_name(
        f".{document_output.name}.partial"
    )
    temporary_manifest = manifest_path.with_name(f".{manifest_path.name}.partial")
    temporary.unlink(missing_ok=True)
    temporary_documents.unlink(missing_ok=True)
    temporary_manifest.unlink(missing_ok=True)
    writer: pq.ParquetWriter | None = None
    row_count = 0
    canonical_rows = sha256()
    document_row_count = 0
    canonical_documents = sha256()
    null_counts = {
        "document_chembl_id": 0,
        "parent_inchikey": 0,
        "target_uniprot": 0,
        "target_sequence_hash": 0,
    }
    conflict_counts = {"document_id": 0, "molecule_id": 0}
    try:
        connection.execute("PRAGMA query_only = ON")
        _install_read_firewall(connection, reads, denied)
        version_rows = _validate_version(connection, spec)
        connection.set_authorizer(None)
        database_certificate = _database_certificate(
            connection, spec, full_integrity=full_integrity
        )
        schemas = _schema_columns(connection)
        schema_digest = sha256(
            json.dumps(schemas, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        _install_read_firewall(connection, reads, denied)
        components = _target_components(connection)
        cursor = connection.execute(PROJECTION_SQL)
        writer = pq.ParquetWriter(temporary, PROJECTION_SCHEMA, compression="zstd")
        for records in _projected_rows(cursor, spec, components, batch_size):
            for record in records:
                canonical_rows.update(
                    (
                        json.dumps(
                            record,
                            ensure_ascii=True,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        + "\n"
                    ).encode("utf-8")
                )
                for name in null_counts:
                    null_counts[name] += record[name] is None
                conflict_counts["document_id"] += not record["document_id_consistent"]
                conflict_counts["molecule_id"] += not record["molecule_id_consistent"]
            table = pa.Table.from_pylist(records, schema=PROJECTION_SCHEMA)
            writer.write_table(table)
            row_count += len(records)
        writer.close()
        writer = None

        cursor = connection.execute(DOCUMENT_SQL)
        writer = pq.ParquetWriter(
            temporary_documents, DOCUMENT_SCHEMA, compression="zstd"
        )
        for records in _document_rows(cursor, spec, batch_size):
            for record in records:
                canonical_documents.update(
                    (
                        json.dumps(
                            record,
                            ensure_ascii=True,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        + "\n"
                    ).encode("utf-8")
                )
            writer.write_table(pa.Table.from_pylist(records, schema=DOCUMENT_SCHEMA))
            document_row_count += len(records)
        writer.close()
        writer = None
    except Exception:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)
        temporary_documents.unlink(missing_ok=True)
        temporary_manifest.unlink(missing_ok=True)
        raise
    finally:
        connection.close()

    forbidden_reads = [
        f"{table}.{column}"
        for table, column in sorted(reads)
        if column in AFFINITY_VALUE_COLUMNS
    ]
    if forbidden_reads:
        temporary.unlink(missing_ok=True)
        temporary_documents.unlink(missing_ok=True)
        raise RuntimeError(f"affinity firewall violation: {forbidden_reads}")

    allowlist_payload = {
        table: sorted(columns) for table, columns in sorted(READ_ALLOWLIST.items())
    }
    projection_bytes = temporary.stat().st_size
    projection_digest = file_sha256(temporary)
    document_bytes = temporary_documents.stat().st_size
    document_digest = file_sha256(temporary_documents)
    manifest = {
        "schema": "a2s-chembl-historical-activity-projection-v1",
        "release": spec.release,
        "index_date": spec.index_date,
        "public_file_date": spec.public_file_date,
        "archive": archive_record,
        "database": {
            "path": str(database),
            "bytes": database.stat().st_size,
            "sha256_before_open": database_digest,
            "schema_columns_sha256": schema_digest,
            "certificate": database_certificate,
            "version_rows": version_rows,
        },
        "projection": {
            "path": str(output),
            "bytes": projection_bytes,
            "sha256": projection_digest,
            "canonical_rows_sha256": canonical_rows.hexdigest(),
            "rows": row_count,
            "columns": PROJECTION_SCHEMA.names,
            "sql_sha256": sha256(PROJECTION_SQL.encode("utf-8")).hexdigest(),
            "null_counts": null_counts,
            "conflict_counts": conflict_counts,
            "filters": [
                "standard_type = Ki",
                "standard_relation = =",
                "target_type = SINGLE PROTEIN",
                "target_organism = Homo sapiens",
            ],
        },
        "document_projection": {
            "path": str(document_output),
            "bytes": document_bytes,
            "sha256": document_digest,
            "canonical_rows_sha256": canonical_documents.hexdigest(),
            "rows": document_row_count,
            "columns": DOCUMENT_SCHEMA.names,
            "sql_sha256": sha256(DOCUMENT_SQL.encode("utf-8")).hexdigest(),
        },
        "firewall": {
            "mode": "main-only SELECT with deny-by-default action and column allowlists",
            "allowlist_sha256": sha256(
                json.dumps(
                    allowlist_payload, sort_keys=True, separators=(",", ":")
                ).encode("utf-8")
            ).hexdigest(),
            "read_columns": [f"{table}.{column}" for table, column in sorted(reads)],
            "numeric_affinity_columns_read": forbidden_reads,
            "denied_actions_during_projection": [list(item) for item in sorted(denied)],
            "affinity_values_materialized": 0,
            "policy_note": (
                "standard relation and units are declared metadata fields in the "
                "H0-S contract; numeric and textual outcome fields remain denied"
            ),
        },
        "identity": {
            "native": "chembl:activity:{activity_id}",
            "secondary": "SHA-256 of source/document/assay/compound/target/endpoint metadata",
            "warning": (
                "cross-release stability and secondary-key collisions must pass an "
                "explicit audit before first_seen_release is assigned"
            ),
        },
        "reproducibility": {
            "implementation_sha256": file_sha256(Path(__file__)),
            "python": platform.python_version(),
            "sqlite_runtime": sqlite3.sqlite_version,
            "pyarrow": pa.__version__,
            "independent_replay_verified": False,
        },
    }
    temporary_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    published_output = False
    published_documents = False
    try:
        temporary.replace(output)
        published_output = True
        temporary_documents.replace(document_output)
        published_documents = True
        temporary_manifest.replace(manifest_path)
    except Exception:
        temporary.unlink(missing_ok=True)
        temporary_documents.unlink(missing_ok=True)
        temporary_manifest.unlink(missing_ok=True)
        if published_output and not manifest_path.exists():
            output.unlink(missing_ok=True)
        if published_documents and not manifest_path.exists():
            document_output.unlink(missing_ok=True)
        raise
    return manifest


def certify_independent_replay(
    database: Path,
    output: Path,
    spec: SnapshotSpec,
    archive_record: dict[str, Any],
) -> dict[str, Any]:
    manifest_path = output.with_suffix(".manifest.json")
    certificate_path = output.with_suffix(".replay.json")
    if not output.exists() or not manifest_path.exists():
        raise FileNotFoundError("primary projection and manifest are required for replay")
    if certificate_path.exists():
        raise FileExistsError(f"refusing to overwrite replay certificate: {certificate_path}")
    primary = json.loads(manifest_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix=".replay-", dir=output.parent) as directory:
        replay_output = Path(directory) / output.name
        replay = project_snapshot(
            database,
            replay_output,
            spec,
            archive_record,
            full_integrity=False,
        )

    comparisons = {
        "database_sha256": (
            primary["database"]["sha256_before_open"]
            == replay["database"]["sha256_before_open"]
        ),
        "primary_integrity_check": (
            primary["database"]["certificate"]["integrity_check"] == "ok"
            and replay["database"]["certificate"]["integrity_check"]
            == "inherited_from_hash-matched_primary_projection"
        ),
        "activity_rows": (
            primary["projection"]["rows"] == replay["projection"]["rows"]
        ),
        "activity_canonical_rows_sha256": (
            primary["projection"]["canonical_rows_sha256"]
            == replay["projection"]["canonical_rows_sha256"]
        ),
        "activity_parquet_sha256": (
            primary["projection"]["sha256"] == replay["projection"]["sha256"]
        ),
        "activity_null_counts": (
            primary["projection"]["null_counts"]
            == replay["projection"]["null_counts"]
        ),
        "activity_conflict_counts": (
            primary["projection"]["conflict_counts"]
            == replay["projection"]["conflict_counts"]
        ),
        "document_rows": (
            primary["document_projection"]["rows"]
            == replay["document_projection"]["rows"]
        ),
        "document_canonical_rows_sha256": (
            primary["document_projection"]["canonical_rows_sha256"]
            == replay["document_projection"]["canonical_rows_sha256"]
        ),
        "document_parquet_sha256": (
            primary["document_projection"]["sha256"]
            == replay["document_projection"]["sha256"]
        ),
        "firewall_reads": primary["firewall"]["read_columns"]
        == replay["firewall"]["read_columns"],
        "firewall_denied_actions": (
            primary["firewall"]["denied_actions_during_projection"]
            == replay["firewall"]["denied_actions_during_projection"]
        ),
    }
    if not all(comparisons.values()):
        raise ValueError(f"independent projection replay mismatch: {comparisons}")
    certificate = {
        "schema": "a2s-chembl-historical-projection-replay-v1",
        "release": spec.release,
        "status": "PASS",
        "primary_manifest_path": str(manifest_path),
        "primary_manifest_sha256": file_sha256(manifest_path),
        "database_sha256": replay["database"]["sha256_before_open"],
        "comparisons": comparisons,
    }
    temporary = certificate_path.with_name(f".{certificate_path.name}.partial")
    temporary.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(certificate_path)
    return certificate


def run_release(
    spec: SnapshotSpec,
    *,
    verify_only: bool = False,
    extract_only: bool = False,
    replay: bool = False,
) -> dict[str, Any]:
    archive = HISTORICAL_ROOT / "archives" / spec.archive_name
    archive_record = verify_archive(archive, spec)
    archive_record["official_url"] = spec.archive_url
    archive_record["publisher_metadata"] = verify_publisher_metadata(spec)
    if verify_only:
        return {"release": spec.release, "archive": archive_record, "verified": True}

    snapshot_dir = HISTORICAL_ROOT / "snapshots" / spec.release
    database = extract_single_sqlite(archive, snapshot_dir)
    if extract_only:
        return {
            "release": spec.release,
            "archive": archive_record,
            "database": {
                "path": str(database),
                "bytes": database.stat().st_size,
                "sha256": file_sha256(database),
                "certified": False,
            },
        }
    output = DEFAULT_OUTPUT_ROOT / spec.release / "activity_identity.parquet"
    if replay and output.exists() and output.with_suffix(".manifest.json").exists():
        return certify_independent_replay(
            database, output, spec, archive_record
        )
    report = project_snapshot(database, output, spec, archive_record)
    if replay:
        report["replay_certificate"] = certify_independent_replay(
            database, output, spec, archive_record
        )
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--release", action="append", choices=sorted(SNAPSHOTS), dest="releases"
    )
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--extract-only", action="store_true")
    parser.add_argument(
        "--replay",
        action="store_true",
        help="repeat an existing/new projection and freeze an equality certificate",
    )
    return parser.parse_args()


def main() -> None:
    from research.shared.dataset_processing import require_run_context

    require_run_context(
        [str(Path.cwd() / "main.py"), "historical-project", *sys.argv[1:]]
    )
    args = parse_args()
    releases = args.releases or sorted(SNAPSHOTS)
    reports = [
        run_release(
            SNAPSHOTS[release],
            verify_only=args.verify_only,
            extract_only=args.extract_only,
            replay=args.replay,
        )
        for release in releases
    ]
    print(json.dumps(reports, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
