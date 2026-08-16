"""Audit stable ChEMBL activity identity and first-observed snapshot presence."""

from __future__ import annotations

import argparse
from array import array
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import groupby
import json
from pathlib import Path
import re
import sqlite3
import sys
import tempfile
from typing import Any, Iterable

import pyarrow as pa
import pyarrow.parquet as pq

from research.a2s.a2s_historical_projection import SNAPSHOTS, file_sha256


PROJECTION_ROOT = Path("dataset/processed/a2s_historical_projection")
DEFAULT_OUTPUT = Path(
    "dataset/processed/a2s_historical_presence/activity_presence.v1.parquet"
)

RELEASES = tuple(SNAPSHOTS)

ACTIVITY_COLUMNS = [
    "activity_id",
    "native_activity_key",
    "measurement_identity_sha256",
    "activity_src_id",
    "activity_source",
    "document_source_id",
    "document_chembl_id",
    "normalized_doi",
    "pubmed_id",
    "patent_id",
    "document_index",
    "assay_chembl_id",
    "assay_src_id",
    "assay_index",
    "assay_source_id",
    "compound_index",
    "record_source_id",
    "source_compound_id",
    "parent_molecule_chembl_id",
    "parent_inchikey",
    "target_chembl_id",
    "target_uniprot",
    "target_component_count",
    "endpoint_type",
    "standard_relation",
    "standard_units",
    "publication_year",
    "document_id_consistent",
    "molecule_id_consistent",
]

PRESENCE_SCHEMA = pa.schema(
    [
        ("identity_key", pa.string()),
        ("measurement_identity_sha256", pa.string()),
        ("match_class", pa.string()),
        ("presence_mask", pa.string()),
        ("first_seen_release", pa.string()),
        ("first_seen_public_file_date", pa.string()),
        ("native_activity_ids_json", pa.string()),
        ("release_activity_ids_json", pa.string()),
        ("persistent", pa.bool_()),
        ("native_collision", pa.bool_()),
        ("release_ambiguous", pa.bool_()),
        ("secondary_ambiguous", pa.bool_()),
        ("lineage_drift", pa.bool_()),
        ("parent_drift", pa.bool_()),
        ("target_drift", pa.bool_()),
        ("relation_unit_drift", pa.bool_()),
        ("document_metadata_drift", pa.bool_()),
        ("document_ids_consistent", pa.bool_()),
        ("molecule_ids_consistent", pa.bool_()),
        ("single_target_component", pa.bool_()),
        ("document_chembl_id", pa.string()),
        ("publication_year", pa.int64()),
        ("normalized_doi", pa.string()),
        ("pubmed_id", pa.int64()),
        ("patent_id", pa.string()),
        ("activity_source", pa.string()),
        ("activity_src_id", pa.int64()),
        ("document_source_id", pa.int64()),
        ("assay_chembl_id", pa.string()),
        ("assay_index", pa.string()),
        ("assay_source_id", pa.string()),
        ("assay_src_id", pa.int64()),
        ("record_source_id", pa.int64()),
        ("compound_index", pa.string()),
        ("source_compound_id", pa.string()),
        ("target_chembl_id", pa.string()),
        ("target_uniprot", pa.string()),
        ("parent_molecule_chembl_id", pa.string()),
        ("parent_inchikey", pa.string()),
        ("endpoint_type", pa.string()),
        ("standard_relation", pa.string()),
        ("standard_units", pa.string()),
        ("first_seen_classification", pa.string()),
        ("stable_identity_eligible", pa.bool_()),
        ("eligibility_reasons_json", pa.string()),
    ]
)


class UnionFind:
    def __init__(self, size: int) -> None:
        if size >= 2**32:
            raise ValueError("presence audit exceeds uint32 occurrence capacity")
        self.parent = array("I", range(size))
        self.rank = array("B", [0]) * size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def _records(path: Path, columns: list[str]) -> Iterable[dict[str, Any]]:
    parquet = pq.ParquetFile(path)
    missing = set(columns).difference(parquet.schema_arrow.names)
    if missing:
        raise ValueError(f"{path} lacks required fields: {sorted(missing)}")
    for batch in parquet.iter_batches(columns=columns, batch_size=50_000):
        yield from batch.to_pylist()


def _lineage_signature(row: dict[str, Any]) -> tuple[Any, ...]:
    return tuple(
        row[name]
        for name in (
            "activity_source",
            "activity_src_id",
            "document_source_id",
            "document_chembl_id",
            "document_index",
            "assay_chembl_id",
            "assay_src_id",
            "assay_index",
            "assay_source_id",
            "compound_index",
            "record_source_id",
            "source_compound_id",
            "endpoint_type",
        )
    )


def _persistent(mask: str) -> bool:
    return re.fullmatch(r"0*1+", mask) is not None


def _one_value(rows: list[dict[str, Any]], name: str) -> tuple[Any, bool]:
    values = {row[name] for row in rows if row[name] is not None}
    return (next(iter(values)) if len(values) == 1 else None, len(values) > 1)


def _first_seen_classification(
    first_index: int,
    row: dict[str, Any],
    document_presence: dict[str, set[str]],
) -> str:
    if first_index == 0:
        return "left_censored_at_chembl_24_1"
    document_tokens = _document_tokens(row)
    earlier_releases = RELEASES[:first_index]
    if document_tokens and any(
        document_tokens.intersection(document_presence[release])
        for release in earlier_releases
    ):
        return "activity_backfill_under_preexisting_document"
    year = row["publication_year"]
    if year is None:
        return "unknown_publication_time"
    previous_year = int(SNAPSHOTS[RELEASES[first_index - 1]].index_date[:4])
    if int(year) <= previous_year:
        return "retrospective_document_ingest"
    return "contemporaneous_candidate"


def _document_tokens(row: dict[str, Any]) -> set[str]:
    tokens: set[str] = set()
    if row.get("document_chembl_id") is not None:
        tokens.add(f"chembl:{str(row['document_chembl_id']).strip().upper()}")
    if row.get("normalized_doi") is not None:
        tokens.add(f"doi:{str(row['normalized_doi']).strip().lower()}")
    if row.get("pubmed_id") is not None:
        tokens.add(f"pmid:{int(row['pubmed_id'])}")
    if row.get("patent_id") is not None:
        patent = re.sub(r"\s+", "", str(row["patent_id"]).upper())
        if patent:
            tokens.add(f"patent:{patent}")
    return tokens


def _identity_record(
    rows: list[dict[str, Any]], document_presence: dict[str, set[str]]
) -> dict[str, Any]:
    rows.sort(key=lambda row: (row["release_index"], int(row["activity_id"])))
    by_release: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_release[row["release"]].append(row)
    mask = "".join("1" if by_release[release] else "0" for release in RELEASES)
    first_index = mask.index("1")
    first = rows[0]
    release_ambiguous = any(len(values) > 1 for values in by_release.values())
    collision = any(row["native_collision"] for row in rows)
    secondary_collision = any(row["secondary_ambiguous"] for row in rows)

    lineage_drift = len({_lineage_signature(row) for row in rows}) > 1
    parent_values = {
        (row["parent_molecule_chembl_id"], row["parent_inchikey"]) for row in rows
    }
    target_values = {
        (row["target_chembl_id"], row["target_uniprot"]) for row in rows
    }
    relation_values = {
        (row["standard_relation"], row["standard_units"]) for row in rows
    }
    document, document_drift = _one_value(rows, "document_chembl_id")
    publication_year, year_drift = _one_value(rows, "publication_year")
    normalized_doi, _ = _one_value(rows, "normalized_doi")
    pubmed_id, _ = _one_value(rows, "pubmed_id")
    patent_id, _ = _one_value(rows, "patent_id")
    activity_source, _ = _one_value(rows, "activity_source")
    activity_src_id, _ = _one_value(rows, "activity_src_id")
    document_source_id, _ = _one_value(rows, "document_source_id")
    assay, _ = _one_value(rows, "assay_chembl_id")
    assay_index, _ = _one_value(rows, "assay_index")
    assay_source_id, _ = _one_value(rows, "assay_source_id")
    assay_src_id, _ = _one_value(rows, "assay_src_id")
    record_source_id, _ = _one_value(rows, "record_source_id")
    compound_index, _ = _one_value(rows, "compound_index")
    source_compound_id, _ = _one_value(rows, "source_compound_id")
    target, _ = _one_value(rows, "target_chembl_id")
    uniprot, _ = _one_value(rows, "target_uniprot")
    parent, _ = _one_value(rows, "parent_molecule_chembl_id")
    inchikey, _ = _one_value(rows, "parent_inchikey")
    endpoint, endpoint_drift = _one_value(rows, "endpoint_type")
    relation, relation_drift = _one_value(rows, "standard_relation")
    units, unit_drift = _one_value(rows, "standard_units")

    native_ids = sorted({int(row["activity_id"]) for row in rows})
    if collision:
        match_class = "native_id_collision"
    elif release_ambiguous:
        match_class = "release_ambiguous"
    elif len(native_ids) > 1:
        match_class = "rekeyed_1to1"
    elif len(rows) > 1:
        match_class = "native_stable"
    else:
        match_class = "singleton"
    identity_key = (
        f"chembl:activity:{native_ids[0]}"
        if len(native_ids) == 1
        else f"chembl:measurement:{first['measurement_identity_sha256']}"
    )

    document_consistent = all(row["document_id_consistent"] for row in rows)
    molecule_consistent = all(row["molecule_id_consistent"] for row in rows)
    single_component = all(row["target_component_count"] == 1 for row in rows)
    conditions = {
        "nonpersistent_presence": not _persistent(mask),
        "native_id_collision": collision,
        "release_ambiguous": release_ambiguous,
        "lineage_drift": lineage_drift,
        "parent_drift": len(parent_values) > 1,
        "target_drift": len(target_values) > 1,
        "relation_unit_drift": len(relation_values) > 1,
        "document_metadata_drift": document_drift or year_drift,
        "inconsistent_document_ids": not document_consistent,
        "inconsistent_molecule_ids": not molecule_consistent,
        "not_single_target_component": not single_component,
        "missing_document": document is None,
        "missing_parent": inchikey is None,
        "missing_target_accession": uniprot is None,
        "endpoint_drift": endpoint_drift,
        "relation_drift": relation_drift,
        "unit_drift": unit_drift,
    }
    reasons = [name for name, failed in conditions.items() if failed]
    return {
        "identity_key": identity_key,
        "measurement_identity_sha256": first["measurement_identity_sha256"],
        "match_class": match_class,
        "presence_mask": mask,
        "first_seen_release": RELEASES[first_index],
        "first_seen_public_file_date": SNAPSHOTS[
            RELEASES[first_index]
        ].public_file_date,
        "native_activity_ids_json": json.dumps(native_ids),
        "release_activity_ids_json": json.dumps(
            {
                release: sorted(int(row["activity_id"]) for row in by_release[release])
                for release in RELEASES
            },
            sort_keys=True,
        ),
        "persistent": _persistent(mask),
        "native_collision": collision,
        "release_ambiguous": release_ambiguous,
        "secondary_ambiguous": secondary_collision,
        "lineage_drift": lineage_drift,
        "parent_drift": len(parent_values) > 1,
        "target_drift": len(target_values) > 1,
        "relation_unit_drift": len(relation_values) > 1,
        "document_metadata_drift": document_drift or year_drift,
        "document_ids_consistent": document_consistent,
        "molecule_ids_consistent": molecule_consistent,
        "single_target_component": single_component,
        "document_chembl_id": document,
        "publication_year": publication_year,
        "normalized_doi": normalized_doi,
        "pubmed_id": pubmed_id,
        "patent_id": patent_id,
        "activity_source": activity_source,
        "activity_src_id": activity_src_id,
        "document_source_id": document_source_id,
        "assay_chembl_id": assay,
        "assay_index": assay_index,
        "assay_source_id": assay_source_id,
        "assay_src_id": assay_src_id,
        "record_source_id": record_source_id,
        "compound_index": compound_index,
        "source_compound_id": source_compound_id,
        "target_chembl_id": target,
        "target_uniprot": uniprot,
        "parent_molecule_chembl_id": parent,
        "parent_inchikey": inchikey,
        "endpoint_type": endpoint,
        "standard_relation": relation,
        "standard_units": units,
        "first_seen_classification": _first_seen_classification(
            first_index, first, document_presence
        ),
        "stable_identity_eligible": not reasons,
        "eligibility_reasons_json": json.dumps(reasons),
    }


def _work_database(
    connection: sqlite3.Connection,
    activity_paths: dict[str, Path],
    document_paths: dict[str, Path],
) -> tuple[int, dict[str, set[str]]]:
    integer_columns = {
        "activity_id",
        "activity_src_id",
        "document_source_id",
        "assay_src_id",
        "record_source_id",
        "target_component_count",
        "publication_year",
        "pubmed_id",
        "document_id_consistent",
        "molecule_id_consistent",
    }
    definitions = ", ".join(
        f'"{name}" {"INTEGER" if name in integer_columns else "TEXT"}'
        for name in ACTIVITY_COLUMNS
    )
    connection.executescript(
        f"""
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = FILE;
        CREATE TABLE occurrence (
          occurrence_id INTEGER PRIMARY KEY,
          release_index INTEGER NOT NULL,
          release TEXT NOT NULL,
          {definitions},
          native_collision INTEGER NOT NULL DEFAULT 0,
          secondary_ambiguous INTEGER NOT NULL DEFAULT 0
        );
        CREATE UNIQUE INDEX occurrence_release_activity
          ON occurrence(release_index, activity_id);
        """
    )
    insert_columns = ["occurrence_id", "release_index", "release", *ACTIVITY_COLUMNS]
    insert_sql = (
        "INSERT INTO occurrence ("
        + ",".join(f'"{name}"' for name in insert_columns)
        + ") VALUES ("
        + ",".join("?" for _ in insert_columns)
        + ")"
    )
    document_presence: dict[str, set[str]] = {}
    occurrence_id = 0
    for release_index, release in enumerate(RELEASES):
        pending: list[tuple[Any, ...]] = []
        for row in _records(activity_paths[release], ACTIVITY_COLUMNS):
            pending.append(
                (
                    occurrence_id,
                    release_index,
                    release,
                    *(row[name] for name in ACTIVITY_COLUMNS),
                )
            )
            occurrence_id += 1
            if len(pending) == 10_000:
                connection.executemany(insert_sql, pending)
                pending.clear()
        if pending:
            connection.executemany(insert_sql, pending)

        documents: set[str] = set()
        for row in _records(
            document_paths[release],
            ["document_chembl_id", "normalized_doi", "pubmed_id", "patent_id"],
        ):
            documents.update(_document_tokens(row))
        document_presence[release] = documents
        connection.commit()
    connection.executescript(
        """
        CREATE INDEX occurrence_activity ON occurrence(activity_id);
        CREATE INDEX occurrence_fingerprint
          ON occurrence(measurement_identity_sha256, release_index);
        """
    )
    return occurrence_id, document_presence


def _link_occurrences(
    connection: sqlite3.Connection, occurrence_count: int
) -> UnionFind:
    union = UnionFind(occurrence_count)
    lineage_columns = [
        "activity_source",
        "activity_src_id",
        "document_source_id",
        "document_chembl_id",
        "document_index",
        "assay_chembl_id",
        "assay_src_id",
        "assay_index",
        "assay_source_id",
        "compound_index",
        "record_source_id",
        "source_compound_id",
        "endpoint_type",
    ]
    native_sql = (
        "SELECT occurrence_id, activity_id, "
        + ",".join(f'"{name}"' for name in lineage_columns)
        + " FROM occurrence ORDER BY activity_id, release_index"
    )
    cursor = connection.execute(native_sql)
    for _, iterator in groupby(cursor, key=lambda row: row[1]):
        group = list(iterator)
        first = int(group[0][0])
        for row in group[1:]:
            union.union(first, int(row[0]))
        if len({tuple(row[2:]) for row in group}) > 1:
            connection.executemany(
                "UPDATE occurrence SET native_collision=1 WHERE occurrence_id=?",
                [(int(row[0]),) for row in group],
            )

    cursor = connection.execute(
        "SELECT occurrence_id, measurement_identity_sha256, release_index "
        "FROM occurrence ORDER BY measurement_identity_sha256, release_index"
    )
    for _, iterator in groupby(cursor, key=lambda row: row[1]):
        group = list(iterator)
        release_counts = Counter(int(row[2]) for row in group)
        if any(count > 1 for count in release_counts.values()):
            connection.executemany(
                "UPDATE occurrence SET secondary_ambiguous=1 WHERE occurrence_id=?",
                [(int(row[0]),) for row in group],
            )
            continue
        first = int(group[0][0])
        for row in group[1:]:
            union.union(first, int(row[0]))
    connection.commit()

    connection.execute(
        "CREATE TABLE membership (occurrence_id INTEGER PRIMARY KEY, component_id INTEGER)"
    )
    for offset in range(0, occurrence_count, 100_000):
        stop = min(offset + 100_000, occurrence_count)
        connection.executemany(
            "INSERT INTO membership VALUES (?,?)",
            ((item, union.find(item)) for item in range(offset, stop)),
        )
    connection.execute("CREATE INDEX membership_component ON membership(component_id)")
    connection.commit()
    return union


def _identity_records(
    connection: sqlite3.Connection,
    document_presence: dict[str, set[str]],
) -> Iterable[dict[str, Any]]:
    cursor = connection.execute(
        "SELECT m.component_id, o.* FROM membership AS m "
        "JOIN occurrence AS o ON o.occurrence_id=m.occurrence_id "
        "ORDER BY m.component_id, o.release_index, o.activity_id"
    )
    columns = [description[0] for description in cursor.description]
    for _, iterator in groupby(cursor, key=lambda row: row[0]):
        rows = []
        for values in iterator:
            row = dict(zip(columns, values))
            row.pop("component_id")
            rows.append(row)
        yield _identity_record(rows, document_presence)


def build_presence_audit(
    activity_paths: dict[str, Path],
    document_paths: dict[str, Path],
    output: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    if set(activity_paths) != set(RELEASES) or set(document_paths) != set(RELEASES):
        raise ValueError(f"all releases are required: {RELEASES}")
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = output.with_suffix(".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite frozen presence audit: {output}")

    canonical = sha256()
    temporary = output.with_name(f".{output.name}.partial")
    temporary.unlink(missing_ok=True)
    work_handle = tempfile.NamedTemporaryFile(
        prefix=".presence-work-", suffix=".sqlite", dir=output.parent, delete=False
    )
    work_path = Path(work_handle.name)
    work_handle.close()
    writer: pq.ParquetWriter | None = None
    counts: Counter[str] = Counter()
    masks: Counter[str] = Counter()
    classifications: Counter[str] = Counter()
    eligible = 0
    identity_count = 0
    occurrence_count = 0
    connection = sqlite3.connect(work_path)
    try:
        occurrence_count, document_presence = _work_database(
            connection, activity_paths, document_paths
        )
        _link_occurrences(connection, occurrence_count)
        writer = pq.ParquetWriter(temporary, PRESENCE_SCHEMA, compression="zstd")
        batch: list[dict[str, Any]] = []
        for record in _identity_records(connection, document_presence):
            batch.append(record)
            identity_count += 1
            counts[record["match_class"]] += 1
            masks[record["presence_mask"]] += 1
            classifications[record["first_seen_classification"]] += 1
            eligible += record["stable_identity_eligible"]
            if len(batch) == 50_000:
                for item in batch:
                    canonical.update(
                        (
                            json.dumps(
                                item,
                                ensure_ascii=True,
                                sort_keys=True,
                                separators=(",", ":"),
                            )
                            + "\n"
                        ).encode("utf-8")
                    )
                writer.write_table(pa.Table.from_pylist(batch, schema=PRESENCE_SCHEMA))
                batch.clear()
        if batch:
            for item in batch:
                canonical.update(
                    (
                        json.dumps(
                            item,
                            ensure_ascii=True,
                            sort_keys=True,
                            separators=(",", ":"),
                        )
                        + "\n"
                    ).encode("utf-8")
                )
            writer.write_table(pa.Table.from_pylist(batch, schema=PRESENCE_SCHEMA))
        writer.close()
        writer = None
    except Exception:
        if writer is not None:
            writer.close()
        temporary.unlink(missing_ok=True)
        raise
    finally:
        connection.close()
        work_path.unlink(missing_ok=True)

    hard_collisions = counts["native_id_collision"]
    status = "H0_IDENTITY_PASS" if hard_collisions == 0 else "DATA_NOT_READY"
    manifest = {
        "schema": "a2s-chembl-historical-presence-v1",
        "status": status,
        "releases": list(RELEASES),
        "inputs": {
            release: {
                "activities": str(activity_paths[release]),
                "activities_sha256": file_sha256(activity_paths[release]),
                "documents": str(document_paths[release]),
                "documents_sha256": file_sha256(document_paths[release]),
            }
            for release in RELEASES
        },
        "output": {
            "path": str(output),
            "rows": identity_count,
            "eligible_stable_identities": eligible,
            "canonical_rows_sha256": canonical.hexdigest(),
            "parquet_sha256": file_sha256(temporary),
        },
        "counts": {
            "occurrences": occurrence_count,
            "match_class": dict(sorted(counts.items())),
            "presence_mask": dict(sorted(masks.items())),
            "first_seen_classification": dict(sorted(classifications.items())),
        },
        "firewall": {
            "input_numeric_affinity_columns": [],
            "affinity_values_read": 0,
            "affinity_values_materialized": 0,
        },
        "admission": {
            "hard_native_collisions": hard_collisions,
            "h1_authorized": status == "H0_IDENTITY_PASS",
            "warning": (
                "first_seen_release is earliest observed database presence, not "
                "measurement or publication date"
            ),
        },
        "implementation_sha256": file_sha256(Path(__file__)),
    }
    temporary_manifest = manifest_path.with_name(f".{manifest_path.name}.partial")
    temporary_manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    published = False
    try:
        temporary.replace(output)
        published = True
        temporary_manifest.replace(manifest_path)
    except Exception:
        temporary.unlink(missing_ok=True)
        temporary_manifest.unlink(missing_ok=True)
        if published and not manifest_path.exists():
            output.unlink(missing_ok=True)
        raise
    return manifest


def certify_presence_replay(
    activity_paths: dict[str, Path],
    document_paths: dict[str, Path],
    output: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    manifest_path = output.with_suffix(".manifest.json")
    certificate_path = output.with_suffix(".replay.json")
    if not output.exists() or not manifest_path.exists():
        raise FileNotFoundError("primary presence audit and manifest are required")
    if certificate_path.exists():
        raise FileExistsError(f"refusing to overwrite replay certificate: {certificate_path}")
    primary = json.loads(manifest_path.read_text(encoding="utf-8"))
    with tempfile.TemporaryDirectory(prefix=".presence-replay-", dir=output.parent) as root:
        replay = build_presence_audit(
            activity_paths, document_paths, Path(root) / output.name
        )
    comparisons = {
        "status": primary["status"] == replay["status"],
        "identity_rows": primary["output"]["rows"] == replay["output"]["rows"],
        "eligible_identities": (
            primary["output"]["eligible_stable_identities"]
            == replay["output"]["eligible_stable_identities"]
        ),
        "canonical_rows_sha256": (
            primary["output"]["canonical_rows_sha256"]
            == replay["output"]["canonical_rows_sha256"]
        ),
        "parquet_sha256": (
            primary["output"]["parquet_sha256"]
            == replay["output"]["parquet_sha256"]
        ),
        "counts": primary["counts"] == replay["counts"],
    }
    if not all(comparisons.values()):
        raise ValueError(f"presence replay mismatch: {comparisons}")
    certificate = {
        "schema": "a2s-chembl-historical-presence-replay-v1",
        "status": "PASS",
        "primary_manifest_path": str(manifest_path),
        "primary_manifest_sha256": file_sha256(manifest_path),
        "comparisons": comparisons,
    }
    temporary = certificate_path.with_name(f".{certificate_path.name}.partial")
    temporary.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(certificate_path)
    return certificate


def _certified_default_inputs() -> tuple[dict[str, Path], dict[str, Path]]:
    activities: dict[str, Path] = {}
    documents: dict[str, Path] = {}
    for release in RELEASES:
        root = PROJECTION_ROOT / release
        manifest_path = root / "activity_identity.manifest.json"
        replay_path = root / "activity_identity.replay.json"
        if not manifest_path.exists() or not replay_path.exists():
            raise FileNotFoundError(f"projection certificate missing for {release}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        replay = json.loads(replay_path.read_text(encoding="utf-8"))
        if (
            manifest.get("release") != release
            or manifest.get("firewall", {}).get("numeric_affinity_columns_read") != []
            or replay.get("status") != "PASS"
            or not all(replay.get("comparisons", {}).values())
            or replay.get("primary_manifest_sha256") != file_sha256(manifest_path)
        ):
            raise ValueError(f"projection is not certified for {release}")
        activities[release] = Path(manifest["projection"]["path"])
        documents[release] = Path(manifest["document_projection"]["path"])
        if (
            file_sha256(activities[release]) != manifest["projection"]["sha256"]
            or file_sha256(documents[release])
            != manifest["document_projection"]["sha256"]
        ):
            raise ValueError(f"projection artifact hash mismatch for {release}")
    return activities, documents


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replay", action="store_true")
    return parser.parse_args()


def main() -> None:
    from research.shared.dataset_processing import require_run_context

    require_run_context(
        [str(Path.cwd() / "main.py"), "historical-presence", *sys.argv[1:]]
    )
    args = parse_args()
    activities, documents = _certified_default_inputs()
    if args.replay and args.output.exists():
        report = certify_presence_replay(activities, documents, args.output)
    else:
        report = build_presence_audit(activities, documents, args.output)
        if args.replay:
            report["replay_certificate"] = certify_presence_replay(
                activities, documents, args.output
            )
    summary = {"status": report["status"]}
    if "admission" in report:
        summary.update(report["admission"])
    else:
        summary["comparisons"] = report["comparisons"]
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
