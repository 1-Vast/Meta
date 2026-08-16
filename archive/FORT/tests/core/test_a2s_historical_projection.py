from __future__ import annotations

from hashlib import md5, sha256
import json
from pathlib import Path
import sqlite3
import tarfile

import pyarrow.parquet as pq
import pytest

from research.a2s import a2s_historical_projection
from research.a2s.a2s_historical_projection import (
    AFFINITY_VALUE_COLUMNS,
    READ_ALLOWLIST,
    SNAPSHOTS,
    SnapshotSpec,
    _install_read_firewall,
    _validate_version,
    certify_independent_replay,
    extract_single_sqlite,
    project_snapshot,
    verify_archive,
)


def test_direct_projection_module_entry_is_blocked(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("A2S_DATASET_RUN_DIR", raising=False)
    monkeypatch.delenv("A2S_DATASET_RUN_ID", raising=False)
    with pytest.raises(SystemExit, match="dataset-run"):
        a2s_historical_projection.main()


def _create_snapshot(path: Path) -> None:
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE version (name TEXT, creation_date TEXT);
        CREATE TABLE activities (
          activity_id INTEGER PRIMARY KEY, assay_id INTEGER, doc_id INTEGER,
          record_id INTEGER, molregno INTEGER, standard_relation TEXT,
          standard_value REAL, standard_units TEXT, standard_type TEXT,
          src_id INTEGER
        );
        CREATE TABLE assays (
          assay_id INTEGER PRIMARY KEY, doc_id INTEGER, description TEXT,
          tid INTEGER, src_id INTEGER, src_assay_id TEXT, chembl_id TEXT,
          aidx TEXT
        );
        CREATE TABLE docs (
          doc_id INTEGER PRIMARY KEY, chembl_id TEXT, year INTEGER,
          pubmed_id INTEGER, doi TEXT, patent_id TEXT, doc_type TEXT,
          src_id INTEGER, ridx TEXT
        );
        CREATE TABLE compound_records (
          record_id INTEGER PRIMARY KEY, molregno INTEGER, doc_id INTEGER,
          src_id INTEGER, src_compound_id TEXT, cidx TEXT
        );
        CREATE TABLE molecule_hierarchy (molregno INTEGER, parent_molregno INTEGER);
        CREATE TABLE molecule_dictionary (molregno INTEGER PRIMARY KEY, chembl_id TEXT);
        CREATE TABLE compound_structures (molregno INTEGER PRIMARY KEY, standard_inchi_key TEXT);
        CREATE TABLE target_dictionary (
          tid INTEGER PRIMARY KEY, target_type TEXT, tax_id INTEGER,
          organism TEXT, chembl_id TEXT
        );
        CREATE TABLE target_components (tid INTEGER, component_id INTEGER, homologue INTEGER);
        CREATE TABLE component_sequences (
          component_id INTEGER PRIMARY KEY, component_type TEXT,
          accession TEXT, sequence_md5sum TEXT, tax_id INTEGER,
          db_source TEXT, db_version TEXT
        );
        CREATE TABLE source (src_id INTEGER PRIMARY KEY, src_short_name TEXT);

        INSERT INTO version VALUES ('ChEMBL 24.1', '2018-06-18');
        INSERT INTO source VALUES (1, 'LITERATURE');
        INSERT INTO docs VALUES (1, 'CHEMBL_DOC_1', 2017, 123, '10.1/EXAMPLE', NULL, 'PUBLICATION', 1, 'RIDX1');
        INSERT INTO target_dictionary VALUES (1, 'SINGLE PROTEIN', 9606, 'Homo sapiens', 'CHEMBL_T1');
        INSERT INTO component_sequences VALUES (1, 'PROTEIN', 'P12345', '973eb56c8acaa2458cd7beae3af41781', 9606, 'SWISS-PROT', '1');
        INSERT INTO target_components VALUES (1, 1, 0);
        INSERT INTO assays VALUES (1, 1, 'binding assay', 1, 1, 'A1', 'CHEMBL_A1', 'AIDX1');
        INSERT INTO compound_records VALUES (1, 10, 1, 1, 'C1', 'CIDX1');
        INSERT INTO molecule_hierarchy VALUES (10, 10);
        INSERT INTO molecule_dictionary VALUES (10, 'CHEMBL_M1');
        INSERT INTO compound_structures VALUES (10, 'ABCDEFGHIJKLMN-AAAAAA-B');
        INSERT INTO activities VALUES (7, 1, 1, 1, 10, '=', 123456789.25, 'nM', 'Ki', 1);
        INSERT INTO activities VALUES (8, 1, 1, 1, 10, '=', 987654321.75, 'nM', 'Kd', 1);
        """
    )
    connection.commit()
    connection.close()


def _spec_for(path: Path) -> SnapshotSpec:
    return SnapshotSpec(
        release="chembl_24_1",
        index_date="2018-12-31",
        public_file_date="2018-01-01",
        archive_bytes=path.stat().st_size,
        archive_sha256=sha256(path.read_bytes()).hexdigest(),
    )


def test_archive_verification_requires_size_and_sha256(tmp_path: Path) -> None:
    archive = tmp_path / "snapshot.tar.gz"
    archive.write_bytes(b"verified bytes")
    spec = _spec_for(archive)

    record = verify_archive(archive, spec)
    assert record["sha256"] == spec.archive_sha256

    wrong = SnapshotSpec(
        spec.release,
        spec.index_date,
        spec.public_file_date,
        spec.archive_bytes,
        "0" * 64,
    )
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        verify_archive(archive, wrong)


def test_chembl_24_1_uses_exact_internal_major_version_name(tmp_path: Path) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    connection = sqlite3.connect(database)
    connection.execute("UPDATE version SET name = 'ChEMBL_24'")
    connection.commit()

    rows = _validate_version(connection, SNAPSHOTS["chembl_24_1"])
    assert rows[0]["name"] == "ChEMBL_24"

    connection.execute("UPDATE version SET name = 'ChEMBL_24_1'")
    connection.commit()
    with pytest.raises(ValueError, match="expected.*ChEMBL_24"):
        _validate_version(connection, SNAPSHOTS["chembl_24_1"])
    connection.close()


def test_extractor_ignores_non_database_members(tmp_path: Path) -> None:
    database = tmp_path / "source.db"
    database.write_bytes(b"sqlite placeholder")
    note = tmp_path / "README.txt"
    note.write_text("not extracted", encoding="ascii")
    archive = tmp_path / "snapshot.tar.gz"
    with tarfile.open(archive, "w:gz") as output:
        output.add(note, arcname="../../README.txt")
        output.add(database, arcname="nested/chembl_test.db")

    extracted = extract_single_sqlite(archive, tmp_path / "snapshot")

    assert extracted.name == "chembl_test.db"
    assert extracted.read_bytes() == b"sqlite placeholder"
    assert not (tmp_path / "snapshot" / "README.txt").exists()


def test_projection_is_outcome_blind_and_deterministic(tmp_path: Path) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    archive = tmp_path / "archive.tar.gz"
    archive.write_bytes(b"archive")
    spec = _spec_for(archive)
    output = tmp_path / "activity_identity.parquet"

    manifest = project_snapshot(
        database,
        output,
        spec,
        verify_archive(archive, spec),
        batch_size=1,
    )
    table = pq.read_table(output)
    row = table.to_pylist()[0]

    assert table.num_rows == 1
    assert row["activity_id"] == 7
    assert row["endpoint_type"] == "Ki"
    assert row["target_uniprot"] == "P12345"
    assert row["target_sequence_hash"] == md5(b"ACDE").hexdigest()
    assert row["target_sequence_hash_algorithm"] == "MD5"
    assert row["connectivity_key"] == "ABCDEFGHIJKLMN"
    assert row["normalized_doi"] == "10.1/example"
    assert row["document_id_consistent"] is True
    assert manifest["firewall"]["numeric_affinity_columns_read"] == []
    assert manifest["firewall"]["affinity_values_materialized"] == 0
    assert all(
        "standard_value" not in name
        for name in manifest["firewall"]["read_columns"]
    )
    assert "standard_value" not in table.column_names


def test_sqlite_firewall_denies_affinity_value_reads(tmp_path: Path) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    connection = sqlite3.connect(database)
    reads: set[tuple[str, str]] = set()
    denied: set[tuple[int, str, str]] = set()
    _install_read_firewall(connection, reads, denied)

    with pytest.raises(sqlite3.DatabaseError, match="prohibited"):
        connection.execute("SELECT standard_value FROM activities").fetchall()

    connection.close()
    assert not AFFINITY_VALUE_COLUMNS.intersection(
        set().union(*READ_ALLOWLIST.values())
    )
    assert ("activities", "standard_value") not in reads
    assert denied


@pytest.mark.parametrize(
    "statement",
    [
        "SELECT activity_id FROM activities WHERE standard_value > 0",
        "SELECT SUM(standard_value) FROM activities",
        "SELECT standard_value + 1 FROM activities",
        "SELECT activity_id FROM activities ORDER BY standard_value",
        "SELECT activity_id FROM activities GROUP BY activity_id HAVING MAX(standard_value) > 0",
        "SELECT activity_id FROM activities WHERE activity_id IN (SELECT activity_id FROM activities WHERE standard_value > 0)",
        "SELECT * FROM activities",
        "UPDATE activities SET standard_value = 0",
        "PRAGMA table_info(activities)",
        "SELECT random()",
    ],
)
def test_sqlite_firewall_blocks_indirect_reads_and_other_actions(
    tmp_path: Path, statement: str
) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    connection = sqlite3.connect(database)
    _install_read_firewall(connection, set(), set())

    with pytest.raises(sqlite3.DatabaseError):
        connection.execute(statement).fetchall()

    connection.close()


def test_sqlite_firewall_blocks_views_attach_and_schema_changes(
    tmp_path: Path,
) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    connection = sqlite3.connect(database)
    connection.execute(
        "CREATE VIEW affinity_view AS SELECT activity_id, standard_value FROM activities"
    )
    _install_read_firewall(connection, set(), set())

    statements = [
        "SELECT * FROM affinity_view",
        f"ATTACH DATABASE '{(tmp_path / 'extra.db').as_posix()}' AS extra",
        "CREATE TABLE extra_table (value REAL)",
    ]
    for statement in statements:
        with pytest.raises(sqlite3.DatabaseError):
            connection.execute(statement).fetchall()

    connection.close()


def test_projection_canonical_hash_replays(tmp_path: Path) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    archive = tmp_path / "archive.tar.gz"
    archive.write_bytes(b"archive")
    spec = _spec_for(archive)
    archive_record = verify_archive(archive, spec)

    first = project_snapshot(
        database, tmp_path / "first.parquet", spec, archive_record, batch_size=1
    )
    second = project_snapshot(
        database, tmp_path / "second.parquet", spec, archive_record, batch_size=2
    )

    assert first["projection"]["canonical_rows_sha256"] == second["projection"][
        "canonical_rows_sha256"
    ]
    assert first["projection"]["sha256"] == second["projection"]["sha256"]


def test_independent_replay_certificate_passes(tmp_path: Path) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    archive = tmp_path / "archive.tar.gz"
    archive.write_bytes(b"archive")
    spec = _spec_for(archive)
    archive_record = verify_archive(archive, spec)
    output = tmp_path / "activity_identity.parquet"
    project_snapshot(database, output, spec, archive_record)

    certificate = certify_independent_replay(
        database, output, spec, archive_record
    )

    assert certificate["status"] == "PASS"
    assert all(certificate["comparisons"].values())
    assert output.with_suffix(".replay.json").exists()


@pytest.mark.parametrize(
    ("section", "field", "changed", "comparison"),
    [
        ("projection", "null_counts", {"parent_inchikey": 999}, "activity_null_counts"),
        (
            "projection",
            "conflict_counts",
            {"document_id": 999, "molecule_id": 0},
            "activity_conflict_counts",
        ),
        (
            "firewall",
            "denied_actions_during_projection",
            [[20, "activities", "standard_value"]],
            "firewall_denied_actions",
        ),
    ],
)
def test_independent_replay_rejects_audit_semantic_mismatch(
    tmp_path: Path,
    section: str,
    field: str,
    changed: object,
    comparison: str,
) -> None:
    database = tmp_path / "chembl_test.db"
    _create_snapshot(database)
    archive = tmp_path / "archive.tar.gz"
    archive.write_bytes(b"archive")
    spec = _spec_for(archive)
    archive_record = verify_archive(archive, spec)
    output = tmp_path / "activity_identity.parquet"
    project_snapshot(database, output, spec, archive_record)

    manifest_path = output.with_suffix(".manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest[section][field] = changed
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    with pytest.raises(ValueError, match=comparison):
        certify_independent_replay(database, output, spec, archive_record)

    assert not output.with_suffix(".replay.json").exists()


def test_projection_refuses_missing_allowlisted_schema(tmp_path: Path) -> None:
    database = tmp_path / "broken.db"
    _create_snapshot(database)
    connection = sqlite3.connect(database)
    connection.execute("ALTER TABLE activities RENAME TO activities_old")
    connection.commit()
    connection.close()
    archive = tmp_path / "archive.tar.gz"
    archive.write_bytes(b"archive")
    spec = _spec_for(archive)

    with pytest.raises(ValueError, match="required ChEMBL table is missing"):
        project_snapshot(
            database,
            tmp_path / "projection.parquet",
            spec,
            verify_archive(archive, spec),
        )
