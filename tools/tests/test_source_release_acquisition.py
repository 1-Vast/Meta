import hashlib
from pathlib import Path
import tarfile

import pytest

from scripts.source_affinity.chembl_static import (
    parse_checksum_listing,
    safe_extract_archive,
    sqlite_schema_sha256,
    verify_archive,
)


def test_parse_checksum_listing_selects_exact_archive():
    digest = "a" * 64
    assert parse_checksum_listing(f"SHA256 file\n{digest}\tchembl_37_sqlite.tar.gz\n") == digest


def test_verify_archive_fails_closed_on_hash_mismatch(tmp_path):
    archive = tmp_path / "archive.tar.gz"
    archive.write_bytes(b"payload")
    with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
        verify_archive(archive, "0" * 64)
    expected = hashlib.sha256(b"payload").hexdigest()
    assert verify_archive(archive, expected)["sha256"] == expected


def test_safe_extract_archive_accepts_regular_sqlite(tmp_path):
    source = tmp_path / "chembl_37.db"
    source.write_bytes(b"SQLite format 3\x00")
    archive = tmp_path / "release.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        bundle.add(source, arcname="chembl_37/chembl_37_sqlite/chembl_37.db")
    output = tmp_path / "extracted"
    files = safe_extract_archive(archive, output)
    assert [path.relative_to(output).as_posix() for path in files] == [
        "chembl_37/chembl_37_sqlite/chembl_37.db"
    ]


def test_safe_extract_archive_rejects_path_traversal(tmp_path):
    source = tmp_path / "source.db"
    source.write_bytes(b"SQLite format 3\x00")
    archive = tmp_path / "release.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        bundle.add(source, arcname="../../escape.db")
    with pytest.raises(RuntimeError, match="unsafe archive member"):
        safe_extract_archive(archive, tmp_path / "extracted")
    assert not (tmp_path / "escape.db").exists()


def test_sqlite_schema_hash_is_deterministic(tmp_path):
    import sqlite3

    database = tmp_path / "test.db"
    connection = sqlite3.connect(database)
    connection.execute("CREATE TABLE rows (id INTEGER PRIMARY KEY, value TEXT)")
    connection.commit()
    connection.close()
    assert sqlite_schema_sha256(database) == sqlite_schema_sha256(database)
