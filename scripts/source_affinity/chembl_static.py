"""Pinned ChEMBL37 SQLite release contract and acquisition helpers."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import shutil
import sqlite3
import tarfile
import time
from typing import Callable

import requests

from scripts.source_affinity.common import sha256_file, write_canonical_json


RELEASE = "37"
ARCHIVE_NAME = "chembl_37_sqlite.tar.gz"
RELEASE_ROOT = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37"
ARCHIVE_URL = f"{RELEASE_ROOT}/{ARCHIVE_NAME}"
CHECKSUM_URL = f"{RELEASE_ROOT}/checksums.txt"
EXPECTED_SHA256 = "33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281"
LICENSE = "CC BY-SA 3.0"
LICENSE_URL = f"{RELEASE_ROOT}/LICENSE"


def parse_checksum_listing(text: str, filename: str = ARCHIVE_NAME) -> str:
    for line in text.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[1] == filename:
            value = fields[0].lower()
            if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError(f"invalid SHA-256 for {filename}")
            return value
    raise ValueError(f"official checksum does not list {filename}")


def fetch_official_checksum(timeout: int = 60) -> tuple[str, bytes]:
    response = requests.get(CHECKSUM_URL, timeout=timeout)
    response.raise_for_status()
    payload = response.content
    return parse_checksum_listing(payload.decode("utf-8")), payload


def download_resumable(
    url: str,
    destination: Path,
    expected_bytes: int | None = None,
    progress: Callable[[int, int | None], None] | None = None,
) -> None:
    partial = destination.with_suffix(destination.suffix + ".partial")
    partial.parent.mkdir(parents=True, exist_ok=True)
    offset = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    with requests.get(url, headers=headers, stream=True, timeout=(30, 120)) as response:
        if offset and response.status_code == 200:
            offset = 0
        elif offset and response.status_code != 206:
            response.raise_for_status()
            raise RuntimeError(f"server rejected byte-range resume: HTTP {response.status_code}")
        else:
            response.raise_for_status()
        mode = "ab" if offset else "wb"
        total = expected_bytes
        if total is None and response.headers.get("Content-Length"):
            total = offset + int(response.headers["Content-Length"])
        written = offset
        with partial.open(mode) as handle:
            for block in response.iter_content(chunk_size=8 * 1024 * 1024):
                if not block:
                    continue
                handle.write(block)
                written += len(block)
                if progress is not None:
                    progress(written, total)
    if expected_bytes is not None and written != expected_bytes:
        raise RuntimeError(f"archive byte mismatch: expected={expected_bytes}, observed={written}")
    partial.replace(destination)


def verify_archive(path: Path, expected_sha256: str = EXPECTED_SHA256) -> dict[str, int | str]:
    observed = sha256_file(path)
    if observed != expected_sha256:
        raise RuntimeError(
            f"archive SHA-256 mismatch: expected={expected_sha256}, observed={observed}"
        )
    return {"path": path.name, "bytes": path.stat().st_size, "sha256": observed}


def safe_extract_archive(archive: Path, destination: Path) -> list[Path]:
    if destination.exists():
        raise FileExistsError(f"extraction destination already exists: {destination}")
    temporary = destination.with_name(destination.name + ".partial")
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True)
    root = temporary.resolve()
    try:
        with tarfile.open(archive, "r:gz") as bundle:
            for member in bundle.getmembers():
                target = (temporary / member.name).resolve()
                if root != target and root not in target.parents:
                    raise RuntimeError(f"unsafe archive member path: {member.name}")
                if not (member.isfile() or member.isdir()):
                    raise RuntimeError(f"unsupported archive member type: {member.name}")
            bundle.extractall(temporary)
        temporary.replace(destination)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise
    return sorted(path for path in destination.rglob("*") if path.is_file())


def sqlite_schema_sha256(path: Path) -> str:
    import hashlib
    import json

    connection = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
    try:
        rows = connection.execute(
            "SELECT type, name, tbl_name, sql FROM sqlite_master "
            "WHERE name NOT LIKE 'sqlite_%' ORDER BY type, name, tbl_name"
        ).fetchall()
    finally:
        connection.close()
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), default=str).encode()
    return hashlib.sha256(payload).hexdigest()


def acquire_release(output_dir: Path, expected_bytes: int | None = None) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    archive = output_dir / ARCHIVE_NAME
    manifest_path = output_dir / "release_manifest.json"
    official_sha, checksum_payload = fetch_official_checksum()
    if official_sha != EXPECTED_SHA256:
        raise RuntimeError(
            f"registered/official checksum disagreement: registered={EXPECTED_SHA256}, "
            f"official={official_sha}"
        )
    checksum_file = output_dir / "checksums.txt"
    checksum_file.write_bytes(checksum_payload)
    license_file = output_dir / "LICENSE"
    if not license_file.is_file():
        license_response = requests.get(LICENSE_URL, timeout=60)
        license_response.raise_for_status()
        license_file.write_bytes(license_response.content)

    if manifest_path.is_file():
        import json

        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        archive_record = verify_archive(archive, official_sha)
        sqlite_path = output_dir / existing["sqlite"]["path"]
        if (
            archive_record != existing["archive"]
            or sha256_file(sqlite_path) != existing["sqlite"]["sha256"]
        ):
            raise RuntimeError("immutable release files disagree with existing manifest")
        if (
            isinstance(existing.get("license"), dict)
            and existing["sqlite"].get("schema_sha256")
            and existing.get("acquisition_code_sha256")
        ):
            return existing

    last_report = 0.0

    def report_progress(written: int, total: int | None) -> None:
        nonlocal last_report
        now = time.monotonic()
        if now - last_report < 30:
            return
        last_report = now
        total_text = str(total) if total is not None else "unknown"
        print(f"download_bytes={written} total_bytes={total_text}", flush=True)

    if not archive.exists():
        download_resumable(ARCHIVE_URL, archive, expected_bytes, report_progress)
    archive_record = verify_archive(archive, official_sha)

    extracted = output_dir / "extracted"
    if not extracted.exists():
        files = safe_extract_archive(archive, extracted)
    else:
        files = sorted(path for path in extracted.rglob("*") if path.is_file())
    sqlite_files = [path for path in files if path.suffix in {".db", ".sqlite", ".sqlite3"}]
    if len(sqlite_files) != 1:
        raise RuntimeError(f"expected one extracted SQLite database, found {len(sqlite_files)}")
    sqlite_path = sqlite_files[0]
    manifest = {
        "schema": "MetaSieve.SourceRelease.v1",
        "stage": "P1R2B-D0-C",
        "source": "ChEMBL",
        "release": RELEASE,
        "format": "sqlite",
        "source_url": ARCHIVE_URL,
        "official_checksum_url": CHECKSUM_URL,
        "official_checksum_verified": True,
        "official_checksum_file_sha256": sha256_file(checksum_file),
        "license": {
            "name": LICENSE,
            "url": LICENSE_URL,
            "path": license_file.name,
            "sha256": sha256_file(license_file),
        },
        "acquisition_code_sha256": sha256_file(Path(__file__)),
        "immutable": True,
        "archive": archive_record,
        "sqlite": {
            "path": sqlite_path.relative_to(output_dir).as_posix(),
            "bytes": sqlite_path.stat().st_size,
            "sha256": sha256_file(sqlite_path),
            "schema_sha256": sqlite_schema_sha256(sqlite_path),
        },
        "acquired_at_utc": datetime.now(timezone.utc).isoformat(),
        "historical_f0r_training_source": False,
        "canonical_corpus_built": False,
        "training_authorized": False,
        "recipient_labels_read": False,
    }
    write_canonical_json(manifest_path, manifest)
    return manifest
