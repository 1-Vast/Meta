from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pytest

from research.a2s import a2s_h0_metadata_inventory
from research.a2s.a2s_h0_metadata_inventory import (
    BINDINGDB_COMPANIONS,
    frozen_bindingdb_companion_metadata,
    historical_admission,
    immutable_write,
    read_bindingdb_header,
)


def test_direct_h0_module_entry_is_blocked(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("A2S_DATASET_RUN_DIR", raising=False)
    monkeypatch.delenv("A2S_DATASET_RUN_ID", raising=False)
    with pytest.raises(SystemExit, match="dataset-run"):
        a2s_h0_metadata_inventory.main()


def test_bindingdb_reader_stops_after_header(tmp_path: Path) -> None:
    archive = tmp_path / "bindingdb.zip"
    header = "BindingDB Reactant_set_id\tKi (nM)\tArticle DOI\r\n"
    row = "1\tSHOULD_NOT_BE_PARSED\t10.1/example\r\n"
    with ZipFile(archive, "w", ZIP_DEFLATED) as output:
        output.writestr("BindingDB.tsv", header + row)

    fields, metadata = read_bindingdb_header(archive)

    assert fields == ["BindingDB Reactant_set_id", "Ki (nM)", "Article DOI"]
    assert metadata["header_bytes"] == len(header.rstrip("\r\n"))
    assert "SHOULD_NOT_BE_PARSED" not in fields


def test_historical_admission_requires_activity_time_or_snapshots() -> None:
    status, blockers = historical_admission(
        dated_documents=True,
        activity_first_seen=False,
        snapshots_present=False,
    )
    assert status == "DATA_NOT_READY"
    assert len(blockers) == 1

    status, blockers = historical_admission(
        dated_documents=True,
        activity_first_seen=False,
        snapshots_present=True,
        presence_certified=False,
    )
    assert status == "DATA_NOT_READY"
    assert "presence audit" in blockers[0]

    status, blockers = historical_admission(
        dated_documents=True,
        activity_first_seen=False,
        snapshots_present=True,
        presence_certified=True,
    )
    assert status == "H0_S_PASS"
    assert blockers == []


def test_immutable_write_rejects_different_content(tmp_path: Path) -> None:
    path = tmp_path / "frozen.json"
    immutable_write(path, b"same")
    immutable_write(path, b"same")

    try:
        immutable_write(path, b"different")
    except FileExistsError:
        pass
    else:
        raise AssertionError("different frozen content was overwritten")


def test_bindingdb_companion_metadata_is_reused_without_network(
    tmp_path: Path,
) -> None:
    path = tmp_path / "h0-v2.json"
    responses = {
        key: {"available": True, "url": url, "content_length": index + 1}
        for index, (key, url) in enumerate(BINDINGDB_COMPANIONS.items())
    }
    path.write_text(
        json.dumps(
            {"sources": {"bindingdb": {"companion_metadata": responses}}}
        ),
        encoding="utf-8",
    )

    frozen = frozen_bindingdb_companion_metadata(path)

    assert frozen["network_requests_in_current_run"] == 0
    assert frozen["responses"] == responses
    assert frozen["evidence_sha256"] == sha256(path.read_bytes()).hexdigest()


def test_bindingdb_companion_metadata_rejects_wrong_url(tmp_path: Path) -> None:
    path = tmp_path / "h0-v2.json"
    responses = {
        key: {"available": True, "url": url}
        for key, url in BINDINGDB_COMPANIONS.items()
    }
    first = next(iter(responses))
    responses[first]["url"] = "https://example.invalid/not-the-frozen-source"
    path.write_text(
        json.dumps(
            {"sources": {"bindingdb": {"companion_metadata": responses}}}
        ),
        encoding="utf-8",
    )

    try:
        frozen_bindingdb_companion_metadata(path)
    except ValueError as error:
        assert "URL mismatch" in str(error)
    else:
        raise AssertionError("wrong companion URL was accepted")
