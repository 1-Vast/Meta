from __future__ import annotations

from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from research.a2s import a2s_historical_presence
from research.a2s.a2s_historical_presence import (
    RELEASES,
    build_presence_audit,
    certify_presence_replay,
)


def test_direct_presence_module_entry_is_blocked(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("A2S_DATASET_RUN_DIR", raising=False)
    monkeypatch.delenv("A2S_DATASET_RUN_ID", raising=False)
    with pytest.raises(SystemExit, match="dataset-run"):
        a2s_historical_presence.main()


def _activity(
    activity_id: int,
    fingerprint: str,
    *,
    document: str = "DOC1",
    year: int = 2017,
) -> dict[str, object]:
    return {
        "activity_id": activity_id,
        "native_activity_key": f"chembl:activity:{activity_id}",
        "measurement_identity_sha256": fingerprint,
        "activity_src_id": 1,
        "activity_source": "LITERATURE",
        "document_source_id": 1,
        "document_chembl_id": document,
        "normalized_doi": f"10.1/{document.lower()}",
        "pubmed_id": activity_id + 1000,
        "patent_id": None,
        "document_index": document,
        "assay_chembl_id": f"ASSAY-{fingerprint}",
        "assay_src_id": 1,
        "assay_index": f"A-{fingerprint}",
        "assay_source_id": None,
        "compound_index": f"C-{fingerprint}",
        "record_source_id": 1,
        "source_compound_id": None,
        "parent_molecule_chembl_id": f"MOL-{fingerprint}",
        "parent_inchikey": "ABCDEFGHIJKLMN-AAAAAA-B",
        "target_chembl_id": "TARGET1",
        "target_uniprot": "P12345",
        "target_component_count": 1,
        "endpoint_type": "Ki",
        "standard_relation": "=",
        "standard_units": "nM",
        "publication_year": year,
        "document_id_consistent": True,
        "molecule_id_consistent": True,
    }


def _write_inputs(
    root: Path, rows: dict[str, list[dict[str, object]]]
) -> tuple[dict[str, Path], dict[str, Path]]:
    activities: dict[str, Path] = {}
    documents: dict[str, Path] = {}
    for release in RELEASES:
        release_root = root / release
        release_root.mkdir(parents=True)
        activity_path = release_root / "activities.parquet"
        document_path = release_root / "documents.parquet"
        if rows[release]:
            activity_table = pa.Table.from_pylist(rows[release])
        else:
            activity_table = pa.Table.from_pylist([_activity(0, "template")]).slice(
                0, 0
            )
        pq.write_table(activity_table, activity_path)
        document_ids = sorted(
            {"DOC1", "DOLD"}.union(
                str(row["document_chembl_id"]) for row in rows[release]
            )
        )
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {
                        "document_chembl_id": document,
                        "normalized_doi": f"10.1/{document.lower()}",
                        "pubmed_id": index + 1000,
                        "patent_id": None,
                    }
                    for index, document in enumerate(document_ids)
                ]
            ),
            document_path,
        )
        activities[release] = activity_path
        documents[release] = document_path
    return activities, documents


def test_presence_masks_rekey_and_backfill_are_explicit(tmp_path: Path) -> None:
    rows = {release: [] for release in RELEASES}
    for release in RELEASES:
        rows[release].append(_activity(1, "f-stable"))
    for release in RELEASES[1:]:
        rows[release].append(
            _activity(2, "f-backfill", document="DOLD", year=2010)
        )
    for release in (RELEASES[0], RELEASES[2], RELEASES[3]):
        rows[release].append(_activity(3, "f-deleted"))
    rows[RELEASES[0]].append(_activity(4, "f-rekey"))
    for release in RELEASES[1:]:
        rows[release].append(_activity(44, "f-rekey"))

    activities, documents = _write_inputs(tmp_path / "inputs", rows)
    output = tmp_path / "presence.parquet"
    report = build_presence_audit(activities, documents, output)
    records = {
        row["measurement_identity_sha256"]: row
        for row in pq.read_table(output).to_pylist()
    }

    assert report["status"] == "H0_IDENTITY_PASS"
    assert records["f-stable"]["presence_mask"] == "1111"
    assert records["f-stable"]["match_class"] == "native_stable"
    assert records["f-backfill"]["presence_mask"] == "0111"
    assert (
        records["f-backfill"]["first_seen_classification"]
        == "activity_backfill_under_preexisting_document"
    )
    assert records["f-deleted"]["presence_mask"] == "1011"
    assert records["f-deleted"]["stable_identity_eligible"] is False
    assert records["f-rekey"]["match_class"] == "rekeyed_1to1"
    assert records["f-rekey"]["presence_mask"] == "1111"


def test_native_id_lineage_collision_stops_identity_admission(tmp_path: Path) -> None:
    rows = {release: [] for release in RELEASES}
    rows[RELEASES[0]].append(_activity(9, "f-old", document="DOC1"))
    rows[RELEASES[1]].append(_activity(9, "f-new", document="DOC2"))
    activities, documents = _write_inputs(tmp_path / "inputs", rows)

    report = build_presence_audit(
        activities, documents, tmp_path / "presence.parquet"
    )
    records = pq.read_table(tmp_path / "presence.parquet").to_pylist()

    assert report["status"] == "DATA_NOT_READY"
    assert report["admission"]["hard_native_collisions"] == 1
    assert records[0]["native_collision"] is True
    assert records[0]["stable_identity_eligible"] is False


def test_presence_audit_replays_with_identical_content_hash(tmp_path: Path) -> None:
    rows = {release: [_activity(1, "f-stable")] for release in RELEASES}
    activities, documents = _write_inputs(tmp_path / "inputs", rows)

    first = build_presence_audit(
        activities, documents, tmp_path / "presence-first.parquet"
    )
    second = build_presence_audit(
        activities, documents, tmp_path / "presence-second.parquet"
    )

    assert first["output"]["canonical_rows_sha256"] == second["output"][
        "canonical_rows_sha256"
    ]
    assert first["output"]["parquet_sha256"] == second["output"][
        "parquet_sha256"
    ]
    assert not list(tmp_path.glob(".presence-work-*.sqlite"))

    certificate = certify_presence_replay(
        activities, documents, tmp_path / "presence-first.parquet"
    )
    assert certificate["status"] == "PASS"
    assert all(certificate["comparisons"].values())
