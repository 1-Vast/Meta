from __future__ import annotations

import pandas as pd

from research.a2s.a2s_natural_tail_audit import build_audit, release_number, tokens


def _document(document: str, release: int, source: int) -> tuple[str, dict[str, object]]:
    return document, {
        "document_chembl_id": document,
        "chembl_release": f"CHEMBL_{release}",
        "src_id": source,
    }


def _closed_target(target: str = "R") -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    rows: list[dict[str, object]] = []
    documents: dict[str, dict[str, object]] = {}
    for index in range(5):
        document, metadata = _document(f"EARLY-{index}", index + 1, 1)
        documents[document] = metadata
        rows.append(
            {
                "source_row": index,
                "target": target,
                "conn": f"EC-{index}",
                "endpoint": "pKi",
                "scaffold": f"ES-{index}",
                "assays": f"EA-{index}",
                "docs": document,
                "hcluster": f"H-{target}",
                "dual_cold_split": "train",
            }
        )
    for index in range(10):
        document, metadata = _document(f"LATE-{index}", 20 + index, 2)
        documents[document] = metadata
        rows.append(
            {
                "source_row": 5 + index,
                "target": target,
                "conn": f"LC-{index}",
                "endpoint": "pKi",
                "scaffold": f"LS-{index}",
                "assays": f"LA-{index}",
                "docs": document,
                "hcluster": f"H-{target}",
                "dual_cold_split": "train",
            }
        )
    return pd.DataFrame(rows), documents


def test_tokens_and_release_parser_are_strict() -> None:
    assert tokens("D2|D1|D2") == frozenset(("D1", "D2"))
    assert release_number("CHEMBL_37") == 37


def test_strict_natural_tail_target_has_nested_closed_episodes() -> None:
    frame, documents = _closed_target()
    report = build_audit(
        frame,
        documents,
        registry_sha256="registry",
        document_metadata_sha256="documents",
    )
    assert report["topology"]["cutflow"]["strict_admitted"] == 1
    record = report["roster"]["records"][0]
    assert record["min_closed_query_depth"] == 10
    for episode in record["episodes"]:
        support = episode["support_by_k"]
        assert support["1"] == support["3"][:1]
        assert support["3"] == support["5"][:3]
        assert set(support["5"]).isdisjoint(episode["query_rows"])


def test_same_source_family_cannot_enter_strict_roster() -> None:
    frame, documents = _closed_target()
    for metadata in documents.values():
        metadata["src_id"] = 1
    report = build_audit(
        frame,
        documents,
        registry_sha256="registry",
        document_metadata_sha256="documents",
    )
    assert report["topology"]["cutflow"]["temporal_envelope"] == 1
    assert report["topology"]["cutflow"]["source_family_envelope"] == 0
    assert report["roster"]["recipients"] == 0
    assert report["decision"]["status"] == "DATA_NOT_READY"
