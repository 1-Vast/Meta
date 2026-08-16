from __future__ import annotations

import pandas as pd

from research.a2s.a2s_cfes_coverage_gate import (
    ALPHAFOLD_API,
    PDBE_API,
    FetchResult,
    parse_alphafold,
    parse_pdbe,
    query_accession,
    summarise,
)


def test_public_metadata_parsers_are_explicit() -> None:
    alphafold = parse_alphafold(
        FetchResult(
            "ok",
            [
                {
                    "latestVersion": 6,
                    "modelCreatedDate": "2025-08-01T00:00:00Z",
                    "cifUrl": "https://example.test/model.cif",
                }
            ],
            None,
        )
    )
    assert alphafold["alphafold_status"] == "ok"
    assert alphafold["alphafold_version"] == 6

    pdbe = parse_pdbe(
        FetchResult(
            "ok",
            {"P12345": [{"pdb_id": "2XYZ"}, {"pdb_id": "1ABC"}, {"pdb_id": "2xyz"}]},
            None,
        ),
        "P12345",
    )
    assert pdbe["pdb_structures"] == 2
    assert pdbe["pdb_ids"] == "1abc|2xyz"


def test_unmapped_target_does_not_call_network() -> None:
    def fail_fetch(_: str) -> FetchResult:
        raise AssertionError("network must not be called for an unmapped target")

    row = query_accession("CHEMBL0", "component", None, fail_fetch)
    assert row["mapped"] is False
    assert row["alphafold_status"] == "unmapped"
    assert row["pdbe_status"] == "unmapped"


def test_query_accession_uses_only_public_structure_endpoints() -> None:
    seen: list[str] = []

    def fetch(url: str) -> FetchResult:
        seen.append(url)
        if "alphafold" in url:
            return FetchResult("ok", [{"latestVersion": 6}], None)
        return FetchResult("ok", {"P12345": [{"pdb_id": "1abc"}]}, None)

    row = query_accession("CHEMBL1", "component", "P12345", fetch)
    assert seen == [
        ALPHAFOLD_API.format(accession="P12345"),
        PDBE_API.format(accession="P12345"),
    ]
    assert row["pdb_structures"] == 1


def test_coverage_summary_applies_all_gates() -> None:
    rows = []
    for index in range(100):
        rows.append(
            {
                "target": f"t{index}",
                "component": f"c{index}",
                "accession": f"p{index}",
                "mapped": True,
                "alphafold_status": "ok" if index < 95 else "not_found",
                "pdbe_status": "ok" if index < 85 else "no_structures",
                "pdb_structures": 2 if index < 85 else 0,
            }
        )
    summary, checks = summarise(pd.DataFrame(rows))
    assert summary["fit_targets"] == 100
    assert summary["components_with_two_pdb_states"] == 85
    assert all(checks.values())

    rows[0]["mapped"] = False
    for index in range(59, 100):
        rows[index]["pdb_structures"] = 0
    _, failed = summarise(pd.DataFrame(rows))
    assert failed["two_pdb_state_fraction"] is False
    assert failed["components_with_two_pdb_states"] is False
