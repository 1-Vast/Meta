from __future__ import annotations

import pandas as pd
import pytest

from scripts.htl_topology import build_audit, token_key, tokens


def _frame() -> pd.DataFrame:
    rows = []
    for index in range(6):
        rows.append(
            {
                "target": "T1",
                "conn": f"C{index}",
                "endpoint": "pKi",
                "n_records": 1,
                "scaffold": f"S{index}",
                "assays": "A2|A1" if index == 0 else "A1",
                "docs": "D2|D1" if index == 0 else "D1",
                "accession": "P1",
                "hcluster": "H1",
                "dual_cold_split": "train",
            }
        )
    rows.extend(
        {
            "target": "T2",
            "conn": "C7",
            "endpoint": "pKi",
            "n_records": 2,
            "scaffold": "S7",
            "assays": "A3",
            "docs": "D3",
            "accession": "P2",
            "hcluster": "H2",
            "dual_cold_split": "train",
        }
        for _ in range(2)
    )
    return pd.DataFrame(rows)


def test_tokens_are_canonical_and_n_eff_is_provenance_deduplicated() -> None:
    assert tokens("D2|D1|D2") == ("D1", "D2")
    assert token_key("A2|A1") == "A1|A2"
    report = build_audit(
        _frame(),
        split="train",
        registry_sha256="registry",
        manifest_sha256="manifest",
    )
    pki = report["global"]["by_endpoint"]["pKi"]
    targets = {row["target"]: row for row in pki["targets"]}
    assert pki["target_count"] == 2
    assert targets["T1"]["n_eff"] == 6
    assert targets["T1"]["unique_documents"] == 2
    assert targets["T1"]["scaffold_closed_query_depth_upper_bound"]["5"] == 1
    assert targets["T2"]["n_eff"] == 1


def test_audit_rejects_mixed_endpoint_or_split() -> None:
    frame = _frame()
    frame.loc[0, "endpoint"] = "IC50"
    with pytest.raises(ValueError, match="unexpected endpoint"):
        build_audit(
            frame,
            split="train",
            registry_sha256="registry",
            manifest_sha256="manifest",
        )

    with pytest.raises(ValueError, match="only accepts train or development"):
        build_audit(
            _frame(),
            split="confirmation",
            registry_sha256="registry",
            manifest_sha256="manifest",
        )
