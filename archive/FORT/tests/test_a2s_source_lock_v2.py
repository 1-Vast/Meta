from __future__ import annotations

import pandas as pd

from research.a2s.a2s_source_lock_v2 import (
    assign_oof_folds,
    assign_roles,
    overlap_audit,
    quarantine_cross_role_rows,
)


def test_cross_role_provenance_rows_are_quarantined() -> None:
    frame = pd.DataFrame(
        [
            (0, "t1", "h1", "fit", "d_shared", "a1"),
            (1, "t1", "h1", "fit", "d_fit", "a2"),
            (2, "t2", "h2", "probe", "d_shared", "a3"),
            (3, "t2", "h2", "probe", "d_probe", "a4"),
            (4, "t3", "h3", "locked", "d_locked", "a5"),
        ],
        columns=["source_row", "target", "hcluster", "role", "docs", "assays"],
    )
    frame["component"] = frame.hcluster
    frame["conn"] = [f"c{row}" for row in range(len(frame))]
    frame["scaffold"] = [f"s{row}" for row in range(len(frame))]

    retained, crossing, reasons = quarantine_cross_role_rows(frame)

    assert crossing["docs"] == {"d_shared"}
    assert set(retained.source_row) == {1, 3, 4}
    assert reasons["docs"] == 2
    assert overlap_audit(retained)["all_mandatory_zero"]


def test_role_and_oof_assignments_are_deterministic_and_balanced() -> None:
    records = [
        {"component": f"h{row}", "rows": 10 + row, "targets": 1}
        for row in range(20)
    ]
    first = assign_roles(records, seed=7)
    assert first == assign_roles(records, seed=7)
    assert set(first.values()) == {"fit", "probe", "locked"}

    rows = []
    for component, role in first.items():
        if role != "fit":
            continue
        size = next(
            int(record["rows"])
            for record in records
            if record["component"] == component
        )
        rows.extend((component, role) for _ in range(size))
    retained = pd.DataFrame(rows, columns=["component", "role"])
    assignment, loads = assign_oof_folds(retained, n_folds=3, seed=7)
    assert set(assignment) == set(retained.component)
    assert max(loads) - min(loads) <= retained.groupby("component").size().max()
