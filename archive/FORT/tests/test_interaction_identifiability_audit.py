from __future__ import annotations

import pandas as pd

from research.shared.interaction_identifiability_audit import aggregate_cells, build_rectangles, summarize


def _rows() -> pd.DataFrame:
    rows = []
    values = {
        ("t1", "A"): (7.0, 7.1),
        ("t1", "B"): (5.0, 5.1),
        ("t2", "A"): (4.0, 4.1),
        ("t2", "B"): (6.0, 6.1),
    }
    for target_ligand, values_pair in values.items():
        target, ligand = target_ligand
        for value in values_pair:
            rows.append(("pKi", "d1", f"a-{target}", target, ligand, value))
    return pd.DataFrame(rows, columns=["endpoint", "doc", "assay", "target", "ligand", "y"])


def test_exact_rectangle_and_reversal():
    cells, replicates = aggregate_cells(_rows())
    rectangles = build_rectangles(cells, {"t1": "h1", "t2": "h2"})
    assert len(rectangles) == 1
    assert rectangles.iloc[0].dd == 4.0
    assert rectangles.iloc[0].order_reversal == 1
    assert len(replicates) == 4


def test_summary_uses_endpoint_specific_and_component_units():
    cells, replicates = aggregate_cells(_rows())
    rectangles = build_rectangles(cells, {"t1": "h1", "t2": "h2"})
    result = summarize(rectangles, replicates, 1729)
    assert result["unit_count"] == 1
    assert result["by_endpoint"]["pKi"]["unit_count"] == 1
    assert result["order_reversal_fraction"] == 1.0


def test_same_homology_pairs_are_excluded_by_strict_default():
    cells, _ = aggregate_cells(_rows())
    rectangles = build_rectangles(cells, {"t1": "h1", "t2": "h1"})
    assert rectangles.empty
    sensitivity = build_rectangles(
        cells, {"t1": "h1", "t2": "h1"}, cross_component_only=False
    )
    assert len(sensitivity) == 1
