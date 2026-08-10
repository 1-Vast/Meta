import sys
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1] / "research" / "crossed_interaction"
sys.path.insert(0, str(MODULE_DIR))

import materialize_x1ar_assays as selection


def _cell(mapping):
    return {"assay_activity_ids": mapping}


def test_common_assay_uses_maximum_minimum_replicates_then_lexical_tie():
    left = _cell({"B": [1, 2], "A": [3, 4], "C": [5]})
    right = _cell({"B": [6, 7], "A": [8, 9], "C": [10, 11, 12]})
    assay, left_ids, right_ids = selection.choose_common_assay(left, right)
    assert assay == "A"
    assert left_ids == [3, 4]
    assert right_ids == [8, 9]


def test_missing_common_assay_fails_rectangle_primary_eligibility():
    cells = [_cell({"A": [1]}), _cell({"B": [2]}),
             _cell({"C": [3]}), _cell({"C": [4]})]
    row = {"rectangle_id": "r", "endpoint": "Ki", "dependency_cluster": "g",
           "panel_id": "p", "selected_at_frozen_cap": True, "cells": cells}
    result = selection.select_rectangle(row)
    assert not result["eligible_primary"]
    assert result["exclusion_reason"] == "protein_a_has_no_common_exact_assay"
