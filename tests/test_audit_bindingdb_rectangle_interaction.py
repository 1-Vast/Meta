import numpy as np

from research.crossed_interaction.audit_bindingdb_rectangle_interaction import (
    additive_cancellation_check,
    panel_rectangles,
    rectangle_value,
)


def test_rectangle_value_cancels_additive_main_effects():
    value = rectangle_value(
        y_ta_la=10.0 + 2.0,
        y_ta_lb=10.0 + 5.0,
        y_tb_la=20.0 + 2.0,
        y_tb_lb=20.0 + 5.0,
    )

    assert value == 0.0
    assert additive_cancellation_check()["pass"] is True


def test_rectangle_value_detects_interaction_term():
    value = rectangle_value(
        y_ta_la=0.0,
        y_ta_lb=1.0,
        y_tb_la=2.0,
        y_tb_lb=5.0,
    )

    assert value == -2.0


def test_panel_rectangles_materializes_complete_2x2():
    cells = [
        {"cell_id": "aa", "target_id": "ta", "ligand_id": "la", "pK": 0.0},
        {"cell_id": "ab", "target_id": "ta", "ligand_id": "lb", "pK": 1.0},
        {"cell_id": "ba", "target_id": "tb", "ligand_id": "la", "pK": 2.0},
        {"cell_id": "bb", "target_id": "tb", "ligand_id": "lb", "pK": 5.0},
    ]
    panel = {
        "panel_id": "p",
        "split": "development",
        "dependency_component": "c",
        "cell_ids": [cell["cell_id"] for cell in cells],
    }

    rows = panel_rectangles(panel, {cell["cell_id"]: cell for cell in cells})

    assert len(rows) == 1
    assert np.isclose(rows[0]["rectangle"], -2.0)
