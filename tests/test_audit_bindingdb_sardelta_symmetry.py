import numpy as np

from research.crossed_interaction.audit_bindingdb_sardelta_symmetry import (
    antisymmetry_error,
    augment_reverse_pairs,
)


def test_augment_reverse_pairs_flips_delta_and_ligand_direction():
    pair = {
        "left_cell_id": "left",
        "right_cell_id": "right",
        "delta_pK": 2.5,
        "ligand_delta": np.asarray([1.0, -3.0], dtype=np.float64),
        "protein": np.asarray([0.2, 0.8], dtype=np.float64),
        "concat": np.asarray([0.2, 0.8, 1.0, -3.0], dtype=np.float64),
        "interaction": np.asarray([0.2, -0.6, 0.8, -2.4], dtype=np.float64),
    }

    forward, reverse = augment_reverse_pairs([pair])

    assert forward["delta_pK"] == 2.5
    assert reverse["left_cell_id"] == "right"
    assert reverse["right_cell_id"] == "left"
    assert reverse["delta_pK"] == -2.5
    assert np.allclose(reverse["ligand_delta"], [-1.0, 3.0])
    assert np.allclose(reverse["interaction"], [-0.2, 0.6, -0.8, 2.4])


def test_antisymmetry_error_pairs_forward_and_reverse_predictions():
    rows = [
        {
            "left_cell_id": "a",
            "right_cell_id": "b",
            "L_prediction": 1.25,
        },
        {
            "left_cell_id": "b",
            "right_cell_id": "a",
            "L_prediction": -1.20,
        },
    ]

    audit = antisymmetry_error(rows, "L")

    assert audit["pairs"] == 1
    assert np.isclose(audit["max_abs_sum"], 0.05)
