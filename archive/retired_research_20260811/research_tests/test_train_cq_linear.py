import numpy as np

from research.crossed_interaction.train_cq_linear import fit_ridge, panel_quotient


def test_panel_quotient_and_linear_witness_recover_interaction():
    rows = [
        {"target_id": "p1", "ligand_id": "a", "pK": 1.0},
        {"target_id": "p1", "ligand_id": "b", "pK": -1.0},
        {"target_id": "p2", "ligand_id": "a", "pK": -1.0},
        {"target_id": "p2", "ligand_id": "b", "pK": 1.0},
    ]
    raw = np.asarray([[1.0], [-1.0], [-1.0], [1.0]])
    response, features, rank = panel_quotient(rows, raw)
    panel = {"X": features, "y": response, "rank": rank}
    weight = fit_ridge([panel], 1e-6)
    assert rank == 1
    assert weight[0] > 0.99
