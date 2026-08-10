import math

from research.crossed_interaction.bindingdb_cq_r1 import (
    panel_residual,
    parse_exact_nm,
)


def test_parse_exact_nm_rejects_censoring():
    assert parse_exact_nm("10") == 10.0
    assert parse_exact_nm("1.5e2") == 150.0
    assert parse_exact_nm("<10") is None
    assert parse_exact_nm("> 10") is None
    assert parse_exact_nm("0") is None


def test_additive_panel_has_zero_quotient():
    rows = []
    for ti, target in enumerate(("t1", "t2", "t3")):
        for li, ligand in enumerate(("l1", "l2", "l3")):
            rows.append({"target_id": target, "ligand_id": ligand, "pK": 5 + ti + 2 * li})
    result = panel_residual(rows)
    assert result is not None
    assert result["retained_rank"] == 4
    assert result["rank_normalized_mse"] < 1e-20
    assert result["orthogonality_error"] < 1e-10


def test_crossed_interaction_survives_quotient():
    rows = []
    values = [[0.0, 1.0], [1.0, 0.0]]
    for ti, target in enumerate(("t1", "t2")):
        for li, ligand in enumerate(("l1", "l2")):
            rows.append({"target_id": target, "ligand_id": ligand, "pK": values[ti][li]})
    result = panel_residual(rows)
    assert result is not None
    assert result["retained_rank"] == 1
    assert math.isclose(result["rank_normalized_mse"], 1.0, rel_tol=1e-10)
