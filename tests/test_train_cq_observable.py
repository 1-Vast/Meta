import numpy as np

from research.crossed_interaction.train_cq_observable import (
    QuotientBlock,
    additive_residual,
    bootstrap_contrast,
    fit_ridge,
    predict,
    score_blocks,
)


def test_additive_residual_removes_target_and_ligand_main_effects():
    targets, ligands, values = [], [], []
    for ti, target in enumerate(("t1", "t2", "t3")):
        for li, ligand in enumerate(("l1", "l2")):
            targets.append(target)
            ligands.append(ligand)
            values.append(10.0 + ti - 2.0 * li)
    residual, rank, orthogonality = additive_residual(
        targets, ligands, np.asarray(values, dtype=np.float64))
    assert rank == 2
    assert np.linalg.norm(residual) < 1e-10
    assert orthogonality < 1e-10


def test_ridge_observable_learns_quotient_interaction_signal():
    targets = ["t1", "t1", "t2", "t2"]
    ligands = ["l1", "l2", "l1", "l2"]
    y_raw = np.asarray([1.0, -1.0, -1.0, 1.0], dtype=np.float64)
    y, retained_rank, _ = additive_residual(targets, ligands, y_raw)
    feature = y[:, None]
    block = QuotientBlock(
        panel_id="p1",
        split="train",
        dependency_component="c1",
        retained_rank=retained_rank,
        y=y,
        features={
            "correct": feature,
            "deranged_protein": np.zeros_like(feature),
            "foreign_ligand": -feature,
        },
    )
    model = fit_ridge([block], "correct", ridge=1e-6)
    prediction = predict(model, block)
    assert np.corrcoef(y, prediction)[0, 1] > 0.99


def test_bootstrap_contrast_passes_only_when_correct_beats_control():
    rows = []
    for component, correct_mse, control_mse in (
            ("c1", 1.0, 2.0), ("c2", 1.5, 3.0), ("c3", 1.2, 2.4)):
        rows.append({
            "dependency_component": component,
            "arm": "correct",
            "rank_normalized_mse": correct_mse,
        })
        rows.append({
            "dependency_component": component,
            "arm": "control",
            "rank_normalized_mse": control_mse,
        })
    result = bootstrap_contrast(
        rows, "correct", "control", draws=999, seed=1)
    assert result["pass"]
    assert result["one_sided_95_lcb"] > 0


def test_zero_additive_control_is_scored_in_quotient_space():
    block = QuotientBlock(
        panel_id="p",
        split="development",
        dependency_component="c",
        retained_rank=1,
        y=np.asarray([0.5, -0.5]),
        features={
            "correct": np.zeros((2, 1)),
            "deranged_protein": np.zeros((2, 1)),
            "foreign_ligand": np.zeros((2, 1)),
        },
    )
    rows, summary = score_blocks([block], {})
    assert {row["arm"] for row in rows} == {"zero_additive"}
    assert summary["zero_additive"]["rank_weighted_mse"] == 0.5
