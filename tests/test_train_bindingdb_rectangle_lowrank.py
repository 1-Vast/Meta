import numpy as np

from research.crossed_interaction.train_bindingdb_rectangle_lowrank import (
    double_center_observed,
    lowrank_reconstruct,
    materialize_matrix,
    score,
    transformation_key,
)


def test_transformation_key_is_orientation_stable():
    assert transformation_key({"ligand_a": "b", "ligand_b": "a"}) == "a>b"


def test_double_center_observed_removes_additive_matrix():
    row = np.asarray([[1.0], [2.0], [3.0]])
    col = np.asarray([[4.0, 5.0, 6.0]])
    matrix = row + col

    residual, _ = double_center_observed(matrix)

    assert np.allclose(residual, 0.0)


def test_lowrank_reconstruct_scores_finite_on_shared_pattern():
    train = np.asarray([
        [1.0, 2.0, 3.0],
        [2.0, 4.0, 6.0],
        [3.0, 6.0, 9.0],
    ])
    dev = train.copy()

    prediction, metadata = lowrank_reconstruct(train, dev, rank=1)
    result = score(dev, prediction)

    assert metadata["rank"] == 1
    assert result["observed_cells"] == 9
    assert np.isfinite(result["prediction_mse"])


def test_materialize_matrix_can_filter_to_empty_shared_transforms():
    rows = [{
        "split": "train",
        "target_a": "ta",
        "target_b": "tb",
        "ligand_a": "la",
        "ligand_b": "lb",
        "rectangle": 1.0,
        "dependency_component": "c",
    }]

    matrix, metadata = materialize_matrix(rows, "train", transforms=["other"], min_targets=1)

    assert matrix.shape == (1, 1)
    assert metadata["observed_cells"] == 0
