import numpy as np

from research.crossed_interaction.audit_bindingdb_sardelta_symmetry import (
    augment_reverse_pairs,
)
from research.crossed_interaction.train_bindingdb_pair_score_difference import (
    score_difference_feature,
)


def test_score_difference_feature_contains_ligand_and_product_differences():
    pair = {
        "protein": np.asarray([2.0], dtype=np.float64),
        "left_ligand": np.asarray([3.0, 5.0], dtype=np.float64),
        "right_ligand": np.asarray([1.0, 7.0], dtype=np.float64),
    }

    feature = score_difference_feature(pair)

    assert np.allclose(feature, [2.0, -2.0, 4.0, -4.0])


def test_score_difference_flips_after_reverse_endpoint_swap():
    pair = {
        "left_cell_id": "left",
        "right_cell_id": "right",
        "delta_pK": 1.0,
        "protein": np.asarray([2.0], dtype=np.float64),
        "ligand_delta": np.asarray([2.0], dtype=np.float64),
        "concat": np.asarray([2.0, 2.0], dtype=np.float64),
        "interaction": np.asarray([4.0], dtype=np.float64),
        "left_ligand": np.asarray([3.0], dtype=np.float64),
        "right_ligand": np.asarray([1.0], dtype=np.float64),
    }
    forward, reverse = augment_reverse_pairs([pair])
    reverse["left_ligand"], reverse["right_ligand"] = (
        reverse["right_ligand"], reverse["left_ligand"])

    assert np.allclose(
        score_difference_feature(reverse),
        -score_difference_feature(forward))
