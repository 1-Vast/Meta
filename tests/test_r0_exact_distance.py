import numpy as np
import pytest

from research.correspondence_router.r0_exact_distance import (
    additive_checkerboard_rps,
    component_macro,
    deterministic_derangement,
    distance_bin_labels,
    negative_log_likelihood,
    paired_component_bootstrap,
    ranked_probability_score,
)


def test_distance_bins_follow_frozen_p1b_boundaries():
    distances = np.array([[0.0, 3.999, 4.0, 5.999, 6.0, 8.0, 10.0, 998.0]])
    assert distance_bin_labels(distances).tolist() == [[0, 0, 1, 1, 2, 3, 4, 4]]
    with pytest.raises(ValueError, match="overflow"):
        distance_bin_labels(np.array([[999.0]]))


def test_ranked_probability_score_is_zero_for_oracle_and_order_sensitive():
    labels = np.array([0, 2, 4])
    oracle = np.eye(5)[labels]
    assert np.array_equal(ranked_probability_score(oracle, labels), np.zeros(3))
    one_bin_away = np.eye(5)[np.array([1, 1, 3])]
    far_away = np.eye(5)[np.array([4, 4, 0])]
    assert ranked_probability_score(one_bin_away, labels).mean() \
        < ranked_probability_score(far_away, labels).mean()
    assert np.array_equal(negative_log_likelihood(oracle, labels), np.zeros(3))


def test_component_macro_does_not_pair_weight_large_systems():
    scores = {"large": 1.0, "small_a": 0.0, "small_b": 0.0}
    components = {"large": "x", "small_a": "y", "small_b": "y"}
    assert component_macro(scores, components) == {"x": 1.0, "y": 0.0}


def test_checkerboard_additive_oracle_recovers_additive_ordered_pattern():
    atom = np.array([-1.5, -0.5, 0.5, 1.5])[:, None]
    residue = np.linspace(-1.0, 1.0, 8)[None, :]
    labels = np.digitize(atom + residue, [-1.0, -0.2, 0.6, 1.4])
    slot = np.repeat(np.arange(4), 2)
    score = additive_checkerboard_rps(labels.astype(np.int64), slot, iterations=60)
    uniform = np.full((*labels.shape, 5), 0.2)
    assert score < ranked_probability_score(uniform, labels).mean()


def test_checkerboard_additive_oracle_rejects_degenerate_grid():
    with pytest.raises(ValueError, match=">=2 atoms"):
        additive_checkerboard_rps(
            np.zeros((1, 3), dtype=np.int64),
            np.arange(3, dtype=np.int64))


def test_derangement_has_no_fixed_points_inside_movable_groups():
    groups = [0, 0, 0, 1, 1, 2]
    mapping = deterministic_derangement(groups, namespace="res-test")
    assert sorted(mapping[:3].tolist()) == [0, 1, 2]
    assert sorted(mapping[3:5].tolist()) == [3, 4]
    assert np.all(mapping[:5] != np.arange(5))
    assert mapping[5] == 5
    assert np.array_equal(
        mapping, deterministic_derangement(groups, namespace="res-test"))


def test_paired_component_bootstrap_is_directional():
    worse = {f"c{i}": float(i + 1) for i in range(8)}
    better = {key: value - 0.5 for key, value in worse.items()}
    result = paired_component_bootstrap(worse, better, seed=7, draws=2000)
    assert result["delta"] == pytest.approx(0.5)
    assert result["lcb95_one_sided"] > 0
