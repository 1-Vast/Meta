"""Mechanical contracts for the heterogeneous-sparsity falsification gate."""

from __future__ import annotations

import numpy as np

from research.a2s import a2s_hotspot_falsification as hotspot


def test_standardisation_uses_only_the_fit_rows():
    values = np.asarray([[0.0], [2.0], [100.0]])
    result = hotspot.standardise(values, np.asarray([True, True, False]))
    assert np.allclose(result[:2].mean(axis=0), 0.0)
    assert np.allclose(result[:2].std(axis=0), 1.0)
    assert result[2, 0] > 50.0


def test_coordinate_truncation_keeps_the_largest_weights():
    head = np.asarray([1.0, -4.0, 3.0, 0.5])
    result = hotspot.truncate_coordinates(head, 2)
    assert np.allclose(result, [0.0, -4.0, 3.0, 0.0])


def test_orthogonal_rotation_preserves_full_linear_predictions():
    rng = np.random.default_rng(1)
    design = rng.normal(size=(30, 8))
    head = rng.normal(size=8)
    rotation, _ = np.linalg.qr(rng.normal(size=(8, 8)))
    assert np.allclose(design @ head, (design @ rotation) @ (rotation.T @ head))


def test_sparse_coordinates_lose_concentration_after_a_generic_rotation():
    rng = np.random.default_rng(2)
    head = np.zeros(30)
    head[[1, 7, 19]] = [3.0, -2.0, 1.0]
    rotation, _ = np.linalg.qr(rng.normal(size=(30, 30)))
    original_mass = np.abs(head)[hotspot.top_coordinates(head, 3)].sum() / np.abs(head).sum()
    rotated = rotation.T @ head
    rotated_mass = np.abs(rotated)[hotspot.top_coordinates(rotated, 3)].sum() / np.abs(rotated).sum()
    assert original_mass == 1.0
    assert rotated_mass < 0.5


def test_support_stability_is_one_for_identical_sets():
    import pandas as pd

    frame = pd.DataFrame(
        {
            "basis": ["x", "x", "x"],
            "target": ["t", "t", "t"],
            "size": [2, 2, 2],
            "coordinates": [(1, 2), (1, 2), (1, 2)],
        }
    )
    result = hotspot.support_stability(frame)
    assert result["mean"] == 1.0
    assert result["pairs"] == 3


def test_decision_requires_rotation_and_low_rank_separation():
    def cell(coord: float, rank: float, fraction: float, effective: int | None) -> dict:
        return {
            "sizes": {
                "s8": {
                    "coord": {"mean": coord, "lower95": coord - 0.01, "retained_fraction": fraction},
                    "rank": {"mean": rank, "lower95": rank - 0.01, "retained_fraction": rank},
                }
            },
            "effective_s60": effective,
            "stability": {"s3": {"mean": 0.6, "median": 0.6, "pairs": 10}},
        }

    summary = {
        "original26": cell(0.04, 0.01, 0.8, 3),
        "original_rot1": cell(0.02, 0.01, 0.3, 8),
    }
    assert hotspot.decide(summary)["verdict"] == "HETEROGENEOUS_SPARSITY_ADMITTED"
    summary["original26"]["sizes"]["s8"]["coord"]["mean"] = 0.005
    assert hotspot.decide(summary)["verdict"] == "HETEROGENEOUS_SPARSITY_NOT_REPRODUCED"
