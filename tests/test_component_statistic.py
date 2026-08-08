import json

import numpy as np
import pytest

from model.component_statistic import component_prediction, support_location
from scripts.evaluate_component_statistic import transform


def test_support_location_is_bounded_and_permutation_invariant():
    residual = np.asarray([0.8, -0.1, 0.4, 0.2, 0.7])
    reference = support_location(residual, ridge=1.0)
    permuted = support_location(residual[[3, 0, 4, 1, 2]], ridge=1.0)
    assert reference == permuted
    assert -0.5 <= reference <= 0.5


def test_component_prediction_keeps_support_location_protein_independent():
    residual = np.asarray([0.2, -0.1, 0.3])
    first, first_location = component_prediction(
        np.asarray([0.1, 0.4]), residual, ridge=2.0,
    )
    second, second_location = component_prediction(
        np.asarray([-0.2, 0.8]), residual, ridge=2.0,
    )
    assert first_location == second_location
    np.testing.assert_allclose(first - first_location, [0.1, 0.4])
    np.testing.assert_allclose(second - second_location, [-0.2, 0.8])


def test_jsonl_transform_rejects_query_labels(tmp_path):
    source = tmp_path / "input.jsonl"
    output = tmp_path / "output.jsonl"
    source.write_text(json.dumps({
        "id": "task-1",
        "biological_surface": [0.1, 0.2],
        "support_residual": [0.2, -0.1],
        "query_labels": [0.4, 0.5],
    }) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="query labels are forbidden"):
        transform(source, output, ridge=1.0, bound=0.5)
