import numpy as np

from research.e0_identifiability.eaff_pilot_contract import (
    assert_paffinity_direction,
    component_macro_contrasts,
    coupling_null,
    fit_pair_ridge,
    task_pair_differences,
)
from research.e0_identifiability.eaff_h0c_contract import centered_interaction
from research.e0_identifiability.x0_contract import (
    crossed_panel_payload,
    rectangle_count,
)


def test_paffinity_direction_is_stronger_larger():
    molar = np.asarray([1e-9, 1e-8, 1e-7])
    assert_paffinity_direction(-np.log10(molar), molar)


def test_coupling_null_preserves_both_marginals():
    rng = np.random.default_rng(4)
    value = rng.uniform(0.1, 2.0, size=(3, 8, 6, 6))
    null = coupling_null(value)
    np.testing.assert_allclose(null.sum(-1), value.sum(-1))
    np.testing.assert_allclose(null.sum((-3, -2)), value.sum((-3, -2)))


def test_centered_interaction_has_zero_marginals_with_negative_coordinates():
    rng = np.random.default_rng(31)
    value = rng.normal(0.01, 0.02, size=(3, 8, 6, 6)) + 0.03
    centered = centered_interaction(value)
    np.testing.assert_allclose(centered.sum(-1), 0.0, atol=1e-12)
    np.testing.assert_allclose(centered.sum((-3, -2)), 0.0, atol=1e-12)


def test_centered_interaction_rejects_nonpositive_mass():
    import pytest

    with pytest.raises(ValueError, match="positive total mass"):
        centered_interaction(np.zeros((1, 8, 6, 6)))


def test_crossed_panel_payload_excludes_target_identity():
    base = {
        "document_chembl_id": "DOC", "endpoint_family": "Ki",
        "assay_organism": "Homo sapiens", "bao_format": "BAO",
        "cell_id": None, "tissue_id": None, "subcellular_fraction": None,
        "relationship_type": "D", "variant_id": None,
        "target_accession": "P1", "assay_chembl_id": "A1",
    }
    changed = dict(base, target_accession="P2", assay_chembl_id="A2")
    assert crossed_panel_payload(base, []) == crossed_panel_payload(changed, [])


def test_rectangle_count_uses_shared_ligand_pairs():
    pairs, rectangles = rectangle_count({
        "p1": {"a", "b", "c"}, "p2": {"a", "b", "c"}, "p3": {"a"},
    })
    assert pairs == 1
    assert rectangles == 3


def test_task_pair_weights_sum_to_one_per_task():
    x = np.arange(24, dtype=float).reshape(6, 4)
    y = np.arange(6, dtype=float)
    tasks = np.asarray(["a"] * 3 + ["b"] * 3)
    _, _, weights = task_pair_differences(x, y, tasks)
    np.testing.assert_allclose(weights[:3].sum(), 1.0)
    np.testing.assert_allclose(weights[3:].sum(), 1.0)


def test_pair_ridge_recovers_a_shared_direction():
    rng = np.random.default_rng(7)
    x = rng.normal(size=(60, 5))
    tasks = np.repeat(np.arange(3), 20)
    truth = np.asarray([0.4, -0.8, 0.2, 0.0, 0.6])
    y = x @ truth + np.repeat([3.0, -2.0, 1.0], 20)
    direction, _ = fit_pair_ridge(x, y, tasks, alpha=1e-8)
    np.testing.assert_allclose(direction, truth, atol=1e-6)


def test_component_macro_does_not_weight_task_rich_component_more():
    rows = [
        {"closure_component_id": "a", "ligand": 0.5, "correct": 0.9,
         "deranged": 0.5, "null": 0.5},
        {"closure_component_id": "a", "ligand": 0.5, "correct": 0.9,
         "deranged": 0.5, "null": 0.5},
        {"closure_component_id": "b", "ligand": 0.5, "correct": 0.5,
         "deranged": 0.5, "null": 0.5},
    ]
    result = component_macro_contrasts(rows)
    np.testing.assert_allclose(result["correct_minus_ligand"], 0.2)
