import numpy as np

from research.meta_fewshot.a0_gauge_audit import (
    controls,
    procrustes_from_support,
    ridge_kernel,
)


def test_ridge_kernel_is_orthogonally_invariant():
    rng = np.random.default_rng(1)
    support, query = rng.normal(size=(5, 3)), rng.normal(size=(7, 3))
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    assert np.allclose(ridge_kernel(support, query, 0.5),
                       ridge_kernel(support @ rotation, query @ rotation, 0.5))


def test_support_procrustes_transfers_exact_rotation_to_query():
    rng = np.random.default_rng(2)
    correct_support, correct_query = rng.normal(size=(5, 3)), rng.normal(size=(7, 3))
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    wrong_support, wrong_query = correct_support @ rotation.T, correct_query @ rotation.T
    fitted = procrustes_from_support(wrong_support, correct_support)
    assert np.allclose(wrong_support @ fitted, correct_support)
    assert np.allclose(wrong_query @ fitted, correct_query)


def test_scale_and_shear_are_negative_controls():
    result = controls()
    assert result["orthogonal_H_max_abs"] < 1e-10
    assert result["scale_0_1_H_relative"] > 0.1
    assert result["scale_10_H_relative"] > 0.1
    assert result["shear_H_relative"] > 0.1
    assert result["rank_one_kernel_finite"]
