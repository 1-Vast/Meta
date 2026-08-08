import numpy as np

from .identifiability_gate import _certificate_features


def test_certificate_features_are_support_order_invariant():
    rng = np.random.default_rng(9)
    correct = rng.normal(size=(13, 3))
    null = rng.normal(size=(13, 3))
    nearest = rng.normal(size=(13, 3))
    support = np.asarray([0, 2, 4, 7, 11])
    residual = rng.normal(size=5)
    reference = _certificate_features(
        correct, null, nearest, support, residual, 100.0, 0.1)
    permutation = np.asarray([3, 0, 4, 1, 2])
    permuted = _certificate_features(
        correct, null, nearest, support[permutation], residual[permutation],
        100.0, 0.1)
    np.testing.assert_allclose(reference, permuted, atol=1e-12)


def test_certificate_has_fixed_finite_dimension():
    rng = np.random.default_rng(10)
    surface = rng.normal(size=(8, 3))
    feature = _certificate_features(
        surface, np.zeros_like(surface), surface[::-1], np.arange(5),
        rng.normal(size=5), 10.0, 1.0)
    assert feature.shape == (14,)
    assert np.all(np.isfinite(feature))
