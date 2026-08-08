import numpy as np

from .law_bridge import DEFAULT_MESH, band_map, operator_law, support_confidence


def test_operator_law_preserves_mass_and_requested_mean():
    for biological, tau, confidence in [(-0.2, 0.1, 0.0), (0.15, 0.05, 0.7), (0.4, 0.4, 1.0)]:
        law = operator_law(biological, tau, confidence)
        expected = np.clip(0.5 + biological + tau, 0.0, 1.0)
        np.testing.assert_allclose(law.f.sum(), 1.0, atol=1e-12)
        np.testing.assert_allclose(law.beta.sum(), 1.0, atol=1e-12)
        np.testing.assert_allclose(law.mean, expected, atol=1e-12)
        assert np.all(law.f >= 0.0) and np.all(law.beta >= 0.0)


def test_band_is_stochastic_banded_and_moment_preserving():
    band = band_map(0.2)
    np.testing.assert_allclose(band.sum(axis=0), 1.0, atol=1e-12)
    np.testing.assert_allclose(DEFAULT_MESH @ band, DEFAULT_MESH, atol=1e-12)
    row, column = np.nonzero(band)
    assert np.all(np.abs(row - column) <= 1)


def test_lower_confidence_increases_dispersion():
    low = operator_law(0.0, 0.0, 0.0)
    high = operator_law(0.0, 0.0, 1.0)
    assert low.variance > high.variance


def test_support_confidence_is_permutation_invariant_and_bounded():
    residual = np.asarray([0.2, -0.1, 0.0, 0.3, 0.1])
    reference = support_confidence(residual)
    assert reference == support_confidence(residual[[3, 0, 4, 1, 2]])
    assert 0.0 <= reference <= 1.0
