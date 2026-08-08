import numpy as np

from model import bands
from model.config import MetaSieveConfig
from model.meta_operator import build_anchors
from research.e0_identifiability.l0_contract import (
    Z_BIO_NAMES,
    assert_bounded_observable,
    band_mean_interval,
    bounded_projection,
    pooled_replicate_sigma,
    z_bio,
)


def test_frozen_anchor_ladder_is_stochastically_ordered():
    cfg = MetaSieveConfig()
    anchors = build_anchors(cfg, device="cpu").numpy()
    ordered = anchors[:6]
    for index in range(len(ordered) - 1):
        low_lower, low_upper = bands.split(ordered[index][None, :])
        high_lower, high_upper = bands.split(ordered[index + 1][None, :])
        assert np.all(high_lower <= low_lower + 1e-12)
        assert np.all(high_upper <= low_upper + 1e-12)


def test_ladder_mean_interval_increases_and_stays_valid():
    cfg = MetaSieveConfig()
    anchors = build_anchors(cfg, device="cpu").numpy()
    intervals = band_mean_interval(anchors[:6], cfg.grid(), cfg.a_max)
    assert np.all(np.diff(intervals[:, 0]) > 0)
    assert np.all(np.diff(intervals[:, 1]) > 0)
    for anchor in anchors:
        assert bands.is_valid(anchor)


def test_z_bio_is_bounded_and_permutation_invariant():
    rng = np.random.default_rng(3)
    phi = rng.gamma(2.0, size=(5, 8, 6, 6))
    heavy = rng.integers(10, 60, size=5).astype(float)
    value = z_bio(phi, heavy)
    assert value.shape == (5, len(Z_BIO_NAMES))
    assert_bounded_observable(value, "z_bio")
    # summing over atoms and residues means a reordering of rows changes nothing
    order = rng.permutation(5)
    assert np.allclose(z_bio(phi[order], heavy[order]), value[order])


def test_bounded_projection_stays_in_unit_cube():
    rng = np.random.default_rng(4)
    train = rng.normal(size=(200, 32))
    value = bounded_projection(train, rng.normal(size=(50, 32)), 7)
    assert value.shape == (50, 7)
    assert value.min() >= 0.0 and value.max() <= 1.0


def test_replicate_sigma_ignores_singleton_cells():
    groups = [np.asarray([1.0]), np.asarray([2.0, 2.0]), np.asarray([0.0, 2.0])]
    sigma, cells, degrees = pooled_replicate_sigma(groups)
    assert cells == 2 and degrees == 2
    assert np.isclose(sigma, 1.0)


def test_narrow_band_cannot_contain_a_step_cdf():
    """The defect that failed L0 closed: step containment is infeasible when narrow."""
    cfg = MetaSieveConfig()
    grid = cfg.grid()
    anchors = build_anchors(cfg, device="cpu").numpy()
    lower, upper = bands.split(anchors[2][None, :])
    step = (grid >= 0.5).astype(float)[None, :]
    assert float(np.mean(upper - lower)) < 1.0 - 1.0 / cfg.n_grid
    assert not (np.all(lower <= step) and np.all(upper >= step))
