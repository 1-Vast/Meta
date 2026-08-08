import numpy as np

from research.e0_identifiability.metrics import concordance
from research.e0_identifiability.run_eaff_x0_feasibility import (
    chi_square_required_n,
    greedy_disjoint_rectangles,
)
from research.e0_identifiability.run_eaff_x0b import (
    breakeven_rho,
    capped_design,
    effective_n,
    pack_cell_disjoint,
    verify_cell_disjoint,
)


def test_within_task_concordance_is_invariant_to_location_and_scale():
    rng = np.random.default_rng(11)
    labels = rng.normal(size=20)
    predictions = rng.normal(size=20)
    base = concordance(labels, predictions)
    assert concordance(labels, predictions + 7.5) == base
    assert concordance(labels, predictions * 3.25) == base
    assert concordance(labels + 7.5, predictions) == base
    assert concordance(labels * 3.25, predictions) == base


def test_constant_within_task_prediction_scores_exactly_chance():
    rng = np.random.default_rng(12)
    labels = rng.normal(size=20)
    assert concordance(labels, np.full(20, 4.2)) == 0.5


def test_frozen_required_n_is_245():
    assert chi_square_required_n(1.25, 0.05, 0.80) == 245


def test_cell_disjoint_packing_is_valid_and_beats_the_conservative_unit():
    targets = {f"t{index}": {f"l{value}" for value in range(8)} for index in range(4)}
    packed = pack_cell_disjoint(targets)
    assert verify_cell_disjoint(packed)
    for left, right, first, second in packed:
        assert left != right and first != second
        for target, ligand in ((left, first), (left, second), (right, first), (right, second)):
            assert ligand in targets[target]
    assert len(packed) > greedy_disjoint_rectangles(targets)


def test_design_effect_reduces_to_cluster_count_at_total_correlation():
    sizes = [500, 200, 100, 40, 10]
    best = max(effective_n(*capped_design(sizes, cap), 1.0)
               for cap in (1, 2, 5, 20, 100, None))
    assert best == float(len(sizes))


def test_effective_n_is_bounded_by_clusters_over_rho():
    sizes = [4000, 1000, 300, 120]
    for rho in (0.05, 0.1, 0.25, 0.5):
        total, influence = capped_design(sizes, None)
        assert effective_n(total, influence, rho) <= len(sizes) / rho + 1e-9


def test_breakeven_rho_is_zero_when_units_are_below_the_requirement():
    assert breakeven_rho(100, 12.0) == 0.0
