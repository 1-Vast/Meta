"""The balanced-pair construction must make symmetry a structural zero.

Stage L's first version sampled `i` and `j` independently, so the orientation
distribution was arbitrary and a symmetric score such as Tanimoto picked up an
incidental non-zero correlation with the signed gap. The report then compared
that number against a learned directional arm as though the two were competing
on the same task.

Under the balanced construction each unordered pair appears once in each
orientation, so the signed target is exactly antisymmetric and **any** symmetric
predictor has identically zero correlation with it — by construction, not by
measurement. These tests pin that.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.research.a2_exact_probe.stage_l2_directional_sar import (  # noqa: E402
    balanced_pairs, correlation, tanimoto_rows,
)


def blocks(sizes):
    out, start = [], 0
    for size in sizes:
        out.append(np.arange(start, start + size))
        start += size
    return out


def test_every_unordered_pair_appears_in_both_orientations():
    left, right, owner = balanced_pairs(blocks([5, 4]), seed=1)
    seen = {}
    for a, b in zip(left, right):
        seen[(int(a), int(b))] = seen.get((int(a), int(b)), 0) + 1
    for (a, b), count in seen.items():
        assert count == 1
        assert (b, a) in seen, f"({a},{b}) has no mirror"
    assert len(left) == len(owner)
    assert len(left) % 2 == 0


def test_the_signed_target_is_exactly_antisymmetric():
    left, right, _ = balanced_pairs(blocks([6, 5]), seed=2)
    y = np.random.default_rng(0).normal(size=11)
    delta = y[left] - y[right]
    assert float(delta.sum()) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_a_symmetric_score_has_identically_zero_signed_correlation(seed):
    """The property the whole construction exists to guarantee."""
    sizes = [7, 6, 5]
    left, right, _ = balanced_pairs(blocks(sizes), seed=seed)
    rng = np.random.default_rng(seed)
    total = sum(sizes)
    fingerprints = (rng.random((total, 64)) > 0.6).astype(np.float32)
    y = rng.normal(size=total)

    similarity = tanimoto_rows(fingerprints[left], fingerprints[right])
    delta = y[left] - y[right]
    # Symmetric predictor, antisymmetric target: the correlation is zero to
    # numerical precision regardless of the data.
    assert abs(correlation(similarity, delta)) < 1e-9
    assert abs(correlation(-similarity, delta)) < 1e-9


def test_an_antisymmetric_score_can_correlate():
    """The construction must not zero out *every* predictor."""
    sizes = [7, 6]
    left, right, _ = balanced_pairs(blocks(sizes), seed=4)
    rng = np.random.default_rng(4)
    total = sum(sizes)
    feature = rng.normal(size=total)
    y = 2.0 * feature + 0.1 * rng.normal(size=total)
    signed = feature[left] - feature[right]
    assert correlation(signed, y[left] - y[right]) > 0.9


def test_the_construction_zeroes_direction_but_not_magnitude():
    """The exact scope of what balancing does.

    Balancing forces a symmetric score's correlation with the **signed** gap to
    machine zero — that is its purpose. It must not also destroy the score's
    relationship with the **magnitude**, which is where Tanimoto is meaningful
    and which the incumbent's transport exploits.

    Asserted as a contrast rather than a threshold: on random fingerprints the
    magnitude correlation is small and seed-dependent, but it is many orders of
    magnitude above the signed correlation, which is identically zero.
    """
    signed, magnitude = [], []
    for seed in (5, 6, 7, 8):
        left, right, _ = balanced_pairs(blocks([20, 18]), seed=seed,
                                        per_target=999)
        rng = np.random.default_rng(seed)
        fingerprints = (rng.random((38, 64)) > 0.6).astype(np.float32)
        y = fingerprints @ rng.normal(size=64)
        similarity = tanimoto_rows(fingerprints[left], fingerprints[right])
        signed.append(abs(correlation(-similarity, y[left] - y[right])))
        magnitude.append(abs(correlation(-similarity, np.abs(y[left] - y[right]))))

    assert max(signed) < 1e-9, signed
    assert float(np.mean(magnitude)) > 1e-3, magnitude
    assert float(np.mean(magnitude)) > 1e5 * max(max(signed), 1e-18)


def test_pairs_never_cross_a_target_boundary():
    sizes = [4, 4, 4]
    left, right, owner = balanced_pairs(blocks(sizes), seed=6)
    for a, b, target in zip(left, right, owner):
        assert a // 4 == b // 4 == target


def test_construction_is_deterministic():
    first = balanced_pairs(blocks([9, 8]), seed=7)
    again = balanced_pairs(blocks([9, 8]), seed=7)
    assert all(np.array_equal(a, b) for a, b in zip(first, again))


def test_the_seed_matters_only_where_sampling_actually_happens():
    """A target with few enough ligands is enumerated exhaustively.

    24 ligands give 276 unordered pairs, above the 48 cap, so the seed selects
    a subset. 9 ligands give 36, below the cap, so every pair is taken and the
    seed is irrelevant — which is the correct behaviour, not a bug.
    """
    small_a = balanced_pairs(blocks([9]), seed=7, per_target=48)
    small_b = balanced_pairs(blocks([9]), seed=8, per_target=48)
    assert all(np.array_equal(a, b) for a, b in zip(small_a, small_b))

    large_a = balanced_pairs(blocks([24]), seed=7, per_target=48)
    large_b = balanced_pairs(blocks([24]), seed=8, per_target=48)
    assert not all(np.array_equal(a, b) for a, b in zip(large_a, large_b))
    assert len(large_a[0]) == 96          # 48 unordered pairs, both orientations


def test_small_targets_are_enumerated_exhaustively():
    """A 4-ligand target has 6 unordered pairs -> 12 oriented rows."""
    left, _, _ = balanced_pairs(blocks([4]), seed=9, per_target=999)
    assert len(left) == 12
