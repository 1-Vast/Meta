"""Contract tests for the S7/L2B Phase 2B residue-residual head.

These pin the numerical and structural guarantees the Phase 2B preregistration
(PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md, sha 5e6688f6...) relies on.
They use tiny synthetic inputs only: no dataset build, no ESM cache, no
checkpoint, no affinity source.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))

p2b = pytest.importorskip("p2b_residue_residual")


# ------------------------------------------------------------- frozen contract
def test_frozen_constants_match_preregistration():
    assert p2b.K == 8
    assert p2b.D_ESM == 1280 and p2b.D_ATOM == 41
    assert p2b.N_PARAMS_EXPECTED == 10568
    assert (p2b.EPOCHS, p2b.LR, p2b.WD, p2b.CLIP) == (6, 1e-3, 1e-4, 5.0)
    assert (p2b.SEED_PARAM, p2b.SEED_SAMPLER, p2b.SEED_BOOT,
            p2b.SEED_CTRL, p2b.SEED_SYNTH) == (20260901, 20260902, 20260903,
                                               20260904, 20260905)
    assert p2b.ORTHO_TOL == 1e-8 and p2b.SYNTH_MIN_AP == 0.50


def test_head_has_exactly_the_registered_parameters_and_no_bias():
    h = p2b.Head()
    assert p2b.n_params(h) == p2b.N_PARAMS_EXPECTED
    names = {n for n, _ in h.named_parameters()}
    assert names == {"U", "V"}
    assert tuple(h.U.shape) == (8, 1280) and tuple(h.V.shape) == (8, 41)


def test_head_initialisation_is_seed_deterministic():
    a, b = p2b.Head(), p2b.Head()
    assert np.array_equal(a.U.detach().numpy(), b.U.detach().numpy())
    assert np.array_equal(a.V.detach().numpy(), b.V.detach().numpy())


# ------------------------------------------------------------- pooling
def test_g_is_invariant_to_atom_permutation():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(23, p2b.D_ATOM))
    assert np.abs(p2b.g_of(A) - p2b.g_of(A[rng.permutation(23)])).max() < 1e-12


# ------------------------------------------------------------- gauge
def test_nuisance_basis_rank_and_degenerate_fallback():
    L = 40
    assert p2b.nuisance_basis(np.linspace(0.0, 1.0, L)).shape[1] == 2
    assert p2b.nuisance_basis(np.full(L, 7.0)).shape[1] == 1      # b is constant
    assert p2b.nuisance_basis(np.zeros(L)).shape[1] == 1          # b is zero


def test_projection_is_orthogonal_to_the_nuisance_basis():
    rng = np.random.default_rng(1)
    for _ in range(5):
        L = int(rng.integers(30, 200))
        Q = p2b.nuisance_basis(rng.normal(size=L))
        d = p2b.project_np(Q, rng.normal(size=L))
        assert p2b.ortho_ratio(Q, d) <= p2b.ORTHO_TOL
        # the constant direction is the first basis column, so the projected
        # residual must sum to zero
        assert abs(float(d.sum())) <= 1e-8 * (1.0 + float(np.linalg.norm(d)))


def test_projection_is_idempotent():
    rng = np.random.default_rng(2)
    Q = p2b.nuisance_basis(rng.normal(size=120))
    d = rng.normal(size=120)
    p1 = p2b.project_np(Q, d)
    assert np.abs(p2b.project_np(Q, p1) - p1).max() < 1e-12


# ---------------------------------------------- prior cancellation & antisymmetry
def test_prior_cancels_and_ligand_order_flips_the_sign():
    rng = np.random.default_rng(3)
    L = 90
    uh = rng.normal(size=(L, p2b.K))
    V = rng.normal(size=(p2b.K, p2b.D_ATOM))
    ga, gb = rng.normal(size=p2b.D_ATOM), rng.normal(size=p2b.D_ATOM)
    bP = rng.normal(size=L)
    Q = p2b.nuisance_basis(bP)
    d_ab = p2b.project_np(Q, uh @ (V @ (ga - gb)))
    d_ba = p2b.project_np(Q, uh @ (V @ (gb - ga)))
    assert np.abs(d_ab + d_ba).max() < 1e-12          # exact antisymmetry
    s_a = bP + p2b.project_np(Q, uh @ (V @ ga))
    s_b = bP + p2b.project_np(Q, uh @ (V @ gb))
    assert np.abs((s_a - s_b) - d_ab).max() < p2b.CANCEL_TOL


# ------------------------------------------------------------- AP estimator
def test_ap_exact_matches_plain_ap_when_there_are_no_ties():
    rng = np.random.default_rng(4)
    N = 400
    s = rng.normal(size=N)
    y = (rng.random(N) < 0.08).astype(np.int8)
    o = np.argsort(-s, kind="stable")
    yy = y[o].astype(float)
    plain = float((np.cumsum(yy) / np.arange(1, N + 1) * yy).sum() / yy.sum())
    assert abs(p2b.ap_exact(s, y) - plain) < 1e-12


def test_ap_exact_matches_the_monte_carlo_expectation_under_ties():
    rng = np.random.default_rng(5)
    N = 300
    s = (rng.random(N) < 0.25).astype(np.float64)      # two large tied blocks
    y = (rng.random(N) < 0.10).astype(np.int8)
    vals = []
    for _ in range(600):
        o = np.lexsort((rng.random(N), -s))
        z = y[o].astype(float)
        vals.append(float((np.cumsum(z) / np.arange(1, N + 1) * z).sum() / z.sum()))
    assert abs(p2b.ap_exact(s, y) - float(np.mean(vals))) < 0.02


def test_ap_exact_returns_none_on_degenerate_labels():
    s = np.arange(10.0)
    assert p2b.ap_exact(s, np.zeros(10, dtype=np.int8)) is None
    assert p2b.ap_exact(s, np.ones(10, dtype=np.int8)) is None


def test_chance_is_the_constant_score_expectation_and_tracks_prevalence():
    y = np.zeros(200, dtype=np.int8)
    y[:20] = 1
    c = p2b.chance_ap(y)
    assert 0.05 < c < 0.20                    # order of the 0.10 prevalence
    assert abs(c - p2b.ap_exact(np.full(200, 3.0), y)) < 1e-12


# ------------------------------------------------------------- pair metrics
def test_pair_metrics_reward_the_correct_direction():
    L = 60
    gain, loss = {1, 2, 3}, {40, 41, 42}
    d = np.zeros(L)
    for r in gain:
        d[r] = 5.0
    for r in loss:
        d[r] = -5.0
    m = p2b.pair_metrics(d, gain, loss, L)
    assert m["ap_gain"] > 0.99 and m["ap_loss"] > 0.99
    assert m["ap_bidir"] > 0.99
    assert m["ap_bidir"] > m["chance_bidir"]
    flipped = p2b.pair_metrics(-d, gain, loss, L)
    assert flipped["ap_bidir"] < m["ap_bidir"]


def test_pair_metrics_is_none_without_a_symmetric_difference():
    assert p2b.pair_metrics(np.zeros(20), set(), set(), 20) is None


# ------------------------------------------------------------- aggregation
def test_aggregation_is_invariant_to_full_pair_duplication():
    pv = {"a": 0.2, "b": 0.8, "c": 0.5}
    pc = {"a": "c1", "b": "c1", "c": "c2"}
    cc = {"c1": "K", "c2": "K"}
    _c1, m1 = p2b.aggregate(pv, pc, cc)
    pv2 = dict(pv, a2=0.2, b2=0.8, c2=0.5)
    pc2 = dict(pc, a2="c1", b2="c1", c2="c2")
    _c2, m2 = p2b.aggregate(pv2, pc2, cc)
    assert abs(m1 - m2) < 1e-12


def test_aggregation_balances_constructs_not_pairs():
    """A construct with many pairs must not outweigh one with few."""
    pv = {f"p{i}": 1.0 for i in range(100)}
    pc = {f"p{i}": "big" for i in range(100)}
    pv["q"] = 0.0
    pc["q"] = "small"
    cc = {"big": "K", "small": "K"}
    _c, m = p2b.aggregate(pv, pc, cc)
    assert abs(m - 0.5) < 1e-12


def test_component_bootstrap_reports_units_and_a_lower_bound():
    a = {f"c{i}": 0.5 for i in range(40)}
    b = {f"c{i}": 0.2 for i in range(40)}
    r = p2b.component_bootstrap(a, b, n_boot=500, seed=1)
    assert r["units"] == 40
    assert abs(r["delta"] - 0.3) < 1e-12
    assert r["lcb95_one_sided"] > 0


# ------------------------------------------------------------- pair building
class _Rec(dict):
    pass


def _rec(key, sk, gk, scaf):
    return {"source_key": key, "seq_key": sk, "graph_key": gk, "scaffold": scaf}


def test_build_pairs_applies_every_frozen_exclusion():
    recs = [_rec("a", "S", "g1", "sc1"), _rec("b", "S", "g2", "sc2"),
            _rec("c", "S", "g1", "sc1"), _rec("d", "S", "g3", "sc1"),
            _rec("e", "S", "g4", "sc4")]
    masks = {"a": frozenset({1, 2}), "b": frozenset({2, 3}),
             "c": frozenset({1, 2}), "d": frozenset({5}),
             "e": frozenset({1, 2})}          # identical mask to a -> zero symdiff
    pairs, excl = p2b.build_pairs(recs, masks)
    got = {(p[1], p[2]) for p in pairs}
    assert ("a", "c") not in got and ("a", "c")[::-1] not in got   # same graph
    assert ("a", "d") not in got                                   # same scaffold
    assert ("a", "e") not in got                                   # zero symdiff
    assert ("a", "b") in got
    assert excl["same_graph"] >= 1
    assert excl["scaffold_not_distinct"] >= 1
    assert excl["zero_symmetric_difference"] >= 1


def test_build_pairs_is_unordered():
    recs = [_rec("a", "S", "g1", "s1"), _rec("b", "S", "g2", "s2")]
    masks = {"a": frozenset({1}), "b": frozenset({2})}
    pairs, _e = p2b.build_pairs(recs, masks)
    assert len(pairs) == 1


# ------------------------------------------------------------- sampler
def test_hierarchical_sampler_visits_every_component_and_respects_caps():
    pairs = []
    for c in range(5):
        for s in range(4):
            for p in range(20):
                pairs.append((f"sk{c}_{s}", f"a{c}_{s}_{p}", f"b{c}_{s}_{p}"))
    cc = {f"sk{c}_{s}": f"K{c}" for c in range(5) for s in range(4)}
    chosen = p2b.hierarchical_sample(pairs, cc, epoch=0)
    comps = {c for c, _sk, _pl in chosen}
    assert comps == {f"K{c}" for c in range(5)}
    per_comp = {}
    for c, sk, pl in chosen:
        per_comp.setdefault(c, []).append(sk)
        assert len(pl) <= p2b.P_MAX
    assert all(len(v) <= p2b.C_MAX for v in per_comp.values())


def test_hierarchical_sampler_is_seed_deterministic_per_epoch():
    pairs = [(f"sk{i%6}", f"a{i}", f"b{i}") for i in range(120)]
    cc = {f"sk{i}": f"K{i%3}" for i in range(6)}
    assert (p2b.hierarchical_sample(pairs, cc, 2)
            == p2b.hierarchical_sample(pairs, cc, 2))
    assert (p2b.hierarchical_sample(pairs, cc, 2)
            != p2b.hierarchical_sample(pairs, cc, 3))


# ------------------------------------------------------------- controls
def test_derangement_has_no_fixed_points():
    rng = np.random.default_rng(7)
    for n in range(2, 30):
        p = p2b.derange(n, rng)
        assert sorted(p) == list(range(n))
        assert all(p[i] != i for i in range(n))
    assert p2b.derange(1, rng) == [0]


def test_context_shuffle_preserves_composition_within_amino_acid_type():
    rng = np.random.default_rng(8)
    seq = "ACDACDACDACD"
    h = rng.normal(size=(len(seq), 5))
    out = p2b.context_shuffle(h, seq, seed=11)
    assert out.shape == h.shape
    for aa in set(seq):
        idx = [i for i, c in enumerate(seq) if c == aa]
        a = np.sort(h[idx].sum(1))
        b = np.sort(out[idx].sum(1))
        assert np.allclose(a, b)              # same multiset of states per type
    # and it is not the identity on a type with several members
    assert not np.allclose(out, h)


# ------------------------------------------------------------- loss
def test_pair_loss_is_none_without_a_change_and_finite_otherwise():
    import torch
    ds = torch.zeros(30, requires_grad=True)
    assert p2b.pair_loss(ds, set(), set(), 30, "cpu") is None
    lo = p2b.pair_loss(ds, {1, 2}, {5}, 30, "cpu")
    assert lo is not None and torch.isfinite(lo)
    lo.backward()
    assert ds.grad is not None and float(ds.grad.abs().sum()) > 0


def test_pair_loss_prefers_the_correct_sign():
    import torch
    good = torch.tensor([3.0] * 3 + [-3.0] * 3 + [0.0] * 24)
    bad = -good
    g = p2b.pair_loss(good, {0, 1, 2}, {3, 4, 5}, 30, "cpu")
    b = p2b.pair_loss(bad, {0, 1, 2}, {3, 4, 5}, 30, "cpu")
    assert float(g) < float(b)
