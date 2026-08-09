"""Contract tests for the single-loss S1R synthetic repair."""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
import s1r_run as R  # noqa: E402


def test_preregistration_is_frozen():
    assert R.sha_file(R.PREREG) == R.PREREG_SHA


def test_pairwise_loss_rewards_correct_bidirectional_order():
    row = {"gi": np.array([0]), "li": np.array([2])}
    perfect = R.rank_loss_np(np.array([2.0, 0.0, -2.0]), row)
    tied = R.rank_loss_np(np.zeros(3), row)
    reversed_ = R.rank_loss_np(np.array([-2.0, 0.0, 2.0]), row)
    assert perfect < tied < reversed_


def test_pairwise_loss_is_antisymmetric_under_ligand_swap():
    d = np.array([1.2, -0.3, 0.1, -0.8])
    ab = {"gi": np.array([0]), "li": np.array([3])}
    ba = {"gi": np.array([3]), "li": np.array([0])}
    assert np.isclose(R.rank_loss_np(d, ab), R.rank_loss_np(-d, ba))


def test_torch_and_numpy_losses_match():
    d = np.array([0.7, -0.2, 0.1, -0.9, 0.3], dtype=np.float64)
    gain, loss = frozenset({0, 4}), frozenset({1, 3})
    row = {"gi": np.array(sorted(gain)), "li": np.array(sorted(loss))}
    got = float(R.rank_loss_torch(torch.tensor(d), gain, loss))
    assert np.isclose(got, R.rank_loss_np(d, row), rtol=1e-12, atol=1e-12)


def test_s1r_changes_only_the_loss_surface():
    source = inspect.getsource(R)
    assert "p2b_run" not in source
    assert "pair_loss(" not in source
    assert "rank_loss_torch" in source
    assert R.AP_THRESHOLD == 0.50
    assert R.PARAM_SEED == 20260901
    assert R.CALIBRATION_SEEDS == (20260921, 20260922, 20260923)
    assert R.SEALED_SEED == 20260998


def test_sealed_teacher_is_constructed_after_persisted_calibration_gate():
    source = inspect.getsource(R.main)
    gate_write = source.index('write_json(OUT / "S1R_CALIBRATION_GATE.json"')
    sealed_construct = source.index("Teacher(ctx, hcache, SEALED_SEED)")
    assert gate_write < sealed_construct

