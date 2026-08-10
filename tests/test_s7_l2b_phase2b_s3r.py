from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest
import torch


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "research" / "s7_l2b_r0r" / "s3r_run.py"
spec = importlib.util.spec_from_file_location("s3r_run_test", MODULE)
s3r = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(s3r)


def test_direct_w_shape_and_unit_norm():
    head = s3r.DirectW()
    assert tuple(head.W.shape) == (1280, 41)
    assert abs(float(head.W.norm()) - 1.0) < 1e-6


def test_direct_w_initialization_is_deterministic():
    left = s3r.DirectW()
    right = s3r.DirectW()
    assert torch.equal(left.W, right.W)


def test_positive_scale_normalization_is_invariant():
    score = torch.linspace(-2.0, 3.0, 97, dtype=torch.float64)
    assert torch.max(torch.abs(s3r.normalized(score) -
                               s3r.normalized(11.0 * score))) < 1e-10


def test_projection_is_orthogonal_and_pair_difference_is_antisymmetric():
    rng = np.random.default_rng(4)
    h = rng.normal(size=(37, 1280))
    W = rng.normal(size=(1280, 41))
    ga = rng.normal(size=41)
    gb = rng.normal(size=41)
    q, _ = np.linalg.qr(rng.normal(size=(37, 2)))
    left = s3r.project_np(q, (h @ W) @ (ga - gb))
    right = s3r.project_np(q, (h @ W) @ (gb - ga))
    assert np.max(np.abs(left + right)) < 1e-10
    assert np.linalg.norm(q.T @ left) / np.linalg.norm(left) < 1e-12


def test_single_ligand_residual_difference_matches_pair_score():
    rng = np.random.default_rng(5)
    h = rng.normal(size=(29, 1280))
    W = rng.normal(size=(1280, 41))
    ga = rng.normal(size=41)
    gb = rng.normal(size=41)
    q, _ = np.linalg.qr(rng.normal(size=(29, 2)))
    da = s3r.project_np(q, (h @ W) @ ga)
    db = s3r.project_np(q, (h @ W) @ gb)
    pair = s3r.project_np(q, (h @ W) @ (ga - gb))
    assert np.max(np.abs((da - db) - pair)) < 1e-10


def test_identical_ligand_has_zero_difference():
    rng = np.random.default_rng(6)
    h = rng.normal(size=(13, 1280))
    W = rng.normal(size=(1280, 41))
    g = rng.normal(size=41)
    assert np.array_equal((h @ W) @ (g - g), np.zeros(13))


def test_common_mask_audit_rejects_silent_intersection():
    metrics = {"x": ({"ap_bidir": 0.1}, "p")}
    arms = {
        "candidate": ({"a": metrics}, {"a": "s"}),
        "control": ({"b": metrics}, {"b": "s"}),
    }
    with pytest.raises(s3r.S3RContractError, match="mask differs"):
        s3r.assert_common_masks(arms, {"s": "c"})


def test_real_stream_has_registered_update_count_on_tiny_shape():
    pairs = []
    component = {}
    for index in range(554):
        sk = f"s{index}"
        component[sk] = f"c{index}"
        pairs.append((sk, f"a{index}", f"b{index}"))
    stream = s3r.build_stream(pairs, component)
    assert len(stream) == 210
    assert [row["update"] for row in stream] == list(range(1, 211))


def test_preregistration_explicitly_defers_r6():
    amendment = (ROOT / "research" / "s7_l2b_r0r" /
                 "PREREG_PHASE2B_S3R_AMENDMENT_01.md").read_text(encoding="utf-8")
    assert "S3R does not compute or inspect the old R6" in amendment
    assert "absolute ligand-feature origin" in amendment


def test_runner_does_not_overwrite_historical_gate():
    source = MODULE.read_text(encoding="utf-8")
    assert '"PHASE2B_S3R_GATE.json"' in source
    assert '"PHASE2B_GATE.json"' not in source


def test_r4_uses_correct_q_without_recomputing_prior():
    source = MODULE.read_text(encoding="utf-8")
    assert "context_shuffle(ctx.h(sk)" in source
    assert "project_np(ctx.Qs[sk], raw)" in source


def test_affinity_sources_are_not_imported():
    source = MODULE.read_text(encoding="utf-8").lower()
    for forbidden in ("chembl", "bindingdb", "davis", "kiba", "recipient"):
        assert forbidden not in source
