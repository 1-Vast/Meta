"""Contract tests for the gauge-free direct-W ordinal witness."""
from __future__ import annotations

import inspect
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
import s2r_run as G  # noqa: E402


def test_preregistration_is_frozen():
    assert G.sha_file(G.PREREG) == G.PREREG_SHA


def test_direct_w_has_exact_parameter_surface_and_unit_norm():
    head = G.DirectW()
    assert tuple(head.W.shape) == (1280, 41)
    assert sum(p.numel() for p in head.parameters()) == 52480
    assert abs(float(head.W.norm()) - 1.0) < 1e-6
    assert not hasattr(head, "U") and not hasattr(head, "V")


def test_direction_normalization_removes_positive_scale():
    d = torch.tensor([0.3, -1.2, 0.8, 0.1], dtype=torch.float64)
    assert torch.max(torch.abs(G.normalized(d) - G.normalized(7.0 * d))) < 1e-10


def test_unit_sphere_projection_is_deterministic():
    a = G.DirectW(3)
    b = G.DirectW(3)
    with torch.no_grad():
        a.W.mul_(13.0)
        a.project_norm()
    assert torch.allclose(a.W, b.W, rtol=1e-6, atol=1e-9)


def test_raw_projection_matches_parent_float64_operation():
    rng = np.random.default_rng(2)
    raw = rng.normal(size=17)
    q, _ = np.linalg.qr(rng.normal(size=(17, 2)))
    parent = torch.from_numpy(raw).double()
    qt = torch.from_numpy(q).double()
    parent = (parent - qt @ (qt.T @ parent)).numpy()
    direct = G.project_np(q, raw)
    assert np.max(np.abs(parent - direct)) < 1e-12


def test_s2r_has_no_factorized_or_parallel_head():
    source = inspect.getsource(G)
    assert "head.U" not in source
    assert "head.V" not in source
    assert "DirectW" in source
    assert G.AP_THRESHOLD == 0.50
    assert G.CALIBRATION_SEEDS == (20260931, 20260932, 20260933)


def test_sealed_seed_is_instantiated_after_calibration_gate_write():
    source = inspect.getsource(G.main)
    gate = source.index('write_json(OUT / "S2R_CALIBRATION_GATE.json"')
    sealed = source.index("run_seed(SEALED_SEED")
    assert gate < sealed
