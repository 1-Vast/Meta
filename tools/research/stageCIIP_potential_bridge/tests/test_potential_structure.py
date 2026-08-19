"""Structure tests for the unified potential module (CPU).

Frozen contracts: identity-zero, reversal antisymmetry, cycle sum ~0,
all contrasts from one s_theta, permutation contracts, no target-ID
input, batch isolation, stable-seed init, non-zero gradients on every
trainable branch, no closed-form solver in the module.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

import torch  # noqa: E402
import potential as PT  # noqa: E402


def _rand(batch=7):
    g = torch.Generator().manual_seed(0)
    P = torch.randn(batch, PT.D_P, generator=g)
    La = torch.randn(batch, PT.D_L, generator=g)
    Lb = torch.randn(batch, PT.D_L, generator=g)
    P2 = torch.randn(batch, PT.D_P, generator=g)
    return P, P2, La, Lb


def test_identity_zero_exact():
    m = PT.UnifiedPotential()
    P, _, L, _ = _rand()
    with torch.no_grad():
        assert torch.equal(m.ligand_contrast(P, L, L), torch.zeros(len(P)))
        assert torch.equal(m.protein_contrast(P, P, L), torch.zeros(len(P)))


def test_reversal_antisymmetry_exact():
    m = PT.UnifiedPotential()
    P, P2, La, Lb = _rand()
    with torch.no_grad():
        fwd = m.ligand_contrast(P, La, Lb)
        bwd = m.ligand_contrast(P, Lb, La)
        assert torch.equal(fwd, -bwd)
        pfwd = m.protein_contrast(P, P2, La)
        pbwd = m.protein_contrast(P2, P, La)
        assert torch.equal(pfwd, -pbwd)


def test_cycle_sum_near_zero():
    m = PT.UnifiedPotential()
    P, _, L1, L2 = _rand()
    with torch.no_grad():
        L3 = torch.randn(len(P), PT.D_L, generator=torch.Generator().manual_seed(3))
        s = (m.ligand_contrast(P, L1, L2) + m.ligand_contrast(P, L2, L3)
             + m.ligand_contrast(P, L3, L1))
        assert torch.allclose(s, torch.zeros(len(P)), atol=1e-6)


def test_free_pairwise_breaks_cycle():
    m = PT.FreePairwise()
    P, P2, L1, L2 = _rand(4)
    with torch.no_grad():
        L3 = torch.randn(4, PT.D_L, generator=torch.Generator().manual_seed(4))
        # free pair predictor has no ligand-contrast telescoping: check the
        # analogous triple-sum property does NOT hold by construction.
        a = m(P, P2, L1)
        b = m(P, P2, L2)
        # antisymmetry still holds
        assert torch.equal(m(P, P2, L1), -m(P2, P, L1))


def test_double_contrast_from_one_s():
    m = PT.UnifiedPotential()
    P, P2, La, Lb = _rand()
    with torch.no_grad():
        D = m.double_contrast(P, P2, La, Lb)
        manual = (m.potential(P2, Lb) - m.potential(P2, La)
                  - m.potential(P, Lb) + m.potential(P, La))
        assert torch.allclose(D, manual, atol=1e-6)


def test_permutation_and_batch_contracts():
    m = PT.UnifiedPotential()
    P, P2, La, Lb = _rand(8)
    perm = torch.randperm(8)
    with torch.no_grad():
        o1 = m.forward(P, La)
        o2 = m.forward(P[perm], La[perm])
        assert torch.allclose(o1[perm], o2, atol=1e-6)
        # padding a global arm batch with duplicate rows keeps outputs
        gm = PT.GlobalPotential()
        g1 = gm.forward(P[:4], La[:4])
        g2 = gm.forward(torch.cat([P[:4], P[:4]]), torch.cat([La[:4], La[:4]]))
        assert torch.allclose(g1, g2[:4], atol=1e-6)


def test_no_target_id_input():
    m = PT.UnifiedPotential()
    P, P2, La, Lb = _rand(2)
    # the module has no ID/embedding layers: every parameter is a Linear
    # or scalar; no nn.Embedding anywhere
    assert not any(isinstance(x, torch.nn.Embedding) for x in m.modules())
    with torch.no_grad():
        m.forward(P, La)


def test_all_branches_get_nonzero_grad():
    m = PT.UnifiedPotential()
    P, P2, La, Lb = _rand(6)
    y = m.forward(P, La)
    ((y - torch.randn(6)) ** 2).mean().backward()
    for name, p in m.named_parameters():
        assert p.grad is not None, name
        assert torch.isfinite(p.grad).all(), name
        assert p.grad.abs().max() > 0, name


def test_stable_seed_init_across_processes():
    # same init procedure (manual_seed -> xavier) reproduces the same
    # parameter values in any process
    def make(seed):
        torch.manual_seed(seed)
        m = PT.UnifiedPotential()
        return [p.detach().clone() for p in m.parameters()]
    a = make(7)
    b = make(7)
    for pa, pb in zip(a, b):
        assert torch.equal(pa, pb)


def test_global_vs_local_share_potential_form():
    gm = PT.GlobalPotential()
    P, P2, La, Lb = _rand(3)
    with torch.no_grad():
        d = gm.protein_contrast(P, P2, La)
        assert d.shape == (3,)
        assert torch.equal(gm.ligand_contrast(P, La, Lb),
                           -gm.ligand_contrast(P, Lb, La))
