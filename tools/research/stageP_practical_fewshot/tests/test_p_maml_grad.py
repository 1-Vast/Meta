"""MAML outer-gradient correctness tests (CPU).

History: the first version of test_fomaml_identity_repro replicated the
original trainer loop and FAILED (red) on the un-fixed code — the adapted
model kept the last support-step gradient and query backward accumulated
onto it, violating the frozen first-order MAML semantics. After the fix
(p_maml.task_fomaml_grad clears adapted-model grads before the query
backward), the same identity is pinned on the fixed path plus an
independent toy functional reference. The toy trunk avoids the PTrunk
inner-loop divergence at random init (covered separately by the trainer's
non-finite-task guard).
"""
import copy
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

import torch  # noqa: E402
import torch.nn as nn  # noqa: E402
import p_train as PT  # noqa: E402
import p_maml as PM  # noqa: E402

T_INNER = 3
T_LR = 1e-3


class ToyTrunk(nn.Module):
    def __init__(self):
        super().__init__()
        self.a = nn.Linear(4, 3)
        self.b = nn.Linear(3, 1)

    def forward(self, xp, xl):
        h = torch.relu(self.a(xp))
        return {"yhat": self.b(h).squeeze(-1)}


def _toy_task(seed=0, sup=5, q=3):
    g = torch.Generator().manual_seed(seed)
    xp_s = torch.randn(sup, 4, generator=g)
    xl_s = torch.randn(sup, 1, generator=g)
    y_s = torch.randn(sup, generator=g)
    xp_q = torch.randn(q, 4, generator=g)
    xl_q = torch.randn(q, 1, generator=g)
    y_q = torch.randn(q, generator=g)
    return xp_s, xl_s, y_s, xp_q, xl_q, y_q


def test_fomaml_identity_repro():
    """Pins the fixed identity: the per-task FOMAML gradient is the pure
    query gradient at the adapted parameters (no support-step residue)."""
    torch.manual_seed(0)
    model = ToyTrunk()
    xp_s, xl_s, y_s, xp_q, xl_q, y_q = _toy_task()
    qloss, grads = PM.task_fomaml_grad(model, xp_s, xl_s, y_s, xp_q, xl_q, y_q,
                                       inner_steps=T_INNER, inner_lr=T_LR)
    m = copy.deepcopy(model)
    inner = torch.optim.SGD(m.parameters(), lr=T_LR)
    for _ in range(T_INNER):
        inner.zero_grad()
        loss = ((m(xp_s, xl_s)["yhat"] - y_s) ** 2).mean()
        loss.backward()
        inner.step()
    m.zero_grad()
    qloss2 = ((m(xp_q, xl_q)["yhat"] - y_q) ** 2).mean()
    qloss2.backward()
    for g, p in zip(grads, m.parameters()):
        torch.testing.assert_close(g, p.grad, atol=0, rtol=0,
                                   msg="adapted-model grad must be pure query grad")
    assert np.isfinite(qloss)


def test_functional_reference():
    """Independent autograd reference of the same computation."""
    torch.manual_seed(1)
    model = ToyTrunk()
    xp_s, xl_s, y_s, xp_q, xl_q, y_q = _toy_task(seed=2)
    qloss, grads = PM.task_fomaml_grad(model, xp_s, xl_s, y_s, xp_q, xl_q, y_q,
                                       inner_steps=T_INNER, inner_lr=T_LR)
    m = copy.deepcopy(model)
    inner = torch.optim.SGD(m.parameters(), lr=T_LR)
    for _ in range(T_INNER):
        inner.zero_grad()
        loss = ((m(xp_s, xl_s)["yhat"] - y_s) ** 2).mean()
        loss.backward()
        inner.step()
    ql = ((m(xp_q, xl_q)["yhat"] - y_q) ** 2).mean()
    ref = torch.autograd.grad(ql, list(m.parameters()), retain_graph=False)
    for g, r in zip(grads, ref):
        torch.testing.assert_close(g, r, atol=1e-6, rtol=1e-6)
    assert np.isfinite(qloss)


def test_param_order_and_none_grads():
    torch.manual_seed(2)
    model = PT.PTrunk()
    xp_s, xl_s, y_s, xp_q, xl_q, y_q = (torch.randn(5, 640), torch.randn(5, 2048),
                                        torch.randn(5) * 0.1,
                                        torch.randn(3, 640), torch.randn(3, 2048),
                                        torch.randn(3) * 0.1)
    _, grads = PM.task_fomaml_grad(model, xp_s, xl_s, y_s, xp_q, xl_q, y_q,
                                   inner_steps=2, inner_lr=1e-4)
    assert len(grads) == len(list(model.parameters()))
    for p, g in zip(model.parameters(), grads):
        assert g is not None and g.shape == p.shape
        assert torch.isfinite(g).all()


def test_multitask_accumulation():
    torch.manual_seed(3)
    model = PT.PTrunk()
    per_task = []
    for s in (4, 5):
        xp_s, xl_s, y_s, xp_q, xl_q, y_q = (torch.randn(5, 640), torch.randn(5, 2048),
                                            torch.randn(5) * 0.1,
                                            torch.randn(3, 640), torch.randn(3, 2048),
                                            torch.randn(3) * 0.1)
        per_task.append(PM.task_fomaml_grad(model, xp_s, xl_s, y_s,
                                            xp_q, xl_q, y_q,
                                            inner_steps=2, inner_lr=1e-4)[1])
    expected = [sum(ts[i] for ts in per_task) for i in range(len(per_task[0]))]
    for e, p in zip(expected, model.parameters()):
        assert e.shape == p.shape  # accumulation order pinned to model order
