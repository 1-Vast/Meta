"""Held-out synthetic admission gates for evidence-locked transport."""
import torch

from model.qpsmp_meta import EvidenceLockedMetaTransport


def _predict(router, batch):
    _, coefficient, reliability, _ = router(*batch[:-1])
    return reliability.unsqueeze(-1).mul(coefficient).mul(batch[5]).sum(-1)


def _fit(router, generator, target_kind, steps=160, tasks=128):
    optimizer = torch.optim.Adam(router.parameters(), lr=1e-2)
    for _ in range(steps):
        batch = _batch(generator, target_kind, tasks)
        prediction = _predict(router, batch)
        loss = (prediction - batch[-1]).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()


def _batch(generator, kind, tasks, queries=8, width=8, mechanisms=4):
    active = torch.randint(mechanisms, (tasks,), generator=generator)
    protein = torch.nn.functional.one_hot(active, width).float()
    support_ligand = torch.randn(tasks, 1, width, generator=generator)
    query_ligand = torch.randn(tasks, queries, width, generator=generator)
    support_phi = torch.zeros(tasks, 1, mechanisms)
    query_phi = torch.zeros(tasks, queries, mechanisms)
    support_sign = torch.sign(torch.randn(tasks, 1, generator=generator))
    query_sign = torch.sign(torch.randn(tasks, queries, generator=generator))
    support_phi[torch.arange(tasks), 0, active] = support_sign[:, 0]
    query_phi[torch.arange(tasks)[:, None], torch.arange(queries), active[:, None]] = query_sign
    amplitude = torch.sign(torch.randn(tasks, generator=generator))
    if kind == "shared":
        residual = amplitude[:, None] * support_sign
        target = amplitude[:, None] * query_sign
    elif kind == "level":
        residual = 0.05 * torch.randn(tasks, 1, generator=generator)
        target = torch.zeros(tasks, queries)
    elif kind == "private":
        residual = amplitude[:, None] * support_sign
        target = torch.sign(torch.randn(tasks, queries, generator=generator))
    else:
        raise ValueError(kind)
    return (protein, support_ligand, support_phi, residual,
            query_ligand, query_phi, target)


def test_oracle_shared_primitive_transfers_to_heldout_tasks_k1():
    torch.manual_seed(41001)
    router = EvidenceLockedMetaTransport(8, 4, hidden_dim=16)
    _fit(router, torch.Generator().manual_seed(41002), "shared")
    heldout = _batch(torch.Generator().manual_seed(73101), "shared", 512)
    prediction = _predict(router, heldout).detach()
    mse = (prediction - heldout[-1]).square().mean()
    assert mse < 0.45 * heldout[-1].square().mean()
    flipped = list(heldout)
    flipped[3] = -flipped[3]
    assert torch.allclose(_predict(router, tuple(flipped)), -prediction,
                          atol=1e-6, rtol=1e-6)


def test_level_only_trains_elmt_to_negligible_correction():
    torch.manual_seed(41003)
    router = EvidenceLockedMetaTransport(8, 4, hidden_dim=16)
    _fit(router, torch.Generator().manual_seed(41004), "level")
    heldout = _batch(torch.Generator().manual_seed(73102), "level", 512)
    correction = _predict(router, heldout).detach()
    assert correction.square().mean() < 5e-4


def test_private_mechanism_does_not_claim_heldout_transfer():
    torch.manual_seed(41005)
    router = EvidenceLockedMetaTransport(8, 4, hidden_dim=16)
    _fit(router, torch.Generator().manual_seed(41006), "private")
    heldout = _batch(torch.Generator().manual_seed(73103), "private", 512)
    prediction = _predict(router, heldout).detach()
    baseline = heldout[-1].square().mean()
    assert (prediction - heldout[-1]).square().mean() >= 0.98 * baseline
    assert prediction.square().mean() < 0.05 * baseline
