"""Held-out synthetic gates for label-locked residual kernels."""
import torch

from model.qpsmp_meta import EvidenceLockedMetaTransport


def _predict(router, batch):
    _, correction, _, _ = router(*batch[:-1])
    return correction


def _fit(router, generator, target_kind, steps=160, tasks=128):
    optimizer = torch.optim.Adam(router.parameters(), lr=1e-2)
    for _ in range(steps):
        batch = _batch(generator, target_kind, tasks)
        prediction = _predict(router, batch)
        loss = (prediction - batch[-1]).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()


def _batch(generator, kind, tasks, queries=8, width=8, supports=4):
    support = torch.zeros(tasks, supports, width)
    support[:, torch.arange(supports), torch.arange(supports)] = 1.0
    active = torch.randint(supports, (tasks, queries), generator=generator)
    query = torch.nn.functional.one_hot(active, width).float()
    support_phi = torch.zeros(tasks, supports, supports)
    query_phi = torch.zeros(tasks, queries, supports)
    amplitude = torch.randn(tasks, supports, generator=generator)
    if kind == "shared":
        residual = amplitude
        target = amplitude.gather(1, active)
    elif kind == "level":
        residual = torch.zeros(tasks, supports)
        target = torch.zeros(tasks, queries)
    elif kind == "private":
        residual = amplitude
        target = torch.randn(tasks, queries, generator=generator)
    else:
        raise ValueError(kind)
    return support, support_phi, residual, query, query_phi, target


def test_shared_interaction_clusters_transfer_to_heldout_tasks():
    torch.manual_seed(41001)
    router = EvidenceLockedMetaTransport(8, 4, hidden_dim=16)
    _fit(router, torch.Generator().manual_seed(41002), "shared")
    heldout = _batch(torch.Generator().manual_seed(73101), "shared", 512)
    prediction = _predict(router, heldout).detach()
    mse = (prediction - heldout[-1]).square().mean()
    assert mse < 0.45 * heldout[-1].square().mean()
    flipped = list(heldout)
    flipped[2] = -flipped[2]
    assert torch.allclose(_predict(router, tuple(flipped)), -prediction,
                          atol=1e-6, rtol=1e-6)


def test_zero_evidence_is_exactly_zero_correction():
    router = EvidenceLockedMetaTransport(8, 4, hidden_dim=16)
    heldout = _batch(torch.Generator().manual_seed(73102), "level", 128)
    assert torch.count_nonzero(_predict(router, heldout)) == 0


def test_private_mechanism_does_not_claim_heldout_transfer():
    torch.manual_seed(41005)
    router = EvidenceLockedMetaTransport(8, 4, hidden_dim=16)
    _fit(router, torch.Generator().manual_seed(41006), "private")
    heldout = _batch(torch.Generator().manual_seed(73103), "private", 512)
    prediction = _predict(router, heldout).detach()
    baseline = heldout[-1].square().mean()
    assert (prediction - heldout[-1]).square().mean() >= 0.98 * baseline
