"""Fast structural falsification gates for TERM before biological training."""
import torch

from model.qpsmp_meta import TriadicEvidenceRouter


def _fit(router, protein, support_ligand, support_phi, residual,
         query_ligand, query_phi, target, steps=120):
    optimizer = torch.optim.Adam(router.parameters(), lr=2e-2)
    for _ in range(steps):
        _, coefficient, reliability, _ = router(
            protein, support_ligand, support_phi, residual,
            query_ligand, query_phi)
        prediction = (reliability.unsqueeze(-1) * coefficient * query_phi).sum(-1)
        loss = (prediction - target).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    return prediction.detach(), reliability.detach()


def test_shared_primitive_oracle_is_one_shot_learnable():
    torch.manual_seed(41001)
    tasks, mechanisms, width, queries = 96, 4, 8, 8
    active = torch.arange(tasks) % mechanisms
    protein = torch.nn.functional.one_hot(active, width).float()
    support_ligand = torch.randn(tasks, 1, width)
    query_ligand = torch.randn(tasks, queries, width)
    support_phi = torch.sign(torch.randn(tasks, 1, mechanisms))
    query_phi = torch.sign(torch.randn(tasks, queries, mechanisms))
    amplitude = torch.where(torch.arange(tasks) % 2 == 0, 1., -1.)
    support_target = amplitude[:, None] * support_phi[
        torch.arange(tasks), :, active]
    query_target = amplitude[:, None] * query_phi[
        torch.arange(tasks), :, active]
    router = TriadicEvidenceRouter(width, mechanisms, hidden_dim=16)
    prediction, _ = _fit(
        router, protein, support_ligand, support_phi, support_target,
        query_ligand, query_phi, query_target)
    baseline = query_target.square().mean()
    assert (prediction - query_target).square().mean() < 0.6 * baseline


def test_private_mechanism_carries_no_query_information():
    """A support-private random bit cannot predict independent query bits."""
    values = []
    for seed in range(41001, 41011):
        generator = torch.Generator().manual_seed(seed)
        support = torch.sign(torch.randn(512, 1, generator=generator))
        query = torch.sign(torch.randn(512, 32, generator=generator))
        copied = support.expand_as(query)
        values.append(float((copied - query).square().mean()
                            - query.square().mean()))
    # Blindly transferring private support evidence must be worse than abstaining.
    assert sum(values) / len(values) > 0.8


def test_level_only_exact_gradient_has_zero_expected_primitive_signal():
    values = []
    for seed in range(41001, 41011):
        generator = torch.Generator().manual_seed(seed)
        residual = torch.randn(512, 5, generator=generator) * 0.1
        primitive = torch.sign(torch.randn(512, 5, 8, generator=generator))
        gradient = TriadicEvidenceRouter.exact_coefficient_gradient(
            residual, primitive).mean((0, 1))
        values.append(float(gradient.square().mean()))
    assert sum(values) / len(values) < 1e-5
