import torch

from research.meta_fewshot.train_v1_development import PairPriorMetaSection


def test_pair_prior_query_api_has_no_query_labels_and_gradients_are_finite():
    torch.manual_seed(2)
    model = PairPriorMetaSection(8, 2, 1.0, 4)
    ligand, pair, y = torch.randn(10, 8), torch.randn(10, 8), torch.randn(10)
    prediction = model.episode(ligand[:5], pair[:5], y[:5], ligand[5:], pair[5:])
    prediction.square().mean().backward()
    assert model.raw_basis.grad is not None
    assert all(parameter.grad is None or torch.isfinite(parameter.grad).all()
               for parameter in model.parameters())


def test_pair_prior_d0_is_support_independent():
    model = PairPriorMetaSection(8, 0, 1.0, 4)
    ligand, pair = torch.randn(10, 8), torch.randn(10, 8)
    first = model.episode(ligand[:5], pair[:5], torch.zeros(5), ligand[5:], pair[5:])
    second = model.episode(ligand[:5], pair[:5], torch.ones(5), ligand[5:], pair[5:])
    assert torch.equal(first, second)


def test_support_row_permutation_is_invariant():
    torch.manual_seed(3)
    model = PairPriorMetaSection(8, 2, 1.0, 4)
    ligand, pair, y = torch.randn(10, 8), torch.randn(10, 8), torch.randn(10)
    expected = model.episode(ligand[:5], pair[:5], y[:5], ligand[5:], pair[5:])
    order = torch.tensor([3, 0, 4, 1, 2])
    actual = model.episode(ligand[order], pair[order], y[order], ligand[5:], pair[5:])
    assert torch.allclose(expected, actual, atol=1e-6)


def test_task_state_rank_cannot_exceed_min_d_k():
    model = PairPriorMetaSection(8, 3, 1.0, 4)
    _, coordinates = model.components(torch.randn(2, 8), torch.randn(2, 8))
    assert torch.linalg.matrix_rank(coordinates) <= min(3, 2)
