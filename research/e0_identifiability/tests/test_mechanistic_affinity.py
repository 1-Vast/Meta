import torch

from research.e0_identifiability.mechanistic_affinity import (
    EndpointCalibration, LocalMechanisticAffinityPotential, pairwise_rank_loss,
)


def _inputs():
    return dict(
        atom_state=torch.randn(2, 3, 128), atom_chemistry=torch.randn(2, 3, 40),
        atom_mask=torch.tensor([[1, 1, 0], [1, 1, 1.]]),
        residue_state=torch.randn(2, 4, 128),
        residue_chemistry=torch.randn(2, 4, 6),
        residue_mask=torch.tensor([[1, 1, 1, 0], [1, 1, 1, 1.]]),
        contact_prob=torch.rand(2, 3, 4), distance_prob=torch.softmax(
            torch.randn(2, 3, 4, 5), dim=-1),
    )


def test_map_is_pair_local_and_permutation_invariant():
    torch.manual_seed(4)
    model = LocalMechanisticAffinityPotential()
    values = _inputs()
    original = model(**values).potential
    permutation = torch.tensor([1, 0, 2])
    permuted = dict(values)
    for key in ("atom_state", "atom_chemistry", "atom_mask"):
        permuted[key] = values[key][:, permutation]
    for key in ("contact_prob", "distance_prob"):
        permuted[key] = values[key][:, permutation]
    assert torch.allclose(original, model(**permuted).potential, atol=1e-6)


def test_endpoint_scales_are_positive_and_rank_loss_rewards_order():
    calibration = EndpointCalibration()
    values = calibration(torch.tensor([1.0, 1.0]), torch.tensor([0, 1]))
    assert torch.isfinite(values).all()
    tasks = torch.zeros(3, dtype=torch.long)
    labels = torch.tensor([0.0, 1.0, 2.0])
    assert pairwise_rank_loss(labels, labels, tasks) < pairwise_rank_loss(-labels, labels, tasks)
