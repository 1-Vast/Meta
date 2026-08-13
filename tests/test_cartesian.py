import torch
import pytest

from model.cartesian import SparseCartesianMechanismEncoder


def _orthogonal(dtype=torch.float64):
    q, _ = torch.linalg.qr(torch.randn(3, 3, dtype=dtype))
    q[:, 0] *= -1  # exercise reflection as well as rotation
    return q


def _edges(batch, nodes):
    edges = []
    for b in range(batch):
        offset = b * nodes
        for i in range(nodes):
            for j in range(nodes):
                if i != j:
                    edges.append((offset + j, offset + i))
    return torch.tensor(edges).T


def test_cartesian_channels_are_o3_equivariant_and_scalar_slots_invariant():
    torch.manual_seed(81)
    model = SparseCartesianMechanismEncoder(
        5, 8, 4, 3, layers=2, rbf_count=8, mechanism_slots=4,
        slot_dim=8, dtype=torch.float64)
    features = torch.randn(2, 5, 5, dtype=torch.float64)
    mask = torch.ones(2, 5, dtype=torch.bool)
    coordinates = torch.randn(2, 5, 3, dtype=torch.float64)
    edges = _edges(2, 5)
    rotation = _orthogonal()
    left = model(features, mask, coordinates, edges)
    right = model(features, mask, coordinates @ rotation.T + 7.0, edges)
    expected_vector = torch.einsum("ij,bncj->bnci", rotation, left.state.vector)
    expected_tensor = torch.einsum(
        "ia,bncad,jd->bncij", rotation, left.state.tensor2, rotation)
    assert torch.allclose(left.state.scalar, right.state.scalar, atol=1e-9, rtol=1e-9)
    assert torch.allclose(expected_vector, right.state.vector, atol=1e-9, rtol=1e-9)
    assert torch.allclose(expected_tensor, right.state.tensor2, atol=1e-9, rtol=1e-9)
    assert torch.allclose(left.mechanism_slots, right.mechanism_slots,
                          atol=1e-9, rtol=1e-9)
    assert torch.allclose(right.state.tensor2, right.state.tensor2.transpose(-1, -2))
    assert torch.allclose(
        right.state.tensor2.diagonal(dim1=-2, dim2=-1).sum(-1),
        torch.zeros_like(right.state.tensor2[..., 0, 0]), atol=1e-9)


def test_cartesian_node_permutation_and_padding_are_safe():
    torch.manual_seed(82)
    model = SparseCartesianMechanismEncoder(
        4, 8, 4, 2, layers=1, mechanism_slots=4, slot_dim=8,
        dtype=torch.float64)
    features = torch.randn(1, 5, 4, dtype=torch.float64)
    coordinates = torch.randn(1, 5, 3, dtype=torch.float64)
    mask = torch.tensor([[1, 1, 1, 1, 0]], dtype=torch.bool)
    edges = _edges(1, 5)
    left = model(features, mask, coordinates, edges)
    order = torch.tensor([2, 0, 3, 1, 4])
    inverse = torch.argsort(order)
    mapping = inverse
    permuted_edges = mapping[edges]
    changed = model(features[:, order], mask[:, order],
                    coordinates[:, order], permuted_edges)
    assert torch.allclose(left.invariant_nodes[:, order], changed.invariant_nodes,
                          atol=1e-9, rtol=1e-9)
    assert torch.allclose(left.mechanism_slots, changed.mechanism_slots,
                          atol=1e-9, rtol=1e-9)
    extreme = coordinates.clone()
    extreme[:, -1] = 1e9
    padded = model(features, mask, extreme, edges)
    assert torch.allclose(left.mechanism_slots, padded.mechanism_slots,
                          atol=1e-9, rtol=1e-9)


def test_missing_geometry_has_zero_non_scalar_ranks_and_finite_gradients():
    torch.manual_seed(83)
    model = SparseCartesianMechanismEncoder(
        4, 8, 4, 2, layers=1, mechanism_slots=4, slot_dim=8)
    features = torch.randn(2, 5, 4, requires_grad=True)
    mask = torch.ones(2, 5, dtype=torch.bool)
    output = model(features, mask)
    assert not output.geometry_available.any()
    assert torch.count_nonzero(output.state.vector) == 0
    assert torch.count_nonzero(output.state.tensor2) == 0
    output.mechanism_slots.square().mean().backward()
    assert features.grad is not None and torch.isfinite(features.grad).all()


def test_cartesian_edges_cannot_cross_flattened_samples():
    model = SparseCartesianMechanismEncoder(
        4, 8, 4, 2, layers=1, mechanism_slots=4, slot_dim=8)
    features = torch.randn(2, 3, 4)
    mask = torch.ones(2, 3, dtype=torch.bool)
    coordinates = torch.randn(2, 3, 3)
    # Node 0 belongs to sample 0 and node 3 belongs to sample 1.
    with pytest.raises(ValueError, match="different flattened samples"):
        model(features, mask, coordinates, torch.tensor([[0], [3]]))
