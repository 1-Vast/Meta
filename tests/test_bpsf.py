import pytest
import torch

from model.bpsf import BipartitePairSectionFormer, QuotientSupportSetOperator
from model.mechanism import GeometrySupervisionHead


def test_pair_section_shapes_masks_and_geometry_head():
    torch.manual_seed(71)
    trunk = BipartitePairSectionFormer(8, 6, blocks=1, latent_count=4, heads=2)
    atoms = torch.randn(3, 5, 8, requires_grad=True)
    residues = torch.randn(3, 7, 8, requires_grad=True)
    atom_mask = torch.tensor([[1, 1, 1, 0, 0], [1] * 5, [1, 1, 0, 0, 0]])
    residue_mask = torch.tensor([[1, 1, 1, 0, 0, 0, 0], [1] * 7, [1, 1, 1, 1, 0, 0, 0]])

    encoded = trunk(atoms, residues, atom_mask, residue_mask, return_pair=True)
    geometry = GeometrySupervisionHead(8)(encoded.pair, encoded.pair_mask)

    assert encoded.endpoint.shape == (3, 8)
    assert encoded.section.shape == (3, 6)
    assert encoded.pair.shape == (3, 5, 7, 8)
    assert geometry.contact_logits.shape == (3, 5, 7)
    assert geometry.distance_logits.shape == (3, 5, 7, 5)
    assert torch.count_nonzero(encoded.pair[0, 3:]) == 0
    (encoded.endpoint.square().mean() + encoded.section.square().mean()
     + geometry.contact_logits.square().mean()).backward()
    assert atoms.grad is not None and residues.grad is not None


def test_quotient_support_operator_invariants_and_gradient():
    torch.manual_seed(72)
    operator = QuotientSupportSetOperator(7, hidden_dim=16, heads=2)
    support = torch.randn(4, 7, requires_grad=True)
    query = torch.randn(3, 7, requires_grad=True)
    residual = torch.randn(4)
    residual = residual - residual.mean()

    state, _, prediction = operator(support, query, residual)
    centered = support - support.mean(0, keepdim=True)
    projection = torch.linalg.pinv(centered) @ centered
    assert torch.allclose(state, projection @ state, atol=1e-5, rtol=1e-5)

    order = torch.tensor([2, 0, 3, 1])
    permuted = operator(support[order], query, residual[order])[2]
    shifted = operator(support, query, residual + 3.0)[2]
    assert torch.allclose(prediction, permuted, atol=1e-6)
    assert torch.allclose(prediction, shifted, atol=1e-6)

    zero = operator(support, query, torch.zeros_like(residual))
    one = operator(support[:1], query, torch.zeros(1))
    assert torch.count_nonzero(zero[0]) == 0
    assert torch.count_nonzero(zero[2]) == 0
    assert torch.count_nonzero(one[0]) == 0
    assert torch.count_nonzero(one[2]) == 0

    prediction.square().mean().backward()
    assert support.grad is not None and query.grad is not None
    assert operator.weight.weight.grad is not None


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA BPSF validation")
def test_pair_section_cuda_amp_is_finite():
    trunk = BipartitePairSectionFormer(
        32, 16, blocks=1, latent_count=8, heads=4, chunk_size=8).cuda()
    atoms = torch.randn(8, 32, 32, device="cuda")
    residues = torch.randn(8, 128, 32, device="cuda")
    atom_mask = torch.ones(8, 32, device="cuda")
    residue_mask = torch.ones(8, 128, device="cuda")
    with torch.autocast("cuda", dtype=torch.float16):
        encoded = trunk(atoms, residues, atom_mask, residue_mask)
        loss = encoded.endpoint.square().mean() + encoded.section.square().mean()
    loss.backward()
    assert torch.isfinite(encoded.endpoint).all()
    assert torch.isfinite(encoded.section).all()
