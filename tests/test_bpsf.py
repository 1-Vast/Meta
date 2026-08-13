import pytest
import torch

from model.bpsf import BipartitePairSectionFormer
from model.geometry_supervision import GeometrySupervisionHead


def test_pair_section_shapes_masks_and_geometry_head():
    torch.manual_seed(71)
    trunk = BipartitePairSectionFormer(8, 6, blocks=1, latent_count=4, heads=2)
    atoms = torch.randn(3, 5, 8, requires_grad=True)
    residues = torch.randn(3, 7, 8, requires_grad=True)
    atom_mask = torch.tensor([[1, 1, 1, 0, 0], [1] * 5, [1, 1, 0, 0, 0]])
    residue_mask = torch.tensor([[1, 1, 1, 0, 0, 0, 0], [1] * 7, [1, 1, 1, 1, 0, 0, 0]])

    encoded = trunk(atoms, residues, atom_mask, residue_mask, return_pair=True)
    geometry = GeometrySupervisionHead(48)(encoded.pair, encoded.pair_mask)

    assert encoded.endpoint.shape == (3, 8)
    assert encoded.section.shape == (3, 6)
    assert encoded.mechanism_slots.shape == (3, 4, 48)
    assert encoded.mechanism_response.shape == (3, 4)
    assert encoded.pair.shape == (3, 5, 7, 48)
    assert geometry.contact_logits.shape == (3, 5, 7)
    assert geometry.distance_logits.shape == (3, 5, 7, 5)
    assert torch.count_nonzero(encoded.pair[0, 3:]) == 0
    (encoded.endpoint.square().mean() + encoded.section.square().mean()
     + geometry.contact_logits.square().mean()).backward()
    assert atoms.grad is not None and residues.grad is not None


def test_protein_projection_and_pair_trunk_are_padding_invariant():
    torch.manual_seed(73)
    trunk = BipartitePairSectionFormer(
        8, 5, pair_dim=16, blocks=1, latent_count=4, heads=2)
    atoms = torch.randn(2, 5, 8)
    residues = torch.randn(2, 7, 8)
    atom_mask = torch.tensor([[1, 1, 1, 0, 0], [1, 1, 0, 0, 0]])
    residue_mask = torch.tensor([[1, 1, 1, 1, 0, 0, 0], [1, 1, 1, 0, 0, 0, 0]])
    adjacency = torch.randn(2, 5, 5)
    left = trunk(atoms, residues, atom_mask, residue_mask, adjacency)

    perturbed_atoms = atoms.clone()
    perturbed_residues = residues.clone()
    perturbed_adjacency = adjacency.clone()
    perturbed_atoms[atom_mask == 0] = 1.0e4
    perturbed_residues[residue_mask == 0] = -1.0e4
    invalid_edges = ~(
        atom_mask[:, :, None].bool() & atom_mask[:, None, :].bool())
    perturbed_adjacency[invalid_edges] = 1.0e4
    right = trunk(
        perturbed_atoms, perturbed_residues, atom_mask, residue_mask,
        perturbed_adjacency)

    assert torch.allclose(left.endpoint, right.endpoint, atol=1e-6, rtol=1e-6)
    assert torch.allclose(left.section, right.section, atol=1e-6, rtol=1e-6)
    assert torch.allclose(
        left.mechanism_slots, right.mechanism_slots, atol=1e-6, rtol=1e-6)


def test_retained_mechanism_slots_exactly_reconstruct_original_readouts():
    torch.manual_seed(74)
    trunk = BipartitePairSectionFormer(
        8, 5, pair_dim=16, blocks=1, latent_count=4, heads=2)
    encoded = trunk(
        torch.randn(2, 5, 8), torch.randn(2, 7, 8),
        torch.ones(2, 5), torch.ones(2, 7))
    endpoint, section = trunk.latent.project_slots(encoded.mechanism_slots)
    assert torch.equal(endpoint, encoded.endpoint)
    assert torch.equal(section, encoded.section)


def test_hypersar_pair_adapter_changes_only_adaptive_path_and_has_gradients():
    torch.manual_seed(75)
    trunk = BipartitePairSectionFormer(
        8, 6, pair_dim=16, blocks=2, latent_count=4, heads=2,
        task_dim=5, adapter_rank=2, adaptive_blocks=1)
    atoms = torch.randn(3, 5, 8, requires_grad=True)
    residues = torch.randn(3, 7, 8, requires_grad=True)
    atom_mask = torch.ones(3, 5)
    residue_mask = torch.ones(3, 7)
    code = torch.randn(3, 5, requires_grad=True)
    base = trunk(atoms, residues, atom_mask, residue_mask)
    adapted = trunk(
        atoms, residues, atom_mask, residue_mask, task_code=code)
    zero = trunk(
        atoms, residues, atom_mask, residue_mask,
        task_code=torch.zeros_like(code))
    assert torch.allclose(base.endpoint, zero.endpoint, atol=1e-7, rtol=1e-7)
    assert not torch.allclose(base.endpoint, adapted.endpoint)
    adapted.endpoint.square().mean().backward()
    assert code.grad is not None
    assert trunk.blocks[-1].adapter_gate.weight.grad is not None


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
