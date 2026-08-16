import torch

from research.meta_fewshot.qpsmp_bpsf_v2 import (
    RelevanceWeightedBPSF, SharedPairLatent,
)
from model.bpsf import BipartitePairSectionFormer


def test_shared_latent_and_relevance_variants_are_trainable_and_masked():
    torch.manual_seed(81)
    base = BipartitePairSectionFormer(
        8, 6, pair_dim=16, blocks=1, latent_count=4, heads=2)
    shared = SharedPairLatent(base.latent)
    relevance = RelevanceWeightedBPSF(
        8, 6, pair_dim=16, blocks=1, latent_count=4, heads=2)
    relevance.latent = shared
    atoms = torch.randn(3, 5, 8, requires_grad=True)
    residues = torch.randn(3, 7, 8, requires_grad=True)
    atom_mask = torch.tensor([[1, 1, 1, 0, 0], [1] * 5, [1, 1, 0, 0, 0]])
    residue_mask = torch.tensor([[1, 1, 1, 0, 0, 0, 0], [1] * 7,
                                 [1, 1, 1, 1, 0, 0, 0]])
    encoded = relevance(atoms, residues, atom_mask, residue_mask)
    loss = encoded.endpoint.square().mean() + encoded.section.square().mean()
    loss.backward()
    assert encoded.endpoint.shape == (3, 8)
    assert encoded.section.shape == (3, 6)
    assert relevance.relevance[-1].weight.grad is not None
    assert shared.section.weight.grad is not None

