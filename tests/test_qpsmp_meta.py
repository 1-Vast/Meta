import torch
import pytest

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.qpsmp_meta import (
    AmortizedTargetConditioner, QPSMPBioModel, SiameseRelativeConditioner,
)


def test_conditioner_is_permutation_invariant_zero_preserving_and_not_support_span_bound():
    torch.manual_seed(31)
    conditioner = AmortizedTargetConditioner(
        6, 6, task_dim=5, hidden_dim=12, blocks=1, dtype=torch.float64)
    interaction = torch.randn(4, 6, dtype=torch.float64, requires_grad=True)
    ligand = torch.randn(4, 6, dtype=torch.float64, requires_grad=True)
    residual = torch.tensor([-1.5, -0.5, 0.5, 1.5], dtype=torch.float64)
    code = conditioner(interaction, ligand, residual)
    order = torch.tensor([2, 0, 3, 1])
    permuted = conditioner(interaction[order], ligand[order], residual[order])
    labels_only_permuted = conditioner(interaction, ligand, residual[order])
    assert code.shape == (5,)
    assert torch.allclose(code, permuted, atol=1e-10, rtol=1e-10)
    assert not torch.allclose(code, labels_only_permuted)
    assert torch.count_nonzero(conditioner(
        interaction, ligand, torch.zeros_like(residual))) == 0
    code.square().mean().backward()
    assert interaction.grad is not None and ligand.grad is not None
    assert conditioner.output[1].weight.grad is not None


def test_siamese_relative_conditioner_is_query_specific_and_zero_preserving():
    torch.manual_seed(34)
    conditioner = SiameseRelativeConditioner(
        6, 6, task_dim=5, hidden_dim=12, dtype=torch.float64)
    reference = torch.randn(5, dtype=torch.float64)
    interaction = torch.randn(4, 6, dtype=torch.float64)
    ligand = torch.randn(4, 6, dtype=torch.float64)
    code = conditioner(reference, interaction, ligand)
    assert code.shape == (4, 5)
    assert not torch.allclose(code[0], code[1])
    assert torch.count_nonzero(conditioner(
        torch.zeros_like(reference), interaction, ligand)) == 0


def test_relative_residual_match_is_jointly_permutation_invariant_and_k1_zero():
    torch.manual_seed(36)
    meta = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2).meta
    support_endpoint = torch.randn(4, 16)
    support_ligand = torch.randn(4, 16)
    residual = torch.tensor([-1.5, -0.5, 0.5, 1.5])
    query_endpoint = torch.randn(3, 16)
    query_ligand = torch.randn(3, 16)
    value = meta.relative_residual_match(
        support_endpoint, support_ligand, residual,
        query_endpoint, query_ligand)
    order = torch.tensor([2, 0, 3, 1])
    permuted = meta.relative_residual_match(
        support_endpoint[order], support_ligand[order], residual[order],
        query_endpoint, query_ligand)
    assert torch.allclose(value, permuted, atol=1e-6, rtol=1e-6)
    reliability = meta.support_match_reliability(
        support_endpoint, support_ligand, residual)
    permuted_reliability = meta.support_match_reliability(
        support_endpoint[order], support_ligand[order], residual[order])
    assert 0 <= reliability <= 1
    assert torch.allclose(reliability, permuted_reliability)
    assert meta.support_match_reliability(
        support_endpoint[:1], support_ligand[:1], torch.zeros(1)) == 0
    assert meta.support_match_loss(
        support_endpoint, support_ligand, residual) > 0
    assert meta.support_match_loss(
        support_endpoint[:1], support_ligand[:1], torch.zeros(1)) == 0
    assert torch.count_nonzero(meta.relative_residual_match(
        support_endpoint[:1], support_ligand[:1], torch.zeros(1),
        query_endpoint, query_ligand)) == 0


def test_active_hypersar_batch_forward_controls_and_gradient():
    torch.manual_seed(32)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=2, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        support_hidden_dim=16, support_blocks=1, adapter_rank=2,
        adaptive_blocks=1)
    batch, k, q, atoms, residues = 2, 3, 4, 5, 6
    pooled = torch.randn(batch, 8)
    tokens = torch.randn(batch, residues, 8)
    protein_mask = torch.ones(batch, residues)
    support_atoms = torch.randn(batch, k, atoms, ATOM_FEAT_DIM)
    query_atoms = torch.randn(batch, q, atoms, ATOM_FEAT_DIM)
    support_bonds = torch.randn(batch, k, atoms, atoms, BOND_FEAT_DIM)
    query_bonds = torch.randn(batch, q, atoms, atoms, BOND_FEAT_DIM)
    support_mask = torch.ones(batch, k, atoms)
    query_mask = torch.ones(batch, q, atoms)
    labels = torch.randn(batch, k)
    output = model(pooled, tokens, protein_mask, support_atoms, support_bonds,
                   support_mask, labels, query_atoms, query_bonds, query_mask)
    zero = model(pooled, tokens, protein_mask, support_atoms, support_bonds,
                 support_mask, labels, query_atoms, query_bonds, query_mask,
                 adapt=False)
    replay = model(
        pooled, tokens, protein_mask, support_atoms, support_bonds,
        support_mask, labels, query_atoms, query_bonds, query_mask,
        task_state_override=output.task_state)
    assert output.prediction.shape == (batch, q)
    assert torch.allclose(output.prediction, replay.prediction)
    assert torch.count_nonzero(zero.task_state) == 0
    assert torch.count_nonzero(zero.sar_adaptation) == 0
    output.prediction.square().mean().backward()
    assert model.meta.conditioner.output[1].weight.grad is not None
    adaptive = model.pair_section.blocks[-1]
    assert adaptive.adapter_gate.weight.grad is not None
    assert adaptive.adapter_down.weight.grad is not None
    assert adaptive.adapter_up.weight.grad is not None
    assert model.ligand_encoder.inp.weight.grad is not None


def test_single_support_has_level_adaptation_but_exact_zero_hypersar_shape():
    torch.manual_seed(33)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2,
        support_hidden_dim=16, support_blocks=1, adapter_rank=2,
        adaptive_blocks=1)
    pooled = torch.randn(8)
    tokens = torch.randn(6, 8)
    protein_mask = torch.ones(6)
    support_atoms = torch.randn(1, 5, ATOM_FEAT_DIM)
    query_atoms = torch.randn(3, 5, ATOM_FEAT_DIM)
    support_bonds = torch.randn(1, 5, 5, BOND_FEAT_DIM)
    query_bonds = torch.randn(3, 5, 5, BOND_FEAT_DIM)
    support_mask = torch.ones(1, 5)
    query_mask = torch.ones(3, 5)
    output = model(
        pooled, tokens, protein_mask, support_atoms, support_bonds,
        support_mask, torch.tensor([2.0]), query_atoms, query_bonds, query_mask)
    assert torch.count_nonzero(output.task_state) == 0
    assert torch.count_nonzero(output.sar_adaptation) == 0


def test_task_state_override_rejects_query_level_or_ambiguous_shapes():
    torch.manual_seed(35)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2,
        support_hidden_dim=16, support_blocks=1, adapter_rank=2,
        adaptive_blocks=1)
    batch, k, q, atoms, residues = 2, 2, 3, 4, 5
    args = (
        torch.randn(batch, 8), torch.randn(batch, residues, 8),
        torch.ones(batch, residues),
        torch.randn(batch, k, atoms, ATOM_FEAT_DIM),
        torch.randn(batch, k, atoms, atoms, BOND_FEAT_DIM),
        torch.ones(batch, k, atoms), torch.randn(batch, k),
        torch.randn(batch, q, atoms, ATOM_FEAT_DIM),
        torch.randn(batch, q, atoms, atoms, BOND_FEAT_DIM),
        torch.ones(batch, q, atoms),
    )
    valid = model(*args).task_state
    assert model(*args, task_state_override=valid).prediction.shape == (batch, q)
    with pytest.raises(ValueError, match="target-level reference code"):
        model(*args, task_state_override=torch.randn(batch, q, 8))
    with pytest.raises(ValueError, match="target-level reference code"):
        model(*args, task_state_override=torch.randn(q, 8))


def test_level_gate_zero_is_exact_support_mean_baseline():
    torch.manual_seed(37)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2)
    with torch.no_grad():
        model.meta.level_gate[-2].weight.zero_()
        model.meta.level_gate[-2].bias.fill_(-100.0)
    support_ligand = torch.randn(3, 16)
    support_endpoint = torch.randn(3, 16)
    query_ligand = torch.randn(4, 16)
    query_endpoint = torch.randn(4, 16)
    labels = torch.tensor([1.0, 2.0, 6.0])
    state = model.meta.infer(
        torch.randn(16), support_ligand, support_endpoint, labels,
        query_ligand, query_endpoint)
    baseline = state.zero_shot + state.level_adjustment
    assert torch.allclose(baseline, labels.mean().expand_as(baseline))
