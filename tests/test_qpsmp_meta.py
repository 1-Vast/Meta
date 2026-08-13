import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.qpsmp_meta import QPSMPBioModel, QPSMPMetaLearner


DTYPE = torch.float64


def meta_inputs(k=3, q=4):
    torch.manual_seed(31)
    protein = torch.randn(6, dtype=DTYPE)
    support_ligand = torch.randn(k, 6, dtype=DTYPE)
    query_ligand = torch.randn(q, 6, dtype=DTYPE)
    support_endpoint = torch.randn(k, 6, dtype=DTYPE)
    query_endpoint = torch.randn(q, 6, dtype=DTYPE)
    support_section = torch.randn(k, 3, dtype=DTYPE)
    query_section = torch.randn(q, 3, dtype=DTYPE)
    return (protein, support_ligand, support_endpoint, support_section,
            query_ligand, query_endpoint, query_section)


def run_meta(model, labels, values=None):
    values = values or meta_inputs(len(labels))
    protein, sl, se, ss, ql, qe, qs = values
    return model(protein, sl, se, ss, labels, ql, qe, qs)


def test_meta_section_is_trainable_permutation_invariant_and_row_span():
    model = QPSMPMetaLearner(6, 3, support_hidden_dim=12,
                             support_blocks=1, dtype=DTYPE)
    labels = torch.randn(3, dtype=DTYPE)
    values = meta_inputs(3)
    left = run_meta(model, labels, values)
    order = torch.tensor([2, 0, 1])
    protein, sl, se, ss, ql, qe, qs = values
    right = model(protein, sl[order], se[order], ss[order], labels[order], ql, qe, qs)
    assert torch.allclose(left.prediction, right.prediction, atol=1e-10)
    centered = ss - ss.mean(0, keepdim=True)
    projection = torch.linalg.pinv(centered) @ centered
    assert torch.allclose(left.task_state, projection @ left.task_state, atol=1e-8)
    left.prediction.square().mean().backward()
    assert model.section_operator.weight.weight.grad is not None
    assert model.cross_head[1].weight.grad is not None


def test_constant_residual_and_single_support_have_exact_zero_sar():
    model = QPSMPMetaLearner(6, 3, support_hidden_dim=12,
                             support_blocks=1, dtype=DTYPE)
    values = meta_inputs(4)
    protein, sl, se, ss, ql, qe, qs = values
    with torch.no_grad():
        support_zero = (model.ligand_baseline(sl).squeeze(-1)
                        + model.protein_level(protein).squeeze(-1)
                        + model.cross_head(se).squeeze(-1))
    output = run_meta(model, support_zero + 2.0, values)
    assert torch.count_nonzero(output.task_state) == 0
    assert torch.count_nonzero(output.sar_adaptation) == 0

    one = meta_inputs(1)
    output = run_meta(model, torch.tensor([1.5], dtype=DTYPE), one)
    assert torch.count_nonzero(output.task_state) == 0
    assert torch.count_nonzero(output.sar_adaptation) == 0


def test_active_bio_model_batch_forward_and_gradient():
    torch.manual_seed(32)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=4, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        support_hidden_dim=16, support_blocks=1, dtype=torch.float32)
    batch, k, q, atoms, residues = 2, 2, 3, 5, 6
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
    assert output.prediction.shape == (batch, q)
    output.prediction.square().mean().backward()
    assert model.pair_section.latent.section.output.weight.grad is not None
    assert model.meta.section_operator.weight.weight.grad is not None
    assert model.ligand_encoder.inp.weight.grad is not None
