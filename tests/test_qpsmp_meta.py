import pytest
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.qpsmp_meta import EvidenceLockedMetaTransport, QPSMPBioModel
from scripts.qpsmp_data import EpisodeBatch, EpisodeSpec
from scripts.train_qpsmp import forward as episode_forward


def _inputs(batch=2, k=3, q=4, atoms=5, residues=6):
    return (
        torch.randn(batch, 8), torch.randn(batch, residues, 8),
        torch.ones(batch, residues),
        torch.randn(batch, k, atoms, ATOM_FEAT_DIM),
        torch.randn(batch, k, atoms, atoms, BOND_FEAT_DIM),
        torch.ones(batch, k, atoms), torch.randn(batch, k),
        torch.randn(batch, q, atoms, ATOM_FEAT_DIM),
        torch.randn(batch, q, atoms, atoms, BOND_FEAT_DIM),
        torch.ones(batch, q, atoms))


def _model():
    return QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        support_blocks=1)


def test_elmt_value_is_label_bound_and_linear():
    torch.manual_seed(31)
    residual = torch.randn(2, 3, dtype=torch.float64)
    primitive = torch.randn(2, 3, 4, dtype=torch.float64)
    value = EvidenceLockedMetaTransport.label_value(residual, primitive)
    assert torch.equal(value, residual.unsqueeze(-1) * primitive)


def test_elmt_is_support_permutation_invariant_k1_active_and_label_bound():
    torch.manual_seed(32)
    term = EvidenceLockedMetaTransport(6, 4, 8, dtype=torch.float64)
    protein = torch.randn(2, 6, dtype=torch.float64)
    support_ligand = torch.randn(2, 3, 6, dtype=torch.float64)
    support_phi = torch.randn(2, 3, 4, dtype=torch.float64)
    residual = torch.tensor([[1., -.5, .25], [-1., .5, .75]], dtype=torch.float64)
    query_ligand = torch.randn(2, 5, 6, dtype=torch.float64)
    query_phi = torch.randn(2, 5, 4, dtype=torch.float64)
    real = term(protein, support_ligand, support_phi, residual,
                query_ligand, query_phi)
    order = torch.tensor([2, 0, 1])
    permuted = term(protein, support_ligand[:, order], support_phi[:, order],
                    residual[:, order], query_ligand, query_phi)
    wrong_binding = term(protein, support_ligand, support_phi,
                         residual.roll(1, 1), query_ligand, query_phi)
    assert real[0].shape == (2, 4, 8)
    assert real[1].shape == (2, 5, 4)
    assert torch.allclose(real[0], permuted[0], atol=1e-10, rtol=1e-10)
    assert torch.allclose(real[1], permuted[1], atol=1e-10, rtol=1e-10)
    assert not torch.allclose(real[1], wrong_binding[1])
    one = term(protein, support_ligand[:, :1], support_phi[:, :1],
               residual[:, :1], query_ligand, query_phi)
    assert torch.count_nonzero(one[3]) > 0
    assert torch.count_nonzero(one[1]) > 0


def test_elmt_coefficients_are_linear_in_values_and_zero_locked():
    torch.manual_seed(321)
    term = EvidenceLockedMetaTransport(6, 4, 8, dtype=torch.float64)
    protein = torch.randn(1, 6, dtype=torch.float64)
    support_ligand = torch.randn(1, 1, 6, dtype=torch.float64)
    support_phi = torch.randn(1, 1, 4, dtype=torch.float64)
    residual = torch.ones(1, 1, dtype=torch.float64)
    query_ligand = torch.randn(1, 2, 6, dtype=torch.float64)
    query_phi = torch.randn(1, 2, 4, dtype=torch.float64)
    left = term(protein, support_ligand, support_phi, residual,
                query_ligand, query_phi)
    right = term(protein, support_ligand, support_phi, 0.25 * residual,
                 query_ligand, query_phi)
    assert torch.allclose(right[1], 0.25 * left[1], atol=1e-12, rtol=1e-12)
    zero = term(protein, support_ligand, support_phi,
                torch.zeros_like(residual), query_ligand, query_phi)
    assert torch.count_nonzero(zero[3]) == 0
    assert torch.count_nonzero(zero[1]) == 0


def test_active_model_k0_k1_and_gradients():
    torch.manual_seed(33)
    model, args = _model(), _inputs(k=1)
    output = model(*args)
    zero = model(*args, adapt=False)
    assert output.prediction.shape == (2, 4)
    assert output.task_state.shape == (2, 4, 8)
    assert output.query_basis.shape == (2, 4, 4)
    assert torch.count_nonzero(output.support_evidence) > 0
    assert torch.count_nonzero(output.sar_adaptation) > 0
    assert torch.equal(zero.prediction, zero.zero_shot)
    assert torch.count_nonzero(zero.task_state) == 0
    assert torch.count_nonzero(zero.sar_adaptation) == 0
    (output.prediction.square().mean()
     + .01 * output.support_match_loss).backward()
    assert model.meta.term.direction[-1].weight.grad is not None
    assert model.meta.term.score[-1].weight.grad is not None
    assert model.residue_query.weight.grad is not None
    assert model.pair_section.latent.interaction.response_weight.grad is not None
    assert model.pair_section.atom_primitive.weight.grad is not None


def test_elmt_rejects_task_state_transplant():
    model, args = _model(), _inputs()
    valid = model(*args).task_state
    with pytest.raises(ValueError, match="does not permit"):
        model(*args, task_state_override=valid)


def test_level_gate_zero_is_exact_support_mean_baseline():
    model, args = _model(), _inputs(batch=2, k=3, q=4)
    with torch.no_grad():
        model.meta.level_gate[-2].weight.zero_()
        model.meta.level_gate[-2].bias.fill_(100.)
    output = model(*args)
    # Level correction is a target-wise scalar and preserves zero-shot order.
    delta_level = output.level_baseline - output.zero_shot
    assert torch.allclose(delta_level, delta_level[:, :1].expand_as(delta_level))
    assert torch.allclose(
        output.level_baseline[:, 1:] - output.level_baseline[:, :-1],
        output.zero_shot[:, 1:] - output.zero_shot[:, :-1])


def test_unbatched_forward_preserves_public_shapes():
    model = _model()
    args = tuple(value.squeeze(0) for value in _inputs(batch=1, k=2, q=3))
    output = model(*args)
    assert output.prediction.shape == (3,)
    assert output.task_state.shape == (4, 8)
    assert output.support_residual_quotient.shape == (2,)


def _edges(samples, nodes):
    return torch.tensor([
        (sample * nodes + source, sample * nodes + target)
        for sample in range(samples) for source in range(nodes)
        for target in range(nodes) if source != target]).T


def test_common_frame_cartesian_bias_is_invariant_and_connected():
    torch.manual_seed(37)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        use_cartesian=True)
    args = _inputs(batch=1, k=2, q=2, atoms=4, residues=3)
    coordinates = torch.randn(1, 4, 7, 3)
    edges = _edges(4, 7)
    available = torch.ones(1, 4, dtype=torch.bool)
    left = model(*args, geometry_coordinates=coordinates,
                 geometry_edge_index=edges, geometry_common_frame=available)
    rotation, _ = torch.linalg.qr(torch.randn(3, 3))
    right = model(*args, geometry_coordinates=coordinates @ rotation.T + 5.,
                  geometry_edge_index=edges, geometry_common_frame=available)
    assert torch.allclose(left.prediction, right.prediction, atol=1e-5, rtol=1e-5)
    left.prediction.square().mean().backward()
    assert model.geometry_scale.grad is not None
    assert model.cartesian_encoder.layers[0].coefficient.weight.grad is not None


def test_training_wrapper_transports_cartesian_fields_and_rejects_frames():
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, use_cartesian=True)
    args = tuple(value.squeeze(0) for value in _inputs(
        batch=1, k=2, q=2, atoms=4, residues=3))
    coordinates = torch.randn(4, 7, 3)
    episode = EpisodeBatch(
        EpisodeSpec("meta_train", "c", "t", (0, 1), (2, 3), "d"),
        *args, torch.randn(2), support_coordinates=coordinates[:2],
        query_coordinates=coordinates[2:], geometry_edge_index=_edges(4, 7),
        geometry_available=torch.ones(4, dtype=torch.bool),
        geometry_common_frame=torch.ones(4, dtype=torch.bool))
    assert episode_forward(model, episode).prediction.shape == (2,)
    batch_args = _inputs(batch=1, k=1, q=1, atoms=3, residues=2)
    coords = torch.randn(1, 2, 5, 3)
    edges = torch.tensor([[0, 1, 5, 6], [1, 0, 6, 5]])
    with pytest.raises(ValueError, match="explicit common-frame"):
        model(*batch_args, geometry_coordinates=coords,
              geometry_edge_index=edges)
