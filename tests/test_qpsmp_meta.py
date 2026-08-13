import torch
import pytest

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.qpsmp_meta import MechanismEvidenceMetaTransformer, QPSMPBioModel
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
        torch.ones(batch, q, atoms),
    )


def _model():
    return QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        support_blocks=1)


def test_mechanism_evidence_is_support_permutation_invariant_and_k1_active():
    torch.manual_seed(31)
    learner = MechanismEvidenceMetaTransformer(
        6, 8, 4, blocks=1, dtype=torch.float64)
    support = torch.randn(2, 3, 4, 6, dtype=torch.float64)
    residual = torch.tensor([[1.0, -0.5, 0.25], [-1.0, 0.5, 0.75]],
                            dtype=torch.float64, requires_grad=True)
    query = torch.randn(2, 5, 4, 6, dtype=torch.float64)
    prompts, gate, evidence, auxiliary = learner(support, residual, query)
    order = torch.tensor([2, 0, 1])
    permuted = learner(support[:, order], residual[:, order], query)
    assert prompts.shape == (2, 4, 8)
    assert gate.shape == (2, 5, 4)
    assert evidence.shape == (2, 3, 4)
    assert torch.allclose(prompts, permuted[0], atol=1e-10, rtol=1e-10)
    assert torch.allclose(gate, permuted[1], atol=1e-10, rtol=1e-10)
    one = learner(support[:, :1], residual[:, :1], query)
    assert torch.count_nonzero(one[2]) > 0
    (gate.square().mean() + auxiliary).backward()
    assert residual.grad is not None
    assert learner.gate[-1].weight.grad is not None


def test_difference_transport_responds_to_label_binding_and_zero_is_exact():
    torch.manual_seed(32)
    learner = MechanismEvidenceMetaTransformer(
        6, 8, 4, blocks=1, dtype=torch.float64)
    support = torch.randn(2, 3, 4, 6, dtype=torch.float64)
    query = torch.randn(2, 5, 4, 6, dtype=torch.float64)
    residual = torch.tensor([[1.0, -0.5, 0.25], [-1.0, 0.5, 0.75]],
                            dtype=torch.float64)
    real = learner(support, residual, query)
    labels_only = learner(support, residual.roll(1, dims=1), query)
    assert not torch.allclose(real[1], labels_only[1])
    zero = learner(support, residual, query, adapt=False)
    assert torch.count_nonzero(zero[0]) == 0
    assert torch.count_nonzero(zero[1]) == 0


def test_active_meta_model_batch_forward_zero_k1_and_gradients():
    torch.manual_seed(33)
    model, args = _model(), _inputs(k=1)
    output = model(*args)
    zero = model(*args, adapt=False)
    assert output.prediction.shape == (2, 4)
    assert output.task_state.shape == (2, 4, 8)
    assert torch.count_nonzero(output.task_state) > 0
    assert torch.count_nonzero(output.support_evidence) > 0
    assert torch.count_nonzero(output.sar_adaptation) > 0
    assert torch.equal(zero.prediction, zero.zero_shot)
    assert torch.count_nonzero(zero.task_state) == 0
    assert torch.count_nonzero(zero.sar_adaptation) == 0
    output.prediction.square().mean().backward()
    assert model.meta.mechanism.gate[-1].weight.grad is not None
    assert model.meta.mechanism.sensitivity[-1].weight.grad is not None
    assert model.meta.reference_delta[-1].weight.grad is not None
    assert model.ligand_encoder.inp.weight.grad is not None


def test_task_state_override_requires_complete_mechanism_prompt_shape():
    torch.manual_seed(34)
    model, args = _model(), _inputs()
    valid = model(*args).task_state
    replay = model(*args, task_state_override=valid)
    assert replay.prediction.shape == (2, 4)
    changed = model(*args, task_state_override=torch.zeros_like(valid))
    assert not torch.allclose(replay.prediction, changed.prediction)
    with pytest.raises(ValueError, match="target mechanism prompts"):
        model(*args, task_state_override=torch.randn(2, 4, 4, 8))
    with pytest.raises(ValueError, match="target mechanism prompts"):
        model(*args, task_state_override=torch.randn(4, 8))


def test_level_gate_zero_is_exact_support_mean_baseline():
    torch.manual_seed(35)
    model = _model()
    with torch.no_grad():
        model.meta.level_gate[-2].weight.zero_()
        model.meta.level_gate[-2].bias.fill_(-100.0)
        model.meta.mechanism.gate[-1].weight.zero_()
    output = model(*_inputs(batch=2, k=3, q=4))
    # Reuse the actual support labels from the generated input tuple.
    args = _inputs(batch=2, k=3, q=4)
    output = model(*args)
    expected = args[6].mean(-1, keepdim=True).expand_as(output.level_baseline)
    assert torch.allclose(output.level_baseline, expected)


def test_unbatched_forward_preserves_public_shapes():
    torch.manual_seed(36)
    model = _model()
    args = _inputs(batch=1, k=2, q=3)
    args = tuple(value.squeeze(0) for value in args)
    output = model(*args)
    assert output.prediction.shape == (3,)
    assert output.task_state.shape == (4, 8)
    assert output.support_residual_quotient.shape == (2,)
    replay = model(*args, task_state_override=output.task_state)
    assert replay.prediction.shape == (3,)


def test_common_frame_cartesian_slots_are_invariant_and_deeply_connected():
    torch.manual_seed(37)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        support_blocks=1, use_cartesian=True)
    args = _inputs(batch=1, k=2, q=2, atoms=4, residues=3)
    coordinates = torch.randn(1, 4, 7, 3)
    # Packed complete graph for four protein--ligand samples.
    edges = []
    for sample in range(4):
        offset = sample * 7
        edges.extend((offset + source, offset + target)
                     for source in range(7) for target in range(7)
                     if source != target)
    edge_index = torch.tensor(edges).T
    left = model(*args, geometry_coordinates=coordinates,
                 geometry_edge_index=edge_index,
                 geometry_common_frame=torch.ones(1, 4, dtype=torch.bool))
    rotation, _ = torch.linalg.qr(torch.randn(3, 3))
    right = model(*args, geometry_coordinates=coordinates @ rotation.T + 5.0,
                  geometry_edge_index=edge_index,
                  geometry_common_frame=torch.ones(1, 4, dtype=torch.bool))
    assert torch.allclose(left.prediction, right.prediction, atol=1e-5, rtol=1e-5)
    left.prediction.square().mean().backward()
    assert model.geometry_scale.grad is not None
    assert model.cartesian_encoder.layers[0].coefficient.weight.grad is not None


def test_training_forward_wrapper_transports_optional_cartesian_episode_fields():
    torch.manual_seed(38)
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, pair_chunk_size=8,
        support_blocks=1, use_cartesian=True)
    args = _inputs(batch=1, k=2, q=2, atoms=4, residues=3)
    args = tuple(value.squeeze(0) for value in args)
    coordinates = torch.randn(4, 7, 3)
    edges = _edges = []
    for sample in range(4):
        offset = sample * 7
        _edges.extend((offset + source, offset + target)
                      for source in range(7) for target in range(7)
                      if source != target)
    edge_index = torch.tensor(_edges).T
    episode = EpisodeBatch(
        EpisodeSpec("meta_train", "component", "target", (0, 1), (2, 3), "donor"),
        args[0], args[1], args[2], args[3], args[4], args[5], args[6],
        args[7], args[8], args[9], torch.randn(2),
        support_coordinates=coordinates[:2], query_coordinates=coordinates[2:],
        geometry_edge_index=edge_index,
        geometry_available=torch.ones(4, dtype=torch.bool),
        geometry_common_frame=torch.ones(4, dtype=torch.bool))
    output = episode_forward(model, episode)
    assert output.prediction.shape == (2,)
    output.prediction.square().mean().backward()
    assert model.geometry_scale.grad is not None


def test_cartesian_interaction_rejects_undeclared_or_independent_frames():
    model = QPSMPBioModel(
        8, hidden_dim=16, task_dim=8, ligand_layers=1, pair_dim=8,
        pair_blocks=1, pair_latents=4, pair_heads=2, support_blocks=1,
        use_cartesian=True)
    args = _inputs(batch=1, k=1, q=1, atoms=3, residues=2)
    coordinates = torch.randn(1, 2, 5, 3)
    edges = torch.tensor([[0, 1, 5, 6], [1, 0, 6, 5]])
    with pytest.raises(ValueError, match="explicit common-frame"):
        model(*args, geometry_coordinates=coordinates,
              geometry_edge_index=edges)
    with pytest.raises(ValueError, match="must share a common frame"):
        model(*args, geometry_coordinates=coordinates,
              geometry_edge_index=edges,
              geometry_available=torch.ones(1, 2, dtype=torch.bool),
              geometry_common_frame=torch.tensor([[True, False]]))
