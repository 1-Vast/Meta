import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.qpsmp_meta import QPSMPBioModel, QPSMPMetaLearner


DTYPE = torch.float64


def episode(model, support_y=None):
    torch.manual_seed(31)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.tensor([True, True, True, True, False])
    support = torch.randn(3, 6, dtype=DTYPE)
    query = torch.randn(4, 6, dtype=DTYPE)
    labels = torch.randn(3, dtype=DTYPE) if support_y is None else support_y
    return model(protein, tokens, mask, support, labels, query)


def test_query_loss_trains_localizer_scalar_basis_and_neural_adapter():
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    output = episode(model)

    output.prediction.square().mean().backward()

    parameters = (
        model.localizer.query.weight,
        model.zero_shot_head.weight,
        model.section_head.weight,
        model.adapter.value.weight,
        model.adapter.update[0].weight,
    )
    assert all(parameter.grad is not None for parameter in parameters)
    assert all(torch.isfinite(parameter.grad).all() for parameter in parameters)


def test_support_permutation_does_not_change_prediction_or_state():
    torch.manual_seed(32)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    support = torch.randn(4, 6, dtype=DTYPE)
    labels = torch.randn(4, dtype=DTYPE)
    query = torch.randn(3, 6, dtype=DTYPE)
    order = torch.tensor([2, 0, 3, 1])

    left = model(protein, tokens, mask, support, labels, query)
    right = model(protein, tokens, mask, support[order], labels[order], query)

    assert torch.allclose(left.task_state, right.task_state, atol=1e-12)
    assert torch.allclose(left.prediction, right.prediction, atol=1e-12)


def test_quotient_null_forces_exact_zero_neural_update():
    torch.manual_seed(33)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    support = torch.randn(4, 6, dtype=DTYPE)
    query = torch.randn(3, 6, dtype=DTYPE)
    with torch.no_grad():
        interaction = model.interaction_features(protein, tokens, mask, support)
        add, _, zero, _ = model.scalar_components(protein, support, interaction)
        labels = add + zero + 2.5

    output = model(protein, tokens, mask, support, labels, query)

    assert torch.allclose(output.support_residual_quotient, torch.zeros(4, dtype=DTYPE))
    assert torch.allclose(output.task_state, torch.zeros(2, dtype=DTYPE), atol=1e-12)


def test_one_support_label_can_change_level_but_not_sar_state():
    torch.manual_seed(34)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    support = torch.randn(1, 6, dtype=DTYPE)
    query = torch.randn(3, 6, dtype=DTYPE)

    output = model(protein, tokens, mask, support, torch.tensor([4.0], dtype=DTYPE), query)

    assert torch.allclose(output.task_state, torch.zeros(2, dtype=DTYPE))
    assert not torch.allclose(output.level_shift, torch.zeros((), dtype=DTYPE))
    assert torch.allclose(output.prediction, output.zero_shot + output.level_adjustment)
    assert 0.0 < output.level_shrinkage < 1.0


def test_zero_support_branch_ignores_support_labels():
    torch.manual_seed(35)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    support = torch.randn(3, 6, dtype=DTYPE)
    query = torch.randn(2, 6, dtype=DTYPE)

    left = model(protein, tokens, mask, support, torch.randn(3, dtype=DTYPE), query, adapt=False)
    right = model(protein, tokens, mask, support, torch.randn(3, dtype=DTYPE), query, adapt=False)

    assert torch.allclose(left.prediction, right.prediction, atol=1e-12)
    assert torch.allclose(left.prediction, left.zero_shot, atol=1e-12)


def test_delta_and_rectangle_use_retained_endpoint_predictions():
    prediction = torch.tensor([1.0, 3.0, -2.0], dtype=DTYPE)
    left = torch.tensor([0])
    right = torch.tensor([1])
    delta = QPSMPMetaLearner.delta(prediction, left, right)

    assert torch.allclose(delta, -QPSMPMetaLearner.delta(prediction, right, left))
    assert torch.allclose(
        QPSMPMetaLearner.rectangle(delta, torch.tensor([0.5], dtype=DTYPE)),
        torch.tensor([1.5], dtype=DTYPE))


def test_interaction_heads_do_not_receive_uncrossed_main_features():
    model = QPSMPMetaLearner(3, 2, dtype=DTYPE)
    with torch.no_grad():
        model.localizer.key.weight.zero_()
        model.localizer.query.weight.zero_()
        model.localizer.value.weight.zero_()
    protein = torch.randn(3, dtype=DTYPE)
    tokens = torch.randn(4, 3, dtype=DTYPE)
    mask = torch.ones(4, dtype=torch.bool)
    ligand = torch.randn(2, 3, dtype=DTYPE)

    crossed = model.interaction_features(protein, tokens, mask, ligand)

    assert torch.allclose(crossed, torch.zeros_like(crossed), atol=1e-12)


def test_cached_float16_protein_bank_is_accepted_by_float32_model():
    model = QPSMPMetaLearner(6, 2, dtype=torch.float32)
    torch.manual_seed(36)
    protein = torch.randn(6, dtype=torch.float16)
    tokens = torch.randn(5, 6, dtype=torch.float16)
    mask = torch.ones(5, dtype=torch.uint8)
    support = torch.randn(2, 6, dtype=torch.float32)
    query = torch.randn(3, 6, dtype=torch.float32)

    output = model(protein, tokens, mask, support, torch.randn(2), query)

    assert output.prediction.dtype == torch.float32
    assert torch.isfinite(output.prediction).all()


def test_query_partition_does_not_change_endpoint_function():
    torch.manual_seed(37)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    support = torch.randn(3, 6, dtype=DTYPE)
    labels = torch.randn(3, dtype=DTYPE)
    query = torch.randn(5, 6, dtype=DTYPE)

    whole = model(protein, tokens, mask, support, labels, query).prediction
    parts = torch.cat([
        model(protein, tokens, mask, support, labels, query[:2]).prediction,
        model(protein, tokens, mask, support, labels, query[2:]).prediction,
    ])

    assert torch.allclose(whole, parts, atol=1e-12)


def test_output_decomposes_into_additive_cross_and_adaptation_channels():
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    output = episode(model)

    assert torch.allclose(
        output.prediction,
        output.zero_shot + output.level_adjustment + output.sar_adaptation,
        atol=1e-12)


def test_zero_sar_state_does_not_gate_zero_shot_backbone_gradients():
    torch.manual_seed(40)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    with torch.no_grad():
        model.adapter.value.weight.zero_()
    output = episode(model)

    output.prediction.square().mean().backward()

    assert torch.allclose(output.evidence_score, torch.zeros((), dtype=DTYPE))
    assert model.zero_shot_head.weight.grad is not None
    assert torch.linalg.vector_norm(model.zero_shot_head.weight.grad) > 0


def test_level_calibration_is_a_learned_shrinkage_of_support_residual_mean():
    torch.manual_seed(41)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    output = episode(model)

    assert torch.allclose(output.level_adjustment, output.level_baseline - output.zero_shot)
    assert 0.0 < output.level_shrinkage < 1.0
    assert torch.allclose(output.shape_scale, torch.tensor(0.1, dtype=DTYPE))
    assert torch.allclose(output.sar_scale, torch.tensor(0.1, dtype=DTYPE))


def test_quotient_null_also_forces_zero_sar_adaptation():
    torch.manual_seed(39)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    support = torch.randn(3, 6, dtype=DTYPE)
    query = torch.randn(2, 6, dtype=DTYPE)
    with torch.no_grad():
        interaction = model.interaction_features(protein, tokens, mask, support)
        add, _, zero, _ = model.scalar_components(protein, support, interaction)
        labels = add + zero + 1.25

    output = model(protein, tokens, mask, support, labels, query)

    assert torch.allclose(output.sar_adaptation, torch.zeros(2, dtype=DTYPE), atol=1e-12)
    assert torch.allclose(output.evidence_score, torch.zeros((), dtype=DTYPE), atol=1e-12)


def test_level_shrinkage_increases_with_support_size():
    torch.manual_seed(43)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    one = episode(model, support_y=torch.randn(3, dtype=DTYPE))
    protein = torch.randn(6, dtype=DTYPE)
    tokens = torch.randn(5, 6, dtype=DTYPE)
    mask = torch.ones(5, dtype=torch.bool)
    query = torch.randn(2, 6, dtype=DTYPE)
    support_one = torch.randn(1, 6, dtype=DTYPE)
    support_five = torch.randn(5, 6, dtype=DTYPE)
    out_one = model(protein, tokens, mask, support_one, torch.randn(1, dtype=DTYPE), query)
    out_five = model(protein, tokens, mask, support_five, torch.randn(5, dtype=DTYPE), query)

    assert out_one.level_shrinkage < one.level_shrinkage < out_five.level_shrinkage


def test_sar_has_first_order_gradient_near_the_quotient_null():
    torch.manual_seed(42)
    model = QPSMPMetaLearner(6, 2, dtype=DTYPE)
    output = episode(model)

    output.sar_adaptation.sum().backward()

    assert model.section_head.weight.grad is not None
    assert torch.linalg.vector_norm(model.section_head.weight.grad) > 0


def test_neural_state_stays_inside_declared_euclidean_ball():
    model = QPSMPMetaLearner(6, 3, state_bound=0.7, dtype=DTYPE)
    output = episode(model)

    assert torch.linalg.vector_norm(output.task_state) < 0.7


def test_bio_model_query_loss_reaches_both_encoders_and_meta_operator():
    torch.manual_seed(38)
    model = QPSMPBioModel(
        8, 6, 2, ligand_layers=1, interaction_mode="atom_residue", dtype=DTYPE)
    protein = torch.randn(8, dtype=torch.float16)
    tokens = torch.randn(5, 8, dtype=torch.float16)
    protein_mask = torch.ones(5, dtype=torch.uint8)
    support_atoms = torch.randn(3, 4, ATOM_FEAT_DIM, dtype=DTYPE)
    query_atoms = torch.randn(2, 4, ATOM_FEAT_DIM, dtype=DTYPE)
    support_bonds = torch.zeros(3, 4, 4, BOND_FEAT_DIM, dtype=DTYPE)
    query_bonds = torch.zeros(2, 4, 4, BOND_FEAT_DIM, dtype=DTYPE)
    support_mask = torch.ones(3, 4, dtype=DTYPE)
    query_mask = torch.ones(2, 4, dtype=DTYPE)

    output = model(
        protein, tokens, protein_mask,
        support_atoms, support_bonds, support_mask,
        torch.randn(3, dtype=DTYPE),
        query_atoms, query_bonds, query_mask)
    output.prediction.square().mean().backward()

    parameters = (
        model.protein_encoder.bank_proj.weight,
        model.ligand_encoder.inp.weight,
        model.atom_residue_field.atom.weight,
        model.meta.meta_posterior.weight[0].weight,
    )
    assert all(parameter.grad is not None for parameter in parameters)
    assert all(torch.isfinite(parameter.grad).all() for parameter in parameters)


def test_bio_model_batched_forward_matches_episode_loop():
    for interaction_mode in ("pooled", "atom_residue"):
        torch.manual_seed(47)
        model = QPSMPBioModel(
            8, 6, 2, ligand_layers=1,
            interaction_mode=interaction_mode, dtype=DTYPE)
        batch_size, support_size, query_size, atoms = 2, 3, 2, 4
        protein = torch.randn(batch_size, 8, dtype=DTYPE)
        tokens = torch.randn(batch_size, 5, 8, dtype=DTYPE)
        protein_mask = torch.ones(batch_size, 5, dtype=torch.bool)
        support_atoms = torch.randn(
            batch_size, support_size, atoms, ATOM_FEAT_DIM, dtype=DTYPE)
        query_atoms = torch.randn(
            batch_size, query_size, atoms, ATOM_FEAT_DIM, dtype=DTYPE)
        support_bonds = torch.zeros(
            batch_size, support_size, atoms, atoms, BOND_FEAT_DIM, dtype=DTYPE)
        query_bonds = torch.zeros(
            batch_size, query_size, atoms, atoms, BOND_FEAT_DIM, dtype=DTYPE)
        support_mask = torch.ones(batch_size, support_size, atoms, dtype=DTYPE)
        query_mask = torch.ones(batch_size, query_size, atoms, dtype=DTYPE)
        support_y = torch.randn(batch_size, support_size, dtype=DTYPE)

        batched = model(
            protein, tokens, protein_mask,
            support_atoms, support_bonds, support_mask, support_y,
            query_atoms, query_bonds, query_mask)
        loop = [model(
            protein[index], tokens[index], protein_mask[index],
            support_atoms[index], support_bonds[index], support_mask[index],
            support_y[index], query_atoms[index], query_bonds[index],
            query_mask[index]) for index in range(batch_size)]

        assert torch.allclose(
            batched.prediction,
            torch.stack([output.prediction for output in loop]),
            atol=1e-10, rtol=1e-8)
        assert torch.allclose(
            batched.task_state,
            torch.stack([output.task_state for output in loop]),
            atol=1e-10, rtol=1e-8)


def test_support_span_state_is_in_centered_support_row_space():
    from model.qpsmp_meta import SupportSpanRidge

    solver = SupportSpanRidge(dtype=DTYPE)
    support = torch.randn(4, 7, dtype=DTYPE)
    query = torch.randn(3, 7, dtype=DTYPE)
    residual = torch.randn(4, dtype=DTYPE)
    residual = residual - residual.mean()

    state, centered_query, prediction = solver(support, query, residual)
    centered_support = support - support.mean(0, keepdim=True)
    projection = torch.linalg.pinv(centered_support) @ centered_support

    assert torch.allclose(state, projection @ state, atol=1e-9, rtol=1e-7)
    assert torch.allclose(prediction, centered_query @ state)


def test_support_span_single_shot_has_exact_zero_sar():
    from model.qpsmp_meta import SupportSpanRidge

    solver = SupportSpanRidge(dtype=DTYPE)
    state, _, prediction = solver(
        torch.randn(1, 5, dtype=DTYPE),
        torch.randn(3, 5, dtype=DTYPE),
        torch.zeros(1, dtype=DTYPE),
    )

    assert torch.count_nonzero(state) == 0
    assert torch.count_nonzero(prediction) == 0


def test_learned_support_span_posterior_is_identifiable_and_trainable():
    from model.qpsmp_meta import LearnedSupportSpanPosterior

    posterior = LearnedSupportSpanPosterior(dtype=DTYPE)
    support = torch.randn(5, 9, dtype=DTYPE)
    query = torch.randn(3, 9, dtype=DTYPE)
    residual = torch.randn(5, dtype=DTYPE)
    residual = residual - residual.mean()

    state, _, prediction = posterior(support, query, residual)
    centered_support = support - support.mean(0, keepdim=True)
    projection = torch.linalg.pinv(centered_support) @ centered_support
    prediction.square().mean().backward()

    assert torch.allclose(state, projection @ state, atol=1e-9, rtol=1e-7)
    assert posterior.weight[0].weight.grad is not None


def test_learned_support_span_has_exact_evidence_null_and_permutation_invariance():
    from model.qpsmp_meta import LearnedSupportSpanPosterior

    posterior = LearnedSupportSpanPosterior(dtype=DTYPE)
    support = torch.randn(4, 6, dtype=DTYPE)
    query = torch.randn(3, 6, dtype=DTYPE)
    zero = torch.zeros(4, dtype=DTYPE)
    state, _, prediction = posterior(support, query, zero)
    assert torch.count_nonzero(state) == 0
    assert torch.count_nonzero(prediction) == 0

    residual = torch.randn(4, dtype=DTYPE)
    residual = residual - residual.mean()
    order = torch.tensor([2, 0, 3, 1])
    original = posterior(support, query, residual)[2]
    permuted = posterior(support[order], query, residual[order])[2]
    assert torch.allclose(original, permuted, atol=1e-10, rtol=1e-8)
