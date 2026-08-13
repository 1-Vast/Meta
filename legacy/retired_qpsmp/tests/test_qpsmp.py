import torch

from model.qpsmp import QPSMPCore


def test_centered_state_is_invariant_to_support_level_shift():
    torch.manual_seed(12)
    model = QPSMPCore(
        feature_dim=3, task_dim=2, ridge=0.4,
        section_radius_bound=1.5, dtype=torch.float64)
    support = torch.randn(4, 3, dtype=torch.float64)
    query = torch.randn(3, 3, dtype=torch.float64)
    support_y = torch.randn(4, dtype=torch.float64)

    left = model(support, support_y, query)
    right = model(support, support_y + 7.0, query)

    assert torch.allclose(left.task_state, right.task_state, atol=1e-12)
    assert torch.allclose(right.level_shift, left.level_shift + 7.0, atol=1e-12)
    assert torch.allclose(right.ridge_prediction, left.ridge_prediction + 7.0, atol=1e-12)


def test_level_statistic_is_invariant_to_support_label_binding():
    torch.manual_seed(17)
    model = QPSMPCore(feature_dim=3, task_dim=2, ridge=0.4, dtype=torch.float64)
    support = torch.randn(4, 3, dtype=torch.float64)
    query = torch.randn(2, 3, dtype=torch.float64)
    support_y = torch.randn(4, dtype=torch.float64)

    left = model(support, support_y, query)
    right = model(support, support_y.flip(0), query)

    assert torch.allclose(left.level_shift, right.level_shift, atol=1e-12)


def test_disabling_support_adaptation_is_independent_of_support_labels():
    torch.manual_seed(18)
    model = QPSMPCore(feature_dim=3, task_dim=2, ridge=0.4, dtype=torch.float64)
    support = torch.randn(4, 3, dtype=torch.float64)
    query = torch.randn(2, 3, dtype=torch.float64)

    left = model(support, torch.randn(4, dtype=torch.float64), query, adapt=False)
    right = model(support, torch.randn(4, dtype=torch.float64), query, adapt=False)

    assert torch.allclose(left.ridge_prediction, right.ridge_prediction, atol=1e-12)
    assert torch.allclose(left.level_shift, torch.zeros((), dtype=torch.float64))
    assert torch.allclose(left.task_state, torch.zeros(2, dtype=torch.float64))


def test_single_support_point_cannot_identify_centered_sar_state():
    torch.manual_seed(13)
    model = QPSMPCore(feature_dim=4, task_dim=3, ridge=1.0, dtype=torch.float64)
    support = torch.randn(1, 4, dtype=torch.float64)
    query = torch.randn(2, 4, dtype=torch.float64)

    output = model(support, torch.tensor([5.0], dtype=torch.float64), query)

    assert torch.allclose(output.task_state, torch.zeros(3, dtype=torch.float64))


def test_delta_and_rectangle_are_derived_antisymmetric_quotients():
    predictions = torch.tensor([1.0, 4.0, -2.0], dtype=torch.float64)
    left = torch.tensor([0])
    right = torch.tensor([1])

    delta_lr = QPSMPCore.delta(predictions, left, right)
    delta_rl = QPSMPCore.delta(predictions, right, left)

    assert torch.allclose(delta_lr, -delta_rl)
    assert torch.allclose(
        QPSMPCore.rectangle(delta_lr, torch.tensor([2.0], dtype=torch.float64)),
        -QPSMPCore.rectangle(torch.tensor([2.0], dtype=torch.float64), delta_lr),
    )


def test_diagnostic_radius_contains_center_row_space_and_external_terms():
    torch.manual_seed(14)
    model = QPSMPCore(
        feature_dim=3, task_dim=2, ridge=0.5,
        section_radius_bound=2.0, dtype=torch.float64)
    support = torch.randn(3, 3, dtype=torch.float64)
    query = torch.randn(4, 3, dtype=torch.float64)
    output = model(
        support, torch.randn(3, dtype=torch.float64), query,
        repr_radius=0.1, trans_radius=0.2, obs_radius=0.3)

    expected = output.center_radius + output.section_radius + 0.6

    assert torch.all(output.section_radius >= 0)
    assert torch.allclose(output.diagnostic_total_radius, expected, atol=1e-12)


def test_zero_task_dimension_keeps_level_calibration_but_has_no_sar_state():
    torch.manual_seed(15)
    model = QPSMPCore(feature_dim=3, task_dim=0, ridge=1.0, dtype=torch.float64)
    support = torch.randn(2, 3, dtype=torch.float64)
    query = torch.randn(2, 3, dtype=torch.float64)

    output = model(support, torch.randn(2, dtype=torch.float64), query)
    baseline, zero_shot = model.scalar_components(query)

    assert output.task_state.numel() == 0
    assert output.query_basis.shape == (2, 0)
    assert torch.allclose(output.section_radius, torch.zeros(2, dtype=torch.float64))
    assert torch.allclose(
        output.section_midpoint, baseline + zero_shot + output.level_shift)


def test_baseline_features_are_separate_from_interaction_features():
    torch.manual_seed(16)
    model = QPSMPCore(
        feature_dim=3, baseline_dim=2, task_dim=1,
        ridge=1.0, dtype=torch.float64)
    support_interaction = torch.randn(3, 3, dtype=torch.float64)
    query_interaction = torch.randn(2, 3, dtype=torch.float64)
    support_baseline = torch.randn(3, 2, dtype=torch.float64)
    query_baseline = torch.randn(2, 2, dtype=torch.float64)

    output = model(
        support_interaction, torch.randn(3, dtype=torch.float64),
        query_interaction,
        support_baseline_features=support_baseline,
        query_baseline_features=query_baseline)

    assert output.query_basis.shape == (2, 1)
    assert torch.allclose(output.query_basis, model.basis(query_interaction))
