from types import SimpleNamespace

import numpy as np
import torch

from scripts.train_qpsmp import (LabelScale, admission_score,
                                 binding_contrastive_loss, centered_task_error,
                                 matched_wrong_labels, pairwise_ranking_loss,
                                 training_label_scale)


def test_label_scale_uses_meta_train_only():
    data = SimpleNamespace(cells=[
        {"split": "meta_train", "pK": 2.0},
        {"split": "meta_train", "pK": 4.0},
        {"split": "meta_val", "pK": 1000.0},
    ])
    scale = training_label_scale(data)
    assert np.isclose(scale.mean, 3.0)
    assert np.isclose(scale.scale, 1.0)


def test_label_scale_reports_raw_pk_squared_error():
    scale = LabelScale(6.0, 2.0)
    assert torch.equal(
        scale.squared_error_pk(torch.tensor([0.0, 1.0]), torch.tensor([1.0, 0.0])),
        torch.tensor([4.0, 4.0]))


def test_pairwise_ranking_loss_rewards_correct_order():
    truth = torch.tensor([1.0, 2.0, 4.0])
    correct = pairwise_ranking_loss(truth, truth, 1.0)
    reversed_order = pairwise_ranking_loss(truth.flip(0), truth, 1.0)
    assert correct < reversed_order


def test_centered_task_error_ignores_level_and_rewards_query_shape():
    truth = torch.tensor([1.0, 2.0, 4.0])
    assert torch.allclose(
        centered_task_error(truth + 9.0, truth), torch.zeros(()))
    assert centered_task_error(truth.flip(0), truth) > 0


def test_binding_contrast_rewards_correct_assignment():
    good = binding_contrastive_loss(
        [torch.tensor(0.1), torch.tensor(1.0)], 0.2)
    bad = binding_contrastive_loss(
        [torch.tensor(1.0), torch.tensor(0.1)], 0.2)
    assert good < bad


def test_one_shot_counterfactual_preserves_residual_magnitude():
    output = SimpleNamespace(
        support_residual_quotient=torch.tensor([0.6]),
        level_adjustment=torch.tensor([0.2, 0.2]))
    label = torch.tensor([1.5])
    wrong = matched_wrong_labels(output, label)
    assert torch.allclose(wrong, torch.tensor([-0.1]))


def test_admission_score_penalizes_harm_and_binding_insensitivity():
    admitted = {"full_mse_pk": 1.0, "sar_cut_mse_pk": 1.2,
                "permuted_mse_pk": 1.1}
    dead = {"full_mse_pk": 1.1, "sar_cut_mse_pk": 1.0,
            "permuted_mse_pk": 1.1}
    assert admission_score(admitted, 0.01) < admission_score(dead, 0.01)
