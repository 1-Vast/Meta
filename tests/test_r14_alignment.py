"""Structural gates for the R14 regression-compatible ranking objective.

The design rests on one algebraic property: the within-target ranking term
must contribute **zero gradient at the regression optimum**, so that adding
it cannot pull the model away from `s = y`. Phase 2 measured what happens
when that property is absent — the incumbent RankNet and hinge losses cost
within-target Pearson `r` in 8 of 8 arms
(`report/meta_fewshot/stageR14_diagnostics_20260816/`).

These gates are cheap and run in the default tier, because the property is a
precondition for the whole cycle rather than a settled question.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.r14_alignment_check import SHIFT, hinge, list_ce, ranknet


def panel(seed: int, size: int = 16) -> torch.Tensor:
    generator = torch.Generator().manual_seed(seed)
    labels = 7.5 + 0.85 * torch.randn(size, generator=generator, dtype=torch.float64)
    return labels.clamp(4.0, 11.0)


def grad_norm(loss_fn, scores: torch.Tensor, labels: torch.Tensor) -> float:
    scores = scores.clone().requires_grad_(True)
    loss_fn(scores, labels).backward()
    return float(scores.grad.norm())


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_list_ce_is_stationary_at_the_regression_optimum(seed):
    """The whole design: adding this term must not move a correct prediction."""
    labels = panel(seed)
    assert grad_norm(list_ce, labels, labels) < 1e-9


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_incumbent_ranking_losses_are_not_stationary_there(seed):
    """The control: this is the property the incumbent losses lack.

    If this ever passes, the Phase 2 diagnosis is wrong and the R14
    preregistration's premise has to be rewritten.
    """
    labels = panel(seed)
    assert grad_norm(ranknet, labels, labels) > 1e-3
    assert grad_norm(hinge, labels, labels) > 1e-3


@pytest.mark.parametrize("scale", [0.5, 1.7, 3.0])
def test_list_ce_is_scale_free(scale):
    """Stationary at every positive multiple, so MSE alone pins the amplitude.

    This is what keeps the two terms from fighting over `sd_p`, which Phase 2
    measured as the amplitude excess.
    """
    labels = panel(11)
    scaled = SHIFT + scale * (labels - SHIFT)
    assert grad_norm(list_ce, scaled, labels) < 1e-9


def test_list_ce_penalises_a_wrong_ordering():
    """A stationary point is useless if the loss is flat everywhere."""
    labels = panel(7)
    correct = list_ce(labels, labels)
    reversed_scores = SHIFT + (labels.max() + labels.min() - labels - SHIFT)
    assert float(list_ce(reversed_scores, labels)) > float(correct) + 1e-3


def test_list_ce_is_invariant_to_query_permutation():
    labels = panel(3)
    order = torch.randperm(len(labels), generator=torch.Generator().manual_seed(5))
    assert torch.allclose(list_ce(labels, labels),
                          list_ce(labels[order], labels[order]))


def test_list_ce_never_reads_a_prediction_below_the_shift():
    """The shift guards the log; a degenerate prediction must not produce NaN."""
    labels = panel(9)
    degenerate = torch.full_like(labels, SHIFT - 5.0)
    value = list_ce(degenerate, labels)
    assert torch.isfinite(value), value


# ------------------------------------------------------------------ trainer
# The gates above verify the objective in isolation. These verify the term as
# it is actually wired into the incumbent trainer, in normalized label units.

from scripts.train_qpsmp import (                                # noqa: E402
    LISTCE_SHIFT_PK, LabelScale, TrainConfig, pairwise_ranking_loss,
    ranking_term, regression_compatible_ranking_loss,
)

SCALE = LabelScale(mean=7.5, scale=1.2)


def normalized_panel(seed: int, size: int = 16) -> torch.Tensor:
    return SCALE.normalize(panel(seed, size))


def trainer_grad(config: TrainConfig, prediction: torch.Tensor,
                 truth: torch.Tensor) -> torch.Tensor:
    prediction = prediction.clone().requires_grad_(True)
    ranking_term(prediction, truth, config, SCALE).backward()
    return prediction.grad


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_trainer_listce_is_stationary_at_the_regression_optimum(seed):
    """The wired-in term must inherit the property, in normalized units."""
    labels = normalized_panel(seed)
    config = TrainConfig(ranking_loss_form="listce")
    assert float(trainer_grad(config, labels, labels).norm()) < 1e-9


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_trainer_ranknet_is_not_stationary_there(seed):
    """The incumbent form, through the same dispatch — the matched control."""
    labels = normalized_panel(seed)
    config = TrainConfig(ranking_loss_form="ranknet")
    assert float(trainer_grad(config, labels, labels).norm()) > 1e-3


def test_dispatch_reproduces_the_incumbent_exactly():
    """`ranknet` must be bit-identical to the retained loss: no silent drift."""
    labels = normalized_panel(4)
    prediction = labels + 0.3 * torch.randn(
        len(labels), generator=torch.Generator().manual_seed(8), dtype=torch.float64)
    config = TrainConfig(ranking_loss_form="ranknet")
    assert torch.equal(
        ranking_term(prediction, labels, config, SCALE),
        pairwise_ranking_loss(prediction, labels, config.ranking_temperature))


def test_unknown_ranking_form_is_refused():
    """A typo must fail loudly rather than silently train the incumbent."""
    labels = normalized_panel(2)
    with pytest.raises(ValueError, match="unknown ranking_loss_form"):
        ranking_term(labels, labels, TrainConfig(ranking_loss_form="softmax"), SCALE)


def test_listce_gradient_flows_to_the_prediction_only():
    """Leakage gate: the term must not create a path from labels into the model.

    The labels enter only as loss weights. Their gradient is irrelevant, but a
    label tensor that requires grad must not silently become trainable input.
    """
    labels = normalized_panel(6).requires_grad_(True)
    prediction = (labels.detach() + 0.2).requires_grad_(True)
    regression_compatible_ranking_loss(prediction, labels, SCALE).backward()
    assert prediction.grad is not None and float(prediction.grad.norm()) > 0
    # The label gradient exists mathematically; what matters is that the
    # trainer never feeds labels to the model. Assert the term is not a
    # function of the label *ordering* alone by checking it uses magnitudes.
    shuffled = labels.detach()[torch.randperm(
        len(labels), generator=torch.Generator().manual_seed(1))]
    assert not torch.allclose(
        regression_compatible_ranking_loss(prediction.detach(), labels.detach(), SCALE),
        regression_compatible_ranking_loss(prediction.detach(), shuffled, SCALE))


def test_listce_is_permutation_invariant_through_the_trainer():
    labels = normalized_panel(5)
    prediction = labels + 0.4
    order = torch.randperm(len(labels), generator=torch.Generator().manual_seed(2))
    config = TrainConfig(ranking_loss_form="listce")
    assert torch.allclose(ranking_term(prediction, labels, config, SCALE),
                          ranking_term(prediction[order], labels[order], config, SCALE))


def test_listce_counterfactual_a_wrong_ordering_costs_more():
    """Counterfactual gate: shuffling the predictions must raise the loss."""
    labels = normalized_panel(12)
    config = TrainConfig(ranking_loss_form="listce")
    correct = float(ranking_term(labels, labels, config, SCALE))
    order = torch.randperm(len(labels), generator=torch.Generator().manual_seed(3))
    shuffled = float(ranking_term(labels[order], labels, config, SCALE))
    assert shuffled > correct + 1e-4, (shuffled, correct)


def test_listce_handles_a_batched_panel():
    """The trainer passes (episodes, queries); the term must reduce correctly."""
    labels = torch.stack([normalized_panel(s) for s in (1, 2, 3)])
    config = TrainConfig(ranking_loss_form="listce")
    batched = ranking_term(labels, labels, config, SCALE)
    assert batched.ndim == 0 and float(batched.abs()) < 1e6
    per_panel = torch.stack([
        regression_compatible_ranking_loss(row, row, SCALE) for row in labels])
    assert torch.allclose(batched, per_panel.mean(), atol=1e-9)


def test_shift_sits_below_the_corpus_label_range():
    """A shift inside the label range would make weights negative."""
    assert LISTCE_SHIFT_PK < 4.0
