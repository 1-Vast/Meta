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
