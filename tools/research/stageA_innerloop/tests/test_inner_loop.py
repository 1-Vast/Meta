"""Structural contract for Stage A. Every property required before training.

Two tests carry most of the weight:

* `test_readout_reproduces_the_model_endpoint` — the inner loop re-evaluates
  the readout outside the model for speed, duplicating three lines of
  `InteractionGrammarModel.encode`. If that duplicate ever drifts, the
  experiment measures a different model than the one it reports. This pins it
  bitwise.
* `test_zero_inner_steps_reproduces_the_production_forward` — `A0` is claimed
  to be matched to the accepted baseline *by construction*, because it is the
  same code path with `inner_steps=0`. This is that claim, as a test.

Synthetic tensors and a small model throughout: these are algebraic properties,
so they need no dataset and must hold in float64 to a tight tolerance.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM      # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel          # noqa: E402
from tools.research.stageA_innerloop.inner_loop import (             # noqa: E402
    ADAPTABLE, AdaptationConfig, adapt, base_weights, encode_task,
    partial_weights, readout, support_query_gradient_cosine,
)
from tools.research.stageA_innerloop.train_meta import (             # noqa: E402
    predict, select_tasks, standardize, task_value,
)

PROTEIN_DIM, SLOTS, ATOMS = 32, 12, 9


def build(seed: int = 0) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12, ligand_layers=2,
        pair_dim=24, pair_latents=6, pair_heads=2,
        use_learned_key=False).double()


def make_parts(support: int, query: int, seed: int = 0) -> dict:
    generator = torch.Generator().manual_seed(seed)

    def graphs(count):
        atoms = torch.rand(1, count, ATOMS, ATOM_FEAT_DIM, generator=generator,
                           dtype=torch.float64)
        bonds = torch.rand(1, count, ATOMS, ATOMS, BOND_FEAT_DIM,
                           generator=generator, dtype=torch.float64)
        bonds = bonds * (bonds > 0.7)
        bonds = 0.5 * (bonds + bonds.transpose(2, 3))
        mask = torch.ones(1, count, ATOMS, dtype=torch.float64)
        fingerprint = (torch.rand(1, count, 64, generator=generator,
                                  dtype=torch.float64) > 0.7).double()
        return atoms, bonds, mask, fingerprint

    support_atoms, support_bonds, support_mask, support_fp = graphs(support)
    query_atoms, query_bonds, query_mask, query_fp = graphs(query)
    return {
        "protein_pooled": torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        "protein_tokens": torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        "protein_mask": torch.ones(1, SLOTS, dtype=torch.float64),
        "protein_chemistry": torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64),
        "support_atoms": support_atoms, "support_bonds": support_bonds,
        "support_mask": support_mask, "support_fingerprint": support_fp,
        "query_atoms": query_atoms, "query_bonds": query_bonds,
        "query_mask": query_mask, "query_fingerprint": query_fp,
        "support_y": torch.randn(1, support, generator=generator, dtype=torch.float64),
        "query_y": torch.randn(1, query, generator=generator, dtype=torch.float64),
    }


def encoded(model, parts):
    return encode_task(
        model, parts["protein_pooled"], parts["protein_tokens"],
        parts["protein_mask"], parts["support_atoms"], parts["support_bonds"],
        parts["support_mask"], parts["query_atoms"], parts["query_bonds"],
        parts["query_mask"], parts["protein_chemistry"])


def model_forward(model, parts, adapt_flag=True):
    return model(
        parts["protein_pooled"], parts["protein_tokens"], parts["protein_mask"],
        parts["support_atoms"], parts["support_bonds"], parts["support_mask"],
        parts["support_y"], parts["query_atoms"], parts["query_bonds"],
        parts["query_mask"], adapt=adapt_flag,
        protein_chemistry=parts["protein_chemistry"],
        support_fingerprint=parts["support_fingerprint"],
        query_fingerprint=parts["query_fingerprint"])


# --- the two load-bearing identities ---------------------------------------

@pytest.mark.parametrize("support,query", [(0, 6), (1, 6), (3, 5), (5, 7)])
def test_readout_reproduces_the_model_endpoint(support, query):
    """The re-implemented readout must equal `encode`'s own endpoint."""
    model = build()
    parts = make_parts(support, query, seed=support + 1)
    task = encoded(model, parts)
    weights = base_weights(model)
    with torch.no_grad():
        mine = readout(model, task.query_hidden, task.query_additive,
                       task.query_occupancy, weights)
        theirs = model_forward(model, parts, adapt_flag=False).zero_shot
    assert torch.allclose(mine, theirs, atol=1e-12, rtol=0.0)


@pytest.mark.parametrize("support,query", [(0, 6), (1, 6), (3, 5), (5, 7)])
def test_zero_inner_steps_reproduces_the_production_forward(support, query):
    """`A0` is matched by construction; this is the construction."""
    model = build(1)
    parts = make_parts(support, query, seed=support + 20)
    config = AdaptationConfig(inner_steps=0)
    with torch.no_grad():
        mine = predict(model, parts, config)["prediction"]
        theirs = model_forward(model, parts).prediction
    assert torch.allclose(mine, theirs, atol=1e-12, rtol=0.0)


# --- k = 0 --------------------------------------------------------------

def test_k0_prediction_is_exactly_the_zero_shot_path():
    model = build(2)
    parts = make_parts(0, 8, seed=3)
    config = AdaptationConfig(inner_steps=3, inner_lr=0.5)
    with torch.no_grad():
        output = predict(model, parts, config)
        frozen = model_forward(model, parts, adapt_flag=False).zero_shot
    assert torch.allclose(output["prediction"], frozen, atol=1e-12, rtol=0.0)
    assert output["inner_trace"] == []
    assert torch.count_nonzero(output["transport"]) == 0


def test_k0_returns_the_base_weights_unmodified():
    model = build(2)
    parts = make_parts(0, 8, seed=4)
    task = encoded(model, parts)
    fast, trace = adapt(model, task, parts["support_y"],
                        AdaptationConfig(inner_steps=5))
    assert trace == []
    for name, value in base_weights(model).items():
        assert torch.equal(fast[name], value)


# --- query labels never enter ---------------------------------------------

@pytest.mark.parametrize("inner_steps", [1, 2, 3])
def test_no_query_label_reaches_the_prediction(inner_steps):
    """Change the query labels arbitrarily; the prediction must not move."""
    model = build(3)
    parts = make_parts(4, 6, seed=5)
    config = AdaptationConfig(inner_steps=inner_steps, inner_lr=0.3)
    with torch.no_grad():
        first = predict(model, parts, config)["prediction"]
        shifted = dict(parts)
        shifted["query_y"] = parts["query_y"] * -7.0 + 3.0
        second = predict(model, shifted, config)["prediction"]
    assert torch.allclose(first, second, atol=0.0, rtol=0.0)


def test_the_inner_loss_reads_only_support_labels():
    """A support-label change must move the inner trace; a query one must not."""
    model = build(3)
    parts = make_parts(4, 6, seed=6)
    config = AdaptationConfig(inner_steps=2, inner_lr=0.3)
    base = predict(model, parts, config)["inner_trace"]
    moved = dict(parts)
    moved["support_y"] = parts["support_y"] + 1.0
    assert predict(model, moved, config)["inner_trace"] != base
    untouched = dict(parts)
    untouched["query_y"] = parts["query_y"] + 1.0
    assert predict(model, untouched, config)["inner_trace"] == base


# --- support labels do something at k = 1 ----------------------------------

def test_k1_has_a_nonzero_support_label_gradient():
    """The support label must reach the prediction differentiably.

    Taken in the second-order mode, where the autograd path through the inner
    gradient exists. First-order adaptation deliberately severs that path (see
    `adapt`), so the substantive check for the arm that actually runs is the
    functional-sensitivity test below.
    """
    model = build(4)
    parts = make_parts(1, 6, seed=7)
    support_y = parts["support_y"].clone().requires_grad_(True)
    parts = {**parts, "support_y": support_y}
    output = predict(model, parts,
                     AdaptationConfig(inner_steps=1, inner_lr=0.3,
                                      first_order=False))
    output["prediction"].sum().backward()
    assert support_y.grad is not None
    assert float(support_y.grad.abs().max()) > 1e-9


@pytest.mark.parametrize("first_order", [True, False])
def test_k1_support_labels_functionally_change_the_prediction(first_order):
    """The property that matters for the arm that runs: labels do something.

    A first-order inner loop has no differentiable path from `support_y`
    through the adaptation, so a gradient test alone would understate the
    mechanism. Perturbing the label and watching the prediction move is the
    direct statement, and it must hold in both modes.
    """
    model = build(4)
    parts = make_parts(1, 6, seed=70)
    config = AdaptationConfig(inner_steps=1, inner_lr=0.5,
                              first_order=first_order)
    with torch.no_grad():
        base = predict(model, parts, config)["prediction"]
        moved = dict(parts)
        moved["support_y"] = parts["support_y"] + 2.0
        shifted = predict(model, moved, config)["prediction"]
    assert float((shifted - base).abs().max()) > 1e-6


def test_k1_adaptation_is_not_a_pure_level_shift():
    """The weight update must reorder queries, not merely offset them.

    A k=1 gain that survives `keep="bias"` is a recalibration. This checks that
    the mechanism is *capable* of shape change; whether it helps is the
    experiment, not the test.
    """
    model = build(4)
    parts = make_parts(1, 8, seed=8)
    config = AdaptationConfig(inner_steps=1, inner_lr=0.5)
    with torch.no_grad():
        base = predict(model, parts, AdaptationConfig(inner_steps=0))["prediction"]
        both = predict(model, parts, config, keep="both")["prediction"]
        bias = predict(model, parts, config, keep="bias")["prediction"]
    shift_both = both - base
    shift_bias = bias - base
    assert float(shift_bias.std()) < 1e-12          # bias alone is pure level
    assert float(shift_both.std()) > 1e-9           # the pair is not


def test_partial_weight_modes_are_exhaustive_and_distinct():
    model = build(4)
    parts = make_parts(3, 6, seed=9)
    task = encoded(model, parts)
    fast, _ = adapt(model, task, parts["support_y"],
                    AdaptationConfig(inner_steps=1, inner_lr=0.5))
    base = base_weights(model)
    assert torch.equal(partial_weights(model, fast, "none")["interaction_head.2.weight"],
                       base["interaction_head.2.weight"])
    assert torch.equal(partial_weights(model, fast, "bias")["interaction_head.2.weight"],
                       base["interaction_head.2.weight"])
    assert torch.equal(partial_weights(model, fast, "weight")["interaction_head.2.bias"],
                       base["interaction_head.2.bias"])
    with pytest.raises(ValueError):
        partial_weights(model, fast, "elsewhere")


# --- scope, isolation, and mutation ----------------------------------------

def test_the_inner_loop_touches_only_the_declared_scope():
    model = build(5)
    parts = make_parts(3, 6, seed=10)
    before = {n: p.detach().clone() for n, p in model.named_parameters()}
    task = encoded(model, parts)
    fast, _ = adapt(model, task, parts["support_y"],
                    AdaptationConfig(inner_steps=3, inner_lr=0.5))
    assert set(fast) == set(ADAPTABLE)
    for name, parameter in model.named_parameters():
        assert torch.equal(parameter.detach(), before[name]), (
            f"the inner loop mutated the persistent parameter {name}")


def test_adaptation_does_not_leak_between_tasks():
    """One task's fast weights must not change another task's prediction."""
    model = build(6)
    first = make_parts(5, 6, seed=11)
    second = make_parts(5, 6, seed=12)
    config = AdaptationConfig(inner_steps=3, inner_lr=0.5)
    with torch.no_grad():
        clean = predict(model, second, config)["prediction"]
        predict(model, first, config)
        after = predict(model, second, config)["prediction"]
    assert torch.allclose(clean, after, atol=0.0, rtol=0.0)


def test_state_dict_is_unchanged_by_a_full_training_style_step():
    model = build(6)
    parts = make_parts(5, 6, seed=13)
    before = {k: v.clone() for k, v in model.state_dict().items()}
    output = predict(model, parts, AdaptationConfig(inner_steps=2, inner_lr=0.4))
    output["prediction"].sum().backward()
    for key, value in model.state_dict().items():
        assert torch.equal(value, before[key])


# --- permutation symmetry ---------------------------------------------------

def test_support_permutation_invariance():
    model = build(7)
    parts = make_parts(5, 6, seed=14)
    order = torch.randperm(5, generator=torch.Generator().manual_seed(1))
    permuted = dict(parts)
    for key in ("support_atoms", "support_bonds", "support_mask",
                "support_fingerprint"):
        permuted[key] = parts[key][:, order]
    permuted["support_y"] = parts["support_y"][:, order]
    config = AdaptationConfig(inner_steps=2, inner_lr=0.3)
    with torch.no_grad():
        assert torch.allclose(predict(model, parts, config)["prediction"],
                              predict(model, permuted, config)["prediction"],
                              atol=1e-11)


def test_query_permutation_equivariance():
    model = build(7)
    parts = make_parts(4, 7, seed=15)
    order = torch.randperm(7, generator=torch.Generator().manual_seed(2))
    permuted = dict(parts)
    for key in ("query_atoms", "query_bonds", "query_mask", "query_fingerprint"):
        permuted[key] = parts[key][:, order]
    config = AdaptationConfig(inner_steps=2, inner_lr=0.3)
    with torch.no_grad():
        straight = predict(model, parts, config)["prediction"]
        shuffled = predict(model, permuted, config)["prediction"]
    assert torch.allclose(straight[:, order], shuffled, atol=1e-11)


# --- meta-gradient ----------------------------------------------------------

def test_first_order_meta_gradient_reaches_the_shared_initialization():
    """The outer loss must train the encoders, not only the adapted readout."""
    model = build(8)
    parts = make_parts(3, 6, seed=16)
    output = predict(model, parts, AdaptationConfig(inner_steps=2, inner_lr=0.3))
    F.mse_loss(output["prediction"], parts["query_y"]).backward()
    named = dict(model.named_parameters())
    for branch in ("ligand_encoder", "protein_encoder", "grammar", "embed",
                   "interaction_head", "section"):
        reached = max(
            float(p.grad.abs().max()) for n, p in named.items()
            if n.startswith(branch) and p.grad is not None)
        assert reached > 0.0, f"no outer gradient reached {branch}"
    for name in ADAPTABLE:
        assert named[name].grad is not None
        assert float(named[name].grad.abs().max()) > 0.0


def test_first_order_detaches_the_inner_gradient():
    """`create_graph=False` is what makes this first order rather than full MAML."""
    model = build(8)
    parts = make_parts(3, 6, seed=17)
    task = encoded(model, parts)
    fast, _ = adapt(model, task, parts["support_y"],
                    AdaptationConfig(inner_steps=1, first_order=True))
    second, _ = adapt(model, task, parts["support_y"],
                      AdaptationConfig(inner_steps=1, first_order=False))
    # Both remain differentiable w.r.t. the base parameters, but only the
    # second-order variant keeps a graph through the inner gradient itself.
    assert fast["interaction_head.2.bias"].requires_grad
    assert second["interaction_head.2.bias"].requires_grad
    assert torch.allclose(fast["interaction_head.2.bias"],
                          second["interaction_head.2.bias"], atol=1e-12)


# --- the selector -----------------------------------------------------------

def test_gradient_cosine_is_finite_and_bounded():
    model = build(9)
    parts = make_parts(4, 6, seed=18)
    task = encoded(model, parts)
    value = support_query_gradient_cosine(
        model, task, parts["support_y"], parts["query_y"])
    assert np.isfinite(value) and -1.0 <= value <= 1.0


def test_gradient_cosine_returns_zero_for_an_empty_side():
    model = build(9)
    parts = make_parts(0, 6, seed=19)
    task = encoded(model, parts)
    assert support_query_gradient_cosine(
        model, task, parts["support_y"], parts["query_y"]) == 0.0


def test_gradient_cosine_is_zero_when_a_gradient_vanishes():
    """A zero-norm gradient carries no direction; 0.0 is the honest encoding."""
    model = build(9)
    parts = make_parts(4, 6, seed=20)
    task = encoded(model, parts)
    with torch.no_grad():
        exact = readout(model, task.support_hidden, task.support_additive,
                        task.support_occupancy, base_weights(model))
    # Support labels equal to the model's own output => zero support gradient.
    assert support_query_gradient_cosine(
        model, task, exact, parts["query_y"]) == 0.0


def test_standardize_handles_a_degenerate_batch():
    assert np.allclose(standardize([2.0, 2.0, 2.0]), 0.0)
    assert np.allclose(standardize([1.0]), 0.0)
    values = standardize([1.0, 2.0, 3.0])
    assert abs(float(values.mean())) < 1e-12
    assert abs(float(values.std()) - 1.0) < 1e-12


def test_selection_prefers_low_loss_and_high_agreement():
    scores = [
        {"post_adaptation_query_loss": 5.0, "gradient_cosine": -0.9},
        {"post_adaptation_query_loss": 0.1, "gradient_cosine": 0.9},
        {"post_adaptation_query_loss": 3.0, "gradient_cosine": 0.0},
    ]
    picked, summary = select_tasks(scores, keep=1)
    assert picked == [1]
    assert summary["candidates"] == 3
    assert 1.0 <= summary["effective_tasks"] <= 3.0


def test_selection_is_stable_when_every_candidate_ties():
    scores = [{"post_adaptation_query_loss": 1.0, "gradient_cosine": 0.5}] * 4
    picked, summary = select_tasks(scores, keep=2)
    assert len(picked) == 2
    assert summary["effective_tasks"] == pytest.approx(4.0, abs=1e-9)


def test_counterfactual_supports_match_the_production_definition():
    """The wrong-label controls must be the accepted recipe's, not new ones.

    `counterfactual_supports` reimplements
    `scripts/train_qpsmp.py::counterfactual_label_assignments` because the
    production helper reads a `QPSMPMetaOutput` the adapted path never builds.
    Reimplementation is only safe while it stays identical, so it is checked
    against the original here.
    """
    from scripts.train_qpsmp import counterfactual_label_assignments
    from tools.research.stageA_innerloop.train_meta import counterfactual_supports

    class FakeOutput:
        def __init__(self, quotient, adjustment):
            self.support_residual_quotient = quotient
            self.level_adjustment = adjustment

    generator = torch.Generator().manual_seed(0)
    for count in (2, 3, 5):
        support_y = torch.randn(1, count, generator=generator, dtype=torch.float64)
        prediction = torch.randn(1, count, generator=generator, dtype=torch.float64)
        mine = counterfactual_supports(support_y, prediction)
        theirs = counterfactual_label_assignments(
            FakeOutput(torch.zeros_like(support_y),
                       torch.zeros(1, 1, dtype=torch.float64)), support_y)
        assert len(mine) == len(theirs) == count - 1
        for left, right in zip(mine, theirs):
            assert torch.allclose(left, right, atol=0.0)

    # k = 1: the equal-magnitude residual flip.
    support_y = torch.randn(1, 1, generator=generator, dtype=torch.float64)
    prediction = torch.randn(1, 1, generator=generator, dtype=torch.float64)
    locked = support_y - prediction
    mine = counterfactual_supports(support_y, prediction)
    theirs = counterfactual_label_assignments(
        FakeOutput(locked, torch.zeros(1, 1, dtype=torch.float64)), support_y)
    assert len(mine) == 1
    assert torch.allclose(mine[0], theirs[0], atol=1e-12)
    # The flip preserves the residual magnitude and reverses its sign.
    assert torch.allclose(mine[0] - prediction, -(support_y - prediction), atol=1e-12)


def test_counterfactual_supports_is_empty_at_k0():
    from tools.research.stageA_innerloop.train_meta import counterfactual_supports
    empty = torch.zeros(1, 0, dtype=torch.float64)
    assert counterfactual_supports(empty, empty) == []


def test_task_value_reports_both_terms():
    model = build(10)
    parts = make_parts(3, 6, seed=21)
    value = task_value(model, parts, AdaptationConfig(inner_steps=1))
    assert set(value) == {"post_adaptation_query_loss", "gradient_cosine"}
    assert np.isfinite(value["post_adaptation_query_loss"])
    assert np.isfinite(value["gradient_cosine"])


# --- checkpoint -------------------------------------------------------------

def test_adaptation_config_round_trips():
    config = AdaptationConfig(inner_steps=2, inner_lr=0.25, first_order=False)
    restored = AdaptationConfig.from_dict(config.to_dict())
    assert restored == config
    assert restored.scope == ADAPTABLE


def test_checkpoint_reloads_strictly_with_its_adaptation_config(tmp_path):
    model = build(11)
    config = AdaptationConfig(inner_steps=2, inner_lr=0.25)
    path = tmp_path / "checkpoint.pt"
    torch.save({"model_state": model.state_dict(),
                "adaptation": config.to_dict()}, path)
    payload = torch.load(path, map_location="cpu", weights_only=False)
    restored = build(12)
    restored.load_state_dict(payload["model_state"], strict=True)
    assert AdaptationConfig.from_dict(payload["adaptation"]) == config
    parts = make_parts(3, 6, seed=22)
    with torch.no_grad():
        assert torch.allclose(predict(model, parts, config)["prediction"],
                              predict(restored, parts, config)["prediction"],
                              atol=0.0)


def test_a_checkpoint_without_adaptation_config_is_rejected():
    with pytest.raises(KeyError):
        AdaptationConfig.from_dict({"inner_lr": 0.1})
