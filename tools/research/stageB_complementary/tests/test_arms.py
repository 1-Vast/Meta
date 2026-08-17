"""Structural contract for the four Stage B arms.

The load-bearing test is `test_C_equals_T_plus_the_meta_correction`: the
candidate must differ from the accepted baseline by exactly one additive term,
or the discriminator does not isolate what it claims to isolate.
"""
from __future__ import annotations

from pathlib import Path
import sys

import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM      # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel          # noqa: E402
from tools.research.stageA_innerloop.train_meta import encode_parts  # noqa: E402
from tools.research.stageB_complementary.arms import (               # noqa: E402
    ADAPTED_BY_MODE, InnerStepSizes, MODES, StageBAdaptation, predict,
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
        fp = (torch.rand(1, count, 64, generator=generator,
                         dtype=torch.float64) > 0.7).double()
        return atoms, bonds, mask, fp

    sa, sb, sm, sf = graphs(support)
    qa, qb, qm, qf = graphs(query)
    return {
        "protein_pooled": torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        "protein_tokens": torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        "protein_mask": torch.ones(1, SLOTS, dtype=torch.float64),
        "protein_chemistry": torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64),
        "support_atoms": sa, "support_bonds": sb, "support_mask": sm,
        "support_fingerprint": sf,
        "query_atoms": qa, "query_bonds": qb, "query_mask": qm,
        "query_fingerprint": qf,
        "support_y": torch.randn(1, support, generator=generator, dtype=torch.float64),
        "query_y": torch.randn(1, query, generator=generator, dtype=torch.float64),
    }


def run(model, parts, mode, **kwargs):
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode=mode, inner_steps=1, inner_lr=0.1)
    return predict(model, parts, task, config, **kwargs)


# --- the identity that makes the discriminator meaningful ------------------

@pytest.mark.parametrize("support", [2, 3, 5])
def test_C_equals_T_plus_the_meta_correction(support):
    model = build()
    parts = make_parts(support, 7, seed=support)
    with torch.no_grad():
        t = run(model, parts, "T")
        c = run(model, parts, "C")
    assert torch.allclose(c["prediction"] - c["meta"], t["prediction"],
                          atol=1e-11), (
        "C must differ from T by exactly the meta-correction")
    assert float(c["meta"].abs().max()) > 1e-9, "the meta term must be active"


@pytest.mark.parametrize("support", [1, 2, 3, 5])
def test_the_split_transport_reproduces_the_incumbent(support):
    """`shrink*(level + transport(centered))` must equal `shrink*sum w r`."""
    from model.similarity_grammar import tanimoto
    model = build(1)
    parts = make_parts(support, 6, seed=20 + support)
    with torch.no_grad():
        t = run(model, parts, "T")
        residual = parts["support_y"] - t["support_zero"]
        shrink = model.transport.shrinkage(support, residual)
        scale = model.transport.similarity_scale
        weight = torch.softmax(scale * tanimoto(
            parts["query_fingerprint"], parts["support_fingerprint"]), -1)
        incumbent = shrink * torch.einsum("bqk,bk->bq", weight, residual)
    assert torch.allclose(t["level"] + t["transport"], incumbent, atol=1e-11)


# --- k = 0 ------------------------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
def test_k0_is_the_ordinary_zero_shot_path(mode):
    model = build(2)
    parts = make_parts(0, 8, seed=3)
    with torch.no_grad():
        output = run(model, parts, mode)
    assert torch.allclose(output["prediction"], output["zero_shot"], atol=0.0)
    for term in ("level", "transport", "meta"):
        assert torch.count_nonzero(output[term]) == 0
    assert output["inner_trace"] == []


@pytest.mark.parametrize("mode", MODES)
def test_k0_does_not_depend_on_the_support_placeholder(mode):
    """Changing an empty support's dtype/shape must not move a k=0 prediction."""
    model = build(2)
    parts = make_parts(0, 8, seed=4)
    with torch.no_grad():
        first = run(model, parts, mode)["prediction"]
        widened = dict(parts)
        widened["support_y"] = torch.zeros(1, 0, dtype=torch.float64)
        second = run(model, widened, mode)["prediction"]
    assert torch.allclose(first, second, atol=0.0)


# --- k = 1: the complementary residual is exactly zero ---------------------

def test_C_has_no_complementary_signal_at_k1():
    """With one support the level absorbs the whole residual, by construction.

    `centered = r - mean(r) = 0` and the leave-one-out transport has no other
    item, so `complementary == 0` exactly and the adapter is inert. This is why
    `C` makes no structural SAR claim at k=1.
    """
    model = build(3)
    parts = make_parts(1, 8, seed=5)
    with torch.no_grad():
        c = run(model, parts, "C")
        t = run(model, parts, "T")
    assert float(c["complementary"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    assert float(c["meta"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    assert torch.allclose(c["prediction"], t["prediction"], atol=1e-11)


def test_H_does_adapt_at_k1_where_C_does_not():
    """The arms must be distinguishable exactly where the hypothesis says."""
    model = build(3)
    parts = make_parts(1, 8, seed=6)
    with torch.no_grad():
        h = run(model, parts, "H")
        c = run(model, parts, "C")
    assert float(h["meta"].abs().max()) > 1e-9
    assert float(c["meta"].abs().max()) < 1e-12


# --- scope ------------------------------------------------------------------

def test_C_never_adapts_the_bias():
    """The level has an explicit term; adapting the bias would fit it twice."""
    assert ADAPTED_BY_MODE["C"] == ("interaction_head.2.weight",)
    assert "interaction_head.2.bias" in ADAPTED_BY_MODE["H"]


@pytest.mark.parametrize("mode", MODES)
def test_no_arm_mutates_the_persistent_model(mode):
    model = build(4)
    parts = make_parts(4, 6, seed=7)
    before = {k: v.clone() for k, v in model.state_dict().items()}
    run(model, parts, mode)
    for key, value in model.state_dict().items():
        assert torch.equal(value, before[key])


@pytest.mark.parametrize("mode", MODES)
def test_no_query_label_reaches_any_arm(mode):
    model = build(4)
    parts = make_parts(4, 6, seed=8)
    with torch.no_grad():
        first = run(model, parts, mode)["prediction"]
        shifted = dict(parts)
        shifted["query_y"] = parts["query_y"] * -5.0 + 2.0
        second = run(model, shifted, mode)["prediction"]
    assert torch.allclose(first, second, atol=0.0)


def test_M_disables_transport_entirely():
    model = build(5)
    parts = make_parts(4, 6, seed=9)
    with torch.no_grad():
        m = run(model, parts, "M")
    assert torch.count_nonzero(m["level"]) == 0
    assert torch.count_nonzero(m["transport"]) == 0
    assert float(m["meta"].abs().max()) > 1e-9


def test_T_runs_no_inner_loop():
    model = build(5)
    parts = make_parts(4, 6, seed=10)
    with torch.no_grad():
        t = run(model, parts, "T")
    assert t["inner_trace"] == []
    assert torch.count_nonzero(t["meta"]) == 0


# --- permutation symmetry ---------------------------------------------------

@pytest.mark.parametrize("mode", MODES)
def test_support_permutation_invariance(mode):
    model = build(6)
    parts = make_parts(5, 6, seed=11)
    order = torch.randperm(5, generator=torch.Generator().manual_seed(1))
    permuted = dict(parts)
    for key in ("support_atoms", "support_bonds", "support_mask",
                "support_fingerprint"):
        permuted[key] = parts[key][:, order]
    permuted["support_y"] = parts["support_y"][:, order]
    with torch.no_grad():
        assert torch.allclose(run(model, parts, mode)["prediction"],
                              run(model, permuted, mode)["prediction"],
                              atol=1e-10)


@pytest.mark.parametrize("mode", MODES)
def test_query_permutation_equivariance(mode):
    model = build(6)
    parts = make_parts(4, 7, seed=12)
    order = torch.randperm(7, generator=torch.Generator().manual_seed(2))
    permuted = dict(parts)
    for key in ("query_atoms", "query_bonds", "query_mask", "query_fingerprint"):
        permuted[key] = parts[key][:, order]
    with torch.no_grad():
        straight = run(model, parts, mode)["prediction"]
        shuffled = run(model, permuted, mode)["prediction"]
    assert torch.allclose(straight[:, order], shuffled, atol=1e-10)


# --- learned step sizes -----------------------------------------------------

def test_bias_only_restriction_is_exactly_a_level_shift():
    """The half of the parameter ablation that really is pure level."""
    model = build(8)
    parts = make_parts(4, 8, seed=40)
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode="H", inner_steps=1, inner_lr=0.3)
    with torch.no_grad():
        full = predict(model, parts, task, config)
        bias = predict(model, parts, task, config, keep="bias")
    shift = bias["meta"]
    assert float(shift.std()) == pytest.approx(0.0, abs=1e-11)
    assert float(full["meta"].std()) > 1e-9


def test_weight_only_restriction_still_carries_a_level_component():
    """Correction 4, on the real model: weight-only is not pure shape."""
    model = build(8)
    parts = make_parts(4, 8, seed=41)
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode="H", inner_steps=1, inner_lr=0.3)
    with torch.no_grad():
        weight = predict(model, parts, task, config, keep="weight")["meta"]
    assert abs(float(weight.mean())) > 1e-9, (
        "a weight-only update carries a non-zero episode mean")


def test_restrict_none_disables_the_correction():
    model = build(8)
    parts = make_parts(4, 8, seed=42)
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode="H", inner_steps=1, inner_lr=0.3)
    with torch.no_grad():
        none = predict(model, parts, task, config, keep="none")
    assert float(none["meta"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError):
        predict(model, parts, task, config, keep="elsewhere")


def test_C_bias_restriction_is_inert_because_C_never_adapts_the_bias():
    model = build(9)
    parts = make_parts(4, 8, seed=43)
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode="C", inner_steps=1, inner_lr=0.3)
    with torch.no_grad():
        bias = predict(model, parts, task, config, keep="bias")
        weight = predict(model, parts, task, config, keep="weight")
        full = predict(model, parts, task, config)
    assert float(bias["meta"].abs().max()) == pytest.approx(0.0, abs=1e-12)
    assert torch.allclose(weight["meta"], full["meta"], atol=1e-12)


def test_learned_step_starts_exactly_at_the_fixed_value():
    steps = InnerStepSizes(initial=0.1, max_step=0.5, dtype=torch.float64)
    assert float(steps.weight_step()) == pytest.approx(0.1, abs=1e-9)
    assert float(steps.bias_step()) == pytest.approx(0.1, abs=1e-9)


def test_learned_step_is_bounded_above_and_cannot_go_negative():
    """The bound that matters is the upper one; the lower limit is a safe zero.

    `sigmoid` underflows to exactly 0.0 at extreme negative logits, so the step
    can reach zero. That is benign — a zero step disables adaptation rather than
    destabilising it — whereas an unbounded step above would diverge, which is
    the failure mode the bound exists to prevent.
    """
    steps = InnerStepSizes(initial=0.1, max_step=0.5, dtype=torch.float64)
    with torch.no_grad():
        steps.raw_weight.fill_(1e6)
        steps.raw_bias.fill_(-1e6)
    assert 0.0 < float(steps.weight_step()) <= 0.5
    assert float(steps.weight_step()) == pytest.approx(0.5, abs=1e-6)
    assert float(steps.bias_step()) >= 0.0
    assert float(steps.bias_step()) < 1e-6


def test_learned_step_rejects_an_initial_outside_the_bound():
    with pytest.raises(ValueError):
        InnerStepSizes(initial=0.6, max_step=0.5)
    with pytest.raises(ValueError):
        InnerStepSizes(initial=0.0, max_step=0.5)


def test_learned_step_routes_weight_and_bias_separately():
    steps = InnerStepSizes(initial=0.1, max_step=0.5, dtype=torch.float64)
    with torch.no_grad():
        steps.raw_bias.fill_(-2.0)
    assert float(steps.for_parameter("interaction_head.2.weight")) == \
        pytest.approx(float(steps.weight_step()), abs=1e-12)
    assert float(steps.for_parameter("interaction_head.2.bias")) == \
        pytest.approx(float(steps.bias_step()), abs=1e-12)
    assert float(steps.bias_step()) != float(steps.weight_step())


def test_learned_step_is_used_by_the_inner_loop():
    model = build(7)
    parts = make_parts(4, 6, seed=13)
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode="H", inner_steps=1, inner_lr=0.1,
                              learned_step=True)
    slow = InnerStepSizes(initial=0.01, max_step=0.5, dtype=torch.float64)
    fast = InnerStepSizes(initial=0.4, max_step=0.5, dtype=torch.float64)
    with torch.no_grad():
        small = predict(model, parts, task, config, steps=slow)["meta"]
        large = predict(model, parts, task, config, steps=fast)["meta"]
    assert float(large.abs().mean()) > float(small.abs().mean())


def test_learned_step_gradient_reaches_its_parameters():
    model = build(7)
    parts = make_parts(4, 6, seed=14)
    task = encode_parts(model, parts)
    config = StageBAdaptation(mode="C", inner_steps=1, learned_step=True)
    steps = InnerStepSizes(initial=0.1, max_step=0.5, dtype=torch.float64)
    output = predict(model, parts, task, config, steps=steps)
    torch.nn.functional.mse_loss(output["prediction"], parts["query_y"]).backward()
    assert steps.raw_weight.grad is not None
    assert abs(float(steps.raw_weight.grad)) > 0.0
