"""Structural contract for the shared discriminator (ladder rung M2).

These pin the claims the harness makes about itself, on the real model classes,
in float64, with no dataset. The load-bearing one is
`test_captured_state_pooling_reproduces_the_models_own_pooling`: the extractor
omits `state_mean`/`state_max` from the cache on the grounds that they are
bitwise identical to `mean_state`/`max_state`, and if that identity were ever
false the cache would be silently missing two representations.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[5]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM       # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel           # noqa: E402
from tools.research.method_ladder._shared import _frozen              # noqa: E402
from tools.research.method_ladder._shared.capture import (            # noqa: E402
    Capture, POOLINGS, extract, pool,
)
from tools.research.method_ladder._shared.panels import (             # noqa: E402
    centered, permuted_protein_assignment, within_target_r,
)
from tools.research.method_ladder._shared.probe import (              # noqa: E402
    LinearProbe, centered_batch_loss, component_folds, evaluate_probe,
    fit_standardizer, train_probe,
)

PROTEIN_DIM, SLOTS, ATOMS, LIGANDS = 32, 12, 9, 7


def build(seed: int = 0) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12, ligand_layers=2,
        pair_dim=24, pair_latents=6, pair_heads=2,
        use_learned_key=False).double()


def protein(seed: int):
    generator = torch.Generator().manual_seed(seed)
    return (torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64),
            torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64),
            torch.ones(1, SLOTS, dtype=torch.float64),
            torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64))


def ligands(count: int, seed: int, padded: int = 0):
    generator = torch.Generator().manual_seed(seed)
    total = ATOMS + padded
    atoms = torch.rand(1, count, total, ATOM_FEAT_DIM, generator=generator,
                       dtype=torch.float64)
    bonds = torch.rand(1, count, total, total, BOND_FEAT_DIM,
                       generator=generator, dtype=torch.float64)
    bonds = bonds * (bonds > 0.7)
    bonds = 0.5 * (bonds + bonds.transpose(2, 3))
    mask = torch.ones(1, count, total, dtype=torch.float64)
    if padded:
        mask[:, :, ATOMS:] = 0.0
    fingerprint = (torch.rand(1, count, 64, generator=generator,
                              dtype=torch.float64) > 0.7).double()
    return atoms, bonds, mask, fingerprint


# --- the identity the cache depends on -------------------------------------

@pytest.mark.parametrize("padded", [0, 5])
def test_captured_state_pooling_reproduces_the_models_own_pooling(padded):
    """`state_mean`/`state_max` are omitted from the cache as exact duplicates.

    `ContactGrammar` masks `state` with a 0/1 atom mask and then pools it, and
    `pool` applies the same gate, so the two must agree to float64 precision.
    If this ever fails the extractor is dropping real information.
    """
    model = build()
    atoms, bonds, mask, fingerprint = ligands(LIGANDS, seed=11, padded=padded)
    capture = Capture(model)
    try:
        with torch.no_grad():
            pooled, tokens, protein_mask, chemistry = protein(1)
            model(pooled, tokens, protein_mask,
                  atoms[:, :0], bonds[:, :0], mask[:, :0],
                  torch.zeros(1, 0, dtype=torch.float64),
                  atoms, bonds, mask, adapt=False, protein_chemistry=chemistry,
                  support_fingerprint=fingerprint[:, :0],
                  query_fingerprint=fingerprint)
        values, atom_mask = dict(capture.values), capture.atom_mask
    finally:
        capture.close()
    assert torch.allclose(pool(values["state"], atom_mask, "mean"),
                          values["mean_state"], atol=1e-12)
    assert torch.allclose(pool(values["state"], atom_mask, "max"),
                          values["max_state"], atol=1e-12)


def test_the_context_slice_is_the_attention_readout():
    """The middle third of the fusion input is `context`, not `atoms`.

    Verified structurally: the fusion input is `cat((a, c, a*c))`, so the last
    third must equal the elementwise product of the first two.
    """
    model = build(2)
    atoms, bonds, mask, fingerprint = ligands(LIGANDS, seed=12)
    seen = {}

    def spy(module, inputs, output):
        seen["fused"] = inputs[0].detach()

    handle = model.grammar.atom_context.register_forward_hook(spy)
    capture = Capture(model)
    try:
        with torch.no_grad():
            pooled, tokens, protein_mask, chemistry = protein(3)
            model(pooled, tokens, protein_mask,
                  atoms[:, :0], bonds[:, :0], mask[:, :0],
                  torch.zeros(1, 0, dtype=torch.float64),
                  atoms, bonds, mask, adapt=False, protein_chemistry=chemistry,
                  support_fingerprint=fingerprint[:, :0],
                  query_fingerprint=fingerprint)
        width = capture.hidden
        fused = seen["fused"]
        assert torch.allclose(fused[..., :width] * fused[..., width:2 * width],
                              fused[..., 2 * width:], atol=1e-12)
        assert torch.allclose(capture.values["context"],
                              fused[..., width:2 * width], atol=1e-12)
    finally:
        capture.close()
        handle.remove()


def test_context_actually_depends_on_the_protein():
    """A positive control: if this were flat, every downstream null is vacuous."""
    model = build(4)
    parts = ligands(LIGANDS, seed=13)

    def run(protein_seed: int):
        capture = Capture(model)
        try:
            with torch.no_grad():
                pooled, tokens, protein_mask, chemistry = protein(protein_seed)
                atoms, bonds, mask, fingerprint = parts
                model(pooled, tokens, protein_mask,
                      atoms[:, :0], bonds[:, :0], mask[:, :0],
                      torch.zeros(1, 0, dtype=torch.float64),
                      atoms, bonds, mask, adapt=False,
                      protein_chemistry=chemistry,
                      support_fingerprint=fingerprint[:, :0],
                      query_fingerprint=fingerprint)
            return capture.values["context"].clone()
        finally:
            capture.close()

    assert (run(1) - run(2)).abs().max() > 1e-6


def test_pooling_ignores_padded_atoms():
    values = torch.arange(2 * 4 * 3, dtype=torch.float64).reshape(2, 4, 3)
    mask = torch.tensor([[1.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]])
    assert torch.allclose(pool(values, mask, "mean")[0],
                          values[0, :2].mean(0), atol=1e-12)
    assert torch.allclose(pool(values, mask, "mean")[1], values[1, 0], atol=1e-12)
    assert torch.allclose(pool(values, mask, "max")[1], values[1, 0], atol=1e-12)
    assert torch.allclose(pool(values, mask, "rms")[1],
                          values[1, 0].abs(), atol=1e-12)


def test_extract_reads_no_label_and_uses_the_k0_path():
    """The features must be exactly the zero-shot path: empty support, adapt off."""
    model = build(5)
    atoms, bonds, mask, fingerprint = ligands(LIGANDS, seed=14)
    features = extract(model, protein(6), atoms, bonds, mask, fingerprint)
    assert features["ligand"].shape == (LIGANDS, 24)
    assert features["context_mean"].shape == (LIGANDS, 24)
    assert features["interaction"].shape == (LIGANDS, 1)
    for name in ("context_mean", "context_max", "context_rms", "state_rms"):
        assert np.isfinite(features[name]).all()
    # The endpoint the model reports at k=0 must equal the sum of its parts.
    total = (features["_ligand_value"] + features["_protein_value"]
             + features["interaction"])
    assert np.allclose(total, features["_endpoint"], atol=1e-9)


def test_the_protein_branch_is_constant_within_a_panel():
    """`protein_value` cannot order ligands; the centered target proves it.

    This is why every ordering number in the ladder is centered: the additive
    protein head is exactly the part centering removes.
    """
    model = build(7)
    atoms, bonds, mask, fingerprint = ligands(LIGANDS, seed=15)
    features = extract(model, protein(8), atoms, bonds, mask, fingerprint)
    # Exactly constant in exact arithmetic — `protein_value` is
    # `protein_head(wide_summary)` and `wide_summary` is one vector expanded
    # over the panel. The residual spread is float32 cancellation in
    # `extract`'s `additive - ligand_value`, ~2e-08 on a value of ~0.24, so the
    # tolerance is set at the float32 noise floor rather than at zero.
    assert features["_protein_value"].std() < 1e-6
    assert abs(centered(features["_protein_value"].ravel())).max() < 1e-6


# --- the statistics --------------------------------------------------------

def test_within_target_r_is_blind_to_level_and_positive_scale():
    truth = np.array([1.0, 2.0, 3.0, 4.5, 2.5, 0.5])
    prediction = np.array([0.4, 1.1, 2.9, 3.2, 1.7, 0.1])
    base = within_target_r(prediction, truth)
    assert within_target_r(prediction + 12.0, truth) == pytest.approx(base, abs=1e-12)
    assert within_target_r(prediction, truth - 7.0) == pytest.approx(base, abs=1e-12)
    assert within_target_r(3.0 * prediction, truth) == pytest.approx(base, abs=1e-12)
    assert within_target_r(-prediction, truth) == pytest.approx(-base, abs=1e-12)


def test_a_constant_prediction_scores_zero_not_nan():
    truth = np.array([1.0, 2.0, 3.0, 4.0])
    assert within_target_r(np.full(4, 2.5), truth) == 0.0


def test_centered_batch_loss_is_level_invariant():
    prediction = torch.tensor([0.3, 1.2, -0.4, 2.0], dtype=torch.float64)
    truth = torch.tensor([1.0, 2.0, 0.5, 3.0], dtype=torch.float64)
    base = float(centered_batch_loss(prediction, truth))
    assert float(centered_batch_loss(prediction + 5.0, truth)) == pytest.approx(base, abs=1e-12)
    assert float(centered_batch_loss(prediction, truth - 3.0)) == pytest.approx(base, abs=1e-12)


def test_a_one_ligand_panel_is_inert_rather_than_nan():
    one = torch.tensor([1.5], dtype=torch.float64)
    assert float(centered_batch_loss(one, one * 2)) == 0.0


def test_component_folds_are_disjoint_and_cover_every_component():
    components = [f"c{i // 3}" for i in range(30)]
    buckets = component_folds(components)
    flat = [c for bucket in buckets for c in bucket]
    assert len(flat) == len(set(flat)) == len(set(components))
    assert set(flat) == set(components)


def test_probe_training_is_bitwise_reproducible():
    rng = np.random.default_rng(0)
    blocks = [(rng.normal(size=(8, 5)), rng.normal(size=8)) for _ in range(64)]
    first = train_probe(blocks, 5, 1e-3, seed=7)
    second = train_probe(blocks, 5, 1e-3, seed=7)
    assert torch.allclose(first.weight.weight, second.weight.weight, atol=0.0)


def test_the_linear_probe_seed_only_moves_the_minibatch_sampler():
    """Documented, because it changes what "three probe seeds" means.

    `LinearProbe` starts from an exact zero init, so the seed enters only
    through `train_probe`'s minibatch draw. With at most `PROBE_BATCH_TARGETS`
    panels the draw is the full set and the seed does nothing — the fit is
    deterministic full-batch descent. Seeds separate runs only once the panel
    count exceeds the batch size, which is the regime the ladder actually uses
    (346 meta_train panels against a batch of 32).
    """
    rng = np.random.default_rng(0)
    small = [(rng.normal(size=(8, 5)), rng.normal(size=8)) for _ in range(6)]
    assert torch.allclose(train_probe(small, 5, 1e-3, seed=7).weight.weight,
                          train_probe(small, 5, 1e-3, seed=8).weight.weight,
                          atol=0.0)
    large = [(rng.normal(size=(8, 5)), rng.normal(size=8)) for _ in range(64)]
    assert not torch.allclose(train_probe(large, 5, 1e-3, seed=7).weight.weight,
                              train_probe(large, 5, 1e-3, seed=8).weight.weight)


def test_a_probe_recovers_a_planted_within_target_direction():
    """Positive control for the discriminator itself.

    If the probe could not find a signal that is present by construction, every
    null it reports would be uninterpretable.
    """
    rng = np.random.default_rng(3)
    direction = rng.normal(size=6)
    blocks = []
    for _ in range(40):
        features = rng.normal(size=(10, 6))
        labels = features @ direction + 5.0 * rng.normal()   # per-panel level
        blocks.append((features, labels))
    probe = train_probe(blocks, 6, 1e-4, seed=1)
    scores = evaluate_probe(probe, blocks)
    assert float(np.mean(scores)) > 0.9


def test_a_probe_on_pure_noise_scores_near_zero():
    rng = np.random.default_rng(4)
    train = [(rng.normal(size=(10, 6)), rng.normal(size=10)) for _ in range(40)]
    held = [(rng.normal(size=(10, 6)), rng.normal(size=10)) for _ in range(40)]
    probe = train_probe(train, 6, 1e-2, seed=1)
    assert abs(float(np.mean(evaluate_probe(probe, held)))) < 0.15


def test_standardizer_uses_only_the_blocks_it_was_given():
    train = [np.array([[0.0, 10.0], [2.0, 20.0]])]
    scaler = fit_standardizer(train)
    assert np.allclose(scaler.mean, [1.0, 15.0])
    assert np.allclose(scaler.apply(np.array([[1.0, 15.0]])), [[0.0, 0.0]])


def test_a_constant_feature_column_does_not_divide_by_zero():
    scaler = fit_standardizer([np.array([[3.0, 1.0], [3.0, 2.0]])])
    out = scaler.apply(np.array([[3.0, 1.5]]))
    assert np.isfinite(out).all() and out[0, 0] == 0.0


def test_protein_permutation_is_a_derangement():
    class Fake:
        def __init__(self, target):
            self.target = target
    panels = tuple(Fake(f"t{i}") for i in range(25))
    mapping = permuted_protein_assignment(panels, _frozen.CONTROL_SEED)
    assert set(mapping) == {p.target for p in panels}
    assert set(mapping.values()) == set(mapping)
    assert all(source != destination for source, destination in mapping.items())


def test_the_frozen_manifest_is_serializable_and_carries_thresholds():
    import json
    manifest = _frozen.frozen_manifest()
    json.dumps(manifest)
    assert manifest["thresholds"]["smallest_effect_of_interest_r"] == 0.05
    assert "no ridge" in manifest["probe"]["solver"]


@pytest.mark.parametrize("interval,expected", [
    ({"mean": 0.20, "lo": 0.10, "hi": 0.30}, "RESOLVED"),
    ({"mean": 0.01, "lo": 0.005, "hi": 0.02}, "RESOLVED_NEGLIGIBLE"),
    ({"mean": 0.00, "lo": -0.02, "hi": 0.02}, "DECISIVE_NULL"),
    ({"mean": 0.05, "lo": -0.30, "hi": 0.40}, "UNDERPOWERED"),
])
def test_verdict_vocabulary(interval, expected):
    assert _frozen.verdict(interval) == expected


def test_poolings_are_the_three_documented_ones():
    assert POOLINGS == ("mean", "max", "rms")
    with pytest.raises(ValueError):
        pool(torch.zeros(1, 2, 3), torch.ones(1, 2), "median")
