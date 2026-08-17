"""`A0` is the accepted baseline, proven rather than asserted.

The whole experiment rests on `A0` being the current accepted recipe with the
inner loop switched off. That claim was nearly false once already: the first
version of `train_meta` silently dropped three auxiliary terms
(`support_match`, the binding contrastive, the protein contrast), which would
have turned any `A1` gain into "the inner loop recovers what we deleted".

This test composes the production episode loss exactly as
`scripts/train_qpsmp.py::train` does and compares it against the Stage A loss
at `inner_steps=0` on real episodes. Needs the corpus, so it is `slow`.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys

import numpy as np
import pytest
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_SLOW") != "1",
    reason="needs the governed corpus; set RUN_SLOW=1")

from scripts.qpsmp_data import QPSMPData                             # noqa: E402
from scripts.train_qpsmp import (                                    # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    batch_counterfactual_episode, binding_contrastive_loss,
    centered_task_error, compact_episode, counterfactual_label_assignments,
    forward, normalized_episode, ranking_term, resolve_architecture,
    training_label_scale, wrong_protein_zero_shot,
)
from tools.research.stageA_innerloop.inner_loop import AdaptationConfig  # noqa: E402
from tools.research.stageA_innerloop.train_meta import (             # noqa: E402
    align_atoms, auxiliary_losses, cache_donor, encode_parts,
    episode_tensors, predict,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"


@pytest.fixture(scope="module")
def setup():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    config = TrainConfig(arch="similarity_only", device="cpu", amp=False,
                         split_directory=str(SPLIT))
    model = resolve_architecture(config.arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks,
        adapter_scale=config.adapter_scale, dtype=torch.float32)
    model.eval()
    return data, config, model, training_label_scale(data)


def production_loss(model, data, episode, config, label_scale) -> dict:
    """The accepted recipe's episode loss, composed exactly as `train` does."""
    full = forward(model, episode)
    query_y = episode.query_y.to(device=full.prediction.device,
                                 dtype=full.prediction.dtype)
    support_size = episode.support_atoms.shape[0]
    loss_full = F.smooth_l1_loss(full.prediction, query_y)
    loss_zero = F.smooth_l1_loss(full.zero_shot, query_y)
    loss_rank = ranking_term(full.prediction, query_y, config, label_scale)
    shape = centered_task_error(full.prediction, query_y)
    binding = loss_full.new_zeros(())
    protein = loss_full.new_zeros(())
    if support_size > 0:
        errors = [(full.prediction - query_y).square().mean()]
        assignments = counterfactual_label_assignments(full, episode.support_y)
        wrong_episode = batch_counterfactual_episode(episode, assignments)
        wrong = forward(model, wrong_episode)
        wrong_truth = wrong_episode.query_y.to(device=wrong.prediction.device,
                                               dtype=wrong.prediction.dtype)
        errors.extend((wrong.prediction - wrong_truth).square().mean(-1).unbind())
        binding = binding_contrastive_loss(errors, config.binding_temperature)
        wrong_zero = wrong_protein_zero_shot(
            model, data, episode, episode.spec.donor_target)
        protein = binding_contrastive_loss(
            [(full.zero_shot - query_y).square().mean(),
             (wrong_zero - query_y).square().mean()],
            config.binding_temperature)
    return {"full": float(loss_full), "zero": float(loss_zero),
            "rank": float(loss_rank), "shape": float(shape),
            "support_match": float(full.support_match_loss),
            "binding": float(binding), "protein": float(protein)}


def stage_a_loss(model, data, episode, config, label_scale) -> dict:
    parts = align_atoms(episode_tensors(model, episode, "cpu", torch.float32))
    parts["_donor_key"] = cache_donor(data, episode.spec.donor_target, "cpu",
                                      torch.float32)
    adaptation = AdaptationConfig(inner_steps=0)
    task = encode_parts(model, parts)
    output = predict(model, parts, adaptation, task=task)
    query_y = parts["query_y"]
    support_match, binding, protein = auxiliary_losses(
        model, parts, task, output, adaptation, config, label_scale)
    return {"full": float(F.smooth_l1_loss(output["prediction"], query_y)),
            "zero": float(F.smooth_l1_loss(output["pre_adaptation_query"],
                                           query_y)),
            "rank": float(ranking_term(output["prediction"], query_y, config,
                                       label_scale)),
            "shape": float(centered_task_error(output["prediction"], query_y)),
            "support_match": float(support_match),
            "binding": float(binding), "protein": float(protein)}


@pytest.mark.parametrize("support_size", [0, 1, 2, 3, 5])
def test_stage_a_at_zero_inner_steps_matches_the_production_loss(
        setup, support_size):
    data, config, model, label_scale = setup
    rng = np.random.default_rng(20260817 + support_size)
    spec = data.draw_episode("meta_train", support_size, 12, rng)
    episode = compact_episode(normalized_episode(
        data.materialize(spec), label_scale))
    with torch.no_grad():
        theirs = production_loss(model, data, episode, config, label_scale)
        mine = stage_a_loss(model, data, episode, config, label_scale)
    for term in theirs:
        assert mine[term] == pytest.approx(theirs[term], abs=2e-5), (
            f"term {term!r} diverges at k={support_size}: "
            f"stage A {mine[term]} vs production {theirs[term]}")


def test_every_production_term_is_present_in_stage_a(setup):
    """A term that is identically zero in both is not evidence of a match."""
    data, config, model, label_scale = setup
    rng = np.random.default_rng(7)
    spec = data.draw_episode("meta_train", 5, 12, rng)
    episode = compact_episode(normalized_episode(
        data.materialize(spec), label_scale))
    with torch.no_grad():
        mine = stage_a_loss(model, data, episode, config, label_scale)
    for term in ("full", "zero", "rank", "shape", "support_match", "binding",
                 "protein"):
        assert abs(mine[term]) > 0.0, f"{term} is inert and proves nothing"
