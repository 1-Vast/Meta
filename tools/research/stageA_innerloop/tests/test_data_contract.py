"""Data-side contract for Stage A: episodes, banks, and selector provenance.

These need the real corpus, so they are marked `slow` and are skipped unless
`RUN_SLOW=1`. The properties they pin are the ones a synthetic tensor cannot
express: that the governed sampler really does give one target per task, that
the evaluation banks really are nested across k, and that the A2 selector
really does read `meta_train` only.
"""
from __future__ import annotations

import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_SLOW") != "1",
    reason="needs the governed corpus; set RUN_SLOW=1")

from scripts.qpsmp_data import QPSMPData                             # noqa: E402
from scripts.train_qpsmp import (                                    # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
EVALUATION_SEED = 73101


@pytest.fixture(scope="module")
def data():
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)


def test_meta_test_is_withheld_and_unreachable(data):
    record = data.seal_record()
    assert record["included"] is False
    assert record["evaluated"] is False
    assert record["sealed_cells_withheld"] > 0
    assert "meta_test" not in data.tasks
    assert all(cell["split"] != "meta_test" for cell in data.cells)


def test_every_episode_has_exactly_one_target_with_disjoint_ligands(data):
    rng = np.random.default_rng(20260817)
    for _ in range(300):
        k = int(rng.integers(0, 6))
        spec = data.draw_episode("meta_train", k, 16, rng)
        rows = [data.cells[i] for i in (*spec.support, *spec.query)]
        assert len({row["target_id"] for row in rows}) == 1
        assert not set(spec.support) & set(spec.query)
        support = [data.cells[i]["ligand_id"] for i in spec.support]
        query = [data.cells[i]["ligand_id"] for i in spec.query]
        assert not set(support) & set(query)
        assert len(set(support)) == len(support)
        assert len(set(query)) == len(query)
        # The donor is a legal wrong protein: a different homology component.
        assert spec.donor_target != spec.target


def test_materialize_refuses_an_overlapping_episode(data):
    from scripts.qpsmp_data import EpisodeSpec
    rng = np.random.default_rng(1)
    spec = data.draw_episode("meta_train", 3, 8, rng)
    broken = EpisodeSpec(spec.split, spec.component, spec.target,
                         spec.support, (spec.support[0], *spec.query),
                         spec.donor_target)
    with pytest.raises(ValueError):
        data.materialize(broken)


def test_evaluation_banks_are_nested_across_support_sizes(data):
    banks = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, 16, 1, EVALUATION_SEED, None)
    index = {k: {(e.target, e.query): e for e in banks[k]} for k in SUPPORT_SIZES}
    keys = set(index[0])
    assert keys, "the nested bank is empty"
    for k in SUPPORT_SIZES:
        assert set(index[k]) == keys, f"the query panel moved at k={k}"
        assert len(banks[k]) == len(banks[0])
    for key in keys:
        prefixes = [index[k][key].support for k in SUPPORT_SIZES]
        for smaller, larger in zip(prefixes, prefixes[1:]):
            assert larger[:len(smaller)] == smaller
        for k in SUPPORT_SIZES:
            episode = index[k][key]
            support = {data.cells[i]["ligand_id"] for i in episode.support}
            query = {data.cells[i]["ligand_id"] for i in episode.query}
            assert not support & query


def test_episode_banks_are_reproducible_under_the_fixed_seed(data):
    first = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, 16, 1, EVALUATION_SEED, None)
    second = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, 16, 1, EVALUATION_SEED, None)
    assert first == second
    other = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, 16, 1, EVALUATION_SEED + 1, None)
    assert other != first


def test_training_draws_and_the_task_selector_read_meta_train_only(data,
                                                                  monkeypatch):
    """The A2 selector must never see a meta_val task.

    Enforced by recording every split the trainer asks the sampler for during a
    short run with selection enabled.
    """
    import torch

    from scripts.train_qpsmp import TrainConfig
    from tools.research.stageA_innerloop import train_meta

    asked: list[str] = []
    original = QPSMPData.draw_episode

    def spy(self, split, *args, **kwargs):
        asked.append(split)
        return original(self, split, *args, **kwargs)

    monkeypatch.setattr(QPSMPData, "draw_episode", spy)
    meta = train_meta.MetaConfig(
        base=TrainConfig(arch="similarity_only", steps=2, seed=1,
                         split_directory=str(SPLIT), device="cpu", amp=False,
                         episodes_per_step=1, val_interval=2,
                         eval_targets_per_component=1, query_size=6,
                         min_query_size=2, hidden_dim=24, task_dim=12,
                         ligand_layers=2, pair_dim=24, pair_latents=6,
                         pair_heads=2),
        inner_steps=1, inner_lr=0.1, task_selection=True)
    with torch.no_grad():
        pass
    train_meta.train(data, meta)
    assert asked, "the trainer drew no episodes"
    assert set(asked) == {"meta_train"}
