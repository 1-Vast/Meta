import numpy as np
import torch

from research.meta_fewshot.train_qpsmp_core import (
    QPSMPTrainConfig,
    level_rows,
    stable_seed,
)


def test_stable_seed_is_deterministic_and_sensitive():
    assert stable_seed("a", 1, "b") == stable_seed("a", 1, "b")
    assert stable_seed("a", 1, "b") != stable_seed("a", 2, "b")


def test_level_rows_use_support_mean_without_query_labels():
    cells = [
        {"cell_id": f"c{i}", "target_id": "t", "pK": float(i)}
        for i in range(8)
    ]
    tensors = {"y": torch.arange(8, dtype=torch.float32)}
    tasks = {"t": np.arange(8)}

    rows = level_rows(
        cells, tensors, tasks, seed=5, k=2, draws=1,
        max_query=3, max_support_k=5)

    assert len(rows) == 3
    assert {row["arm"] for row in rows} == {"level"}
    assert {row["target_id"] for row in rows} == {"t"}
    assert all("pK" not in row for row in rows)


def test_qpsmp_config_defaults_are_small_smoke_budget():
    config = QPSMPTrainConfig()

    assert config.steps == 300
    assert config.test_draws == 3
    assert config.ridge > 0
    assert config.task_dim <= max(config.support_sizes)
