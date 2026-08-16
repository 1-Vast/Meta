from __future__ import annotations

import numpy as np
import pandas as pd

from research.a2s import a2s_baseline
from research.a2s.a2s_baseline import _metrics, _units, build_episodes, role_targets


def _rows() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for index in range(8):
        rows.append(
            {
                "target": "R",
                "conn": f"C{index}",
                "endpoint": "pKi",
                "affinity": 5.0 + index / 10,
                "scaffold": f"S{index}",
                "assays": f"A{index}",
                "docs": f"D{index}",
                "accession": "P-R",
                "hcluster": "H-R",
                "dual_cold_split": "train",
            }
        )
    rows.append({**rows[0], "affinity": 9.0})
    return pd.DataFrame(rows)


def test_units_deduplicate_and_episode_has_disjoint_support_query() -> None:
    units = _units(_rows())
    assert len(units) == 8
    source, recipient, counts = role_targets(units)
    assert source == set()
    assert recipient == {"R"}
    assert counts["R"] == 8
    episodes = build_episodes(units, recipient, k=3, min_query=5)
    assert len(episodes) == 1
    episode = episodes[0]
    assert len(episode.support) == 3
    assert len(episode.query) == 5
    assert set(episode.support).isdisjoint(episode.query)


def test_units_preserves_global_registry_row_ids_after_filtering() -> None:
    frame = _rows().reset_index(names="source_row")
    filtered = frame.loc[frame.source_row >= 3].copy()
    units = _units(filtered)
    assert units.source_row.min() == 3
    assert units.source_row.max() == 8


def test_load_units_uses_global_registry_index(monkeypatch) -> None:
    source = _rows().set_index(pd.Index(range(100, 109)))
    indexed = source.reset_index(names="source_row")
    archive = {
        "feat": np.zeros((len(indexed), 2), dtype=np.float32),
        "conn_sha": np.asarray(a2s_baseline.connection_hash(indexed)),
    }
    monkeypatch.setattr(a2s_baseline.pd, "read_parquet", lambda *args, **kwargs: source)
    monkeypatch.setattr(a2s_baseline.np, "load", lambda *args, **kwargs: archive)
    units = a2s_baseline.load_units("pKi")
    assert units.source_row.min() == 100
    assert units.source_row.max() == 107


def test_metrics_reports_perfect_fit() -> None:
    result = _metrics(
        label=pd.Series([1.0, 2.0, 3.0]).to_numpy(),
        prediction=pd.Series([1.0, 2.0, 3.0]).to_numpy(),
    )
    assert result["rmse"] == 0.0
    assert result["mae"] == 0.0
    assert result["spearman"] == 1.0
    assert result["pairwise_accuracy"] == 1.0
