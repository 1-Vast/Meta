from __future__ import annotations

import numpy as np
import pandas as pd

from research.shared.residual_bilinear_probe import _bootstrap_gain, _strict_episodes


def test_strict_episode_closes_metadata_axes() -> None:
    frame = pd.DataFrame(
        [
            ("t1", "s1", "c1", "d1", "a1", 1, "h1"),
            ("t1", "s2", "c2", "d2", "a2", 2, "h1"),
            ("t1", "s3", "c3", "d3", "a3", 3, "h1"),
            ("t1", "s4", "c4", "d4", "a4", 4, "h1"),
            ("t1", "s5", "c5", "d5", "a5", 5, "h1"),
            ("t1", "s6", "c6", "d6", "a6", 6, "h1"),
            ("t1", "s7", "c7", "d7", "a7", 7, "h1"),
        ],
        columns=["target", "scaffold", "conn", "docs", "assays", "source_row", "hcluster"],
    )
    episodes, skipped = _strict_episodes(frame, support_size=5, query_cap=8)
    assert skipped == 0
    assert len(episodes) == 1
    episode = episodes[0]
    assert len(episode.support) == 5
    assert episode.query == (6, 7)


def test_component_bootstrap_gain_is_reproducible() -> None:
    arm = {"t1": {"rmse": 1.0}, "t2": {"rmse": 2.0}}
    reference = {"t1": {"rmse": 1.5}, "t2": {"rmse": 2.5}}
    components = {"t1": "h1", "t2": "h2"}
    first = _bootstrap_gain(arm, reference, components, "rmse", False, 1729, draws=64)
    second = _bootstrap_gain(arm, reference, components, "rmse", False, 1729, draws=64)
    assert first == second
    assert np.isclose(first["mean"], 0.5)
