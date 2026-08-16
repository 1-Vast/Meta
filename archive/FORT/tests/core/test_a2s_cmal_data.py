from __future__ import annotations

import json

import numpy as np
import pandas as pd

from research.a2s.a2s_cmal_data import (
    SUPPORT_K,
    _first_observation,
    assign_source_splits,
    attach_counterfactuals,
    audit_package,
    build_source_episodes,
)


def test_assign_source_splits_keeps_components_intact() -> None:
    sources = pd.DataFrame({
        "target_uid": [f"t{row}" for row in range(24)],
        "component_id": [row // 2 for row in range(24)],
    })
    split = assign_source_splits(sources, seed=7)
    assert set(split.meta_split) == {"meta_train", "meta_validation", "meta_test"}
    assert split.groupby("component_id").meta_split.nunique().max() == 1


def test_first_observation_uses_target_compound_composite_key() -> None:
    frame = pd.DataFrame({
        "target_uid": ["t0", "t1", "t0"],
        "compound_parent_uid": ["shared", "shared", "shared"],
        "document_uid": ["d0", "d1", "d2"],
        "document_year": [2018, 2019, 2020],
    })
    first = _first_observation(frame)
    assert set(zip(first.target_uid, first.compound_parent_uid)) == {
        ("t0", "shared"), ("t1", "shared")
    }
    assert first.loc[first.target_uid == "t0", "document_year"].item() == 2018


def _metadata() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    targets = []
    for target_row, target in enumerate(("t0", "t1", "t2")):
        targets.append({
            "target_uid": target,
            "component_id": target_row,
            "role": "source",
            "meta_split": "meta_train",
            "homology_warm": False,
        })
        for parent_row in range(12):
            year = 2018 if parent_row < 6 else 2020
            rows.append({
                "target_uid": target,
                "compound_parent_uid": f"p{target_row}_{parent_row}",
                "measurement_uid": f"m{target_row}_{parent_row}",
                "assay_context_uid": f"a{target_row}_{parent_row}",
                "document_uid": f"doc_{target_row}_{parent_row}",
                "document_year": year,
            })
    return pd.DataFrame(rows), pd.DataFrame(targets)


def test_source_episodes_are_nested_label_blind_and_time_ordered() -> None:
    metadata, targets = _metadata()
    episodes = build_source_episodes(
        metadata,
        targets,
        seed=11,
        draws=2,
        protocols=("ordered",),
        query_min=5,
        query_max=6,
    )
    assert set(episodes.k) == set(SUPPORT_K)
    assert (episodes.support_max_year < episodes.query_min_year).all()
    for row in episodes.itertuples():
        assert len(json.loads(row.support_measurement_uids)) == row.k
        assert len(json.loads(row.query_measurement_uids)) == len(
            json.loads(row.query_parent_uids)
        )
    for _, frame in episodes.groupby(["target_uid", "draw_id"]):
        support = {
            int(row.k): set(json.loads(row.support_parent_uids))
            for row in frame.itertuples()
        }
        assert support[1] <= support[3] <= support[5]
        for row in frame.itertuples():
            assert not set(json.loads(row.support_parent_uids)) & set(
                json.loads(row.query_parent_uids)
            )


def test_counterfactual_mappings_are_target_mismatched_and_query_blind() -> None:
    metadata, targets = _metadata()
    episodes = build_source_episodes(
        metadata,
        targets,
        seed=13,
        draws=2,
        protocols=("ordered",),
        query_min=5,
        query_max=6,
    )
    parents = sorted(metadata.compound_parent_uid.unique())
    parent_index = {parent: row for row, parent in enumerate(parents)}
    ecfp = np.zeros((len(parents), 8), dtype=np.float32)
    for parent, row in parent_index.items():
        target = int(parent.split("_")[0][1:])
        ecfp[row, target] = 1.0
        ecfp[row, 3 + (row % 5)] = 1.0
    protein = np.eye(3, dtype=np.float32)
    scaffolds = {parent: f"s{row % 2}" for row, parent in enumerate(parents)}
    mapped = attach_counterfactuals(
        episodes,
        parent_uids=parents,
        ecfp4=ecfp,
        scaffolds=scaffolds,
        target_uids=list(targets.target_uid),
        target_features=protein,
        seed=17,
        device="cpu",
    )
    id_target = dict(zip(mapped.episode_id, mapped.target_uid))
    for column in (
        "random_negative_episode_id",
        "protein_hard_negative_episode_id",
        "chemical_match_negative_episode_id",
    ):
        assert (mapped[column].map(id_target) != mapped.target_uid).all()

    by_id = mapped.set_index("episode_id")
    for row in mapped.itertuples():
        own = {scaffolds[parent] for parent in json.loads(row.support_parent_uids)}
        selected = by_id.loc[row.chemical_match_negative_episode_id]
        selected_set = {
            scaffolds[parent] for parent in json.loads(selected.support_parent_uids)
        }
        selected_jaccard = len(own & selected_set) / len(own | selected_set)
        candidates = mapped[
            (mapped.role == row.role)
            & (mapped.meta_split == row.meta_split)
            & (mapped.protocol == row.protocol)
            & (mapped.k == row.k)
            & (mapped.target_uid != row.target_uid)
        ]
        best = max(
            len(own & {scaffolds[parent] for parent in json.loads(value)})
            / len(own | {scaffolds[parent] for parent in json.loads(value)})
            for value in candidates.support_parent_uids
        )
        assert selected_jaccard == best
        assert row.chemical_match_scaffold_jaccard == best

    audit = audit_package(mapped, targets)
    assert audit["label_columns_read"] == []
    assert audit["violations"]["support_query_parent"] == 0
    assert audit["violations"]["ordered_time"] == 0
    assert not any(audit["violations"]["negative_same_target"].values())
