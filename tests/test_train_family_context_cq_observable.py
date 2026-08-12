import numpy as np

from research.crossed_interaction.train_family_context_cq_observable import (
    build_family_context,
    family_context_descriptor,
    pooled_slot_descriptor,
    train_targets_from_panels,
)


def test_pooled_slot_descriptor_uses_only_masked_slots():
    blocks = np.asarray([[1.0, 3.0], [100.0, 100.0], [5.0, 7.0]])
    mask = np.asarray([True, False, True])
    assert np.array_equal(pooled_slot_descriptor(blocks, mask), np.asarray([3.0, 5.0]))


def test_train_targets_from_panels_ignores_development_panels():
    cells = {
        "a": {"target_id": "train_target"},
        "b": {"target_id": "development_target"},
    }
    panels = [
        {"split": "train", "cell_ids": ["a"]},
        {"split": "development", "cell_ids": ["b"]},
    ]
    assert train_targets_from_panels(cells, panels) == {"train_target"}


def test_build_family_context_uses_train_only_centroids_and_fallback():
    cells = {
        "a": {"target_id": "t1", "protein_group_40": "g1"},
        "b": {"target_id": "t2", "protein_group_40": "g1"},
        "c": {"target_id": "t3", "protein_group_40": "g1"},
        "d": {"target_id": "t4", "protein_group_40": "g2"},
    }
    panels = [
        {"split": "train", "cell_ids": ["a", "b"]},
        {"split": "development", "cell_ids": ["c", "d"]},
    ]
    protein_descriptors = {
        "t1": np.asarray([1.0, 1.0]),
        "t2": np.asarray([3.0, 3.0]),
        "t3": np.asarray([100.0, 100.0]),
        "t4": np.asarray([5.0, 7.0]),
    }
    context, metadata = build_family_context(
        cells, panels, protein_descriptors, min_train_family_size=2)
    assert np.array_equal(context.family_centroids["g1"], np.asarray([2.0, 2.0]))
    assert np.array_equal(context.global_centroid, np.asarray([2.0, 2.0]))
    assert metadata["groups_with_train_centroid"] == 1
    assert metadata["targets_using_global_centroid"] == 1


def test_family_context_descriptor_concatenates_contrast_and_centroid():
    cells = {
        "a": {"target_id": "t1", "protein_group_40": "g1"},
        "b": {"target_id": "t2", "protein_group_40": "g1"},
    }
    panels = [{"split": "train", "cell_ids": ["a", "b"]}]
    protein_descriptors = {
        "t1": np.asarray([1.0, 1.0]),
        "t2": np.asarray([3.0, 3.0]),
    }
    context, _ = build_family_context(cells, panels, protein_descriptors)
    descriptor = family_context_descriptor("t1", protein_descriptors, context)
    assert np.array_equal(descriptor, np.asarray([-1.0, -1.0, 2.0, 2.0]))
