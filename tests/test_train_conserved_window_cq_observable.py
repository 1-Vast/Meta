import numpy as np

from research.crossed_interaction.train_conserved_window_cq_observable import (
    build_conserved_window_context,
    conserved_window_descriptor,
    selected_window_descriptor,
)


def test_selected_window_descriptor_concatenates_contrast_and_centroid():
    slots = np.asarray([[1.0], [5.0], [9.0]])
    selected = np.asarray([0, 2])
    centroid = np.asarray([[2.0], [7.0]])
    descriptor = selected_window_descriptor(
        slots, selected, centroid, mode="family_conserved")
    assert np.array_equal(descriptor, np.asarray([-1.0, 2.0, 2.0, 7.0]))


def test_selected_window_descriptor_raw_mode_uses_windows_only():
    slots = np.asarray([[1.0], [5.0], [9.0]])
    selected = np.asarray([0, 2])
    centroid = np.asarray([[2.0], [7.0]])
    assert np.array_equal(
        selected_window_descriptor(slots, selected, centroid, mode="raw_windows"),
        np.asarray([1.0, 9.0]))


def test_build_conserved_window_context_uses_train_only_family_windows():
    cells = {
        "a": {"target_id": "t1", "protein_group_40": "g1"},
        "b": {"target_id": "t2", "protein_group_40": "g1"},
        "c": {"target_id": "t3", "protein_group_40": "g1"},
    }
    panels = [
        {"split": "train", "cell_ids": ["a", "b"]},
        {"split": "development", "cell_ids": ["c"]},
    ]
    protein_slots = {
        "t1": np.asarray([[1.0], [10.0], [100.0]]),
        "t2": np.asarray([[1.0], [12.0], [300.0]]),
        "t3": np.asarray([[99.0], [99.0], [99.0]]),
    }
    context, metadata = build_conserved_window_context(
        cells, panels, protein_slots, top_windows=1)
    assert context.family_slots["g1"].tolist() == [0]
    assert np.array_equal(context.family_centroids["g1"], np.asarray([[1.0]]))
    assert metadata["groups_with_conserved_windows"] == 1


def test_conserved_window_descriptor_falls_back_to_global_windows():
    cells = {
        "a": {"target_id": "t1", "protein_group_40": "g1"},
        "b": {"target_id": "t2", "protein_group_40": "g2"},
    }
    panels = [{"split": "train", "cell_ids": ["a"]}]
    protein_slots = {
        "t1": np.asarray([[1.0], [3.0]]),
        "t2": np.asarray([[5.0], [7.0]]),
    }
    context, _ = build_conserved_window_context(
        cells, panels, protein_slots, top_windows=1, min_train_family_size=2)
    descriptor = conserved_window_descriptor("t2", protein_slots, context)
    assert descriptor.shape == (2,)
