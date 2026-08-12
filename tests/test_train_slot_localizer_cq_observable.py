import numpy as np

from research.crossed_interaction.train_slot_localizer_cq_observable import (
    learn_slot_localizer,
    protein_slot_blocks,
)


def test_protein_slot_blocks_preserves_slot_axis_and_blocks_hidden_dim():
    residues = np.arange(4 * 6, dtype=np.float64).reshape(4, 6)
    mask = np.asarray([1, 0, 1, 1], dtype=bool)
    blocks = protein_slot_blocks(residues, mask, hidden_blocks=3)
    assert blocks.shape == (4, 3)
    assert np.array_equal(blocks[1], np.zeros(3))
    assert np.isfinite(blocks).all()


def test_learn_slot_localizer_selects_train_quotient_interaction_slot():
    cells = {
        "c11": {"target_id": "t1", "ligand_id": "l1", "pK": 1.0},
        "c12": {"target_id": "t1", "ligand_id": "l2", "pK": -1.0},
        "c21": {"target_id": "t2", "ligand_id": "l1", "pK": -1.0},
        "c22": {"target_id": "t2", "ligand_id": "l2", "pK": 1.0},
        "dev": {"target_id": "t1", "ligand_id": "l1", "pK": 0.0},
    }
    panels = [
        {
            "panel_id": "train_panel",
            "split": "train",
            "cell_ids": ["c11", "c12", "c21", "c22"],
        },
        {
            "panel_id": "dev_panel",
            "split": "development",
            "cell_ids": ["dev"],
        },
    ]
    proteins = {
        "t1": np.asarray([[1.0], [0.0], [0.0]]),
        "t2": np.asarray([[-1.0], [0.0], [0.0]]),
    }
    ligands = {
        "l1": np.asarray([1.0]),
        "l2": np.asarray([-1.0]),
    }

    localizer = learn_slot_localizer(
        cells, panels, proteins, ligands, top_slots=1, seed=1)
    assert localizer["mode"] == "supervised"
    assert localizer["selected_slots"] == [0]
    assert localizer["train_panels_used"] == 1


def test_uniform_localizer_is_deterministic_without_training_panels():
    proteins = {"t": np.zeros((8, 1), dtype=np.float64)}
    localizer = learn_slot_localizer(
        {}, [], proteins, {}, top_slots=4, mode="uniform")
    assert localizer["mode"] == "uniform"
    assert localizer["selected_slots"] == [0, 2, 4, 7]
    assert localizer["train_panels_used"] == 0


def test_learn_slot_localizer_rejects_too_many_slots():
    cells = {}
    panels = []
    proteins = {"t": np.zeros((2, 1), dtype=np.float64)}
    ligands = {"l": np.ones(1, dtype=np.float64)}
    try:
        learn_slot_localizer(cells, panels, proteins, ligands, top_slots=3)
    except ValueError as exc:
        assert "top_slots" in str(exc)
    else:
        raise AssertionError("expected top_slots validation failure")
