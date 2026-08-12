import numpy as np
import torch

from research.crossed_interaction.train_attention_localizer_cq_observable import (
    SlotAttentionObservable,
    TrainBlock,
    projection_residual_matrix,
    train_attention_localizer,
)


def test_projection_residual_matrix_removes_additive_effects():
    target_ids = ["t1", "t1", "t2", "t2"]
    ligand_ids = ["l1", "l2", "l1", "l2"]
    projection = projection_residual_matrix(target_ids, ligand_ids)
    additive = np.asarray([3.0, 5.0, 7.0, 9.0])
    assert np.allclose(projection @ additive, np.zeros(4))
    interaction = np.asarray([1.0, -1.0, -1.0, 1.0])
    assert np.allclose(projection @ interaction, interaction)


def test_attention_weights_respect_mask():
    model = SlotAttentionObservable(slot_dim=2, ligand_dim=3, attention_dim=4)
    slots = torch.randn(2, 3, 2)
    mask = torch.tensor([[True, False, True], [False, True, True]])
    ligand = torch.randn(2, 3)
    weights = model.attention_weights(slots, mask, ligand)
    assert weights.shape == (2, 3)
    assert torch.allclose(weights.sum(dim=1), torch.ones(2))
    assert torch.equal(weights[~mask], torch.zeros_like(weights[~mask]))


def test_train_attention_localizer_learns_train_quotient_proxy():
    proteins = {
        "t1": np.asarray([[1.0], [0.0]], dtype=np.float64),
        "t2": np.asarray([[-1.0], [0.0]], dtype=np.float64),
    }
    masks = {key: np.ones(2, dtype=bool) for key in proteins}
    ligands = {
        "l1": np.asarray([1.0], dtype=np.float64),
        "l2": np.asarray([-1.0], dtype=np.float64),
    }
    y = np.asarray([1.0, -1.0, -1.0, 1.0], dtype=np.float64)
    block = TrainBlock(
        panel_id="train",
        protein_keys=["t1", "t1", "t2", "t2"],
        ligand_keys=["l1", "l2", "l1", "l2"],
        y=y,
        residual_matrix=projection_residual_matrix(
            ["t1", "t1", "t2", "t2"], ["l1", "l2", "l1", "l2"]),
    )
    model, metadata = train_attention_localizer(
        [block], proteins, masks, ligands, mode="attention",
        attention_dim=2, epochs=30, learning_rate=0.05,
        weight_decay=0.0, seed=7, device="cpu")
    assert metadata["mode"] == "attention"
    assert metadata["loss_history"][-1]["train_row_mse"] < 0.5
    assert isinstance(model, SlotAttentionObservable)
