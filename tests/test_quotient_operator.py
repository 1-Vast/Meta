import numpy as np
import pytest
import torch

from research.crossed_interaction.quotient_operator import (
    nuisance_basis,
    panel_design,
    quotient_residual,
    rank_normalized_loss,
)


TARGETS = ["p1", "p1", "p2", "p2"]
LIGANDS = ["a", "b", "a", "b"]


def _panel(device="cpu"):
    design = panel_design(TARGETS, LIGANDS)
    basis, rank = nuisance_basis(design)
    return design, torch.tensor(basis, dtype=torch.float64, device=device), rank


def test_additive_null_is_annihilated_and_rank_is_one():
    design, basis, retained_rank = _panel()
    additive = torch.tensor(design @ np.arange(design.shape[1]), dtype=torch.float64)
    assert retained_rank == 1
    assert quotient_residual(additive, basis).abs().max().item() < 1e-10


def test_disconnected_graph_rank_and_permutation_invariance():
    targets = ["p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4"]
    ligands = ["a", "b", "a", "b", "c", "d", "c", "d"]
    design = panel_design(targets, ligands)
    _, retained_rank = nuisance_basis(design)
    assert retained_rank == 2
    values = np.asarray([1.0, 2.0, 4.0, 8.0, 2.0, 3.0, 5.0, 9.0])
    basis, _ = nuisance_basis(design)
    reference = np.linalg.norm(values - basis @ (basis.T @ values))
    order = np.asarray([7, 2, 4, 0, 6, 3, 1, 5])
    permuted = panel_design(
        [targets[index] for index in order], [ligands[index] for index in order]
    )
    permuted_basis, _ = nuisance_basis(permuted)
    observed = np.linalg.norm(values[order] - permuted_basis @ (permuted_basis.T @ values[order]))
    assert observed == pytest.approx(reference, abs=1e-12)


def test_loss_has_nonzero_gradient():
    _, basis, retained_rank = _panel()
    response = torch.tensor([0.0, 1.0, 1.0, 0.0], dtype=torch.float64)
    prediction = torch.zeros(4, dtype=torch.float64, requires_grad=True)
    loss = rank_normalized_loss(response, prediction, basis, retained_rank)
    loss.backward()
    assert prediction.grad is not None
    assert prediction.grad.abs().max().item() > 0


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA unavailable")
def test_cpu_gpu_projection_agree():
    _, cpu_basis, retained_rank = _panel()
    response = torch.tensor([0.0, 1.0, 1.0, 0.0], dtype=torch.float64)
    cpu = rank_normalized_loss(response, torch.zeros_like(response), cpu_basis, retained_rank)
    gpu = rank_normalized_loss(
        response.cuda(), torch.zeros_like(response).cuda(), cpu_basis.cuda(), retained_rank
    )
    assert gpu.cpu().item() == pytest.approx(cpu.item(), abs=1e-12)
