"""Canonical complete-panel quotient projection for affinity interactions."""

from __future__ import annotations

import numpy as np
import torch


def panel_design(target_ids: list[str], ligand_ids: list[str]) -> np.ndarray:
    if len(target_ids) != len(ligand_ids) or not target_ids:
        raise ValueError("a panel requires equally sized, non-empty identity lists")
    targets = {value: index for index, value in enumerate(sorted(set(target_ids)))}
    ligands = {value: index for index, value in enumerate(sorted(set(ligand_ids)))}
    design = np.zeros(
        (len(target_ids), 1 + len(targets) + len(ligands)), dtype=np.float64
    )
    design[:, 0] = 1.0
    for row, (target, ligand) in enumerate(zip(target_ids, ligand_ids)):
        design[row, 1 + targets[target]] = 1.0
        design[row, 1 + len(targets) + ligands[ligand]] = 1.0
    return design


def nuisance_basis(design: np.ndarray) -> tuple[np.ndarray, int]:
    matrix = np.asarray(design, dtype=np.float64)
    if matrix.ndim != 2 or matrix.shape[0] == 0:
        raise ValueError("design must be a non-empty matrix")
    left, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    tolerance = singular.max(initial=0.0) * max(matrix.shape) * np.finfo(np.float64).eps
    rank = int(np.count_nonzero(singular > tolerance))
    retained_rank = matrix.shape[0] - rank
    if retained_rank <= 0:
        raise ValueError("panel has no interaction quotient")
    return left[:, :rank], retained_rank


def quotient_residual(values: torch.Tensor, basis: torch.Tensor) -> torch.Tensor:
    if values.ndim not in (1, 2):
        raise ValueError("values must have shape [edges] or [edges, features]")
    if basis.ndim != 2 or basis.shape[0] != values.shape[0]:
        raise ValueError("basis and values disagree on panel edges")
    if values.dtype != torch.float64 or basis.dtype != torch.float64:
        raise ValueError("quotient projection requires float64")
    return values - basis @ (basis.transpose(0, 1) @ values)


def rank_normalized_loss(
    response: torch.Tensor,
    prediction: torch.Tensor,
    basis: torch.Tensor,
    retained_rank: int,
) -> torch.Tensor:
    if retained_rank <= 0:
        raise ValueError("retained_rank must be positive")
    residual = quotient_residual(response - prediction, basis)
    return residual.square().sum() / retained_rank
