"""Gradient-based few-shot adaptation baseline."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Literal

import torch
from torch import nn
from torch.nn import functional as F

from .interaction import InteractionEncoder


Backbone = Literal["transformer", "mamba", "hybrid"]


@dataclass(frozen=True)
class AdaptState:
    intercept: torch.Tensor
    scale: torch.Tensor
    code: torch.Tensor
    supportligand: torch.Tensor
    supportbase: torch.Tensor

    def zerocode(self) -> "AdaptState":
        return replace(self, code=torch.zeros_like(self.code))


class SupportEncoder(nn.Module):
    def __init__(self, liganddim: int, taskdim: int) -> None:
        super().__init__()
        hidden = max(64, taskdim * 8)
        self.element = nn.Sequential(
            nn.Linear(liganddim + 2, hidden),
            nn.GELU(),
            nn.Linear(hidden, taskdim),
        )

    def forward(
        self,
        ligand: torch.Tensor,
        base: torch.Tensor,
        residual: torch.Tensor,
    ) -> torch.Tensor:
        values = torch.cat((ligand, base[:, None], residual[:, None]), dim=1)
        return self.element(values).mean(dim=0)


class GradientAdapter(nn.Module):
    """MAML-style baseline with an exact zero-code interaction null."""

    def __init__(
        self,
        *,
        proteindim: int,
        liganddim: int,
        dmodel: int = 256,
        taskdim: int = 8,
        stages: int = 2,
        landmarks: int = 32,
        backbone: Backbone = "hybrid",
    ) -> None:
        super().__init__()
        self.taskdim = taskdim
        self.interaction = InteractionEncoder(
            proteindim,
            liganddim,
            dmodel,
            taskdim,
            stages,
            landmarks,
            backbone,
        )
        self.support = SupportEncoder(liganddim, taskdim)
        self.readout = nn.Parameter(torch.empty(taskdim, taskdim))
        nn.init.orthogonal_(self.readout, gain=1e-2)

    def pairfeatures(
        self,
        proteintokens: torch.Tensor,
        ligand: torch.Tensor,
        code: torch.Tensor,
    ) -> torch.Tensor:
        return self.interaction(proteintokens, ligand, code)


def buildadapter(
    *,
    protein_dim: int,
    ligand_dim: int,
    d_model: int = 256,
    task_dim: int = 8,
    hybrid_stages: int = 2,
    landmarks: int = 32,
    backbone: Backbone = "hybrid",
) -> GradientAdapter:
    return GradientAdapter(
        proteindim=protein_dim,
        liganddim=ligand_dim,
        dmodel=d_model,
        taskdim=task_dim,
        stages=hybrid_stages,
        landmarks=landmarks,
        backbone=backbone,
    )


def buildmatched(
    *,
    protein_dim: int,
    ligand_dim: int,
    d_model: int = 256,
    task_dim: int = 8,
    hybrid_stages: int = 2,
    landmarks: int = 32,
) -> dict[Backbone, GradientAdapter]:
    """Build comparison arms without artificial parameter padding."""

    return {
        kind: buildadapter(
            protein_dim=protein_dim,
            ligand_dim=ligand_dim,
            d_model=d_model,
            task_dim=task_dim,
            hybrid_stages=hybrid_stages,
            landmarks=landmarks,
            backbone=kind,
        )
        for kind in ("transformer", "mamba", "hybrid")
    }


def countparams(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)


def fitcalibration(
    supportbase: torch.Tensor,
    supportlabel: torch.Tensor,
    ridge: float = 1e-4,
) -> tuple[torch.Tensor, torch.Tensor]:
    design = torch.stack((torch.ones_like(supportbase), supportbase), dim=1)
    gram = design.T @ design + ridge * torch.eye(2, device=design.device, dtype=design.dtype)
    coefficients = torch.linalg.solve(gram, design.T @ supportlabel)
    return coefficients[0], coefficients[1]


def orthogonalize(
    supportfeature: torch.Tensor,
    supportbase: torch.Tensor,
    queryfeature: torch.Tensor,
    querybase: torch.Tensor,
    ridge: float = 1e-4,
) -> torch.Tensor:
    supportdesign = torch.stack((torch.ones_like(supportbase), supportbase), dim=1)
    querydesign = torch.stack((torch.ones_like(querybase), querybase), dim=1)
    gram = supportdesign.T @ supportdesign
    gram = gram + ridge * torch.eye(2, device=gram.device, dtype=gram.dtype)
    coefficients = torch.linalg.solve(gram, supportdesign.T @ supportfeature)
    return queryfeature - querydesign @ coefficients


def predict(
    model: GradientAdapter,
    state: AdaptState,
    proteintokens: torch.Tensor,
    queryligand: torch.Tensor,
    querybase: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    supportfeature = model.pairfeatures(proteintokens, state.supportligand, state.code)
    queryfeature = model.pairfeatures(proteintokens, queryligand, state.code)
    feature = orthogonalize(
        supportfeature,
        state.supportbase,
        queryfeature,
        querybase,
    )
    residual = feature @ model.readout @ state.code
    calibration = state.intercept + state.scale * querybase
    return calibration + residual, residual


def adapttarget(
    model: GradientAdapter,
    *,
    protein_tokens: torch.Tensor,
    support_ligand: torch.Tensor,
    support_y: torch.Tensor,
    support_b0: torch.Tensor,
    inner_steps: int = 2,
    inner_lr: float = 0.05,
) -> AdaptState:
    if not 1 <= inner_steps <= 3:
        raise ValueError("inner_steps must be in [1, 3]")
    if support_ligand.shape[0] != support_y.numel() or support_y.numel() != support_b0.numel():
        raise ValueError("support tensors must share their first dimension")

    intercept, scale = fitcalibration(support_b0.float(), support_y.float())
    residual = support_y.float() - (intercept + scale * support_b0.float())
    code = model.support(support_ligand.float(), support_b0.float(), residual)
    intercept = intercept.detach().requires_grad_(True)
    scale = scale.detach().requires_grad_(True)
    code.requires_grad_(True)

    for _ in range(inner_steps):
        state = AdaptState(intercept, scale, code, support_ligand, support_b0)
        prediction, _ = predict(model, state, protein_tokens, support_ligand, support_b0)
        loss = F.huber_loss(prediction, support_y.float())
        gradintercept, gradscale, gradcode = torch.autograd.grad(
            loss, (intercept, scale, code), create_graph=False, retain_graph=True
        )
        intercept = (intercept - inner_lr * gradintercept).detach().requires_grad_(True)
        scale = (scale - inner_lr * gradscale).detach().requires_grad_(True)
        code = code - inner_lr * gradcode

    return AdaptState(intercept.detach(), scale.detach(), code, support_ligand, support_b0)


def predictquery(
    model: GradientAdapter,
    state: AdaptState,
    *,
    protein_tokens: torch.Tensor,
    query_ligand: torch.Tensor,
    query_b0: torch.Tensor,
) -> dict[str, torch.Tensor]:
    prediction, residual = predict(model, state, protein_tokens, query_ligand, query_b0)
    calibration = state.intercept + state.scale * query_b0
    return {
        "prediction": prediction,
        "calibration": calibration,
        "residual": residual,
        "task_code": state.code,
    }


__all__ = [
    "GradientAdapter",
    "adapttarget",
    "buildadapter",
    "buildmatched",
    "countparams",
    "predictquery",
]
