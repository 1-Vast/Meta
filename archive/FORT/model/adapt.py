"""Calibration-null contrast posterior for few-shot target adaptation."""

from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import nn

from .ligandbase import LigandBaseline
from .reorder import ReorderingModel, ReorderingState


@dataclass(frozen=True)
class TargetState:
    """Label-safe features and the single support residual vector."""

    supportligand: torch.Tensor
    supportlabel: torch.Tensor
    supportb0: torch.Tensor
    supportprediction: torch.Tensor
    ligandfeature: torch.Tensor
    protein: ReorderingState


class TargetAdapter(nn.Module):
    """Adapt a frozen base with identifiable support-label contrasts."""

    def __init__(
        self,
        ligand: LigandBaseline,
        protein: ReorderingModel,
        *,
        inclusion: float = 0.05,
        mode: str = "joint",
    ) -> None:
        super().__init__()
        if not 0.0 < inclusion < 1.0:
            raise ValueError("protein inclusion prior must lie in (0, 1)")
        if mode not in {"joint", "soft", "exact"}:
            raise ValueError("posterior mode must be joint, soft, or exact")
        self.ligand = ligand
        self.protein = protein
        self.inclusion = inclusion
        self.mode = mode
        self.ligand.requires_grad_(False)
        self.ligand.eval()

    def train(self, mode: bool = True) -> TargetAdapter:
        super().train(mode)
        self.ligand.eval()
        return self

    @staticmethod
    def validatesupport(
        supportligand: torch.Tensor,
        supportlabel: torch.Tensor,
        supportb0: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if supportligand.ndim != 2:
            raise ValueError("support ligand features must have shape [rows, ligand_dim]")
        label = supportlabel.float().reshape(-1)
        base = supportb0.float().reshape(-1)
        if supportligand.shape[0] != label.numel() or label.numel() != base.numel():
            raise ValueError("support tensors must share their first dimension")
        return label, base

    def adapt(
        self,
        *,
        proteintokens: torch.Tensor,
        supportligand: torch.Tensor,
        supportlabel: torch.Tensor,
        supportb0: torch.Tensor,
    ) -> TargetState:
        label, base = self.validatesupport(
            supportligand, supportlabel, supportb0
        )
        with torch.no_grad():
            representation = self.ligand.ligand(supportligand)
            ligandfeature = self.ligand.bayes.features(representation)
            supportprediction = self._leaveoneoutcalibration(label, base)
        protein = self.protein.adapt(
            proteintokens=proteintokens,
            supportligand=supportligand,
        )
        return TargetState(
            supportligand=supportligand,
            supportlabel=label,
            supportb0=base,
            supportprediction=supportprediction,
            ligandfeature=ligandfeature,
            protein=protein,
        )

    @torch.no_grad()
    def _leaveoneoutcalibration(
        self,
        supportlabel: torch.Tensor,
        supportb0: torch.Tensor,
    ) -> torch.Tensor:
        """Predict each support label from an affine calibration fit without it."""

        rows = supportlabel.numel()
        if rows == 0:
            return supportlabel.new_zeros(0)
        predictions = []
        for heldout in range(rows):
            keep = torch.arange(rows, device=supportlabel.device) != heldout
            output = self.protein.posterior.calibration(
                supportlabel[keep] - supportb0[keep],
                supportb0[keep],
                supportb0[heldout : heldout + 1],
            )
            prediction = (
                supportb0[heldout]
                + output["weight"] * output["mean"].reshape(())
            )
            predictions.append(prediction)
        return torch.stack(predictions)

    @torch.no_grad()
    def basepredict(
        self,
        state: TargetState,
        queryligand: torch.Tensor,
        queryb0: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        self.ligand.eval()
        return self.ligand(
            queryligand,
            state.supportligand,
            state.supportlabel,
            query_base=queryb0,
            support_base=state.supportb0,
        )

    @staticmethod
    def _fallback(base: dict[str, torch.Tensor]) -> dict[str, torch.Tensor]:
        zero = torch.zeros_like(base["pred"])
        scalar = base["pred"].new_zeros(())
        return {
            "prediction": base["pred"],
            "pred": base["pred"],
            "ligandprediction": base["pred"],
            "jointprediction": base["pred"],
            "appliedprotein": zero,
            "ligandvariance": base["epistemic"],
            "jointvariance": base["epistemic"],
            "latentvariance": base["epistemic"],
            "observationvariance": base["aleatoric"],
            "totalvariance": base["total_variance"],
            "proteinprobability": scalar,
            "proteinweight": scalar,
            "proteinlogbf": scalar,
        }

    def predict(
        self,
        state: TargetState,
        *,
        queryligand: torch.Tensor,
        queryb0: torch.Tensor,
        useprotein: bool = True,
        mode: str | None = None,
    ) -> dict[str, torch.Tensor]:
        querybase = queryb0.float().reshape(-1)
        if queryligand.ndim != 2 or queryligand.shape[0] != querybase.numel():
            raise ValueError("query tensors must share their first dimension")
        selected = self.mode if mode is None else mode
        if selected not in {"joint", "soft", "exact"}:
            raise ValueError("posterior mode must be joint, soft, or exact")

        base = self.basepredict(state, queryligand, querybase)
        if not useprotein or state.supportlabel.numel() <= 1:
            return self._fallback(base)

        calibration = self.protein.posterior.calibration(
            state.supportlabel - state.supportb0,
            state.supportb0,
            querybase,
        )
        calibrated = querybase + calibration["weight"] * calibration["mean"]
        calibrationvariance = (
            calibration["weight"] * calibration["variance"]
            + calibration["weight"]
            * (1.0 - calibration["weight"])
            * calibration["mean"].square()
        )
        residual = state.supportlabel - state.supportprediction
        if self.protein.proteinconditioned:
            reordered = self.protein.predict(
                state.protein,
                queryligand=queryligand,
                supportresidual=residual,
                supportbase=state.supportb0,
                querybase=calibrated,
                querycalibrationbase=querybase,
                basevariance=calibrationvariance,
                observationvariance=base["aleatoric"],
                rankingmode=selected,
                applycalibration=False,
            )
        else:
            with torch.no_grad():
                queryrepresentation = self.ligand.ligand(queryligand)
                queryfeature = self.ligand.bayes.features(queryrepresentation)[:, 1:]
            supportfeature = state.ligandfeature[:, 1:]
            reordered = self.protein.posterior(
                queryfeature,
                supportfeature,
                residual,
                calibrated,
                state.protein.basis,
                state.protein.priorcov,
                basevariance=calibrationvariance,
                observationvariance=base["aleatoric"],
                supportbase=state.supportb0,
                querycalibrationbase=querybase,
                rankingmode=selected,
                applycalibration=False,
            )
        prediction = reordered["prediction"]
        return {
            "prediction": prediction,
            "pred": prediction,
            "ligandprediction": base["pred"],
            "calibrationprediction": calibrated,
            "calibrationvariance": calibrationvariance,
            "jointprediction": calibrated + reordered["rawrankingmean"],
            "appliedprotein": reordered["appliedranking"],
            "ligandvariance": base["epistemic"],
            "jointvariance": calibrationvariance + reordered["rawrankingvariance"],
            "latentvariance": reordered["latentvariance"],
            "observationvariance": base["aleatoric"],
            "totalvariance": reordered["totalvariance"],
            "proteinprobability": reordered["rankprobability"],
            "proteinweight": reordered["rankweight"],
            "proteinlogbf": reordered["rankinglogbf"],
        }


__all__ = ["TargetAdapter", "TargetState"]
