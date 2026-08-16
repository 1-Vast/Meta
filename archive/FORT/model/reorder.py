"""Protein-conditioned finite-rank Bayesian ligand reordering."""

from __future__ import annotations

from dataclasses import dataclass
import math

import torch
from torch import nn

from .interaction import InteractionEncoder


JITTER = 1e-6


def stiefelqr(raw: torch.Tensor) -> torch.Tensor:
    """Map an ambient matrix to sign-stable orthonormal columns."""

    q, r = torch.linalg.qr(raw, mode="reduced")
    sign = torch.sign(torch.diagonal(r, dim1=-2, dim2=-1))
    sign = torch.where(sign == 0, torch.ones_like(sign), sign)
    return q * sign.unsqueeze(-2)


def helmert(k: int, *, device=None, dtype=torch.float32) -> torch.Tensor:
    """Return an orthonormal basis for the support contrast space."""

    if k <= 1:
        return torch.zeros(0, max(k, 0), device=device, dtype=dtype)
    contrast = torch.zeros(k - 1, k, device=device, dtype=dtype)
    for index in range(1, k):
        scale = 1.0 / math.sqrt(index * (index + 1))
        contrast[index - 1, :index] = scale
        contrast[index - 1, index] = -index * scale
    return contrast


def calibrationgeometry(
    supportbase: torch.Tensor,
    querybase: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Return orthogonal calibration and contrast designs for ``[1, base]``."""

    supportbase = supportbase.float().reshape(-1)
    querybase = querybase.float().reshape(-1)
    k = supportbase.numel()
    if k == 0:
        return (
            supportbase.new_zeros((0, 0)),
            supportbase.new_zeros((querybase.numel(), 0)),
            supportbase.new_zeros((0, 0)),
        )

    constant = torch.ones(k, 1, device=supportbase.device, dtype=supportbase.dtype)
    columns = [constant / math.sqrt(k)]
    querycolumns = [
        torch.ones(
            querybase.numel(), 1, device=querybase.device, dtype=querybase.dtype
        )
        / math.sqrt(k)
    ]
    centered = supportbase - supportbase.mean()
    norm = torch.linalg.vector_norm(centered)
    tolerance = torch.finfo(supportbase.dtype).eps * max(k, 1) * supportbase.abs().max().clamp_min(1.0)
    if bool(norm > tolerance):
        columns.append((centered / norm)[:, None])
        querycolumns.append(((querybase - supportbase.mean()) / norm)[:, None])

    calibration = torch.cat(columns, dim=1)
    querycalibration = torch.cat(querycolumns, dim=1)
    complete = torch.linalg.qr(calibration, mode="complete").Q
    contrast = complete[:, calibration.shape[1] :].T
    return calibration, querycalibration, contrast


class ProteinSubspace(nn.Module):
    """Map a label-free protein representation to an orthonormal basis."""

    def __init__(self, targetdim: int, ambientdim: int, rank: int, hidden: int = 64) -> None:
        super().__init__()
        self.ambientdim = ambientdim
        self.rank = rank
        self.net = nn.Sequential(
            nn.Linear(targetdim, hidden),
            nn.Tanh(),
            nn.Linear(hidden, ambientdim * rank),
        )

    def forward(self, target: torch.Tensor) -> torch.Tensor:
        raw = self.net(target.float()).reshape(self.ambientdim, self.rank)
        return stiefelqr(raw)


def modelprobability(logbf: torch.Tensor, prior: float) -> torch.Tensor:
    if not 0.0 < prior < 1.0:
        raise ValueError("model inclusion prior must lie in (0, 1)")
    logodds = math.log(prior / (1.0 - prior))
    return torch.sigmoid(logbf + logodds)


class ReorderingPosterior(nn.Module):
    """Separate calibration and ranking posteriors with coherent BMA variance."""

    def __init__(
        self,
        ambientdim: int,
        rank: int,
        calibrationprior: float = 1.0,
        noiseinit: float = 1.0,
        rankinginclusion: float = 0.05,
        calibrationinclusion: float = 0.5,
        gatethreshold: float = 0.0,
        rankinghardgate: bool = True,
        calibrationhardgate: bool = False,
    ) -> None:
        super().__init__()
        if not 1 <= rank <= ambientdim:
            raise ValueError("rank must be between one and the ambient dimension")
        self.ambientdim = ambientdim
        self.rank = rank
        self.rankinginclusion = rankinginclusion
        self.calibrationinclusion = calibrationinclusion
        self.rankinghardgate = rankinghardgate
        self.calibrationhardgate = calibrationhardgate
        self.logcalibrationprior = nn.Parameter(torch.tensor(math.log(calibrationprior)))
        self.lognoise = nn.Parameter(torch.tensor(math.log(noiseinit)))
        self.register_buffer("gatethreshold", torch.tensor(float(gatethreshold)))

    def noisevariance(self) -> torch.Tensor:
        return torch.exp(self.lognoise.clamp(-6.0, 6.0)).square()

    def gate(
        self,
        logbf: torch.Tensor,
        prior: float,
        hardgate: bool,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        probability = modelprobability(logbf, prior)
        if self.training or not hardgate:
            return probability, probability
        active = (logbf > self.gatethreshold).to(probability.dtype)
        return probability, probability * active

    def ranking(
        self,
        query: torch.Tensor,
        support: torch.Tensor,
        residual: torch.Tensor,
        priorcov: torch.Tensor,
        supportbase: torch.Tensor | None = None,
        querybase: torch.Tensor | None = None,
        mode: str | None = None,
    ) -> dict[str, torch.Tensor]:
        k = support.shape[0]
        zero = torch.zeros(query.shape[0], device=query.device, dtype=query.dtype)
        if supportbase is None:
            supportbase = support.new_zeros(k)
        if querybase is None:
            querybase = query.new_zeros(query.shape[0])
        calibration, querycalibration, contrast = calibrationgeometry(
            supportbase, querybase
        )
        if contrast.shape[0] == 0:
            scalar = support.new_zeros(())
            return {
                "mean": zero,
                "variance": zero,
                "logbf": scalar,
                "probability": scalar,
                "weight": scalar,
            }

        centeredresidual = contrast @ residual
        centeredfeature = contrast @ support
        calibrationcoefficient = calibration.T @ support
        centeredquery = query - querycalibration @ calibrationcoefficient
        noise = self.noisevariance().to(support.dtype)
        eye = torch.eye(self.rank, device=support.device, dtype=support.dtype)
        priorchol = torch.linalg.cholesky(priorcov + JITTER * eye)
        priorprecision = torch.cholesky_solve(eye, priorchol)
        precision = priorprecision + centeredfeature.T @ centeredfeature / noise
        chol = torch.linalg.cholesky(precision + JITTER * eye)
        mean = torch.cholesky_solve(
            (centeredfeature.T @ centeredresidual / noise)[:, None], chol
        ).squeeze(1)
        querymean = centeredquery @ mean
        solved = torch.linalg.solve_triangular(chol, centeredquery.T, upper=False)
        queryvariance = solved.square().sum(0)

        samples = contrast.shape[0]
        sampleeye = torch.eye(samples, device=support.device, dtype=support.dtype)
        alternative = noise * sampleeye + centeredfeature @ priorcov @ centeredfeature.T
        alternativechol = torch.linalg.cholesky(alternative + JITTER * sampleeye)
        alternativequad = centeredresidual @ torch.cholesky_solve(
            centeredresidual[:, None], alternativechol
        ).squeeze(1)
        alternativelogdet = 2.0 * torch.log(torch.diagonal(alternativechol)).sum()
        nullquad = centeredresidual.square().sum() / noise
        nulllogdet = samples * torch.log(noise)
        logbf = (
            -0.5 * (alternativequad + alternativelogdet)
            + 0.5 * (nullquad + nulllogdet)
        ).clamp(-20.0, 20.0)
        probability = modelprobability(logbf, self.rankinginclusion)
        selected = mode
        if selected is None:
            selected = "soft" if self.training or not self.rankinghardgate else "exact"
        if selected == "joint":
            weight = torch.ones_like(probability)
        elif selected == "soft":
            weight = probability
        elif selected == "exact":
            weight = (logbf > self.gatethreshold).to(probability.dtype)
        else:
            raise ValueError("ranking mode must be joint, soft, or exact")
        return {
            "mean": querymean,
            "variance": queryvariance,
            "logbf": logbf,
            "probability": probability,
            "weight": weight,
            "posteriormean": mean,
            "precisionchol": chol,
        }

    def calibration(
        self,
        residual: torch.Tensor,
        supportbase: torch.Tensor | None = None,
        querybase: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        k = residual.shape[0]
        if k == 0:
            scalar = residual.new_zeros(())
            queryrows = 0 if querybase is None else querybase.numel()
            zero = residual.new_zeros(queryrows)
            return {
                "mean": zero,
                "variance": zero,
                "logbf": scalar,
                "probability": scalar,
                "weight": scalar,
            }
        if supportbase is None:
            supportbase = residual.new_zeros(k)
        if querybase is None:
            querybase = residual.new_zeros(1)
        design, querydesign, _ = calibrationgeometry(supportbase, querybase)
        noise = self.noisevariance().to(residual.dtype)
        priorvariance = torch.exp(self.logcalibrationprior.clamp(-8.0, 8.0)).to(residual.dtype)
        dimension = design.shape[1]
        eye = torch.eye(dimension, device=residual.device, dtype=residual.dtype)
        precision = eye / priorvariance + design.T @ design / noise
        chol = torch.linalg.cholesky(precision + JITTER * eye)
        posteriormean = torch.cholesky_solve(
            (design.T @ residual / noise)[:, None], chol
        ).squeeze(1)
        querymean = querydesign @ posteriormean
        solved = torch.linalg.solve_triangular(chol, querydesign.T, upper=False)
        queryvariance = solved.square().sum(0)
        rows = torch.eye(k, device=residual.device, dtype=residual.dtype)
        alternative = noise * rows + priorvariance * design @ design.T
        alternativechol = torch.linalg.cholesky(alternative + JITTER * rows)
        alternativequad = residual @ torch.cholesky_solve(
            residual[:, None], alternativechol
        ).squeeze(1)
        nullquad = residual.square().sum() / noise
        logbf = (
            -0.5
            * (
                alternativequad
                + 2.0 * torch.log(torch.diagonal(alternativechol)).sum()
            )
            + 0.5 * (nullquad + k * torch.log(noise))
        ).clamp(-20.0, 20.0)
        probability, weight = self.gate(
            logbf, self.calibrationinclusion, self.calibrationhardgate
        )
        return {
            "mean": querymean,
            "variance": queryvariance,
            "logbf": logbf,
            "probability": probability,
            "weight": weight,
        }

    def forward(
        self,
        queryfeature: torch.Tensor,
        supportfeature: torch.Tensor | None,
        supportresidual: torch.Tensor | None,
        querybase: torch.Tensor,
        basis: torch.Tensor,
        priorcov: torch.Tensor,
        basevariance: torch.Tensor | float | None = None,
        observationvariance: torch.Tensor | float | None = None,
        supportbase: torch.Tensor | None = None,
        querycalibrationbase: torch.Tensor | None = None,
        rankingmode: str | None = None,
        applycalibration: bool = True,
    ) -> dict[str, torch.Tensor]:
        queryrank = queryfeature.float() @ basis
        if basevariance is None:
            basevariance = torch.zeros_like(querybase)
        else:
            basevariance = torch.as_tensor(
                basevariance, device=querybase.device, dtype=querybase.dtype
            ).expand_as(querybase)
        if observationvariance is None:
            observationvariance = self.noisevariance().to(querybase.dtype).expand_as(querybase)
        else:
            observationvariance = torch.as_tensor(
                observationvariance, device=querybase.device, dtype=querybase.dtype
            ).expand_as(querybase)

        if supportfeature is None or supportfeature.shape[0] == 0:
            residual = querybase.new_zeros(0)
            supportrank = queryrank.new_zeros((0, self.rank))
        else:
            if supportresidual is None or supportresidual.numel() != supportfeature.shape[0]:
                raise ValueError("support residuals must match support features")
            residual = supportresidual.float().reshape(-1)
            supportrank = supportfeature.float() @ basis

        if supportbase is None:
            supportbase = residual.new_zeros(residual.shape[0])
        if querycalibrationbase is None:
            querycalibrationbase = querybase
        ranking = self.ranking(
            queryrank,
            supportrank,
            residual,
            priorcov,
            supportbase,
            querycalibrationbase,
            rankingmode,
        )
        calibration = self.calibration(
            residual, supportbase, querycalibrationbase
        )
        appliedranking = ranking["weight"] * ranking["mean"]
        appliedcalibration = calibration["weight"] * calibration["mean"]
        rankingvariance = (
            ranking["weight"] * ranking["variance"]
            + ranking["weight"] * (1.0 - ranking["weight"]) * ranking["mean"].square()
        )
        calibrationvariance = (
            calibration["weight"] * calibration["variance"]
            + calibration["weight"]
            * (1.0 - calibration["weight"])
            * calibration["mean"].square()
        )
        if not applycalibration:
            appliedcalibration = torch.zeros_like(querybase)
            calibrationvariance = torch.zeros_like(querybase)
        latentvariance = basevariance + rankingvariance + calibrationvariance
        prediction = querybase + appliedcalibration + appliedranking
        return {
            "prediction": prediction,
            "pred": prediction,
            "calibrationprediction": querybase + appliedcalibration,
            "rawrankingmean": ranking["mean"],
            "rawrankingvariance": ranking["variance"],
            "appliedranking": appliedranking,
            "rawcalibrationmean": calibration["mean"],
            "rawcalibrationvariance": calibration["variance"],
            "appliedcalibration": appliedcalibration,
            "rankingvariance": rankingvariance,
            "calibrationvariance": calibrationvariance,
            "latentvariance": latentvariance,
            "observationvariance": observationvariance,
            "totalvariance": latentvariance + observationvariance,
            "rankprobability": ranking["probability"],
            "rankweight": ranking["weight"],
            "calibrationprobability": calibration["probability"],
            "calibrationweight": calibration["weight"],
            "rankinglogbf": ranking["logbf"],
            "calibrationlogbf": calibration["logbf"],
        }


@dataclass(frozen=True)
class ReorderingState:
    proteinfeature: torch.Tensor
    supportfeature: torch.Tensor
    basis: torch.Tensor
    priorcov: torch.Tensor


class ReorderingModel(nn.Module):
    """Label-safe encoder and finite-rank posterior used by the primary path."""

    def __init__(
        self,
        *,
        proteindim: int,
        liganddim: int,
        dmodel: int = 256,
        ambientdim: int = 8,
        rank: int = 2,
        primaryk: int = 5,
        stages: int = 2,
        landmarks: int = 32,
        backbone: str = "hybrid",
        proteinconditioned: bool = True,
        interactiononly: bool = True,
        rankinghardgate: bool = True,
        calibrationhardgate: bool = False,
    ) -> None:
        super().__init__()
        if rank > primaryk - 1:
            raise ValueError("posterior rank cannot exceed primary k minus one")
        self.ambientdim = ambientdim
        self.rank = rank
        self.proteinconditioned = proteinconditioned
        if proteinconditioned:
            self.interaction = InteractionEncoder(
                proteindim,
                liganddim,
                dmodel,
                ambientdim,
                stages,
                landmarks,
                backbone,
                conditioned=False,
                interactiononly=interactiononly,
            )
            self.subspace = ProteinSubspace(dmodel, ambientdim, rank)
            self.register_parameter("globalbasis", None)
        else:
            self.interaction = None
            self.subspace = None
            self.globalbasis = nn.Parameter(torch.eye(ambientdim, rank))
        self.priorraw = nn.Parameter(torch.eye(rank))
        self.posterior = ReorderingPosterior(
            ambientdim,
            rank,
            rankinghardgate=rankinghardgate,
            calibrationhardgate=calibrationhardgate,
        )

    def priorcovariance(self) -> torch.Tensor:
        diagonal = torch.nn.functional.softplus(torch.diagonal(self.priorraw)) + 1e-3
        factor = torch.tril(self.priorraw, diagonal=-1) + torch.diag(diagonal)
        return factor @ factor.T

    def makebasis(self, proteinfeature: torch.Tensor) -> torch.Tensor:
        if self.subspace is not None:
            return self.subspace(proteinfeature.mean(dim=0))
        return stiefelqr(self.globalbasis)

    def adapt(
        self,
        *,
        proteintokens: torch.Tensor,
        supportligand: torch.Tensor,
    ) -> ReorderingState:
        if supportligand.ndim != 2:
            raise ValueError("support ligand features must be a matrix")
        if self.interaction is None:
            proteinfeature = supportligand.new_zeros((0, self.ambientdim))
            supportfeature = supportligand.new_zeros(
                (supportligand.shape[0], self.ambientdim)
            )
            basis = stiefelqr(self.globalbasis)
        else:
            proteinfeature = self.interaction.encodeprotein(proteintokens)
            supportfeature = self.interaction.pairfromprotein(
                proteinfeature, supportligand
            )
            basis = self.makebasis(proteinfeature)
        return ReorderingState(
            proteinfeature=proteinfeature,
            supportfeature=supportfeature,
            basis=basis,
            priorcov=self.priorcovariance(),
        )

    def rankfeatures(
        self,
        state: ReorderingState,
        queryligand: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Return support-centred protein interaction coordinates."""

        if self.interaction is None:
            raise RuntimeError("protein rank features require a protein-conditioned model")

        support = state.supportfeature @ state.basis
        query = self.interaction.pairfromprotein(
            state.proteinfeature, queryligand
        ) @ state.basis
        center = support.mean(dim=0, keepdim=True)
        return support - center, query - center

    def predict(
        self,
        state: ReorderingState,
        *,
        queryligand: torch.Tensor,
        supportresidual: torch.Tensor,
        supportbase: torch.Tensor,
        querybase: torch.Tensor,
        basevariance: torch.Tensor,
        observationvariance: torch.Tensor,
        rankingmode: str,
        applycalibration: bool = False,
        querycalibrationbase: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        if self.interaction is None:
            raise RuntimeError("protein prediction requires a protein-conditioned model")
        queryfeature = self.interaction.pairfromprotein(
            state.proteinfeature, queryligand
        )
        return self.posterior(
            queryfeature,
            state.supportfeature,
            supportresidual,
            querybase,
            state.basis,
            state.priorcov,
            basevariance=basevariance,
            observationvariance=observationvariance,
            supportbase=supportbase,
            querycalibrationbase=querycalibrationbase,
            rankingmode=rankingmode,
            applycalibration=applycalibration,
        )

__all__ = [
    "ProteinSubspace",
    "ReorderingModel",
    "ReorderingPosterior",
    "ReorderingState",
    "calibrationgeometry",
    "helmert",
    "modelprobability",
    "stiefelqr",
]
