"""V2 core candidate: calibration-separated section trained by crossed contrast.

Motivation, derived from the recorded failures rather than from architecture
taste:

1.  The v0 gain is predominantly target-level calibration (R2-E0).
2.  Replacing support **and** query protein with the same wrong protein restores
    v0 performance (V1 2x2 factorial).  The section therefore uses a
    self-consistent coordinate system, not partner identity.
3.  Fact 2 is not an accident.  Write the pair coordinate as an affine function
    of a ligand embedding, ``m(P, L) = A_P phi(L) + b_P``.  Under the centred
    section the intercept ``b_P`` cancels **exactly** and the prediction depends
    on the protein only through the Gram metric ``G_P = A_P^T A_P``.  If ``G_P``
    is constant in ``P`` then correct/correct equals wrong/wrong identically.

    Consequence: making the section calibration-orthogonal is necessary to stop
    mis-attributing the offset, but it *removes the protein's additive channel*
    and therefore cannot on its own create partner specificity.  The protein
    must modulate the within-task ligand **metric**, and something in training
    must force it to.

This module implements the two coupled pieces.

``AnisotropicPairCoordinate``
    ``m(P, L) = W(P)^T phi(L)`` with ``W(P) = W_0 + sum_j g_j(P) W_j`` and
    ``||g(P)||_2 <= 1``.  The protein selects a member of a small learned
    family of interaction bases; it does not get a free embedding, and no new
    pooled pair head or wider MLP is introduced.  The reported kernel
    ``W(P) W(P)^T`` is invariant to the latent gauge ``W(P) -> W(P) O``.

``crossed_contrast_loss``  (X-CON)
    On measured rectangles ``{P1, P2} x {L1, L2}`` with all four cells measured
    under the same panel/assay, supervise the model's *own few-shot* second
    difference against the measured one:

        [yhat(P1,L1) - yhat(P1,L2)] - [yhat(P2,L1) - yhat(P2,L2)]
            vs.  y11 - y12 - y21 + y22

    Exactly cancelled from this quantity: the ligand population ``mu_L``, the
    per-task intercept ``b_t``, and every protein main effect.  What survives is
    ``G_{P1} - G_{P2}``.  The objective therefore has zero gradient through the
    calibration channel and the ligand prior, and non-zero gradient only where
    the representation makes the within-task ligand metric protein-specific.
    Only measured cells are used; no unmeasured pair is labelled a non-binder.

Nothing here is authorized for ``model/`` or ``scripts/``.  Precedents that
must be cited when reporting: within-assay relative ranking (MBP), direct
target-difference selectivity regression, PCM/IMC/Macau side-information
factorization, R2-D2/ALPaCA closed-form adaptation, and the meta-learning
memorization literature.  The claimed novelty is the *second-order* crossed
contrast taken through a support-identifiable centred section, not any single
component.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import Tensor, nn


# --------------------------------------------------------------------------- #
# protein-anisotropic pair coordinate
# --------------------------------------------------------------------------- #
class AnisotropicPairCoordinate(nn.Module):
    """``m(P, L) = W(P)^T phi(L)``, ``W(P) = W_0 + sum_j g_j(P) W_j``.

    ``family_size`` is the number of learned interaction bases the protein may
    select between.  It is deliberately small: the protein contributes a
    bounded, low-dimensional selection, never an unconstrained embedding.
    """

    def __init__(self, ligand_dim: int, section_dim: int, protein_dim: int,
                 family_size: int = 8) -> None:
        super().__init__()
        if not 1 <= section_dim <= 5:
            raise ValueError("section_dim must satisfy 1 <= d <= 5 (support budget)")
        if family_size < 1:
            raise ValueError("family_size must be positive")
        self.section_dim = int(section_dim)
        self.family_size = int(family_size)
        scale = ligand_dim ** -0.5
        self.base = nn.Parameter(torch.randn(ligand_dim, section_dim) * scale)
        self.family = nn.Parameter(
            torch.randn(family_size, ligand_dim, section_dim) * scale)
        self.selector = nn.Linear(protein_dim, family_size)

    def gate(self, protein_features: Tensor) -> Tensor:
        """Bounded family selector: ``||g||_2 <= 1`` by construction."""
        raw = torch.tanh(self.selector(protein_features))
        norm = raw.norm(dim=-1, keepdim=True).clamp(min=1.0)
        return raw / norm

    def maps(self, protein_features: Tensor) -> Tensor:
        """Return ``W(P)`` with shape ``[n, ligand_dim, section_dim]``."""
        gate = self.gate(protein_features)
        return self.base + torch.einsum("nj,jhd->nhd", gate, self.family)

    def forward(self, protein_features: Tensor, ligand_features: Tensor) -> Tensor:
        return torch.einsum("nh,nhd->nd", ligand_features,
                            self.maps(protein_features))

    def metric(self, protein_features: Tensor) -> Tensor:
        """Gauge-invariant ligand metric ``G_P = W(P) W(P)^T``.

        Invariant under ``W(P) -> W(P) O`` for orthogonal ``O``, so it is a
        legitimate reported quantity under the frozen gauge-invariance
        constraint.  Returned for diagnostics only.
        """
        maps = self.maps(protein_features)
        return torch.einsum("nhd,ngd->nhg", maps, maps)


# --------------------------------------------------------------------------- #
# calibration-separated section
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SectionOutput:
    prediction: Tensor
    intercept: Tensor
    correction: Tensor
    rank: int


def centered_section(support_coordinates: Tensor, support_residual: Tensor,
                     query_coordinates: Tensor, ridge: float) -> SectionOutput:
    """Explicit intercept plus centred positive-ridge dual section.

    Label-dependent degrees of freedom per task: one intercept plus
    ``rank(M_c) <= k - 1``, so at most ``k`` in total.  This respects the frozen
    support-capacity ceiling (at most ``k`` continuous dimensions of member
    identity from ``k`` observations).
    """
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    if support_coordinates.ndim != 2 or query_coordinates.ndim != 2:
        raise ValueError("coordinates must be [n, d] matrices")
    if support_residual.ndim != 1 or len(support_residual) != len(support_coordinates):
        raise ValueError("support residual must be [k] and match the coordinates")
    k = len(support_residual)
    intercept = support_residual.mean()
    mean = support_coordinates.mean(dim=0, keepdim=True)
    centered_support = support_coordinates - mean
    centered_query = query_coordinates - mean
    identity = torch.eye(k, dtype=support_coordinates.dtype,
                         device=support_coordinates.device)
    gram = centered_support @ centered_support.T
    dual = torch.linalg.solve(gram + ridge * identity, support_residual - intercept)
    correction = centered_query @ (centered_support.T @ dual)
    return SectionOutput(
        prediction=intercept + correction,
        intercept=intercept,
        correction=correction,
        rank=int(torch.linalg.matrix_rank(centered_support).item()),
    )


class CrossedContrastSection(nn.Module):
    """Ligand population + anisotropic coordinate + calibration-separated section."""

    def __init__(self, ligand_dim: int, protein_dim: int, section_dim: int = 2,
                 ridge: float = 1.0, family_size: int = 8) -> None:
        super().__init__()
        self.ligand_population = nn.Linear(ligand_dim, 1)
        self.coordinate = AnisotropicPairCoordinate(
            ligand_dim, section_dim, protein_dim, family_size)
        self.ridge = float(ridge)

    def population(self, ligand_features: Tensor) -> Tensor:
        return self.ligand_population(ligand_features).squeeze(-1)

    def predict(self, support_protein: Tensor, support_ligand: Tensor,
                support_y: Tensor, query_protein: Tensor,
                query_ligand: Tensor) -> Tensor:
        support_mu = self.population(support_ligand)
        query_mu = self.population(query_ligand)
        support_coordinates = self.coordinate(support_protein, support_ligand)
        query_coordinates = self.coordinate(query_protein, query_ligand)
        section = centered_section(support_coordinates, support_y - support_mu,
                                   query_coordinates, self.ridge)
        return query_mu + section.prediction


# --------------------------------------------------------------------------- #
# X-CON objective
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Rectangle:
    """One measured 2x2 rectangle plus the two legitimate support episodes.

    ``protein_a``/``protein_b``: ``[protein_dim]`` feature rows.
    ``ligand_1``/``ligand_2``:   ``[ligand_dim]`` feature rows.
    ``values``:                  ``[2, 2]`` measured affinities ordered
                                 ``[[y_a1, y_a2], [y_b1, y_b2]]``.
    ``support_*``:               support ligand features ``[k, ligand_dim]`` and
                                 labels ``[k]``; the support must exclude both
                                 rectangle ligands and every scaffold in them.
    """

    protein_a: Tensor
    protein_b: Tensor
    ligand_1: Tensor
    ligand_2: Tensor
    values: Tensor
    support_a_ligand: Tensor
    support_a_y: Tensor
    support_b_ligand: Tensor
    support_b_y: Tensor


def measured_crossed_difference(values: Tensor) -> Tensor:
    """``y_a1 - y_a2 - y_b1 + y_b2``: removes every protein and ligand main effect."""
    if values.shape[-2:] != (2, 2):
        raise ValueError("values must end in a measured 2x2 rectangle")
    return (values[..., 0, 0] - values[..., 0, 1]
            - values[..., 1, 0] + values[..., 1, 1])


def predicted_crossed_difference(model: CrossedContrastSection,
                                 rectangle: Rectangle) -> Tensor:
    """Second difference of the model's own few-shot predictions.

    The per-task intercept cancels inside each protein and the ligand
    population cancels between the two proteins, so this quantity is
    identically free of both channels.
    """
    query_ligand = torch.stack([rectangle.ligand_1, rectangle.ligand_2])
    predictions = []
    for protein, support_ligand, support_y in (
        (rectangle.protein_a, rectangle.support_a_ligand, rectangle.support_a_y),
        (rectangle.protein_b, rectangle.support_b_ligand, rectangle.support_b_y),
    ):
        support_protein = protein.expand(len(support_ligand), -1)
        query_protein = protein.expand(len(query_ligand), -1)
        predictions.append(model.predict(
            support_protein, support_ligand, support_y,
            query_protein, query_ligand))
    first, second = predictions
    return (first[0] - first[1]) - (second[0] - second[1])


def crossed_contrast_loss(model: CrossedContrastSection,
                          rectangles: list[Rectangle]) -> Tensor:
    """X-CON: mean squared error on the measured second difference."""
    if not rectangles:
        raise ValueError("at least one measured rectangle is required")
    terms = []
    for rectangle in rectangles:
        predicted = predicted_crossed_difference(model, rectangle)
        measured = measured_crossed_difference(rectangle.values)
        terms.append((predicted - measured).square())
    return torch.stack(terms).mean()


def episodic_loss(model: CrossedContrastSection, support_protein: Tensor,
                  support_ligand: Tensor, support_y: Tensor,
                  query_protein: Tensor, query_ligand: Tensor,
                  query_y: Tensor) -> Tensor:
    """Ordinary outer few-shot query loss; retains the absolute scale."""
    prediction = model.predict(support_protein, support_ligand, support_y,
                               query_protein, query_ligand)
    return (prediction - query_y).square().mean()


def combined_objective(model: CrossedContrastSection, episode_terms: list[Tensor],
                       rectangles: list[Rectangle], weight: float) -> Tensor:
    """Absolute episodic loss plus the weighted crossed contrast.

    The absolute term is required because the crossed difference does not
    identify any offset; the crossed term is required because the absolute
    term is satisfiable by calibration alone.  Neither alone is sufficient and
    the weight is a declared hyperparameter selected on development data only.
    """
    if weight < 0:
        raise ValueError("weight must be non-negative")
    absolute = torch.stack(episode_terms).mean()
    return absolute + weight * crossed_contrast_loss(model, rectangles)


# --------------------------------------------------------------------------- #
# falsifiers
# --------------------------------------------------------------------------- #
def wrong_wrong_invariance(model: CrossedContrastSection, support_ligand: Tensor,
                           support_y: Tensor, query_ligand: Tensor,
                           correct_protein: Tensor,
                           wrong_protein: Tensor) -> Tensor:
    """Absolute prediction change when both support and query protein are swapped.

    The V1 factorial showed this is ~0 for the incumbent.  Any candidate that
    keeps it at ~0 has not acquired partner identity, whatever its MSE.
    """
    def run(protein: Tensor) -> Tensor:
        return model.predict(
            protein.expand(len(support_ligand), -1), support_ligand, support_y,
            protein.expand(len(query_ligand), -1), query_ligand)

    return (run(correct_protein) - run(wrong_protein)).abs()


def calibration_leakage(model: CrossedContrastSection, support_protein: Tensor,
                        support_ligand: Tensor, support_y: Tensor,
                        query_protein: Tensor, query_ligand: Tensor,
                        shift: float) -> Tensor:
    """Prediction change under a constant shift of every support label.

    For a calibration-separated section this must equal ``shift`` exactly: the
    constant enters the explicit intercept and nothing else.  Any deviation is
    calibration leaking into the section.
    """
    base = model.predict(support_protein, support_ligand, support_y,
                         query_protein, query_ligand)
    shifted = model.predict(support_protein, support_ligand, support_y + shift,
                            query_protein, query_ligand)
    return shifted - base
