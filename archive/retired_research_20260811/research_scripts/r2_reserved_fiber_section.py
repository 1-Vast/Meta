"""R2-M1: Reserved-Fiber Meta-Section (RFMS).

Motivation
----------
In v0 the section coordinates are a function of BOTH partner and ligand,

    m(P, L) = B^T phi(P, L),        yhat_q = pop(L_q) + m_q . w,
    w = argmin_w ||M w - rho||^2 + lam ||w||^2 .

For a *fixed support residual*, the unregularized ridge kernel is invariant
under a common episode-wise G in GL(d), and positive ridge preserves exact
orthogonal invariance.  The complete v0 predictor has an additional learned
population-coordinate term, however, and a wrong-protein substitution need not
be such a common transform.  R2-E1 in fact falsified the registered near-gauge
hypothesis, so this algebra is a boundary condition rather than an explanation
of the observed wrong/wrong result.

RFMS changes the two enabling conditions, but does not by itself prove partner
identification.

  (1) The frame is protein-independent.  psi(L) in R^d is a function of the
      LIGAND ONLY, shared across every target. There is no protein-dependent
      basis, hence no protein-dependent gauge to cancel.

  (2) The adaptation is deliberately rank-deficient.  R^d = V_S (+) V_C with
      dim V_S = d_s <= k and dim V_C = d_c = d - d_s >= 1. The support cuts only
      V_S. The V_C component is supplied solely by the partner map c0(P) and is
      not reachable by any support fit.

Prediction and certificate:

    yhat(L) = psi(L) . ( c0(P) + Pi_S delta ),   delta = ridge fit on support
    realizability rho = || y_s - Psi_s (c0 + Pi_S delta) ||_inf

`rho` below is the residual of the ridge solution. It is only an upper bound on
the true L-infinity fiber feasibility residual and must not be called an
empty-fiber certificate.

Candidate constraints
---------------------
  * d <= 5 and d_s <= k retain the R2 research budget; these are not claimed
    as consequences of the authoritative frozen CSMO operator.
  * c0(P) in R^{d_c}, d_c <= 4: this is NOT a high-dimensional target embedding.
  * Support enters only through the section cut; permutation symmetry holds.
  * Only block-preserving orthogonal basis changes preserve the hard split and
    ridge penalty without further transformation.

Falsifier
---------
Wrong/wrong can still recover when the reserved partner difference has zero
query exposure or is compensated through the adaptable span. A future audit
must therefore report the reserved quotient exposure, not only between-target
variance of c0.

Usage
-----
    python -m research.meta_fewshot.r2_reserved_fiber_section \
        --protein-features dataset/processed/multipanel/esm2_t30_bdb.npz \
        --d 5 --d-support 2 --k 5
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import torch
from torch import Tensor, nn

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "report/meta_fewshot/r2_reserved_fiber.json"
SCHEMA = "MetaSieve.R2ReservedFiber.v1"


def reserved_exposure_from_frames(adaptable_s: Tensor, reserved_s: Tensor,
                                  adaptable_q: Tensor, reserved_q: Tensor,
                                  ridge: float) -> Tensor:
    """Return Xi = Bq - Aq(As^T As + ridge I)^-1 As^T Bs."""
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    gram = adaptable_s.T @ adaptable_s
    identity = torch.eye(gram.shape[0], dtype=gram.dtype, device=gram.device)
    transfer = torch.linalg.solve(gram + ridge * identity,
                                  adaptable_s.T @ reserved_s)
    return reserved_q - adaptable_q @ transfer


class ReservedFiberSection(nn.Module):
    """Protein-independent ligand frame + partner-pinned reserved subspace."""

    def __init__(self, ligand_dim: int, protein_dim: int, d: int, d_support: int,
                 ridge: float = 1.0) -> None:
        super().__init__()
        if not 1 <= d <= 5:
            raise ValueError("d must satisfy 1 <= d <= 5 (support capacity ceiling)")
        if not 1 <= d_support < d:
            raise ValueError("need 1 <= d_support < d so that d_c >= 1")
        if ridge <= 0:
            raise ValueError("ridge must be strictly positive")
        self.d, self.d_support, self.ridge = d, d_support, float(ridge)
        # frame: ligand only. no protein argument anywhere in this module.
        self.frame = nn.Sequential(
            nn.Linear(ligand_dim, 128), nn.GELU(), nn.Linear(128, d))
        self.offset = nn.Linear(ligand_dim, 1)
        # partner map: low capacity, writes only into the reserved subspace.
        self.partner = nn.Sequential(
            nn.Linear(protein_dim, 64), nn.GELU(), nn.Linear(64, d - d_support))

    def psi(self, ligand: Tensor) -> Tensor:
        return self.frame(ligand)

    def c0(self, protein: Tensor) -> Tensor:
        """Coefficient vector with the adaptable block pinned to zero."""
        reserved = self.partner(protein)
        pad = torch.zeros(*reserved.shape[:-1], self.d_support,
                          dtype=reserved.dtype, device=reserved.device)
        return torch.cat([pad, reserved], dim=-1)

    def single_c0(self, protein: Tensor) -> Tensor:
        """Normalize one partner represented as `(p,)` or `(1, p)` to `(d,)`."""
        coefficient = self.c0(protein)
        if coefficient.ndim == 1:
            return coefficient
        if coefficient.ndim == 2 and coefficient.shape[0] == 1:
            return coefficient.squeeze(0)
        raise ValueError("an episode requires exactly one support/query protein")

    def episode(self, support_ligand: Tensor, support_y: Tensor, protein: Tensor,
                query_ligand: Tensor, query_protein: Tensor | None = None
                ) -> tuple[Tensor, Tensor]:
        if len(support_ligand) < self.d_support:
            raise ValueError("support size must be at least d_support")
        psi_s = self.psi(support_ligand)
        psi_q = self.psi(query_ligand)
        base_s = self.offset(support_ligand).squeeze(-1)
        base_q = self.offset(query_ligand).squeeze(-1)
        support_coefficient = self.single_c0(protein)
        query_coefficient = self.single_c0(
            query_protein if query_protein is not None else protein)
        residual = support_y - base_s - psi_s @ support_coefficient
        # the support may only move the first d_support coordinates
        adaptable = psi_s[:, : self.d_support]
        gram = adaptable.T @ adaptable
        identity = torch.eye(self.d_support, dtype=gram.dtype, device=gram.device)
        delta = torch.linalg.solve(gram + self.ridge * identity, adaptable.T @ residual)
        support_full = support_coefficient + torch.cat(
            [delta, torch.zeros(self.d - self.d_support, dtype=delta.dtype,
                                device=delta.device)])
        query_full = query_coefficient + torch.cat(
            [delta, torch.zeros(self.d - self.d_support, dtype=delta.dtype,
                                device=delta.device)])
        prediction = base_q + psi_q @ query_full
        ridge_residual_upper_bound = (support_y - base_s - psi_s @ support_full).abs().max()
        return prediction, ridge_residual_upper_bound

    def reserved_exposure(self, support_ligand: Tensor, query_ligand: Tensor) -> Tensor:
        """Xi mapping a shared reserved coefficient change to query predictions.

        Xi = b - a (A^T A + lambda I)^-1 A^T C.  Nonzero c0 variance is
        insufficient if Xi annihilates the observed partner difference.
        """
        psi_s, psi_q = self.psi(support_ligand), self.psi(query_ligand)
        adaptable_s, reserved_s = psi_s[:, :self.d_support], psi_s[:, self.d_support:]
        adaptable_q, reserved_q = psi_q[:, :self.d_support], psi_q[:, self.d_support:]
        return reserved_exposure_from_frames(
            adaptable_s, reserved_s, adaptable_q, reserved_q, self.ridge)


def episode_loss(model: ReservedFiberSection, batch):
    """Query MSE only; wrong-protein controls are never training targets."""
    total = 0.0
    for support_ligand, support_y, protein, query_ligand, query_y in batch:
        prediction, rho_correct = model.episode(support_ligand, support_y, protein, query_ligand)
        total = total + ((prediction - query_y) ** 2).mean()
    return total / len(batch)


ARMS = ("correct", "zero_support", "foreign_support", "permuted_support",
        "wrong_protein_support_only", "wrong_protein_query_only",
        "wrong_protein_both", "ligand_only")


def describe() -> dict:
    return {
        "schema": SCHEMA,
        "status": "RFMS_TRAINING_BLOCKED_BY_R2_E1_FALSIFIER",
        "mechanism": "RESERVED_FIBER_META_SECTION",
        "invariance_boundary": (
            "the hard V_S/V_C split removes arbitrary block-mixing gauge; "
            "block-preserving orthogonal changes remain reparameterizations"),
        "wrong_wrong_boundary": (
            "a break requires nonzero Xi@(c_correct-c_wrong); c0 variance alone "
            "does not guarantee a prediction difference"),
        "required_controls": list(ARMS),
        "gates": {
            "meta_effect": "RFMS < d0",
            "support_specificity": "correct < zero, foreign, permuted",
            "partner_specificity": "correct < wrong_protein_both  (REQUIRED)",
            "ligand_specificity": "Pearson/Spearman/CI improve, not only MSE",
            "independence": "cluster/component-macro, not pair bootstrap",
            "no_collapse": "c0_between_target_variance > 0",
            "reserved_exposure": "Xi has nonzero singular values on observed queries",
            "functional_partner_gap": "Xi@(c_correct-c_wrong) is nonzero and measured-pair aligned",
        },
        "forbidden_training_signal": "random wrong-protein margin or non-binder label",
        "fiber_certificate_required": (
            "solve inf_delta ||y-base-Cc-A delta||_inf against a frozen noise epsilon; "
            "ridge residual is only an upper bound"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-features", type=str, default=None)
    parser.add_argument("--d", type=int, default=5)
    parser.add_argument("--d-support", type=int, default=2)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--describe", action="store_true")
    args = parser.parse_args()
    if args.d_support > args.k:
        parser.error("d-support cannot exceed k")
    if args.describe or args.protein_features is None:
        payload = describe() | {"note": "pass --protein-features to train"}
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    raise SystemExit(
        "RFMS training is blocked: R2 E1 falsified H0-regime and the candidate "
        "requires a new preregistration with quotient-exposure and exact-fiber controls")


if __name__ == "__main__":
    raise SystemExit(main())
