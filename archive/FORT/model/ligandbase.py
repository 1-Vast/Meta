"""Ligand-only predictor with a target-level residual posterior.

The target enters ONLY through the B1 support likelihood -- there is no target-identity feature path,
which is the purest form of the dual-cold contract. At k=0 the prediction is the ligand-descriptor base;
at k>0 the exact target-level Bayesian residual posterior (Cholesky solve, epistemic variance) adapts it
from globally-cold support labels. Support permutation invariance, the strict k=0 fallback and Cholesky
stability are inherited from BayesianResidualPosterior.
"""
from __future__ import annotations

import torch
from torch import nn

from .ligand import FingerprintEncoder
from .posterior import BayesianResidualPosterior


class LigandBaseline(nn.Module):
    # NOTE (audit 2026-07-25): the real ligand tensor is 1034 = Morgan-1024 + 10 standardized
    # physicochemical descriptors, concatenated into ONE vector (not projected separately). The old
    # 1024 default silently dropped the descriptor block for any caller that did not override it;
    # `bits` now defaults to the audited feature width.
    def __init__(self, bits: int = 1034, rep: int = 128, pdim: int = 64) -> None:
        super().__init__()
        self.ligand = FingerprintEncoder(bits, rep)
        self.base_head = nn.Linear(rep, 1)
        self.logvar_head = nn.Linear(rep, 1)
        self.bayes = BayesianResidualPosterior(rep=rep, pdim=pdim)

    def forward(
        self,
        query_fp: torch.Tensor,
        support_fp: torch.Tensor | None = None,
        support_y: torch.Tensor | None = None,
        query_base: torch.Tensor | None = None,
        support_base: torch.Tensor | None = None,
    ) -> dict:
        query_rep = self.ligand(query_fp)
        base_q = (
            self.base_head(query_rep).squeeze(-1)
            if query_base is None
            else query_base.float()
        )
        logvar = self.logvar_head(query_rep).squeeze(-1)
        if support_fp is None or support_fp.shape[0] == 0:
            out = self.bayes(query_rep, None, None, base_q)
        else:
            support_rep = self.ligand(support_fp)
            base_s = (
                self.base_head(support_rep).squeeze(-1)
                if support_base is None
                else support_base.float()
            )
            residual_s = support_y.reshape(-1) - base_s
            out = self.bayes(query_rep, support_rep, residual_s, base_q)
        aleatoric = torch.exp(logvar.clamp(-8.0, 8.0))
        return {
            "pred": out["pred"], "base": base_q, "logvar": logvar,
            "epistemic": out["epistemic"], "aleatoric": aleatoric,
            "total_variance": aleatoric + out["epistemic"], "effective_k": out["effective_k"],
        }


__all__ = ["LigandBaseline"]
