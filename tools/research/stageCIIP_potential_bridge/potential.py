"""Stage CIIP-1 unified potential module (single deployable function).

f_theta(P, L) = b_P(P) + b_L(L) + s_theta(P, L),  s = alpha(P)^T psi(L).

ALL contrasts come from the SAME s_theta:
- ligand contrast   g(P, A->B) = s(P, B) - s(P, A)
- protein contrast  g(A->B, L) = s(B, L) - s(A, L)
- CIIP double contrast D = [s(Pa,B)-s(Pa,A)] - [s(Pb,B)-s(Pb,A)]
Identity-zero, antisymmetry and cycle-zero hold BY CONSTRUCTION for the
integrable potential; the free pairwise arms (diagnostic only) do not.

No target ID, no closed-form solver, no test-time gradients anywhere in
this module. Stable SHA-256 keyed rng (x0_common.stable_rng).
"""
from __future__ import annotations

import torch
import torch.nn as nn

D_P = 1700   # KLIFS pocket one-hot (85 positions x 20 aa), frozen X0 asset
D_L = 2048   # ECFP4, frozen X0 asset
HID = 64
RANK = 8


class UnifiedPotential(nn.Module):
    """f = b_P + b_L + s; s = alpha(P)^T psi(L)."""

    def __init__(self, d_p=D_P, d_l=D_L, hid=HID, rank=RANK):
        super().__init__()
        self.p_enc = nn.Linear(d_p, hid)
        self.l_enc = nn.Linear(d_l, hid)
        self.b_P = nn.Linear(hid, 1)
        self.b_L = nn.Linear(hid, 1)
        self.alpha = nn.Linear(hid, rank)
        self.psi = nn.Linear(hid, rank)
        self.mu = nn.Parameter(torch.zeros(1))
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                nn.init.zeros_(m.bias)

    def _enc(self, P, L):
        return torch.relu(self.p_enc(P)), torch.relu(self.l_enc(L))

    def potential(self, P, L):
        ep, el = self._enc(P, L)
        return (self.alpha(ep) * self.psi(el)).sum(-1)

    def forward(self, P, L):
        ep, el = self._enc(P, L)
        return (self.mu + self.b_P(ep).squeeze(-1)
                + self.b_L(el).squeeze(-1) + self.potential(P, L))

    # --- all contrasts from the same s_theta ---
    def ligand_contrast(self, P, La, Lb):
        """g(P, A->B) = s(P,B) - s(P,A)."""
        return self.potential(P, Lb) - self.potential(P, La)

    def protein_contrast(self, Pa, Pb, L):
        """g(A->B, L) = s(B,L) - s(A,L)."""
        return self.potential(Pb, L) - self.potential(Pa, L)

    def double_contrast(self, Pa, Pb, La, Lb):
        """CIIP crossed double difference."""
        return self.ligand_contrast(Pb, La, Lb) - self.ligand_contrast(Pa, La, Lb)

    def centered_mutation_effect(self, P_wt, P_v, L_mat):
        """c_hat(v) per ligand = d_hat - mean_l(d_hat), d_hat =
        s(P_v,L) - s(P_wt,L). Returns (n_ligands,)."""
        d = self.protein_contrast(P_wt, P_v, L_mat)
        return d - d.mean()


class GlobalPotential(nn.Module):
    """Global-compression diagnostic: P -> mean over pocket positions
    (20-dim aa composition) before the same alpha/psi machinery. Same
    capacity class, different protein representation."""

    def __init__(self, hid=HID, rank=RANK):
        super().__init__()
        self.p_enc = nn.Linear(20, hid)
        self.l_enc = nn.Linear(D_L, hid)
        self.b_P = nn.Linear(hid, 1)
        self.b_L = nn.Linear(hid, 1)
        self.alpha = nn.Linear(hid, rank)
        self.psi = nn.Linear(hid, rank)
        self.mu = nn.Parameter(torch.zeros(1))
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                nn.init.zeros_(m.bias)

    @staticmethod
    def reduce(P):
        """(n, 1700) one-hot -> (n, 20) pocket aa composition."""
        return P.view(-1, 85, 20).sum(1)

    def potential(self, P, L):
        ep = torch.relu(self.p_enc(self.reduce(P)))
        el = torch.relu(self.l_enc(L))
        return (self.alpha(ep) * self.psi(el)).sum(-1)

    def forward(self, P, L):
        ep = torch.relu(self.p_enc(self.reduce(P)))
        el = torch.relu(self.l_enc(L))
        return (self.mu + self.b_P(ep).squeeze(-1)
                + self.b_L(el).squeeze(-1) + self.potential(P, L))

    def ligand_contrast(self, P, La, Lb):
        return self.potential(P, Lb) - self.potential(P, La)

    def protein_contrast(self, Pa, Pb, L):
        return self.potential(Pb, L) - self.potential(Pa, L)

    def centered_mutation_effect(self, P_wt, P_v, L_mat):
        d = self.protein_contrast(P_wt, P_v, L_mat)
        return d - d.mean()


class FreePairwise(nn.Module):
    """DIAGNOSTIC ONLY: free protein-pair predictor for the mutation
    effect. g_free = h([P_wt,P_v,L]) - h([P_v,P_wt,L]). Antisymmetric
    by construction but NOT integrable (cycle test fails); never a
    production mechanism."""

    def __init__(self, d_p=D_P, d_l=D_L, hid=HID):
        super().__init__()
        self.h = nn.Sequential(nn.Linear(2 * d_p + d_l, hid), nn.ReLU(),
                               nn.Linear(hid, hid), nn.ReLU(),
                               nn.Linear(hid, 1))
        for m in [self.h[0], self.h[2]]:
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
        nn.init.xavier_uniform_(self.h[4].weight)
        nn.init.zeros_(self.h[4].bias)

    def __call__(self, P_wt, P_v, L):
        fwd = self.h(torch.cat([P_wt, P_v, L], -1))
        bwd = self.h(torch.cat([P_v, P_wt, L], -1))
        return (fwd - bwd).squeeze(-1)


class FreeLigandPair(nn.Module):
    """DIAGNOSTIC ONLY: free ligand-pair predictor (ligand contrast)."""

    def __init__(self, d_p=D_P, d_l=D_L, hid=HID):
        super().__init__()
        self.h = nn.Sequential(nn.Linear(d_p + 2 * d_l, hid), nn.ReLU(),
                               nn.Linear(hid, hid), nn.ReLU(),
                               nn.Linear(hid, 1))
        for m in [self.h[0], self.h[2]]:
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
        nn.init.xavier_uniform_(self.h[4].weight)
        nn.init.zeros_(self.h[4].bias)

    def __call__(self, P, La, Lb):
        fwd = self.h(torch.cat([P, La, Lb], -1))
        bwd = self.h(torch.cat([P, Lb, La], -1))
        return (fwd - bwd).squeeze(-1)
