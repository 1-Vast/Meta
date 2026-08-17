"""Stage W W1 local interaction operator.

    pocket_states = cross_attention(learned_queries, ordered_protein_slots)
    ligand_tokens = pharmacophore_atom_tokens(core, R_a, R_b)
    interaction   = cross_attention(ligand_tokens, pocket_states)
    R(tau, p)     = MLP(attention(learned_query, interaction))
    D_hat(tau,p1,p2) = R(tau,p1) - R(tau,p2)

The difference form gives identity, antisymmetry and cycle consistency for
every parameter setting. No pooled protein, target id, component id or assay
id enters the local operator.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import torch
import torch.nn as nn

PROTEIN_DIM = 640
SLOTS = 128
DIM = 128
POCKET_QUERIES = 8
ATOM_DIM = 21
CAT_DIM = 6
TOKEN_DIM = ATOM_DIM + CAT_DIM
MAX_ATOMS = 24
ELEMENTS = ("C", "N", "O", "S", "P", "F", "Cl", "Br", "I")


def atom_features(mol):
    out = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            continue
        symbol = atom.GetSymbol()
        element = [0.0] * (len(ELEMENTS) + 1)
        try:
            element[ELEMENTS.index(symbol)] = 1.0
        except ValueError:
            element[-1] = 1.0
        hyb = str(atom.GetHybridization())
        hyb_vec = [0.0] * 5
        hyb_vec[{"SP": 0, "SP2": 1, "SP3": 2}.get(hyb, 3 if hyb else 4)] = 1.0
        charge = atom.GetFormalCharge()
        charge_vec = [1.0 if charge < 0 else 0.0,
                      1.0 if charge == 0 else 0.0,
                      1.0 if charge > 0 else 0.0]
        out.append(element + [
            float(atom.GetIsAromatic()), float(atom.IsInRing()),
            float(atom.GetDegree()) / 4.0] + charge_vec + hyb_vec)
    return out


def ligand_tokens(murcko_core: str, category_a, category_b):
    """Pharmacophore atom tokens for the soft-family context."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(murcko_core)
    features = atom_features(mol) if mol is not None else []
    if len(features) > MAX_ATOMS:
        features = features[:MAX_ATOMS]
    n_atoms = len(features)
    tokens = []
    for vector in features:
        tokens.append(vector + [0.0] * (TOKEN_DIM - len(vector)))
    def cat_vector(category):
        values = list(category)
        if len(values) < 6:
            values = values + [0] * (6 - len(values))
        return [float(values[0]) / 4.0, float(values[1]), float(values[2]),
                float(values[3]) / 2.0, float(values[4]) / 2.0,
                float(values[5] - 1) if values[5] in (0, 1, 2) else 0.0]
    tokens.append([0.0] * ATOM_DIM + cat_vector(category_a))
    tokens.append([0.0] * ATOM_DIM + cat_vector(category_b))
    tensor = torch.tensor(np.asarray(tokens, dtype=np.float32))
    mask = torch.ones(tensor.shape[0], dtype=torch.bool)
    return tensor, mask


def sinusoidal_positions(length: int, dim: int) -> torch.Tensor:
    position = torch.arange(length, dtype=torch.float32).unsqueeze(1)
    div = torch.exp(torch.arange(0, dim, 2, dtype=torch.float32)
                    * (-math.log(10000.0) / dim))
    table = torch.zeros(length, dim)
    table[:, 0::2] = torch.sin(position * div)
    table[:, 1::2] = torch.cos(position * div)
    return table


class _FFN(nn.Module):
    def __init__(self, dim: int, hidden: int = 256, dropout: float = 0.0):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, hidden), nn.SiLU(), nn.Dropout(dropout),
            nn.Linear(hidden, dim), nn.Dropout(dropout))

    def forward(self, x):
        return self.net(x)


class LocalResponse(nn.Module):
    """R(tau, p): ligand pharmacophore tokens x multiple latent pocket states."""

    def __init__(self, protein_dim: int = PROTEIN_DIM, dim: int = DIM,
                 atom_dim: int = TOKEN_DIM, slots: int = SLOTS,
                 pocket_queries: int = POCKET_QUERIES, dropout: float = 0.0):
        super().__init__()
        self.dim = dim
        self.protein_proj = nn.Linear(protein_dim, dim)
        self.protein_norm = nn.LayerNorm(dim)
        self.register_buffer("protein_pos",
                             sinusoidal_positions(slots, dim).unsqueeze(0))
        self.ligand_proj = nn.Linear(atom_dim, dim)
        ligand_layer = nn.TransformerEncoderLayer(
            d_model=dim, nhead=4, dim_feedforward=256, dropout=dropout,
            batch_first=True, activation="gelu", norm_first=True)
        self.ligand_encoder = nn.TransformerEncoder(ligand_layer, num_layers=2)
        self.pocket_queries = nn.Parameter(torch.randn(pocket_queries, dim) * 0.02)
        pocket_layer = nn.TransformerEncoderLayer(
            d_model=dim, nhead=4, dim_feedforward=256, dropout=dropout,
            batch_first=True, activation="gelu", norm_first=True)
        self.pocket_encoder = nn.TransformerEncoder(pocket_layer, num_layers=2)
        self.interaction_attn = nn.MultiheadAttention(
            dim, 4, dropout=dropout, batch_first=True)
        self.interaction_norm = nn.LayerNorm(dim)
        self.interaction_query = nn.Parameter(torch.zeros(1, dim))
        self.readout_attn = nn.MultiheadAttention(
            dim, 4, dropout=dropout, batch_first=True)
        self.head = nn.Sequential(
            nn.Linear(dim, 128), nn.SiLU(), nn.Linear(128, 1, bias=False))

    def protein_tokens(self, slots, mask):
        tokens = self.protein_proj(slots) + self.protein_pos
        tokens = self.protein_norm(tokens)
        return tokens

    def ligand_encode(self, ligand, ligand_mask):
        tokens = self.ligand_proj(ligand)
        return self.ligand_encoder(tokens, src_key_padding_mask=~ligand_mask)

    def pocket_states(self, protein_tokens, protein_mask):
        queries = self.pocket_queries.unsqueeze(0).expand(
            protein_tokens.shape[0], -1, -1)
        # pocket queries attend over protein slot tokens (keys/values).
        attended, _ = self.interaction_attn(
            queries, protein_tokens, protein_tokens,
            key_padding_mask=~protein_mask)
        states = self.interaction_norm(queries + attended)
        # refine pocket states with self-attention.
        states = self.pocket_encoder(states)
        return states

    def forward(self, ligand, ligand_mask, protein_slots, protein_mask):
        ligand_encoded = self.ligand_encode(ligand, ligand_mask)
        protein = self.protein_tokens(protein_slots, protein_mask)
        pockets = self.pocket_states(protein, protein_mask)
        # ligand tokens query pocket states.
        interacted, _ = self.interaction_attn(
            ligand_encoded, pockets, pockets)
        interacted = self.interaction_norm(ligand_encoded + interacted)
        query = self.interaction_query.unsqueeze(0).expand(
            interacted.shape[0], -1, -1)
        context, _ = self.readout_attn(query, interacted, interacted,
                                       key_padding_mask=~ligand_mask)
        return self.head(context.squeeze(1)).squeeze(-1)


class GlobalResponse(nn.Module):
    """Arm B: edit summary + global ESM pooled protein summary."""

    def __init__(self, protein_dim: int = PROTEIN_DIM, dim: int = DIM,
                 atom_dim: int = TOKEN_DIM):
        super().__init__()
        self.ligand_proj = nn.Linear(atom_dim, dim)
        self.protein_proj = nn.Linear(protein_dim, dim)
        self.head = nn.Sequential(
            nn.Linear(2 * dim, 128), nn.SiLU(), nn.Linear(128, 1, bias=False))

    def forward(self, ligand, ligand_mask, protein_slots, protein_mask):
        ligand_summary = (self.ligand_proj(ligand) * ligand_mask.unsqueeze(-1)
                          ).sum(1) / ligand_mask.sum(1, keepdim=True).clamp(min=1)
        weights = protein_mask.float().unsqueeze(-1)
        protein_summary = (protein_slots * weights).sum(1) / weights.sum(
            1).clamp(min=1)
        protein_summary = self.protein_proj(protein_summary)
        return self.head(torch.cat((ligand_summary, protein_summary), dim=-1)
                         ).squeeze(-1)


class ZeroResponse(nn.Module):
    def forward(self, ligand, ligand_mask, protein_slots, protein_mask):
        return torch.zeros(ligand.shape[0], device=ligand.device)


class DoubleDifferenceModel(nn.Module):
    def __init__(self, mode: str = "local"):
        super().__init__()
        self.mode = mode
        if mode == "zero":
            self.response = ZeroResponse()
        elif mode == "global":
            self.response = GlobalResponse()
        elif mode == "local":
            self.response = LocalResponse()
        else:
            raise ValueError(mode)

    def forward(self, ligand, ligand_mask, protein_left, mask_left,
                protein_right, mask_right):
        left = self.response(ligand, ligand_mask, protein_left, mask_left)
        right = self.response(ligand, ligand_mask, protein_right, mask_right)
        return left - right
