"""Small O(3)-equivariant Cartesian mechanism encoder.

This module is optional: coordinates are never synthesized from the active
sequence/2D DTA banks.  It encodes one coordinate frame at a time; independent
protein and ligand frames must be reduced to invariant slots before fusion.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
from torch import Tensor


@dataclass(frozen=True)
class CartesianTensorState:
    scalar: Tensor
    vector: Tensor
    tensor2: Tensor
    node_mask: Tensor


@dataclass(frozen=True)
class GeometryMechanismEncoding:
    state: CartesianTensorState
    invariant_nodes: Tensor
    mechanism_slots: Tensor
    geometry_available: Tensor


class SparseCartesianLayer(nn.Module):
    """Sparse scalar/vector/symmetric-traceless rank-2 message passing."""

    def __init__(self, scalar_dim: int, vector_channels: int,
                 tensor_channels: int, rbf_count: int = 16,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.vector_channels = int(vector_channels)
        self.tensor_channels = int(tensor_channels)
        self.register_buffer(
            "rbf_centers", torch.linspace(0.0, 12.0, rbf_count, dtype=dtype))
        self.radial = nn.Sequential(
            nn.Linear(rbf_count, scalar_dim, dtype=dtype), nn.SiLU())
        self.coefficient = nn.Linear(
            2 * scalar_dim,
            scalar_dim + 2 * vector_channels + 2 * tensor_channels,
            dtype=dtype)
        invariant_dim = scalar_dim + vector_channels + tensor_channels
        self.scalar_update = nn.Sequential(
            nn.LayerNorm(scalar_dim + invariant_dim, dtype=dtype),
            nn.Linear(scalar_dim + invariant_dim, 2 * scalar_dim, dtype=dtype),
            nn.SiLU(), nn.Linear(2 * scalar_dim, scalar_dim, dtype=dtype))
        self.vector_mix = nn.Linear(
            vector_channels, vector_channels, bias=False, dtype=dtype)
        self.tensor_mix = nn.Linear(
            tensor_channels, tensor_channels, bias=False, dtype=dtype)

    def forward(self, state: CartesianTensorState, coordinates: Tensor,
                edge_index: Tensor) -> CartesianTensorState:
        scalar, vector, tensor2, mask = (
            state.scalar, state.vector, state.tensor2, state.node_mask)
        batch, nodes = scalar.shape[:2]
        if coordinates.shape != (batch, nodes, 3):
            raise ValueError("coordinates do not match Cartesian node state")
        if edge_index.ndim != 2 or edge_index.shape[0] != 2:
            raise ValueError("edge_index must have shape [2,E]")
        flat_scalar = scalar.flatten(0, 1)
        flat_vector = vector.flatten(0, 1)
        flat_tensor = tensor2.flatten(0, 1)
        flat_coord = coordinates.flatten(0, 1)
        flat_mask = mask.flatten().bool()
        source, target = edge_index.long()
        if source.numel() and (int(source.min()) < 0 or
                               int(torch.maximum(source.max(), target.max())) >= batch * nodes):
            raise ValueError("Cartesian edge index is out of range")
        if source.numel() and not torch.equal(
                torch.div(source, nodes, rounding_mode="floor"),
                torch.div(target, nodes, rounding_mode="floor")):
            raise ValueError(
                "Cartesian edges cannot connect different flattened samples")
        valid = flat_mask[source] & flat_mask[target]
        source, target = source[valid], target[valid]
        displacement = flat_coord[source] - flat_coord[target]
        distance = torch.linalg.vector_norm(displacement, dim=-1)
        unit = displacement / distance.clamp_min(1e-8).unsqueeze(-1)
        unit = torch.where(distance[:, None] > 1e-8, unit, torch.zeros_like(unit))
        width = (self.rbf_centers[1] - self.rbf_centers[0]).clamp_min(1e-6)
        rbf = torch.exp(-((distance[:, None] - self.rbf_centers) / width).square())
        coefficient = self.coefficient(torch.cat((
            flat_scalar[source], self.radial(rbf)), -1))
        sizes = (scalar.shape[-1], self.vector_channels, self.vector_channels,
                 self.tensor_channels, self.tensor_channels)
        c0, c1, g1, c2, g2 = torch.split(coefficient, sizes, -1)
        identity = torch.eye(3, device=scalar.device, dtype=scalar.dtype)
        quadrupole = (unit[:, :, None] * unit[:, None, :]
                      - identity.unsqueeze(0) / 3.0)
        message0 = c0
        message1 = c1[..., None] * unit[:, None] + g1[..., None] * flat_vector[source]
        message2 = c2[..., None, None] * quadrupole[:, None] \
            + g2[..., None, None] * flat_tensor[source]
        aggregate0 = torch.zeros_like(flat_scalar).index_add(0, target, message0)
        aggregate1 = torch.zeros_like(flat_vector).index_add(0, target, message1)
        aggregate2 = torch.zeros_like(flat_tensor).index_add(0, target, message2)
        degree = torch.zeros(batch * nodes, device=scalar.device, dtype=scalar.dtype)
        degree.index_add_(0, target, torch.ones_like(distance))
        degree = degree.clamp_min(1.0)
        aggregate0 = aggregate0 / degree[:, None]
        aggregate1 = aggregate1 / degree[:, None, None]
        aggregate2 = aggregate2 / degree[:, None, None, None]
        invariants = torch.cat((
            aggregate0,
            torch.linalg.vector_norm(aggregate1, dim=-1),
            torch.linalg.matrix_norm(aggregate2, dim=(-2, -1))), -1)
        flat_scalar = flat_scalar + self.scalar_update(torch.cat((
            flat_scalar, invariants), -1))
        flat_vector = flat_vector + self.vector_mix(
            aggregate1.transpose(-1, -2)).transpose(-1, -2)
        flat_tensor = flat_tensor + self.tensor_mix(
            aggregate2.permute(0, 2, 3, 1)).permute(0, 3, 1, 2)
        gate = flat_mask.to(scalar.dtype)
        return CartesianTensorState(
            (flat_scalar * gate[:, None]).reshape_as(scalar),
            (flat_vector * gate[:, None, None]).reshape_as(vector),
            (flat_tensor * gate[:, None, None, None]).reshape_as(tensor2),
            mask)


class SparseCartesianMechanismEncoder(nn.Module):
    """Optional coordinate encoder with an exact coordinate-free fallback."""

    def __init__(self, input_dim: int, scalar_dim: int = 32,
                 vector_channels: int = 16, tensor_channels: int = 8,
                 layers: int = 2, rbf_count: int = 16,
                 mechanism_slots: int = 8, slot_dim: int = 32,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.vector_channels = int(vector_channels)
        self.tensor_channels = int(tensor_channels)
        self.scalar = nn.Sequential(
            nn.Linear(input_dim, scalar_dim, dtype=dtype), nn.SiLU(),
            nn.LayerNorm(scalar_dim, dtype=dtype))
        self.layers = nn.ModuleList(
            SparseCartesianLayer(
                scalar_dim, vector_channels, tensor_channels, rbf_count, dtype)
            for _ in range(layers))
        invariant_dim = scalar_dim + vector_channels + tensor_channels
        self.slot_seed = nn.Parameter(
            torch.empty(mechanism_slots, slot_dim, dtype=dtype))
        nn.init.normal_(self.slot_seed, std=slot_dim ** -0.5)
        self.invariant_projection = nn.Linear(
            invariant_dim, slot_dim, dtype=dtype)
        heads = next(h for h in (8, 4, 2, 1) if slot_dim % h == 0)
        self.slot_attention = nn.MultiheadAttention(
            slot_dim, heads, batch_first=True, dtype=dtype)

    def forward(self, features: Tensor, node_mask: Tensor,
                coordinates: Tensor | None = None,
                edge_index: Tensor | None = None) -> GeometryMechanismEncoding:
        if features.ndim != 3 or node_mask.shape != features.shape[:2]:
            raise ValueError("Cartesian features/mask have incompatible shapes")
        scalar = self.scalar(features) * node_mask.unsqueeze(-1)
        vector = scalar.new_zeros(
            *scalar.shape[:2], self.vector_channels, 3)
        tensor2 = scalar.new_zeros(
            *scalar.shape[:2], self.tensor_channels, 3, 3)
        state = CartesianTensorState(scalar, vector, tensor2, node_mask)
        available = torch.zeros(
            features.shape[0], dtype=torch.bool, device=features.device)
        if coordinates is not None or edge_index is not None:
            if coordinates is None or edge_index is None:
                raise ValueError("coordinates and Cartesian edges must be provided together")
            available = node_mask.bool().any(-1)
            for layer in self.layers:
                state = layer(state, coordinates, edge_index)
        invariant = torch.cat((
            state.scalar,
            torch.linalg.vector_norm(state.vector, dim=-1),
            torch.linalg.matrix_norm(state.tensor2, dim=(-2, -1))), -1)
        values = self.invariant_projection(invariant)
        safe_mask = node_mask.bool().clone()
        missing = ~safe_mask.any(-1)
        if bool(missing.any()):
            safe_mask[missing, 0] = True
            values = values.clone()
            values[missing, 0] = 0
        seed = self.slot_seed.unsqueeze(0).expand(features.shape[0], -1, -1)
        slots, _ = self.slot_attention(
            seed, values, values, key_padding_mask=~safe_mask,
            need_weights=False)
        slots = slots * available[:, None, None]
        return GeometryMechanismEncoding(state, invariant, slots, available)
