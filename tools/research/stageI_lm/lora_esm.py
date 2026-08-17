"""Stage I: live ESM-2 150M protein encoder with LoRA adapters.

External-representation + training lane: the frozen 150M protein bank is
replaced by a live ESM-2 150M encoder whose attention projections carry
low-rank adapters (r=8, alpha=16). The base LM stays frozen; only the
adapters and the DTA trunk train, in ONE single-stage optimization run.
Chunked encoding (<=1022 residues per chunk, mean/slot pooled) matches the
governed bank policy; gradients flow within chunks. Reported as external
data (model snapshot facebook/esm2_t30_150M_UR50D revision
a695f6045e2e32885fa60af20c13cb35398ce30c, local cache).
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class LoRALinear(nn.Module):
    """Linear layer plus a low-rank adapter delta (alpha/r scaled)."""

    def __init__(self, base: nn.Linear, rank: int = 8, alpha: float = 16.0):
        super().__init__()
        self.base_weight = nn.Parameter(base.weight.detach().clone())
        self.base_bias = (nn.Parameter(base.bias.detach().clone())
                          if base.bias is not None else None)
        self.base_weight.requires_grad_(False)
        if self.base_bias is not None:
            self.base_bias.requires_grad_(False)
        self.rank = int(rank)
        self.scale = float(alpha) / float(rank)
        self.lora_a = nn.Parameter(torch.zeros(base.in_features, rank))
        self.lora_b = nn.Parameter(torch.zeros(rank, base.out_features))
        nn.init.kaiming_uniform_(self.lora_a, a=5 ** 0.5)
        nn.init.zeros_(self.lora_b)

    def forward(self, x):
        delta = (x @ self.lora_a) @ self.lora_b * self.scale
        return F.linear(x, self.base_weight, self.base_bias) + delta


class LiveESMProteinEncoder(nn.Module):
    """Chunked ESM-2 encoding to (pooled, 128-slot residues, mask)."""

    def __init__(self, model_dir: str, device, rank: int = 8, alpha: float = 16.0,
                 slots: int = 128, max_chunk: int = 1022):
        super().__init__()
        from transformers import AutoTokenizer, EsmModel
        self.tokenizer = AutoTokenizer.from_pretrained(model_dir,
                                                       local_files_only=True)
        self.esm = EsmModel.from_pretrained(model_dir, local_files_only=True,
                                            add_pooling_layer=False)
        for parameter in self.esm.parameters():
            parameter.requires_grad_(False)
        for layer in self.esm.encoder.layer:
            self_attn = layer.attention.self
            for name in ("query", "key", "value"):
                setattr(self_attn, name,
                        LoRALinear(getattr(self_attn, name), rank, alpha))
            layer.attention.output.dense = LoRALinear(
                layer.attention.output.dense, rank, alpha)
        self.esm = self.esm.to(device)
        self.device = device
        self.slots = int(slots)
        self.max_chunk = int(max_chunk)
        self.hidden = int(self.esm.config.hidden_size)

    def lora_parameters(self):
        out = []
        for module in self.modules():
            if isinstance(module, LoRALinear):
                out.append(module.lora_a)
                out.append(module.lora_b)
        return out

    def encode(self, sequence, max_chunks: int | None = None):
        """Encode with optional gradient bounding: chunks beyond max_chunks
        are encoded without gradient, keeping GPU memory bounded during
        training while all residues still contribute features."""
        pieces = []
        for index, start in enumerate(range(0, len(sequence), self.max_chunk)):
            chunk = sequence[start:start + self.max_chunk]
            tokens = self.tokenizer(chunk, return_tensors="pt",
                                    add_special_tokens=True)
            tokens = {k: v.to(self.device) for k, v in tokens.items()}
            from contextlib import nullcontext
            context = (nullcontext() if max_chunks is None
                       or index < max_chunks else torch.no_grad())
            with context, torch.autocast(
                    device_type="cuda", dtype=torch.bfloat16,
                    enabled=self.device.startswith("cuda")):
                hidden = self.esm(
                    **tokens).last_hidden_state[0, 1:len(chunk) + 1]
            pieces.append(hidden.float())
        residues = torch.cat(pieces, dim=0)
        length = residues.shape[0]
        slot_index = torch.div(torch.arange(length, device=self.device)
                               * self.slots, length, rounding_mode="floor")
        slots = torch.zeros(self.slots, self.hidden, device=self.device,
                            dtype=torch.float32)
        slots.index_add_(0, slot_index, residues)
        counts = torch.bincount(slot_index, minlength=self.slots)
        mask = counts.gt(0)
        slots[mask] /= counts[mask].unsqueeze(-1)
        pooled = residues.mean(0)
        return pooled, slots, mask

    def lora_state(self):
        return {name: value.detach().cpu() for name, value in
                self.named_parameters() if "lora" in name}

    def load_lora_state(self, state):
        for name, value in state.items():
            for module_name, module in self.named_modules():
                if isinstance(module, LoRALinear):
                    target = module_name + ".lora_a"
                    if name == target:
                        module.lora_a.data.copy_(value)
                    target = module_name + ".lora_b"
                    if name == target:
                        module.lora_b.data.copy_(value)
