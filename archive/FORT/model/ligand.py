"""Ligand fingerprint encoder.

The input may include Morgan bits and appended physicochemical descriptors.
"""
from torch import nn
import torch


class FingerprintEncoder(nn.Module):
    def __init__(self, bits: int = 1034, rep: int = 128):
        super().__init__()
        self.inputdim = bits
        self.net = nn.Sequential(
            nn.Linear(bits, 256), nn.ReLU(),
            nn.Linear(256, rep),
        )

    def forward(self, fp: torch.Tensor) -> torch.Tensor:
        if fp.shape[-1] != self.inputdim:
            raise ValueError(
                f"expected {self.inputdim} ligand features, received {fp.shape[-1]}"
            )
        return self.net(fp)


__all__ = ["FingerprintEncoder"]
