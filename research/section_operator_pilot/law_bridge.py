"""Moment-preserving law-valued bridge for the F6I component statistic."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


DEFAULT_MESH = np.linspace(0.0, 1.0, 7, dtype=np.float64)


@dataclass(frozen=True)
class OperatorLaw:
    mesh: np.ndarray
    z: np.ndarray
    f: np.ndarray
    band: np.ndarray
    beta: np.ndarray

    @property
    def mean(self) -> float:
        return float(self.mesh @ self.beta)

    @property
    def variance(self) -> float:
        return float(np.square(self.mesh - self.mean) @ self.beta)


def _validate_mesh(mesh):
    mesh = np.asarray(mesh, dtype=np.float64)
    if mesh.ndim != 1 or len(mesh) < 2:
        raise ValueError("mesh must be a one-dimensional array with at least two points")
    if not np.all(np.isfinite(mesh)) or not np.all(np.diff(mesh) > 0):
        raise ValueError("mesh must be finite and strictly increasing")
    if not np.allclose(np.diff(mesh), np.diff(mesh)[0], atol=1e-12, rtol=1e-12):
        raise ValueError("the moment-preserving tridiagonal band requires a uniform mesh")
    return mesh


def _mean_at_eta(mesh, target, precision, eta):
    log_weight = eta * mesh - precision * np.square(mesh - target)
    log_weight -= float(np.max(log_weight))
    weight = np.exp(log_weight)
    probability = weight / float(weight.sum())
    return float(mesh @ probability), probability


def f_map(mean, confidence, mesh=DEFAULT_MESH):
    """Return a strictly normalized discrete law with the requested barycentre."""
    mesh = _validate_mesh(mesh)
    mean = float(np.clip(mean, mesh[0], mesh[-1]))
    confidence = float(np.clip(confidence, 0.0, 1.0))
    if np.isclose(mean, mesh[0], atol=1e-14):
        out = np.zeros(len(mesh)); out[0] = 1.0
        return out
    if np.isclose(mean, mesh[-1], atol=1e-14):
        out = np.zeros(len(mesh)); out[-1] = 1.0
        return out
    precision = 2.0 + 48.0 * confidence
    low, high = -1.0, 1.0
    while _mean_at_eta(mesh, mean, precision, low)[0] > mean:
        low *= 2.0
    while _mean_at_eta(mesh, mean, precision, high)[0] < mean:
        high *= 2.0
    probability = None
    for _ in range(100):
        midpoint = 0.5 * (low + high)
        value, probability = _mean_at_eta(mesh, mean, precision, midpoint)
        if value < mean:
            low = midpoint
        else:
            high = midpoint
    probability = np.asarray(probability, dtype=np.float64)
    probability /= probability.sum()
    return probability


def band_map(confidence, mesh=DEFAULT_MESH):
    """Column-stochastic nearest-neighbour diffusion preserving the first moment."""
    mesh = _validate_mesh(mesh)
    confidence = float(np.clip(confidence, 0.0, 1.0))
    diffusion = 0.45 * (1.0 - confidence)
    band = np.eye(len(mesh), dtype=np.float64)
    for column in range(1, len(mesh) - 1):
        band[column, column] = 1.0 - diffusion
        band[column - 1, column] = diffusion / 2.0
        band[column + 1, column] = diffusion / 2.0
    return band


def operator_law(biological, tau, confidence, base=0.5, mesh=DEFAULT_MESH):
    """Evaluate ``z -> F(z) -> B(z)F(z) -> K(beta)`` on a fixed mesh."""
    mesh = _validate_mesh(mesh)
    biological = float(np.clip(biological, -0.5, 0.5))
    tau = float(np.clip(tau, -0.5, 0.5))
    confidence = float(np.clip(confidence, 0.0, 1.0))
    mean = float(np.clip(float(base) + biological + tau, mesh[0], mesh[-1]))
    f = f_map(mean, confidence, mesh)
    band = band_map(confidence, mesh)
    beta = band @ f
    beta = np.clip(beta, 0.0, None)
    beta /= beta.sum()
    z = np.asarray([biological, tau, confidence], dtype=np.float64)
    return OperatorLaw(mesh=mesh.copy(), z=z, f=f, band=band, beta=beta)


def support_confidence(residual, scale=0.10):
    """Bounded, permutation-invariant reliability from robust support dispersion."""
    residual = np.asarray(residual, dtype=np.float64)
    if residual.ndim != 1 or not len(residual) or not np.all(np.isfinite(residual)):
        raise ValueError("residual must be a non-empty finite vector")
    median = float(np.median(residual))
    mad = float(np.median(np.abs(residual - median)))
    effective_noise = 1.4826 * mad / np.sqrt(len(residual))
    return float(np.clip(1.0 / (1.0 + effective_noise / float(scale)), 0.0, 1.0))
