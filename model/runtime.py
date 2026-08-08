"""CUDA runtime contract for the executable meta-operator workflow.

The frozen mathematics is device agnostic.  This project instantiation is not:
training and inference are required to run on CUDA and fail closed otherwise.
Host conversions are reserved for reporting and NumPy-only mathematical audits.
"""

from __future__ import annotations

import os
import platform
import random

# Required by deterministic cuBLAS GEMM on CUDA 10.2 and newer.  Set this before
# importing torch so the first CUDA context observes it.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np
import torch


DEFAULT_DTYPE = torch.float64


def require_cuda(device=None) -> torch.device:
    """Resolve a CUDA device and reject CPU execution explicitly."""
    resolved = torch.device("cuda:0" if device is None else device)
    if resolved.type != "cuda":
        raise RuntimeError(f"CUDA execution is required; received device={resolved}")
    if not torch.cuda.is_available():
        raise RuntimeError(
            "CUDA execution is required, but torch.cuda.is_available() is False"
        )
    if resolved.index is not None and resolved.index >= torch.cuda.device_count():
        raise RuntimeError(
            f"CUDA device {resolved.index} is unavailable; "
            f"device_count={torch.cuda.device_count()}"
        )
    return resolved


def configure_cuda(seed: int, deterministic: bool = True, device=None) -> torch.device:
    """Configure the process before constructing CUDA tensors or modules."""
    resolved = require_cuda(device)
    random.seed(seed)
    np.random.seed(seed % (2**32))
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = not deterministic
    torch.backends.cudnn.deterministic = deterministic
    torch.use_deterministic_algorithms(deterministic)
    return resolved


def module_device(module: torch.nn.Module) -> torch.device:
    """Return the unique device used by a module's parameters and buffers."""
    devices = {x.device for x in (*module.parameters(), *module.buffers())}
    if not devices:
        raise RuntimeError(f"{type(module).__name__} has no parameters or buffers")
    if len(devices) != 1:
        raise RuntimeError(f"module spans multiple devices: {sorted(map(str, devices))}")
    return next(iter(devices))


def assert_cuda_tensor(value: torch.Tensor, name: str) -> None:
    if not isinstance(value, torch.Tensor) or not value.is_cuda:
        device = getattr(value, "device", None)
        raise RuntimeError(f"{name} must be a CUDA tensor; received device={device}")


def assert_finite(value: torch.Tensor, name: str) -> None:
    if not bool(torch.isfinite(value).all().item()):
        raise FloatingPointError(f"{name} contains NaN or infinity")


def to_numpy(value: torch.Tensor) -> np.ndarray:
    """Explicit reporting/audit boundary from CUDA tensors to NumPy."""
    return value.detach().cpu().numpy()


def runtime_report(device=None, dtype=DEFAULT_DTYPE) -> dict:
    resolved = require_cuda(device)
    return {
        "python": platform.python_version(),
        "python_executable": os.path.realpath(os.sys.executable),
        "torch": torch.__version__,
        "cuda_runtime": torch.version.cuda,
        "device": str(resolved),
        "device_name": torch.cuda.get_device_name(resolved),
        "compute_capability": ".".join(map(str, torch.cuda.get_device_capability(resolved))),
        "dtype": str(dtype).removeprefix("torch."),
        "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
        "cudnn_deterministic": torch.backends.cudnn.deterministic,
        "cudnn_benchmark": torch.backends.cudnn.benchmark,
    }
