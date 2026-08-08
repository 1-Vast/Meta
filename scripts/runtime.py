"""Performance helpers for model workflows and CUDA batch handling."""
from __future__ import annotations

from collections.abc import Mapping

import torch


def seeded_generator(seed: int) -> torch.Generator:
    return torch.Generator().manual_seed(seed)


def optimized_loader_options(workers: int) -> dict:
    if workers < 1:
        raise ValueError("optimized loader requires at least one worker")
    return {"num_workers": workers, "pin_memory": True,
            "persistent_workers": True, "prefetch_factor": 2}


def move_to_device(value, device: torch.device, *, non_blocking: bool):
    """Move nested tensor batches or PyG Batch objects without changing values."""
    if isinstance(value, torch.Tensor):
        return value.to(device, non_blocking=non_blocking)
    if isinstance(value, Mapping):
        return type(value)((key, move_to_device(item, device, non_blocking=non_blocking))
                           for key, item in value.items())
    if isinstance(value, tuple):
        return tuple(move_to_device(item, device, non_blocking=non_blocking) for item in value)
    if isinstance(value, list):
        return [move_to_device(item, device, non_blocking=non_blocking) for item in value]
    if hasattr(value, "to"):
        return value.to(device, non_blocking=non_blocking)
    return value


def scalar_history_to_list(values: list[torch.Tensor]) -> list[float]:
    if not values:
        return []
    if any(value.numel() != 1 for value in values):
        raise ValueError("training history values must be scalar tensors")
    return torch.stack([value.detach().reshape(()) for value in values]).cpu().tolist()


def assert_nested_equal(left, right, path: str = "batch") -> None:
    if isinstance(left, torch.Tensor) and isinstance(right, torch.Tensor):
        if not torch.equal(left.cpu(), right.cpu()):
            raise AssertionError(f"{path} tensor values differ")
        return
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        if list(left) != list(right):
            raise AssertionError(f"{path} mapping keys differ")
        for key in left:
            assert_nested_equal(left[key], right[key], f"{path}.{key}")
        return
    if isinstance(left, (tuple, list)) and isinstance(right, type(left)):
        if len(left) != len(right):
            raise AssertionError(f"{path} sequence lengths differ")
        for index, (left_item, right_item) in enumerate(zip(left, right)):
            assert_nested_equal(left_item, right_item, f"{path}[{index}]")
        return
    if left != right:
        raise AssertionError(f"{path} values differ")


def assert_nested_cuda(value, path: str = "batch") -> None:
    if isinstance(value, torch.Tensor):
        if not value.is_cuda:
            raise RuntimeError(f"{path} must be CUDA resident; got {value.device}")
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            assert_nested_cuda(item, f"{path}.{key}")
        return
    if isinstance(value, (tuple, list)):
        for index, item in enumerate(value):
            assert_nested_cuda(item, f"{path}[{index}]")


def assert_optimizer_cuda(optimizer: torch.optim.Optimizer) -> None:
    for parameter, state in optimizer.state.items():
        if not parameter.is_cuda:
            raise RuntimeError("optimizer owns a non-CUDA parameter")
        assert_nested_cuda(state, "optimizer.state")
