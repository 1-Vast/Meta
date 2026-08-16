"""CUDA contracts for the isolated target-function memory gate."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import torch

from research.shared.memory import (
    ExpertBank,
    FeatureMap,
    FunctionMemory,
    MemoryConfig,
    fitexperts,
    fitfeatures,
    strictsplit,
    wrongtargets,
)


pytestmark = pytest.mark.skipif(
    not torch.cuda.is_available(), reason="function-memory tests require CUDA"
)


def bank() -> ExpertBank:
    coefficient = torch.tensor([[1.0, 0.0], [-1.0, 0.0]], device="cuda")
    protein = torch.tensor([[1.0, 0.0], [-1.0, 0.0]], device="cuda")
    return ExpertBank(
        targets=("A", "B"),
        coefficient=coefficient,
        protein=protein,
        proteincenter=torch.zeros(2, device="cuda"),
        mapalpha=torch.zeros_like(coefficient),
        fitnoise=0.5,
    )


def tensors() -> tuple[torch.Tensor, ...]:
    support = torch.zeros(5, 2, device="cuda")
    label = torch.zeros(5, device="cuda")
    base = torch.zeros(5, device="cuda")
    query = torch.tensor([[2.0, 0.0]], device="cuda")
    querybase = torch.zeros(1, device="cuda")
    return support, label, base, query, querybase


def testproteinpriorchangesfunction() -> None:
    config = MemoryConfig(
        noise=1.0,
        proteinweight=20.0,
        functionweight=0.0,
        topk=1,
        dynamicbias=-30.0,
        nullbias=-30.0,
    )
    memory = FunctionMemory(bank(), config)
    support, label, base, query, querybase = tensors()
    correct = memory.predict(
        support, label, base, query, querybase, torch.tensor([1.0, 0.0], device="cuda")
    )
    wrong = memory.predict(
        support, label, base, query, querybase, torch.tensor([-1.0, 0.0], device="cuda")
    )
    free = memory.predict(support, label, base, query, querybase, None)
    assert correct["prediction"].item() > 1.9
    assert wrong["prediction"].item() < -1.9
    assert abs(free["prediction"].item()) < 1e-5
    assert correct["prediction"].is_cuda


def testnullexpertisexplicit() -> None:
    memory = FunctionMemory(
        bank(),
        MemoryConfig(
            noise=1.0,
            proteinweight=0.0,
            topk=2,
            dynamicbias=-30.0,
            nullbias=5.0,
        ),
    )
    output = memory.predict(*tensors(), torch.tensor([1.0, 0.0], device="cuda"))
    assert output["nullprobability"].item() > 0.95
    assert output["weight"].shape == (4,)


def testexactfivesupport() -> None:
    memory = FunctionMemory(bank(), MemoryConfig(topk=2))
    support, label, base, query, querybase = tensors()
    with pytest.raises(ValueError, match="exactly five"):
        memory.predict(
            support[:4], label[:4], base[:4], query, querybase, None
        )


def testoffsetmarginalisfinite() -> None:
    memory = FunctionMemory(bank(), MemoryConfig(noise=0.7, topk=2))
    score = torch.randn(5, 3, device="cuda")
    residual = torch.randn(5, device="cuda")
    likelihood, offset = memory._evidence(score, residual)
    assert likelihood.shape == (3,)
    assert offset.shape == (3,)
    assert likelihood.is_cuda and torch.isfinite(likelihood).all()
    assert offset.is_cuda and torch.isfinite(offset).all()


def teststrictcomponentclosure() -> None:
    frame = pd.DataFrame(
        {
            "target": [f"T{index}" for index in range(12)],
            "component": [f"P{index // 2}" for index in range(12)],
        }
    )
    fit, gate = strictsplit(frame, 0.25, 17)
    assert set(fit.component).isdisjoint(gate.component)
    assert len(fit) + len(gate) == len(frame)


def testwrongproteiniscrosscomponent() -> None:
    frame = pd.DataFrame(
        {
            "target": ["A", "B", "C", "D"],
            "component": ["P1", "P1", "P2", "P3"],
        }
    )
    mapping = wrongtargets(frame, 29)
    components = frame.set_index("target").component.to_dict()
    assert mapping == wrongtargets(frame, 29)
    assert all(components[target] != components[other] for target, other in mapping.items())


def testfeatureselectionuseswithintargetsignal() -> None:
    feature = np.zeros((12, 1034), dtype=np.float32)
    feature[:, 0] = np.tile([0.0, 1.0], 6)
    feature[:, 1] = np.repeat([0.0, 1.0], 6)
    frame = pd.DataFrame(
        {
            "target": ["A"] * 6 + ["B"] * 6,
            "source": np.arange(12),
            "affinity": feature[:, 0] * 2.0 + np.repeat([3.0, 8.0], 6),
        }
    )
    mapping = fitfeatures(frame, feature, np.zeros(12, dtype=np.float32), latent=1)
    assert int(mapping.columns[0]) == 0
    assert mapping.center.is_cuda and mapping.scale.is_cuda
    assert len(mapping.columns) == 11


def testexpertfitstaysongpu() -> None:
    feature = np.zeros((16, 1034), dtype=np.float32)
    feature[:, 0] = np.tile([0.0, 1.0], 8)
    frame = pd.DataFrame(
        {
            "target": ["A"] * 8 + ["B"] * 8,
            "source": np.arange(16),
            "affinity": np.r_[feature[:8, 0] * 2.0, feature[8:, 0] * -2.0],
        }
    )
    mapping = FeatureMap(
        columns=torch.tensor([0], device="cuda"),
        center=torch.tensor([0.5], device="cuda"),
        scale=torch.tensor([0.5], device="cuda"),
    )
    pooled = {
        "A": torch.tensor([1.0, 0.0], device="cuda"),
        "B": torch.tensor([0.0, 1.0], device="cuda"),
    }
    fitted = fitexperts(
        frame,
        feature,
        np.zeros(16, dtype=np.float32),
        mapping,
        pooled,
        ridge=0.1,
        mapridge=1.0,
    )
    assert fitted.coefficient.is_cuda and fitted.mapalpha.is_cuda
    assert fitted.coefficient[0, 0] > 0
    assert fitted.coefficient[1, 0] < 0
    assert fitted.targets == ("A", "B")
