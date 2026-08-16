"""CUDA contracts for the isolated sequence-conditioned pair-prior gate."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import torch

from research.shared.pairprior import PairPrior, evidencecontrast, gaussiannll
from research.shared.priorgate import balancedorder, splitcomponents, wrongtargets
from research.shared.priorgate import applybase, ridgebase
from research.shared.rankgate import anchored


pytestmark = pytest.mark.skipif(
    not torch.cuda.is_available(), reason="pair-prior tests require CUDA"
)


def testkzeroprior() -> None:
    torch.manual_seed(3)
    prior = PairPrior(proteindim=16, liganddim=12, width=32, heads=4).cuda().eval()
    protein = torch.randn(5, 4, 16, device="cuda")
    ligand = torch.randn(5, 12, device="cuda")
    direct = prior.predict(protein, ligand)
    explicit = prior.predict(
        protein, ligand, supportlabel=torch.empty(0, device="cuda")
    )

    assert direct["prediction"].is_cuda
    assert torch.equal(direct["prediction"], explicit["prediction"])
    assert torch.equal(direct["logvariance"], explicit["logvariance"])
    with pytest.raises(ValueError, match="k=0"):
        prior.predict(protein, ligand, supportlabel=torch.ones(1, device="cuda"))


def testproteinsensitivity() -> None:
    torch.manual_seed(5)
    prior = PairPrior(proteindim=16, liganddim=12, width=32, heads=4).cuda().eval()
    protein = torch.randn(6, 4, 16, device="cuda", requires_grad=True)
    ligand = torch.randn(6, 12, device="cuda")
    prediction = prior.predict(protein, ligand)["prediction"]
    permuted = prior.predict(protein.roll(1, dims=0), ligand)["prediction"]
    prediction.sum().backward()

    assert not torch.allclose(prediction, permuted)
    assert protein.grad is not None
    assert torch.count_nonzero(protein.grad) > 0


def testcontrastgradient() -> None:
    torch.manual_seed(7)
    prior = PairPrior(proteindim=16, liganddim=12, width=32, heads=4).cuda()
    protein = torch.randn(8, 4, 16, device="cuda")
    ligand = torch.randn(8, 12, device="cuda")
    label = torch.randn(8, device="cuda")
    correct = prior.predict(protein, ligand)
    wrong = prior.predict(protein.roll(1, dims=0), ligand)
    loss = gaussiannll(correct, label).mean() + evidencecontrast(
        correct, wrong, label
    )
    loss.backward()

    assert loss.is_cuda and torch.isfinite(loss)
    assert all(
        parameter.grad is None or torch.isfinite(parameter.grad).all()
        for parameter in prior.parameters()
    )


def testcomponentclosure() -> None:
    frame = pd.DataFrame(
        {
            "target": [f"T{index}" for index in range(10) for _ in range(2)],
            "hcluster": [f"H{index // 2}" for index in range(10) for _ in range(2)],
        }
    )
    fit, gate = splitcomponents(frame, holdout=0.2, seed=11)

    assert not set(fit.hcluster).intersection(gate.hcluster)
    assert len(fit) + len(gate) == len(frame)


def testbalancedcontrol() -> None:
    frame = pd.DataFrame(
        {
            "target": ["A"] * 5 + ["B"] * 3 + ["C"] * 2,
            "hcluster": ["H1"] * 5 + ["H2"] * 3 + ["H3"] * 2,
        }
    )
    order = balancedorder(frame, seed=13)
    first = frame.loc[order[:3], "target"].tolist()
    mapping = wrongtargets(frame, seed=17)

    assert np.array_equal(np.sort(order), np.arange(len(frame)))
    assert len(set(first)) == 3
    assert mapping == wrongtargets(frame, seed=17)
    component = frame.groupby("target").hcluster.first().to_dict()
    assert all(component[target] != component[wrong] for target, wrong in mapping.items())


def testridgebasegpu() -> None:
    feature = np.arange(30, dtype=np.float32).reshape(10, 3) / 10.0
    label = 2.0 + feature @ np.array([0.5, -0.2, 0.3], dtype=np.float32)
    fit = pd.DataFrame(
        {"source": np.arange(8), "affinity": label[:8]}
    )
    gate = pd.DataFrame(
        {"source": np.arange(8, 10), "affinity": label[8:]}
    )

    base, info = ridgebase(fit, gate, feature, ridge=1e-4, batchsize=2)

    assert np.sqrt(np.mean((base[:8] - label[:8]) ** 2)) < 1e-3
    assert np.sqrt(np.mean((base[8:] - label[8:]) ** 2)) < 1e-3
    assert info["kind"] == "ligand ridge"


def testapplybasegpu() -> None:
    correction = torch.tensor([0.2, -0.3], device="cuda")
    baseline = torch.tensor([1.0, 2.0], device="cuda")
    output = {
        "prediction": correction,
        "logvariance": torch.zeros(2, device="cuda"),
    }

    result = applybase(output, baseline)

    assert torch.equal(result["prediction"], baseline + correction)
    assert torch.equal(result["logvariance"], output["logvariance"])


def testanchoredshiftgpu() -> None:
    torch.manual_seed(19)
    model = PairPrior(proteindim=16, liganddim=12, width=32, heads=4).cuda()
    protein = torch.randn(4, 16, device="cuda")
    support = torch.randn(5, 12, device="cuda")
    query = torch.randn(7, 12, device="cuda")
    supportlabel = torch.randn(5, device="cuda")
    supportbase = torch.randn(5, device="cuda")
    querybase = torch.randn(7, device="cuda")
    first = anchored(
        model,
        protein,
        support,
        query,
        supportlabel,
        supportbase,
        querybase,
    )
    with torch.no_grad():
        model.head[-1].bias[0].add_(7.0)
    shifted = anchored(
        model,
        protein,
        support,
        query,
        supportlabel,
        supportbase,
        querybase,
    )

    assert torch.allclose(first, shifted, atol=2e-6, rtol=0.0)


def testanchorednullgpu() -> None:
    model = PairPrior(proteindim=16, liganddim=12, width=32, heads=4).cuda()
    with torch.no_grad():
        model.head[-1].weight[0].zero_()
        model.head[-1].bias[0].fill_(3.0)
    protein = torch.randn(4, 16, device="cuda")
    support = torch.randn(5, 12, device="cuda")
    query = torch.randn(7, 12, device="cuda")
    supportlabel = torch.randn(5, device="cuda")
    supportbase = torch.randn(5, device="cuda")
    querybase = torch.randn(7, device="cuda")

    prediction = anchored(
        model,
        protein,
        support,
        query,
        supportlabel,
        supportbase,
        querybase,
    )

    expected = querybase + (supportlabel - supportbase).mean()
    assert torch.equal(prediction, expected)
