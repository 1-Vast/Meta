import json
from pathlib import Path

import pytest
import torch

from model.config import DEFAULT
from model.meta_operator import (FixedGate, Gate, query_only_model,
                                 scale_only_model,
                                 support_mean_calibration_model,
                                 uniform_mixture_model)


DTYPE = torch.float64
ROOT = Path(__file__).resolve().parents[1]


def test_unknown_learned_gate_parameterization_is_rejected():
    if not torch.cuda.is_available():
        pytest.skip("model construction requires CUDA")
    with pytest.raises(ValueError, match="unknown gate parameterization"):
        Gate(DEFAULT.d_z, DEFAULT.n_views, 8, 2, param="uniform", device="cuda")


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CSMO requires CUDA")
def test_uniform_routing_is_exact_and_input_independent():
    model = uniform_mixture_model(DEFAULT, device="cuda")
    z = torch.rand(7, DEFAULT.d_z, dtype=DTYPE, device="cuda")
    weights = model.gate_weights(z)
    assert isinstance(model.gate, FixedGate)
    assert torch.equal(weights, torch.full_like(weights, 1.0 / DEFAULT.n_views))
    assert not tuple(model.gate.parameters())


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CSMO requires CUDA")
def test_operator_baselines_follow_declared_statistic_coordinates():
    assert support_mean_calibration_model(DEFAULT, device="cuda").views == [(12,)]
    assert scale_only_model(DEFAULT, device="cuda").views == [(13,)]
    assert query_only_model(DEFAULT, device="cuda").views == [(0, 1), (2, 3), (26, 27)]


def test_serialized_deployment_declares_six_contexts():
    deployment = json.loads((ROOT / "config" / "default.json").read_text())
    context = deployment["statistic_interface"]["context_map"]
    assert context["bins"] == [3, 2, 1]
    assert context["n_context"] == DEFAULT.n_context
