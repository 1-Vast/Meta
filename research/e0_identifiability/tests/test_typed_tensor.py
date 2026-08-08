import numpy as np
import torch

from research.e0_identifiability.typed_energy_tensor_contract import (
    CP_TENSOR_PARAMETERS, FULL_TENSOR_PARAMETERS,
)
from research.e0_identifiability.run_typed_tensor import (
    CPTypedTensor, FullTypedTensor, _analytic_cp_tensor, _typed_statistics,
)


def test_registered_parameter_counts_and_cp_exact_realization():
    full = FullTypedTensor()
    cp = CPTypedTensor()
    assert sum(value.numel() for value in full.parameters()) == FULL_TENSOR_PARAMETERS
    assert sum(value.numel() for value in cp.parameters()) == CP_TENSOR_PARAMETERS
    rng = np.random.default_rng(3)
    weights = rng.normal(size=(8, 6)).astype(np.float32)
    distance = np.asarray([1.0, 0.7, 0.2, -0.2, -0.6], dtype=np.float32)
    tensor, error = _analytic_cp_tensor(weights, distance)
    assert tensor.shape == (8, 6, 5)
    assert error < 1e-6


def test_typed_statistics_are_atom_permutation_invariant():
    rows = [{"ligand_state_key": "l", "active_protein_key": "p", "example_id": 0}]
    atom = np.zeros((3, 40), dtype=np.float32)
    atom[np.arange(3), 32 + np.arange(3)] = 1
    residue = np.eye(6, dtype=np.float32)[:2]
    contact = np.arange(6, dtype=np.float32).reshape(3, 2) + 1
    distance = np.zeros((3, 2, 5), dtype=np.float32)
    distance[..., 2] = 1
    proteins = {"p": {"chemistry": residue}}
    ligands = {"l": {"chemistry": atom}}
    original = _typed_statistics(
        rows, proteins, ligands, {0: {"contact": contact, "distance": distance}})
    reversed_value = _typed_statistics(
        rows, proteins, {"l": {"chemistry": atom[::-1].copy()}},
        {0: {"contact": contact[::-1].copy(), "distance": distance[::-1].copy()}})
    assert np.allclose(original, reversed_value, atol=1e-7)


def test_full_tensor_computes_registered_inner_product():
    model = FullTypedTensor()
    with torch.no_grad():
        model.energy.copy_(torch.arange(240).reshape(8, 6, 5))
    features = torch.ones(2, 8, 6, 5)
    assert torch.allclose(model(features), torch.full((2,), float(sum(range(240)))))
