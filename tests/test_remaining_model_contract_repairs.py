from dataclasses import replace

import pytest
import torch

from model import bands
from model.config import DEFAULT, profile
from model.mathematical import FrozenHN
from model.meta_operator import (CSMO, assert_exact_hn_node_budget,
                                 build_band_operator, build_population_bands, context_index,
                                 deployment_manifest,
                                 validate_deployment_manifest)


DTYPE = torch.float64


def _artifact_hashes():
    return {
        "state_schema_hash": "state",
        "view_registry_hash": "views",
        "context_registry_hash": "context",
        "mechanism_schema_hash": "mechanism",
        "protein_bank_hash": "proteins",
        "ligand_bank_hash": "ligands",
        "pair_bank_hash": "pairs",
        "archetype_manifest_hash": "archetypes",
    }


def test_population_fallback_stays_within_reachable_context_partition():
    cfg = profile("local", M=4, dkw_n_min=3, dkw_eps_min=0.01)
    context = torch.tensor([0, 0, 0, 1, 1, 1])
    labels = torch.tensor([0.10, 0.20, 0.30, 0.75, 0.85, 0.95], dtype=DTYPE)
    populations = build_population_bands(cfg, context, labels, dtype=DTYPE)
    broad = torch.as_tensor(
        bands.join(torch.zeros(cfg.n_grid), torch.ones(cfg.n_grid)), dtype=DTYPE
    )
    assert not torch.allclose(populations[0], populations[1])
    assert not torch.allclose(populations[0], broad)


def test_deployment_manifest_rejects_changed_table():
    cfg = profile("local", M=4)
    table = torch.arange(
        cfg.n_context * cfg.band_dim * cfg.n_coef, dtype=DTYPE
    ).reshape(cfg.n_context, cfg.band_dim, cfg.n_coef)
    manifest = deployment_manifest(
        cfg, table, frontend_hash="frontend", source_manifest_hash="source",
        artifact_hashes=_artifact_hashes(),
    )
    validate_deployment_manifest(
        manifest, cfg, table, frontend_hash="frontend", source_manifest_hash="source",
        artifact_hashes=_artifact_hashes(),
    )
    changed = table.clone()
    changed[0, 0, 0] += 1.0
    with pytest.raises(ValueError, match="contract mismatch"):
        validate_deployment_manifest(
            manifest, cfg, changed,
            frontend_hash="frontend", source_manifest_hash="source",
            artifact_hashes=_artifact_hashes(),
        )


def test_deployment_manifest_requires_all_semantic_artifacts():
    cfg = profile("local", M=4)
    table = torch.zeros(cfg.n_context, cfg.band_dim, cfg.n_coef, dtype=DTYPE)
    with pytest.raises(ValueError, match="artifact hashes are incomplete"):
        deployment_manifest(
            cfg, table, frontend_hash="frontend", source_manifest_hash="source"
        )
    with pytest.raises(ValueError, match="frontend and source"):
        deployment_manifest(
            cfg, table, frontend_hash="", source_manifest_hash="source",
            artifact_hashes=_artifact_hashes(),
        )


def test_statistic_domain_fails_closed_instead_of_clamping():
    cfg = profile("local")
    z = torch.zeros(1, cfg.d_z, dtype=DTYPE)
    z[0, 0] = 1.01
    with pytest.raises(ValueError, match=r"outside the declared \[0,1\] domain"):
        context_index(z, cfg)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="model contract requires CUDA")
def test_band_operator_factory_uses_the_local_context_contract():
    cfg = profile("local", M=4, dkw_n_min=2)
    z = torch.zeros(4, cfg.d_z, dtype=DTYPE)
    labels = torch.tensor([0.1, 0.2, 0.8, 0.9], dtype=DTYPE)
    operator = build_band_operator(cfg, z, labels, device="cuda", dtype=DTYPE)
    assert operator.B_table.shape == (cfg.n_context, cfg.band_dim, cfg.n_coef)


@pytest.mark.skipif(not torch.cuda.is_available(), reason="model contract requires CUDA")
def test_frozen_hn_exactness_guards_are_fail_closed():
    small = replace(DEFAULT, d_z=3)
    duplicate = FrozenHN(small, res_N=1, coords=[0, 0, 0], device="cuda")
    assert not duplicate.full_Z
    softmax = FrozenHN(small, res_N=1, coords=[0, 1], node_param="softmax", device="cuda")
    witness = torch.full((softmax.nu_N, small.n_coef), 1.0 / small.n_coef).numpy()
    with pytest.raises(ValueError, match="projected-node"):
        softmax.set_witness(witness)
    with pytest.raises(ValueError, match="projected-node"):
        CSMO.embed_frozen_hn(softmax)
    with pytest.raises(ValueError, match="safety budget"):
        assert_exact_hn_node_budget(DEFAULT, 1, range(DEFAULT.d_z))
