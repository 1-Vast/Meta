from dataclasses import replace

import pytest
import torch

from model import bands
from model.config import DEFAULT, profile
from model.mathematical import FrozenHN
from model.meta_operator import (CSMO, assert_exact_hn_node_budget,
                                 build_population_bands, deployment_manifest,
                                 validate_deployment_manifest)


DTYPE = torch.float64


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
        cfg, table, frontend_hash="frontend", source_manifest_hash="source"
    )
    validate_deployment_manifest(
        manifest, cfg, table, frontend_hash="frontend", source_manifest_hash="source"
    )
    changed = table.clone()
    changed[0, 0, 0] += 1.0
    with pytest.raises(ValueError, match="contract mismatch"):
        validate_deployment_manifest(
            manifest, cfg, changed,
            frontend_hash="frontend", source_manifest_hash="source",
        )


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
