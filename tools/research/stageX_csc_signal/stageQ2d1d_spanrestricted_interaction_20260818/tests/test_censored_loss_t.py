"""Equivalence: runner censored_loss_t vs q2.censored_loss (frozen speed fix)."""
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
X0C = STAGE.parent / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(STAGE))

from q2 import censored_loss  # noqa: E402
import runner_d as R  # noqa: E402


def test_censored_loss_t_equivalence():
    rng = np.random.default_rng(7)
    n = 500
    z_obs = rng.normal(0, 2, size=n)
    det = rng.random(n) > 0.3
    z_obs[~det] = np.nan
    blo = rng.normal(-5, 1, size=n)
    bhi = rng.normal(5, 1, size=n)
    blo[det] = 0.0
    bhi[det] = 0.0
    blo[~det][rng.random((~det).sum()) < 0.5] = -np.inf
    bhi[~det][rng.random((~det).sum()) < 0.5] = np.inf
    yhat = torch.from_numpy(rng.normal(0, 3, size=n)).float()
    ref = censored_loss({"yhat": yhat}, z_obs, det, blo, bhi, "cpu")
    zt = torch.from_numpy(z_obs.astype(np.float32))
    dt = torch.from_numpy(det)
    lo = torch.from_numpy(blo.astype(np.float32))
    hi = torch.from_numpy(bhi.astype(np.float32))
    got = R.censored_loss_t({"yhat": yhat}, zt, dt, lo, hi)
    assert float(abs(ref - got)) < 1e-6, f"ref {float(ref):.6f} got {float(got):.6f}"


def test_censored_loss_t_cuda_if_available():
    if not torch.cuda.is_available():
        pytest.skip("no cuda")
    rng = np.random.default_rng(7)
    n = 200
    z_obs = rng.normal(0, 2, size=n)
    det = rng.random(n) > 0.3
    z_obs[~det] = np.nan
    blo = np.zeros(n)
    bhi = np.zeros(n)
    blo[~det] = -np.inf
    bhi[~det] = 5.0
    yhat = torch.from_numpy(rng.normal(0, 3, size=n)).float().cuda()
    ref = censored_loss({"yhat": yhat.cpu()}, z_obs, det, blo, bhi, "cpu")
    got = R.censored_loss_t({"yhat": yhat}, torch.from_numpy(z_obs.astype(np.float32)).cuda(),
                            torch.from_numpy(det).cuda(),
                            torch.from_numpy(blo.astype(np.float32)).cuda(),
                            torch.from_numpy(bhi.astype(np.float32)).cuda())
    assert float(abs(ref - got.cpu())) < 1e-5
