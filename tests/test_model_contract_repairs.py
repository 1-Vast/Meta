import numpy as np
import pytest
import torch

from model import bands
from model.mathematical import (midband_cdf, readout_midband_mean,
                                readout_midband_mean_torch,
                                readout_midband_median)


DTYPE = torch.float64


def _valid_beta():
    lower = np.array([0.10, 0.40, 0.80, 0.95])
    upper = np.array([0.30, 0.60, 0.90, 1.00])
    beta = bands.join(lower, upper)
    bands.assert_valid(beta)
    return beta, np.array([0.0, 0.25, 0.75, 1.0])


def test_midpoint_readouts_stay_inside_the_band_without_normalization():
    beta, grid = _valid_beta()
    lo, up = bands.split(beta)
    mid = midband_cdf(beta)
    assert np.all(mid >= lo)
    assert np.all(mid <= up)
    assert np.all(np.diff(mid) >= 0)
    assert mid[-1] == 1.0
    assert np.isfinite(readout_midband_mean(beta, grid))
    assert np.isfinite(readout_midband_median(beta, grid))
    beta_t = torch.as_tensor(beta, dtype=DTYPE).unsqueeze(0)
    expected = torch.tensor([readout_midband_mean(beta, grid)], dtype=DTYPE)
    assert torch.allclose(
        readout_midband_mean_torch(beta_t, torch.as_tensor(grid, dtype=DTYPE)),
        expected,
    )


def test_midpoint_readouts_fail_closed_for_an_invalid_band():
    beta, grid = _valid_beta()
    beta[-1] = 0.9
    with pytest.raises(AssertionError, match="midband readout input"):
        readout_midband_mean(beta, grid)
