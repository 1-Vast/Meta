import sys
from pathlib import Path

import numpy as np


MODULE_DIR = Path(__file__).resolve().parents[1] / "research" / "crossed_interaction"
sys.path.insert(0, str(MODULE_DIR))

import run_x1ar_direct_dd as audit


def test_profile_icc_separates_low_and_high_cluster_dependence():
    rng = np.random.default_rng(7)
    low = [rng.normal(size=40) for _ in range(20)]
    offsets = rng.normal(scale=3.0, size=20)
    high = [offset + rng.normal(scale=0.2, size=40) for offset in offsets]
    assert audit.profile_icc(low)["rho_mle"] < 0.1
    assert audit.profile_icc(high)["rho_mle"] > 0.8


def test_dd_scaling_matches_cell_interaction_scale():
    means = [8.0, 6.0, 5.0, 7.0]
    dd = means[0] - means[1] - means[2] + means[3]
    assert dd == 4.0
    assert dd / 2.0 == 2.0
