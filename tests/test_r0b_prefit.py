import numpy as np
import pytest

from research.correspondence_router.run_r0b_prefit import bootstrap_mde80


def test_bootstrap_mde80_is_shift_invariant_and_positive():
    values = np.asarray([0.1, 0.2, 0.4, 0.8], dtype=np.float64)
    first = bootstrap_mde80(values, seed=1700, draws=1000)
    shifted = bootstrap_mde80(values + 17.0, seed=1700, draws=1000)
    assert first > 0
    assert shifted == pytest.approx(first, abs=1e-15)


def test_bootstrap_mde80_rejects_invalid_input():
    with pytest.raises(ValueError, match="at least two"):
        bootstrap_mde80(np.asarray([0.1]), seed=1700)
    with pytest.raises(ValueError, match="positive"):
        bootstrap_mde80(np.asarray([0.1, 0.2]), seed=1700, draws=0)
