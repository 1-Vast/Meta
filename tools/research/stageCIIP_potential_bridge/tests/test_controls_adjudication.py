import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
sys.path.insert(0, str(STAGE))

from adjudicate_controls import load_and_adjudicate  # noqa: E402


def test_frozen_control_result_adjudicates_not_supported():
    out = load_and_adjudicate(STAGE / "CONTROL_RESULT.json")
    assert out["verdict"] == "ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED"
    assert not out["checks"]["correct_nonconstant_coverage_exceeds_random"]
    assert not out["checks"]["correct_beats_random_window_r2_and_sign"]
    assert out["authorization"]["ciip_1b"] == "NOT_AUTHORIZED"
