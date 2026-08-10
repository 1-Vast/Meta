import sys
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1] / "research" / "crossed_interaction"
sys.path.insert(0, str(MODULE_DIR))

import run_x1b_interaction as audit


def test_cluster_t_uses_components_not_rows():
    result = audit.cluster_t_inference({"a": [2.0] * 20, "b": [1.0] * 10, "c": [3.0] * 10})
    assert result["components"] == 3
    assert result["rectangles"] == 40
    assert result["estimate_tau2"] == 2.0
    assert len(result["leave_one_component_out"]) == 3


def test_exact_sign_enumeration_is_one_sided_and_deterministic():
    values = {"a": [1.0], "b": [1.0]}
    assert audit.exact_rademacher_pvalue(values) == 0.25


def test_endpoint_verdict_distinguishes_margin_from_detection():
    rows = [{"endpoint": "Ki", "dependency_cluster": key, "Z": value}
            for key, value in (("a", 0.01), ("b", 0.011), ("c", 0.009), ("d", 0.01))]
    result = audit.adjudicate_endpoint(rows, {"pass": True}, {"sigma_ucb95": 1.0}, "Ki")
    assert result["gates"]["B3_positive_lcb"]
    assert not result["gates"]["B4_design_margin"]
    assert result["terminal_verdict"] == "X1B_INTERACTION_PRESENT_BELOW_DESIGN_MARGIN"
