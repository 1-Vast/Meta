import sys
import json
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1] / "research" / "crossed_interaction"
sys.path.insert(0, str(MODULE_DIR))

import cycle_quotient_census as audit


def test_tree_has_zero_cycle_dimension():
    result = audit.graph_stats({("p1", "l1"), ("p1", "l2"), ("p2", "l2")})
    assert result["cycle_dimension"] == 0


def test_rectangle_is_one_dimensional_cycle_space():
    edges = {("p1", "l1"), ("p1", "l2"), ("p2", "l1"), ("p2", "l2")}
    assert audit.graph_stats(edges)["cycle_dimension"] == 1


def test_six_cycle_is_one_dimensional_without_rectangle():
    edges = {("p1", "l1"), ("p2", "l1"), ("p2", "l2"),
             ("p3", "l2"), ("p3", "l3"), ("p1", "l3")}
    assert audit.graph_stats(edges)["cycle_dimension"] == 1


def test_frozen_census_records_algebra_without_claiming_independence():
    path = (Path(__file__).resolve().parents[1] / "report" / "crossed_interaction" /
            "cycle_quotient_feasibility" / "census.json")
    report = json.loads(path.read_text(encoding="utf-8"))
    assert report["affinity_value_reads"] == 0
    assert not report["training_authorized"]
    assert report["endpoints"]["Ki"]["raw_panel_cycle_dimension"] == 29677
    assert report["endpoints"]["Kd"]["raw_panel_cycle_dimension"] == 3279
    assert report["endpoints"]["Ki"]["exact_assay_cycle_dimension"] == 0
    assert report["endpoints"]["Kd"]["exact_assay_cycle_dimension"] == 0
    assert report["endpoints"]["Ki"]["largest_dependency_component_share"] > 0.25
    assert report["endpoints"]["Kd"]["largest_dependency_component_share"] > 0.25
