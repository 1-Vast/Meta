import sys
import json
from pathlib import Path


MODULE_DIR = Path(__file__).resolve().parents[1] / "research" / "crossed_interaction"
sys.path.insert(0, str(MODULE_DIR))

import materialize_x0b_rectangles as manifest


def test_packing_is_cell_disjoint_and_deterministic():
    observed = {
        "p1": {"a", "b", "c", "d"},
        "p2": {"a", "b", "c", "d"},
        "p3": {"a", "b"},
    }
    first = manifest.pack_cell_disjoint(observed)
    assert first == manifest.pack_cell_disjoint(observed)
    used = set()
    for left, right, ligand_a, ligand_b in first:
        cells = {(left, ligand_a), (left, ligand_b),
                 (right, ligand_a), (right, ligand_b)}
        assert len(cells) == 4
        assert not used & cells
        used |= cells


def test_rectangle_id_is_order_sensitive_but_reproducible():
    values = ("p1", "p2", "l1", "l2")
    assert manifest.rectangle_id("Ki", "panel", values) == \
        manifest.rectangle_id("Ki", "panel", values)
    assert manifest.rectangle_id("Ki", "panel", values) != \
        manifest.rectangle_id("Ki", "panel", ("p2", "p1", "l1", "l2"))


def test_materialized_manifest_matches_frozen_x0b_counts_and_hash():
    artifact = manifest.OUTPUT
    recorded = json.loads((artifact / "manifest.json").read_text(encoding="utf-8"))
    assert recorded["affinity_value_fields_selected"] == 0
    assert recorded["counts"]["Ki"]["rectangles"] == 11168
    assert recorded["counts"]["Ki"]["selected_rectangles"] == 827
    assert recorded["counts"]["Kd"]["rectangles"] == 1041
    assert recorded["counts"]["Kd"]["selected_rectangles"] == 605
    assert manifest.sha256_file(artifact / "rectangles.jsonl") == \
        recorded["rectangles_sha256"]
    assert manifest.sha256_file(manifest.CONTRACT) == recorded["contract_sha256"]
    assert manifest.sha256_file(Path(manifest.__file__)) == recorded["runner_sha256"]
