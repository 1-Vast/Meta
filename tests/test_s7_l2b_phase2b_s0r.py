"""Contract tests for the label-blind complete-panel S0R replay."""
from __future__ import annotations

import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
import s0r_run as S  # noqa: E402


def test_preregistration_and_metadata_inputs_are_frozen():
    assert S.sha_file(S.PREREG) == S.PREREG_SHA
    for name, expected in S.INPUT_HASHES.items():
        assert S.sha_file(S.META_DIR / name) == expected


def test_metadata_only_records_have_no_structural_label_fields():
    forbidden = {"edges", "positive_binary_edges", "mask", "gain_set",
                 "loss_set", "symmetric_difference"}
    rows = S.load_jsonl(S.META_DIR / "metadata_only_records.jsonl")
    assert len(rows) == 13477
    assert all(not (forbidden & set(row)) for row in rows)


def test_runner_does_not_import_real_phase2b_pair_builder():
    source = inspect.getsource(S)
    assert "import p2b_run" not in source
    assert "P2B.prepare" not in source
    assert "build_pairs(" not in source


def test_frozen_pair_files_have_registered_counts_and_unique_rows():
    train = S.load_jsonl(S.META_DIR / "train_pairs.jsonl")
    held = S.load_jsonl(S.META_DIR / "heldoutA_pairs.jsonl")
    assert len(train) == S.EXPECTED["train_universe"]
    assert len(held) == S.EXPECTED["heldout"]
    assert len({(x["seq_key"], x["a"], x["b"]) for x in train}) == len(train)
    assert len({(x["seq_key"], x["a"], x["b"]) for x in held}) == len(held)


def test_stream_serialization_is_semantically_exact(tmp_path):
    stream = [
        (0, 0, [("c1", [("p1", [("p1", "a", "b"),
                                  ("p1", "c", "d")])])]),
        (0, 1, [("c2", [("p2", [("p2", "e", "f")])])]),
    ]
    path = tmp_path / "stream.jsonl"
    S.serialize_stream(path, stream)
    loaded = S.deserialize_stream(path)
    assert loaded == stream
    assert S.stream_hash(loaded) == S.stream_hash(stream)


def test_balanced_ray_scaling_preserves_product_without_gauge_distortion():
    rng = np.random.default_rng(7)
    U = rng.normal(size=(8, 13))
    V = rng.normal(size=(8, 5))
    scale = 20.1977
    root = math.sqrt(scale)
    Us, Vs = root * U, root * V
    assert np.allclose(Us.T @ Vs, scale * (U.T @ V), rtol=1e-13, atol=1e-13)
    assert np.isclose(np.linalg.norm(Us) / np.linalg.norm(Vs),
                      np.linalg.norm(U) / np.linalg.norm(V), rtol=1e-14)


def test_registered_trajectory_decision_uses_component_ucb(monkeypatch):
    comp0 = {f"c{i}": 0.9 for i in range(81)}
    comp100 = {f"c{i}": 0.7 for i in range(81)}
    trajectory = {"checkpoints": [
        {"update": 0, "train": {"component_macro_bce": 0.6},
         "heldout": {"component_macro_ap_bidir": 0.9, "component_ap": comp0}},
        {"update": 100, "train": {"component_macro_bce": 0.4},
         "heldout": {"component_macro_ap_bidir": 0.7, "component_ap": comp100}},
    ]}
    out = S.trajectory_decision(trajectory)
    assert out["MISALIGNED"] is True
    assert out["component_inference"]["units"] == 81
    assert out["component_inference"]["ucb95_one_sided"] < 0


def test_machine_outputs_are_written_by_runner_not_posthoc():
    source = inspect.getsource(S.main)
    assert "PHASE2B_S0R_VERDICT.json" in source
    assert "PHASE2B_POST_S0R_NOT_RUN.json" in source
    assert "SYNTHETIC_CONTRACT_INVALID" in source


def test_census_declares_zero_cross_split_component_overlap():
    census = json.loads(
        (S.OUT / "S0R_METADATA_ONLY_CENSUS.json").read_text(encoding="utf-8")
    )
    assert census["components"]["overlap"] == 0
    assert census["pairs"]["heldoutA"] == S.EXPECTED["heldout"]
    assert census["pairs"]["hash_stratified_train_panel"] == S.EXPECTED["train_panel"]

