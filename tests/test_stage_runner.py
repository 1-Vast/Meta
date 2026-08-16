"""Tests for the robust stage runner (engineering contract 2026-08-16)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.stage_smoke import check_artifacts


def _fake_output(root: Path) -> Path:
    output = root / "run"
    output.mkdir(parents=True)
    (output / "checkpoint.pt").write_bytes(b"ckpt")
    (output / "progress.jsonl").write_text("{}\n", encoding="utf-8")
    (output / "PREDICTIONS_meta_val.jsonl").write_text(
        json.dumps({"k": 0, "target": "t", "component": "c"}) + "\n",
        encoding="utf-8")
    payload = {
        "config": {"seed": 1},
        "split_assignment_sha256": "abc",
        "checkpoint_sha256": "def",
        "donors": {"evaluation_wrong_protein_pool": "meta_val",
                   "whitening_pool": "meta_train"},
        "meta_test": {"included": False, "evaluated": False},
    }
    (output / "RESULT.json").write_text(
        json.dumps(payload), encoding="utf-8")
    return output


def test_check_artifacts_accepts_a_complete_run(tmp_path):
    output = _fake_output(tmp_path)
    result = check_artifacts(output)
    assert result["ok"], result["problems"]


def test_check_artifacts_rejects_missing_files(tmp_path):
    output = _fake_output(tmp_path)
    (output / "checkpoint.pt").unlink()
    result = check_artifacts(output)
    assert not result["ok"]
    assert any("checkpoint.pt" in problem for problem in result["problems"])


def test_check_artifacts_rejects_an_unsealed_meta_test(tmp_path):
    output = _fake_output(tmp_path)
    payload = json.loads((output / "RESULT.json").read_text(encoding="utf-8"))
    payload["meta_test"] = {"included": True, "evaluated": True}
    (output / "RESULT.json").write_text(
        json.dumps(payload), encoding="utf-8")
    result = check_artifacts(output)
    assert not result["ok"]
    assert any("meta_test" in problem for problem in result["problems"])


def test_check_artifacts_rejects_malformed_predictions(tmp_path):
    output = _fake_output(tmp_path)
    (output / "PREDICTIONS_meta_val.jsonl").write_text("not json\n",
                                                       encoding="utf-8")
    result = check_artifacts(output)
    assert not result["ok"]


def test_run_stage_treats_partial_dirs_as_incomplete(tmp_path):
    """A leftover partial output directory must not count as completed."""
    import runpy
    import scripts.run_stage as run_stage

    partial = tmp_path / "partial"
    partial.mkdir()
    (partial / "progress.jsonl").write_text("{}\n", encoding="utf-8")
    assert not run_stage._completed(partial)
    complete = _fake_output(tmp_path)
    assert run_stage._completed(complete)
