"""Contract tests for the 2026-08-16 Stage 0 fixes.

Every test here pins one rule of the frozen evaluation contract:

1. the sealed meta_test split is dropped physically unless authorized;
2. wrong-protein donors at evaluation come from the evaluation split itself;
3. the whitening mean/covariance is fitted on meta_train only, always;
4. training-time counterfactual donors come from meta_train (both proteins seen);
5. no default training path evaluates meta_test.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_level_shape import matched_donors
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
PYTHON = Path(sys.executable)


@pytest.fixture(scope="module")
def data():
    """One shared sealed double-cold dataset for every test in this module."""
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT, include_meta_test=False)


def test_meta_test_is_physically_sealed_by_default():
    sealed = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                       split_directory=SPLIT, include_meta_test=False)
    assert "meta_test" not in sealed.tasks
    assert all(cell["split"] != "meta_test" for cell in sealed.cells)
    assert "meta_val" in sealed.tasks and "meta_train" in sealed.tasks


def test_meta_test_is_present_only_when_authorized():
    unsealed = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                         split_directory=SPLIT, include_meta_test=True)
    assert "meta_test" in unsealed.tasks
    assert len(unsealed.tasks["meta_test"]) == 22


def test_evaluation_donors_come_from_the_evaluation_split(data):
    donors = matched_donors(data, "meta_val", donor_pool="meta_val",
                            whitening_pool="meta_train")
    component = {cell["target_id"]: cell["protein_group_40"]
                 for cell in data.cells}
    meta_val = set(data.tasks["meta_val"])
    assert set(donors) == meta_val
    for target, donor in donors.items():
        assert donor in meta_val
        assert donor != target
        assert component[donor] != component[target]


def test_training_donors_come_from_meta_train(data):
    donors = matched_donors(data, "meta_train", donor_pool="meta_train")
    component = {cell["target_id"]: cell["protein_group_40"]
                 for cell in data.cells}
    meta_train = set(data.tasks["meta_train"])
    assert set(donors) == meta_train
    for target, donor in donors.items():
        assert donor in meta_train
        assert component[donor] != component[target]


def test_whitening_is_fitted_on_meta_train_only(data):
    """The fit pool is a real, exercised choice, not a cosmetic parameter."""
    train_fit = matched_donors(data, "meta_val", donor_pool="meta_val",
                               whitening_pool="meta_train")
    eval_fit = matched_donors(data, "meta_val", donor_pool="meta_val",
                              whitening_pool="meta_val")
    differences = [target for target in train_fit
                   if train_fit[target] != eval_fit[target]]
    assert differences, "whitening fit pool did not change any donor selection"


def test_matched_donors_reproduces_the_meta_train_transform(data):
    """The returned donor is the cross-component argmax under the transform
    whose mean and covariance come from meta_train and nothing else."""
    train_targets = sorted(data.tasks["meta_train"])
    eval_targets = sorted(data.tasks["meta_val"])
    pooled = {target: np.asarray(data.protein_for_target(target)[0],
                                 dtype=np.float32)
              for target in train_targets + eval_targets}
    matrix = np.stack([pooled[target] for target in train_targets])
    center = matrix.mean(0, keepdims=True)
    deviation = matrix - center
    covariance = deviation.T @ deviation / max(len(matrix) - 1, 1)
    values, vectors = np.linalg.eigh(covariance.astype(np.float64))
    whiten = ((vectors / np.sqrt(np.maximum(values, 1e-3)))
              @ vectors.T).astype(np.float32)

    def transform(vector):
        out = (vector - center[0]) @ whiten.T
        return out / max(float(np.linalg.norm(out)), 1e-9)

    component = {cell["target_id"]: cell["protein_group_40"]
                 for cell in data.cells}
    bank = np.stack([transform(pooled[target]) for target in eval_targets])
    donors = matched_donors(data, "meta_val", donor_pool="meta_val",
                            whitening_pool="meta_train")
    for target in eval_targets:
        similarity = bank @ transform(pooled[target])
        chosen = None
        for index in np.argsort(-similarity):
            candidate = eval_targets[int(index)]
            if component[candidate] != component[target]:
                chosen = candidate
                break
        assert donors[target] == chosen


def test_cli_seal_flags_exist_and_default_off():
    """The opt-in seal flags are present on every default evaluation path."""
    for module in ("scripts.train_qpsmp", "scripts.evaluate_qpsmp",
                   "scripts.evaluate_checkpoint_nested"):
        completed = subprocess.run(
            [str(PYTHON), "-m", module, "--help"], cwd=ROOT,
            capture_output=True, text=True, timeout=300,
            env=dict(os.environ, PYTHONIOENCODING="utf-8"))
        assert completed.returncode == 0, completed.stderr[-2000:]
        assert "--include-meta-test" in completed.stdout, module
    help_text = subprocess.run(
        [str(PYTHON), "-m", "scripts.evaluate_qpsmp", "--help"], cwd=ROOT,
        capture_output=True, text=True, timeout=300,
        env=dict(os.environ, PYTHONIOENCODING="utf-8")).stdout
    assert "meta_val" in help_text        # default split is the development split


@pytest.mark.skipif(os.environ.get("RUN_SLOW") != "1",
                    reason="slow: loads the full corpus and trains 2 steps")
def test_training_does_not_evaluate_meta_test_by_default(tmp_path):
    """A default train_qpsmp run must record meta_test as unevaluated."""
    output = tmp_path / "run"
    command = [
        str(PYTHON), "-m", "scripts.train_qpsmp",
        "--split-directory", str(SPLIT),
        "--output", str(output),
        "--steps", "2", "--episodes-per-step", "1",
        "--val-interval", "1", "--arch", "similarity_only",
    ]
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    completed = subprocess.run(command, cwd=ROOT, env=env,
                               capture_output=True, text=True, timeout=1800)
    assert completed.returncode == 0, completed.stderr[-3000:]
    payload = json.loads((output / "RESULT.json").read_text(encoding="utf-8"))
    assert payload["meta_test"]["evaluated"] is False
    assert payload["meta_test"]["included"] is False
    assert payload["test"] == {}
    assert payload["checkpoint_sha256"]


@pytest.mark.skipif(os.environ.get("RUN_SLOW") != "1",
                    reason="slow: loads the full corpus and trains 2 steps")
def test_train_level_shape_saves_contract_artifacts(tmp_path):
    """RESULT.json carries donor pools, checkpoint hash and activation stats."""
    output = tmp_path / "run"
    command = [
        str(PYTHON), "-m", "scripts.train_level_shape",
        "--split-directory", str(SPLIT),
        "--output", str(output),
        "--steps", "2", "--episodes-per-step", "1",
        "--val-interval", "1", "--device", "cpu",
    ]
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    completed = subprocess.run(command, cwd=ROOT, env=env,
                               capture_output=True, text=True, timeout=3600)
    assert completed.returncode == 0, completed.stderr[-3000:]
    payload = json.loads((output / "RESULT.json").read_text(encoding="utf-8"))
    assert payload["donors"]["evaluation_wrong_protein_pool"] == "meta_val"
    assert payload["donors"]["whitening_pool"] == "meta_train"
    assert payload["meta_test"]["evaluated"] is False
    assert payload["checkpoint_sha256"]
    rows = [json.loads(line) for line in
            (output / "PREDICTIONS_meta_val.jsonl").read_text(
                encoding="utf-8").splitlines()]
    assert rows and "centered_mean_pk" in rows[0]
    assert "transport_abs_mean_pk" in rows[0]


@pytest.mark.skipif(os.environ.get("RUN_SLOW") != "1",
                    reason="slow: loads the full corpus and trains 2 steps")
def test_train_reltransport_saves_contract_artifacts(tmp_path):
    """The relative-transport runner honors the same contract."""
    output = tmp_path / "run"
    command = [
        str(PYTHON), "-m", "scripts.train_reltransport",
        "--split-directory", str(SPLIT),
        "--output", str(output),
        "--steps", "2", "--episodes-per-step", "1",
        "--val-interval", "1", "--device", "cpu",
    ]
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    completed = subprocess.run(command, cwd=ROOT, env=env,
                               capture_output=True, text=True, timeout=3600)
    assert completed.returncode == 0, completed.stderr[-3000:]
    payload = json.loads((output / "RESULT.json").read_text(encoding="utf-8"))
    assert payload["donors"]["evaluation_wrong_protein_pool"] == "meta_val"
    assert payload["donors"]["whitening_pool"] == "meta_train"
    assert payload["meta_test"]["evaluated"] is False
    assert payload["checkpoint_sha256"]
    assert "gradient_summary" in payload["training"]
    assert "gradient_coverage" in payload["training"]
