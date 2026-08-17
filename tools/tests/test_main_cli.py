import importlib
import json
import sys
from types import SimpleNamespace

import pytest

import main as cli
from scripts import train_qpsmp
from scripts import preprocess_dataset, project_status, verify_dataset


def test_root_dispatch_exposes_only_promoted_script_modules():
    assert cli.COMMANDS
    assert all(command.module.startswith("scripts.")
               for command in cli.COMMANDS.values())
    assert all("research" not in command.module
               for command in cli.COMMANDS.values())


def test_root_and_group_help_expose_only_retained_commands(capsys):
    assert cli.main(["--help"]) == 0
    help_text = capsys.readouterr().out
    for name in ("status", "archive status", "verify tests",
                 "qpsmp train", "qpsmp evaluate", "geometry pretrain",
                 "data prepare", "data verify"):
        assert name in help_text
    for rejected in ("r0c", "tbasis", "scheduler", "crossed"):
        assert rejected not in help_text.lower()
    assert cli.main(["qpsmp", "--help"]) == 0
    assert "train" in capsys.readouterr().out


def test_dispatch_is_lazy_forwards_argv_and_restores_process_state(monkeypatch):
    calls = []

    class Module:
        @staticmethod
        def main():
            calls.append(list(sys.argv))
            return 7

    monkeypatch.setattr(importlib, "import_module", lambda name: Module)
    before = list(sys.argv)
    assert cli.main(["archive", "status", "--json"]) == 7
    assert calls == [["scripts.project_status", "--archive", "--json"]]
    assert sys.argv == before


def test_cuda_commands_default_to_cuda0_and_reject_other_devices(monkeypatch, capsys):
    calls = []
    module = SimpleNamespace(main=lambda: calls.append(list(sys.argv)) or 0)
    monkeypatch.setattr(importlib, "import_module", lambda name: module)
    assert cli.main(["qpsmp", "train"]) == 0
    assert calls[0][-2:] == ["--device", "cuda:0"]
    assert cli.main([
        "qpsmp", "train", "--output", "report/meta_fewshot/new",
        "--device", "cpu",
    ]) == 2
    assert "require --device cuda:0" in capsys.readouterr().err
    assert len(calls) == 1


def test_leaf_help_is_forwarded_without_cuda_initialization(monkeypatch):
    calls = []
    module = SimpleNamespace(main=lambda: calls.append(list(sys.argv)) or 0)
    monkeypatch.setattr(importlib, "import_module", lambda name: module)
    assert cli.main(["qpsmp", "train", "--help"]) == 0
    assert calls == [["scripts.train_qpsmp", "--help"]]


def test_status_and_archive_views_are_read_only_contracts():
    status = project_status.load_status()
    assert status["core_task"]["target_is_meta_task"] is True
    assert "current_stage" in status and "unresolved" in status
    archive = project_status.load_status(archive_only=True)
    assert "not active entry points" in archive["policy"]


def test_dataset_verifier_reports_invalid_without_writing(tmp_path, monkeypatch, capsys):
    compiled = tmp_path / "compiled"
    compiled.mkdir()
    (compiled / "rows.jsonl").write_text("", encoding="utf-8")
    monkeypatch.setattr(sys, "argv", ["verify_dataset", str(compiled)])
    assert verify_dataset.main() == 1
    assert json.loads(capsys.readouterr().out)["valid"] is False
    assert sorted(path.name for path in compiled.iterdir()) == ["rows.jsonl"]


def test_compile_dataset_rejects_existing_output_before_writing(tmp_path):
    output = tmp_path / "existing"
    output.mkdir()
    marker = output / "marker"
    marker.write_text("keep", encoding="utf-8")
    with pytest.raises(FileExistsError, match="already exists"):
        preprocess_dataset.compile_dataset(tmp_path / "missing.json", output)
    assert marker.read_text(encoding="utf-8") == "keep"


def test_qpsmp_requires_safe_new_output_before_data_loading(tmp_path, monkeypatch):
    existing = tmp_path / "existing"
    existing.mkdir()
    monkeypatch.setattr(sys, "argv", [
        "train_qpsmp", "--output", str(existing),
    ])
    with pytest.raises(FileExistsError, match="already exists"):
        train_qpsmp.main()


def test_data_prepare_cli_rejects_output_outside_processed_root(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "preprocess_dataset", str(tmp_path / "missing.json"),
        str(tmp_path / "outside"),
    ])
    with pytest.raises(ValueError, match="must be a new child"):
        preprocess_dataset.main()
