from __future__ import annotations

from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shutil
import sys

import pytest

import research.shared.dataset_processing as dataset_processing
from research.shared.dataset_processing import (
    DRUG_PYTHON,
    DatasetProcessPolicyError,
    _H0_INVENTORY_OUTPUT,
    _automatic_bundles,
    _presence_outputs,
    _process_identity_matches,
    _projection_outputs,
    _started_identity_complete,
    _validate_command_policy,
    execute_run,
    reconcile_orphans,
    record_retrospective,
    relocate_file_run,
    require_run_context,
    terminate_active_run,
    verify_ledger,
)


def _allow_test_command(*_arguments: object) -> None:
    return None


def _json_hash(value: object) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _write_single_ledger_entry(path: Path, entry: dict[str, object]) -> None:
    unsigned = {key: value for key, value in entry.items() if key != "entry_sha256"}
    entry["entry_sha256"] = _json_hash(unsigned)
    path.write_text(
        json.dumps(entry, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def test_successful_run_freezes_command_logs_artifacts_and_chain(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    source = workspace / "source.bin"
    source.write_bytes(b"source bytes")
    output = workspace / "result.bin"
    code = (
        "from pathlib import Path; import sys; "
        "Path(sys.argv[1]).write_bytes(b'result bytes'); "
        "print('captured stdout'); print('captured stderr', file=sys.stderr)"
    )
    expected = sha256(b"result bytes").hexdigest()

    result = execute_run(
        operation="test-success",
        command=[sys.executable, "-c", code, str(output)],
        cwd=workspace,
        root=workspace / "history",
        inputs=[source],
        outputs=[output],
        code_inputs=[source],
        expected_output_bytes={output: len(b"result bytes")},
        expected_output_sha256={output: expected},
        _policy_validator=_allow_test_command,
    )

    run_dir = Path(result["run_path"])
    prepared = json.loads((run_dir / "prepared.json").read_text(encoding="utf-8"))
    frozen = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert result["status"] == "SUCCESS"
    assert prepared["command_argv"][2] == code
    assert prepared["cwd"] == str(workspace)
    assert prepared["inputs_before"][0]["sha256"] == sha256(b"source bytes").hexdigest()
    assert frozen["process_status"] == "completed"
    assert frozen["validation_status"] == "pass"
    assert frozen["exit_code"] == 0
    started = json.loads((run_dir / "started.json").read_text(encoding="utf-8"))
    assert _started_identity_complete(started)
    assert frozen["validated_outputs_before_publish"][0]["sha256"] == expected
    assert "captured stdout" in (run_dir / "stdout.log").read_text(encoding="utf-8")
    assert "captured stderr" in (run_dir / "stderr.log").read_text(encoding="utf-8")
    assert (run_dir / "intent.json").exists()
    assert (run_dir / "started.json").exists()
    assert (run_dir / "seal.json").exists()
    assert verify_ledger(workspace / "history", workspace)["entries"] == 1


def test_nonzero_exit_and_missing_output_are_indexed_as_failure(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    output = workspace / "missing.bin"

    result = execute_run(
        operation="test-failure",
        command=[sys.executable, "-c", "import sys; sys.exit(7)"],
        cwd=workspace,
        root=workspace / "history",
        outputs=[output],
        _policy_validator=_allow_test_command,
    )

    assert result["status"] == "FAILED"
    assert result["process_status"] == "nonzero_exit"
    assert result["validation_status"] == "fail"
    assert result["exit_code"] == 7
    assert verify_ledger(workspace / "history", workspace)["status"] == "PASS"


def test_exit_zero_cannot_override_expected_hash_failure(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    output = workspace / "result.bin"
    code = "from pathlib import Path; import sys; Path(sys.argv[1]).write_bytes(b'wrong')"

    result = execute_run(
        operation="test-hash-failure",
        command=[sys.executable, "-c", code, str(output)],
        cwd=workspace,
        root=workspace / "history",
        outputs=[output],
        expected_output_sha256={output: "0" * 64},
        _policy_validator=_allow_test_command,
    )

    assert result["exit_code"] == 0
    assert result["process_status"] == "completed"
    assert result["validation_status"] == "fail"
    assert result["status"] == "FAILED"


def test_preexisting_output_is_never_overwritten(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    output = workspace / "result.bin"
    output.write_bytes(b"keep")

    result = execute_run(
        operation="test-no-overwrite",
        command=[sys.executable, "-c", "raise SystemExit('must not run')"],
        cwd=workspace,
        root=workspace / "history",
        outputs=[output],
        _policy_validator=_allow_test_command,
    )

    assert result["process_status"] == "preparation_failed"
    assert output.read_bytes() == b"keep"
    assert Path(result["run_path"], "run.json").exists()


def test_tampering_is_detected(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    output = workspace / "result.bin"
    code = "from pathlib import Path; import sys; Path(sys.argv[1]).write_bytes(b'ok')"
    result = execute_run(
        operation="test-tamper",
        command=[sys.executable, "-c", code, str(output)],
        cwd=workspace,
        root=workspace / "history",
        outputs=[output],
        _policy_validator=_allow_test_command,
    )
    run_dir = Path(result["run_path"])
    (run_dir / "stdout.log").write_text("tampered", encoding="utf-8")

    with pytest.raises(ValueError, match="sealed file mismatch"):
        verify_ledger(workspace / "history", workspace)


def test_ledger_rejects_hash_valid_run_path_outside_history(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    result = execute_run(
        operation="test-external-run-path",
        command=[sys.executable, "-c", "print('ok')"],
        cwd=workspace,
        root=workspace / "history",
        _policy_validator=_allow_test_command,
    )
    external = workspace / "external-run"
    shutil.copytree(Path(result["run_path"]), external)
    ledger = workspace / "history" / "ledger.jsonl"
    entry = json.loads(ledger.read_text(encoding="utf-8"))
    entry["run_path"] = str(external)
    _write_single_ledger_entry(ledger, entry)

    with pytest.raises(ValueError, match="run_path escapes"):
        verify_ledger(workspace / "history", workspace)


def test_ledger_rejects_hash_valid_external_sealed_file(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    result = execute_run(
        operation="test-external-sealed-file",
        command=[sys.executable, "-c", "print('ok')"],
        cwd=workspace,
        root=workspace / "history",
        _policy_validator=_allow_test_command,
    )
    external = workspace / "external.bin"
    external.write_bytes(b"outside processing history")
    run_dir = Path(result["run_path"])
    seal_path = run_dir / "seal.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    seal["files"][str(external)] = sha256(external.read_bytes()).hexdigest()
    seal_unsigned = {
        key: value for key, value in seal.items() if key != "content_sha256"
    }
    seal["content_sha256"] = _json_hash(seal_unsigned)
    seal_path.write_text(
        json.dumps(seal, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    ledger = workspace / "history" / "ledger.jsonl"
    entry = json.loads(ledger.read_text(encoding="utf-8"))
    entry["seal_sha256"] = sha256(seal_path.read_bytes()).hexdigest()
    _write_single_ledger_entry(ledger, entry)

    with pytest.raises(ValueError, match="invalid sealed file name"):
        verify_ledger(workspace / "history", workspace)


def test_retrospective_record_never_claims_contemporaneous_logs(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    artifact = workspace / "observed.bin"
    artifact.write_bytes(b"observed now")
    event_file = workspace / "event.json"
    event_file.write_text(
        json.dumps(
            {
                "schema": "a2s-dataset-retrospective-event-v1",
                "operation": "past-download",
                "historical_event_status": "FAILED",
                "evidence": {"known_sha256": "f" * 64},
                "unknown_fields": ["pid", "exit_code", "stdout", "stderr"],
                "limitations": ["No contemporaneous stdout, PID, or exit code survived."],
            }
        ),
        encoding="utf-8",
    )

    result = record_retrospective(
        event_file=event_file,
        artifacts=[artifact],
        cwd=workspace,
        root=workspace / "history",
    )

    prepared = json.loads(
        Path(result["run_path"], "prepared.json").read_text(encoding="utf-8")
    )
    assert result["status"] == "RECORDED_RETROSPECTIVE"
    assert prepared["record_type"] == "retrospective_reconstruction"
    assert prepared["contemporaneous_stdout_available"] is False
    assert prepared["artifacts_observed_at_reconstruction"][0]["sha256"] == sha256(
        b"observed now"
    ).hexdigest()
    assert verify_ledger(workspace / "history", workspace)["entries"] == 1


def test_policy_blocks_training_and_non_drug_python(tmp_path: Path) -> None:
    with pytest.raises(DatasetProcessPolicyError, match="drug"):
        _validate_command_policy(
            [r"C:\Windows\System32\cmd.exe", "main.py", "train"], tmp_path
        )


def test_metadata_policy_freezes_projection_primary_and_replay_outputs(
    tmp_path: Path,
) -> None:
    workspace = tmp_path.resolve()
    main = workspace / "main.py"
    release = "chembl_24_1"
    projection, documents, manifest, certificate = _projection_outputs(
        workspace, release
    )
    primary = [
        str(DRUG_PYTHON),
        str(main),
        "historical-project",
        "--release",
        release,
    ]

    _validate_command_policy(
        primary, workspace, [projection, documents, manifest]
    )
    with pytest.raises(DatasetProcessPolicyError, match="requires all frozen"):
        _validate_command_policy([*primary, "--replay"], workspace, [certificate])
    with pytest.raises(DatasetProcessPolicyError, match="exact committed output"):
        _validate_command_policy(primary, workspace, [projection, manifest])
    with pytest.raises(DatasetProcessPolicyError, match="exact committed output"):
        _validate_command_policy(
            primary, workspace, [projection, documents, manifest, manifest]
        )
    with pytest.raises(DatasetProcessPolicyError, match="exact committed output"):
        _validate_command_policy(
            primary, workspace, [projection, documents, manifest, certificate]
        )
    with pytest.raises(DatasetProcessPolicyError, match="requires exactly"):
        _validate_command_policy(
            [*primary, "--extract-only"],
            workspace,
            [projection, documents, manifest],
        )
    for forbidden in (
        [*primary, "--verify-only"],
        [*primary, "--replay", "--extract-only"],
        [*primary, "--release", "chembl_27"],
        [*primary, "--unknown"],
    ):
        with pytest.raises(DatasetProcessPolicyError, match="requires exactly"):
            _validate_command_policy(
                forbidden, workspace, [projection, documents, manifest]
            )

    for path in (projection, documents, manifest):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"primary")
    replay = [*primary, "--replay"]
    _validate_command_policy(replay, workspace, [certificate])
    with pytest.raises(DatasetProcessPolicyError, match="exact committed output"):
        _validate_command_policy(replay, workspace, [projection, certificate])


def test_metadata_policy_freezes_presence_and_inventory_modes(tmp_path: Path) -> None:
    workspace = tmp_path.resolve()
    main = workspace / "main.py"
    presence, manifest, certificate = _presence_outputs(workspace)
    presence_command = [str(DRUG_PYTHON), str(main), "historical-presence"]

    _validate_command_policy(
        presence_command, workspace, [presence, manifest]
    )
    with pytest.raises(DatasetProcessPolicyError, match="requires the frozen"):
        _validate_command_policy(
            [*presence_command, "--replay"], workspace, [certificate]
        )
    with pytest.raises(DatasetProcessPolicyError, match="permits only"):
        _validate_command_policy(
            [*presence_command, "--output", str(presence)],
            workspace,
            [presence, manifest],
        )
    presence.parent.mkdir(parents=True, exist_ok=True)
    presence.write_bytes(b"presence")
    manifest.write_bytes(b"manifest")
    _validate_command_policy(
        [*presence_command, "--replay"], workspace, [certificate]
    )

    inventory = (workspace / _H0_INVENTORY_OUTPUT).resolve()
    inventory_command = [str(DRUG_PYTHON), str(main), "h0-inventory"]
    _validate_command_policy(inventory_command, workspace, [inventory])
    _validate_command_policy(
        [*inventory_command, "--output", str(inventory)], workspace, [inventory]
    )
    with pytest.raises(DatasetProcessPolicyError, match="local v3"):
        _validate_command_policy(
            [*inventory_command, "--acquire-remote-metadata"],
            workspace,
            [inventory],
        )
    with pytest.raises(DatasetProcessPolicyError, match="publish mappings"):
        _validate_command_policy(
            inventory_command,
            workspace,
            [inventory],
            publish_outputs={str(inventory): str(inventory.with_suffix(".published"))},
        )


def test_task_markdown_is_an_automatic_config_binding(tmp_path: Path) -> None:
    workspace = tmp_path.resolve()
    (workspace / "main.py").write_text("", encoding="utf-8")
    task = workspace / "task.md"
    task.write_text("frozen task\n", encoding="utf-8")

    _, config = _automatic_bundles(
        [str(DRUG_PYTHON), str(workspace / "main.py"), "h0-inventory"],
        workspace,
    )

    assert task in config


def test_process_identity_rejects_reused_pid() -> None:
    expected = {
        "pid": 123,
        "creation_time_filetime": 456,
        "executable_image_path": r"D:\anaconda\envs\drug\python.exe",
        "command_sha256": "a" * 64,
    }
    observed = {
        "pid": 123,
        "creation_time_filetime": 456,
        "executable_image_path": r"d:/ANACONDA/envs/drug/python.exe",
    }

    assert _started_identity_complete(expected)
    assert _process_identity_matches(expected, observed)
    assert not _process_identity_matches(
        expected, {**observed, "creation_time_filetime": 457}
    )


def test_terminate_never_signals_a_reused_pid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace = tmp_path / "workspace"
    target_run_id = "target-run"
    target_dir = workspace / "history" / "runs" / target_run_id
    target_dir.mkdir(parents=True)
    started = {
        "pid": 4321,
        "creation_time_filetime": 100,
        "executable_image_path": r"D:\anaconda\envs\drug\python.exe",
        "command_sha256": "a" * 64,
    }
    (target_dir / "started.json").write_text(json.dumps(started), encoding="utf-8")
    (target_dir / "prepared.json").write_text("{}", encoding="utf-8")

    class FakeKernel:
        terminate_calls = 0

        def TerminateProcess(self, _handle: object, _code: int) -> bool:
            self.terminate_calls += 1
            return True

        def CloseHandle(self, _handle: object) -> bool:
            return True

    kernel = FakeKernel()
    monkeypatch.setattr(
        dataset_processing,
        "_open_windows_process",
        lambda _pid, _access: (kernel, object(), 0),
    )
    monkeypatch.setattr(
        dataset_processing,
        "_windows_process_details",
        lambda _kernel, _handle, pid: {
            "pid": pid,
            "creation_time_filetime": 101,
            "executable_image_path": started["executable_image_path"],
            "running": True,
            "observed_exit_code": 259,
        },
    )

    with pytest.raises(DatasetProcessPolicyError, match="identity mismatch"):
        terminate_active_run(
            target_run_id=target_run_id,
            reason="test only",
            cwd=workspace,
            root=workspace / "history",
        )

    assert kernel.terminate_calls == 0


def test_legacy_started_identity_is_never_reconciled_or_terminated(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "workspace"
    run_dir = workspace / "history" / "runs" / "legacy-active"
    run_dir.mkdir(parents=True)
    (run_dir / "intent.json").write_text(
        json.dumps(
            {
                "schema": "a2s-dataset-processing-intent-v1",
                "run_id": run_dir.name,
                "record_type": "contemporaneous_execution",
                "operation": "legacy-active",
                "declared_outputs": [],
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "prepared.json").write_text("{}", encoding="utf-8")
    (run_dir / "started.json").write_text(
        json.dumps({"pid": 999_999, "started_utc": "legacy"}), encoding="utf-8"
    )
    (run_dir / "events.jsonl").write_text("", encoding="utf-8")

    report = reconcile_orphans(workspace / "history", workspace)

    assert report["status"] == "IDENTITY_UNKNOWN_OPERATIONS_PRESENT"
    assert report["identity_unknown"] == [run_dir.name]
    assert not (run_dir / "run.json").exists()
    with pytest.raises(DatasetProcessPolicyError, match="IDENTITY_UNKNOWN"):
        terminate_active_run(
            target_run_id=run_dir.name,
            reason="must not signal",
            cwd=workspace,
            root=workspace / "history",
        )


def test_policy_failure_never_records_secret_or_secret_hash(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    stage = workspace / "chembl_37_sqlite.tar.gz.partial"
    canonical = workspace / "chembl_37_sqlite.tar.gz"
    secret = "Authorization: bearer test-secret-value"
    command = [
        r"C:\Windows\System32\curl.exe",
        "--fail",
        "--location",
        "--continue-at",
        "-",
        "--header",
        secret,
        "--output",
        str(stage),
        "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/chembl_37_sqlite.tar.gz",
    ]
    result = execute_run(
        operation="test-secret-rejection",
        command=command,
        cwd=workspace,
        root=workspace / "history",
        outputs=[stage],
        expected_output_bytes={stage: 1},
        expected_output_sha256={stage: "0" * 64},
        publish_outputs={stage: canonical},
    )
    recorded = b"\n".join(
        path.read_bytes()
        for path in Path(result["run_path"]).iterdir()
        if path.is_file()
    )

    assert result["status"] == "FAILED"
    assert secret.encode("utf-8") not in recorded
    assert sha256(secret.encode("utf-8")).hexdigest().encode("ascii") not in recorded


def test_curl_policy_requires_official_staging_and_publish(tmp_path: Path) -> None:
    stage = (tmp_path / "chembl_37_sqlite.tar.gz.partial").resolve()
    canonical = (tmp_path / "chembl_37_sqlite.tar.gz").resolve()
    url = (
        "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
        "chembl_37/chembl_37_sqlite.tar.gz"
    )
    command = [
        r"C:\Windows\System32\curl.exe",
        "--fail",
        "--location",
        "--continue-at",
        "-",
        "--output",
        str(stage),
        url,
    ]

    _validate_command_policy(
        command,
        tmp_path,
        [stage],
        {str(stage): 123},
        {str(stage): "0" * 64},
        {str(stage): str(canonical)},
    )
    with pytest.raises(DatasetProcessPolicyError, match="credentials"):
        _validate_command_policy(
            [*command, "--header", "Authorization: secret"],
            tmp_path,
            [stage],
            {str(stage): 123},
            {str(stage): "0" * 64},
            {str(stage): str(canonical)},
        )


def test_validated_staging_output_is_atomically_published(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    stage = workspace / "artifact.bin.partial"
    canonical = workspace / "artifact.bin"
    code = "from pathlib import Path; import sys; Path(sys.argv[1]).write_bytes(b'valid')"
    result = execute_run(
        operation="test-publish",
        command=[sys.executable, "-c", code, str(stage)],
        cwd=workspace,
        root=workspace / "history",
        outputs=[stage],
        expected_output_bytes={stage: 5},
        expected_output_sha256={stage: sha256(b"valid").hexdigest()},
        publish_outputs={stage: canonical},
        _policy_validator=_allow_test_command,
    )

    assert result["status"] == "SUCCESS"
    assert result["publish_status"] == "published"
    assert canonical.read_bytes() == b"valid"
    assert not stage.exists()
    run_dir = Path(result["run_path"])
    claim = json.loads((run_dir / "output_claim.json").read_text(encoding="utf-8"))
    assert claim["run_id"] == result["run_id"]
    assert {item["path"] for item in claim["paths"]} == {
        str(stage.resolve()),
        str(canonical.resolve()),
    }
    seal = json.loads((run_dir / "seal.json").read_text(encoding="utf-8"))
    assert "output_claim.json" in seal["files"]


def test_direct_protected_command_requires_matching_run_context(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("A2S_DATASET_RUN_DIR", raising=False)
    monkeypatch.delenv("A2S_DATASET_RUN_ID", raising=False)
    with pytest.raises(SystemExit, match="dataset-run"):
        require_run_context(["main.py", "historical-project"], tmp_path)


def test_concurrent_runs_form_one_ledger_chain(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    history = workspace / "history"

    def run(index: int) -> str:
        output = workspace / f"result-{index}.bin"
        code = "from pathlib import Path; import sys; Path(sys.argv[1]).write_bytes(b'x')"
        return execute_run(
            operation=f"concurrent-{index}",
            command=[sys.executable, "-c", code, str(output)],
            cwd=workspace,
            root=history,
            outputs=[output],
            _policy_validator=_allow_test_command,
        )["status"]

    with ThreadPoolExecutor(max_workers=4) as pool:
        statuses = list(pool.map(run, range(4)))

    assert statuses == ["SUCCESS"] * 4
    assert verify_ledger(history, workspace)["entries"] == 4


def test_concurrent_runs_cannot_claim_the_same_output(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    history = workspace / "history"
    output = workspace / "shared.bin"
    code = (
        "from pathlib import Path; import sys, time; "
        "time.sleep(0.75); Path(sys.argv[1]).write_bytes(b'one owner')"
    )

    def run(index: int) -> dict[str, object]:
        return execute_run(
            operation=f"exclusive-{index}",
            command=[sys.executable, "-c", code, str(output)],
            cwd=workspace,
            root=history,
            outputs=[output],
            _policy_validator=_allow_test_command,
        )

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(run, range(2)))

    assert sorted(result["status"] for result in results) == ["FAILED", "SUCCESS"]
    failed = next(result for result in results if result["status"] == "FAILED")
    assert failed["process_status"] == "preparation_failed"
    assert "already reserved" in failed["error"]["message"]
    assert output.read_bytes() == b"one owner"
    assert verify_ledger(history, workspace)["entries"] == 2


def test_orphan_reconciliation_preserves_unknown_exit_fields(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    run_dir = workspace / "history" / "runs" / "orphan-run"
    run_dir.mkdir(parents=True)
    (run_dir / "intent.json").write_text(
        json.dumps(
            {
                "schema": "a2s-dataset-processing-intent-v1",
                "run_id": "orphan-run",
                "record_type": "contemporaneous_execution",
                "operation": "orphan-test",
                "declared_outputs": [],
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "events.jsonl").write_text("", encoding="utf-8")
    (run_dir / "stdout.log").write_bytes(b"partial output")
    (run_dir / "stderr.log").write_bytes(b"")

    report = reconcile_orphans(workspace / "history", workspace)
    recovered = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))

    assert report["recovered"] == ["orphan-run"]
    assert recovered["process_status"] == "orphaned_unknown"
    assert recovered["exit_code"] is None
    assert recovered["child_ended_utc"] is None
    assert verify_ledger(workspace / "history", workspace)["status"] == "PASS"


def test_orphan_after_publish_records_staging_canonical_and_partial_artifacts(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "workspace"
    run_dir = workspace / "history" / "runs" / "publish-orphan"
    run_dir.mkdir(parents=True)
    stage = workspace / "archive.tar.gz.partial"
    canonical = workspace / "archive.tar.gz"
    projection = workspace / "activity_identity.parquet"
    projection_partial = workspace / ".activity_identity.parquet.partial"
    canonical.write_bytes(b"published bytes")
    projection_partial.write_bytes(b"projection partial")
    (run_dir / "intent.json").write_text(
        json.dumps(
            {
                "schema": "a2s-dataset-processing-intent-v1",
                "run_id": run_dir.name,
                "record_type": "contemporaneous_execution",
                "operation": "publish-orphan",
                "declared_outputs": [str(stage), str(projection)],
                "publish_outputs": {str(stage): str(canonical)},
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "events.jsonl").write_text("", encoding="utf-8")

    report = reconcile_orphans(workspace / "history", workspace)
    recovered = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    artifacts = recovered["artifacts_observed_after_recovery"]

    assert report["recovered"] == [run_dir.name]
    assert recovered["status"] == "FAILED"
    assert recovered["publish_success_inferred"] is False
    assert artifacts["publish_staging_sources"][0]["kind"] == "missing"
    assert artifacts["publish_canonical_destinations"][0]["sha256"] == sha256(
        b"published bytes"
    ).hexdigest()
    assert artifacts["dot_partial_candidates"][0]["sha256"] == sha256(
        b"projection partial"
    ).hexdigest()


def test_orphan_relocation_records_both_endpoints_without_inferring_success(
    tmp_path: Path,
) -> None:
    workspace = tmp_path / "workspace"
    run_dir = workspace / "history" / "runs" / "relocation-orphan"
    run_dir.mkdir(parents=True)
    source = workspace / "legacy.bin"
    destination = workspace / "legacy.bin.prelogging-partial"
    destination.write_bytes(b"relocated bytes")
    (run_dir / "intent.json").write_text(
        json.dumps(
            {
                "schema": "a2s-dataset-processing-intent-v1",
                "run_id": run_dir.name,
                "record_type": "contemporaneous_internal_relocation",
                "operation": "relocation-orphan",
                "source": str(source),
                "destination": str(destination),
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "events.jsonl").write_text("", encoding="utf-8")

    reconcile_orphans(workspace / "history", workspace)
    recovered = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    artifacts = recovered["artifacts_observed_after_recovery"]

    assert recovered["status"] == "FAILED"
    assert artifacts["relocation_source"]["kind"] == "missing"
    assert artifacts["relocation_destination"]["sha256"] == sha256(
        b"relocated bytes"
    ).hexdigest()


def test_verify_reports_a_valid_unindexed_active_run_as_non_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    workspace = tmp_path / "workspace"
    run_dir = workspace / "history" / "runs" / "active-run"
    run_dir.mkdir(parents=True)
    command = [str(DRUG_PYTHON), str(workspace / "main.py"), "h0-inventory"]
    (run_dir / "intent.json").write_text(
        json.dumps(
            {
                "run_id": run_dir.name,
                "record_type": "contemporaneous_execution",
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "prepared.json").write_text(
        json.dumps(
            {
                "run_id": run_dir.name,
                "record_type": "contemporaneous_execution",
                "command_argv": command,
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "started.json").write_text(
        json.dumps(
            {
                "run_id": run_dir.name,
                "pid": 123,
                "creation_time_filetime": 456,
                "executable_image_path": r"D:\anaconda\envs\drug\python.exe",
                "command_sha256": _json_hash(command),
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(
        dataset_processing,
        "_inspect_started_process",
        lambda _started: ("MATCH", {"pid": 123}),
    )

    report = verify_ledger(workspace / "history", workspace)

    assert report["status"] == "RUNNING_OPERATIONS_PRESENT"
    assert report["running_unindexed_runs"] == [run_dir.name]


def test_network_request_log_is_sealed_when_present(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    output = workspace / "result.bin"
    code = (
        "from pathlib import Path; import os, sys; "
        "Path(os.environ['A2S_DATASET_RUN_DIR'], 'network_requests.jsonl')"
        ".write_text('{\\\"url\\\":\\\"https://example.invalid\\\"}\\n'); "
        "Path(sys.argv[1]).write_bytes(b'ok')"
    )
    result = execute_run(
        operation="test-network-seal",
        command=[sys.executable, "-c", code, str(output)],
        cwd=workspace,
        root=workspace / "history",
        outputs=[output],
        _policy_validator=_allow_test_command,
    )
    run_dir = Path(result["run_path"])
    seal = json.loads((run_dir / "seal.json").read_text(encoding="utf-8"))

    assert result["status"] == "SUCCESS"
    assert "network_requests.jsonl" in seal["files"]
    (run_dir / "network_requests.jsonl").write_text("tampered\n", encoding="utf-8")
    with pytest.raises(ValueError, match="sealed file mismatch"):
        verify_ledger(workspace / "history", workspace)


def test_legacy_relocation_preserves_bytes_and_hash(tmp_path: Path) -> None:
    workspace = tmp_path / "workspace"
    archives = workspace / "dataset" / "public" / "chembl_historical" / "archives"
    archives.mkdir(parents=True)
    source = archives / "legacy.tar.gz"
    destination = archives / "legacy.tar.gz.prelogging-partial"
    source.write_bytes(b"partial bytes")

    result = relocate_file_run(
        operation="relocate-legacy",
        source=source,
        destination=destination,
        expected_bytes=len(b"partial bytes"),
        cwd=workspace,
        root=workspace / "history",
    )

    assert result["status"] == "SUCCESS"
    assert not source.exists()
    assert destination.read_bytes() == b"partial bytes"
    assert verify_ledger(workspace / "history", workspace)["entries"] == 1
