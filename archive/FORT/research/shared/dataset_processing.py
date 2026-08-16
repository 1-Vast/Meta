"""Run dataset operations with immutable, tamper-evident audit records."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import traceback
from typing import Any, Iterable, Mapping, Sequence
from uuid import uuid4
from urllib.parse import urlsplit


DEFAULT_ROOT = Path("dataset/processing_history/v1")
DRUG_PYTHON = Path(r"D:\anaconda\envs\drug\python.exe")
_OPERATION = re.compile(r"^[a-z0-9][a-z0-9._-]{1,79}$")
_ENVIRONMENT_KEYS = (
    "COMSPEC",
    "NUMBER_OF_PROCESSORS",
    "OS",
    "PATH",
    "PATHEXT",
    "PROCESSOR_ARCHITECTURE",
    "SystemDrive",
    "SystemRoot",
    "TEMP",
    "TMP",
    "WINDIR",
)
_HISTORICAL_RELEASES = {"chembl_24_1", "chembl_27", "chembl_31", "chembl_37"}
_PROJECTION_ROOT = Path("dataset/processed/a2s_historical_projection")
_PRESENCE_OUTPUT = Path(
    "dataset/processed/a2s_historical_presence/activity_presence.v1.parquet"
)
_H0_INVENTORY_OUTPUT = Path("dataset/processed/a2s_h0_metadata_inventory.v3.json")
_A2S_VALIDATION_OUTPUT = Path("dataset/processed/a2s_validation_small.v1")

_PROCESS_TERMINATE = 0x0001
_PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
_SYNCHRONIZE = 0x00100000
_STILL_ACTIVE = 259
_WAIT_OBJECT_0 = 0x00000000
_WAIT_TIMEOUT = 0x00000102
_ERROR_INVALID_PARAMETER = 87
_WINDOWS_KERNEL32: Any | None = None
_SEALED_RUN_FILES = (
    "intent.json",
    "prepared.json",
    "started.json",
    "recovery.json",
    "runner_error.log",
    "run.json",
    "events.jsonl",
    "stdout.log",
    "stderr.log",
    "network_requests.jsonl",
    "output_claim.json",
)


class DatasetProcessPolicyError(RuntimeError):
    """Raised when a command is outside the frozen metadata-processing policy."""


_SENSITIVE_ARGUMENTS = {
    "--anyauth",
    "--basic",
    "--config",
    "--cookie",
    "--cookie-jar",
    "--digest",
    "--header",
    "--netrc",
    "--netrc-file",
    "--oauth2-bearer",
    "--proxy-header",
    "--proxy-user",
    "--user",
    "-b",
    "-h",
    "-k",
    "-u",
}
_SENSITIVE_SHORT_PREFIXES = ("-b", "-h", "-k", "-u")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace(
        "+00:00", "Z"
    )


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def _json_hash(value: Any) -> str:
    return sha256(_canonical_bytes(value)).hexdigest()


def _file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json_once(path: Path, value: Any) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, ensure_ascii=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def _append_event(path: Path, event: str, **details: Any) -> None:
    record = {"utc": _utc_now(), "event": event, **details}
    line = _canonical_bytes(record) + b"\n"
    with path.open("ab") as stream:
        stream.write(line)
        stream.flush()
        os.fsync(stream.fileno())


def _resolve_path(path: str | Path, cwd: Path) -> Path:
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = cwd / candidate
    return candidate.resolve(strict=False)


def _snapshot_file(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "path": str(path),
        "kind": "file",
        "bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": _file_hash(path),
    }


def _snapshot_path(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"path": str(path), "kind": "missing", "exists": False}
    if path.is_file():
        return _snapshot_file(path)
    if not path.is_dir():
        raise ValueError(f"unsupported dataset path type: {path}")
    files = [_snapshot_file(item) for item in sorted(path.rglob("*")) if item.is_file()]
    return {
        "path": str(path),
        "kind": "directory",
        "files": len(files),
        "bytes": sum(item["bytes"] for item in files),
        "tree_sha256": _json_hash(files),
        "entries": files,
    }


def _snapshots(paths: Iterable[str | Path], cwd: Path) -> list[dict[str, Any]]:
    resolved = sorted({_resolve_path(path, cwd) for path in paths}, key=str)
    return [_snapshot_path(path) for path in resolved]


def _bundle(snapshot: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    return {"files": list(snapshot), "bundle_sha256": _json_hash(snapshot)}


def _snapshot_identity(record: Mapping[str, Any]) -> tuple[Any, ...]:
    return (
        record.get("kind"),
        record.get("bytes"),
        record.get("sha256"),
        record.get("tree_sha256"),
    )


def _controlled_environment(overrides: Mapping[str, str] | None = None) -> dict[str, str]:
    environment = {
        key: os.environ[key]
        for key in _ENVIRONMENT_KEYS
        if key in os.environ
    }
    environment.update({"PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1"})
    if overrides:
        raise DatasetProcessPolicyError(
            "environment overrides are forbidden for dataset-processing commands"
        )
    return dict(sorted(environment.items(), key=lambda item: item[0].lower()))


def _git_output(cwd: Path, arguments: Sequence[str]) -> bytes | None:
    result = subprocess.run(
        ["git", *arguments],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.stdout if result.returncode == 0 else None


def _git_state(cwd: Path) -> dict[str, Any]:
    commit = _git_output(cwd, ["rev-parse", "HEAD"])
    diff = _git_output(cwd, ["diff", "--binary", "--no-ext-diff", "HEAD"])
    status = _git_output(cwd, ["status", "--porcelain=v1", "--untracked-files=all"])
    return {
        "commit": commit.decode("ascii").strip() if commit else None,
        "tracked_diff_sha256": sha256(diff).hexdigest() if diff is not None else None,
        "tracked_diff_bytes": len(diff) if diff is not None else None,
        "status_sha256": sha256(status).hexdigest() if status is not None else None,
        "status_entries": status.count(b"\n") if status is not None else None,
    }


def _resolve_command(command: Sequence[str], environment: Mapping[str, str]) -> list[str]:
    if not command:
        raise ValueError("dataset operation command is empty")
    executable = command[0]
    candidate = Path(executable)
    if candidate.is_absolute():
        resolved = candidate.resolve(strict=True)
    else:
        located = shutil.which(executable, path=environment.get("PATH"))
        if located is None:
            raise FileNotFoundError(f"cannot resolve executable: {executable}")
        resolved = Path(located).resolve(strict=True)
    if not resolved.is_file():
        raise FileNotFoundError(f"executable is not a file: {resolved}")
    return [str(resolved), *command[1:]]


def _sensitive_argument(argument: str) -> tuple[str, bool] | None:
    lowered = argument.lower()
    flag = lowered.split("=", 1)[0]
    if flag in _SENSITIVE_ARGUMENTS:
        has_value = "=" in argument or flag not in {
            "--anyauth",
            "--basic",
            "--digest",
            "--netrc",
        }
        return flag, has_value
    for prefix in _SENSITIVE_SHORT_PREFIXES:
        if lowered.startswith(prefix) and len(argument) > len(prefix):
            return prefix, True
    return None


def _redacted_command(command: Sequence[str]) -> list[str]:
    redacted: list[str] = []
    hide_next = False
    for argument in command:
        lowered = argument.lower()
        if hide_next:
            redacted.append("<redacted>")
            hide_next = False
            continue
        sensitive = _sensitive_argument(argument)
        if sensitive is not None:
            flag, has_value = sensitive
            redacted.append(flag)
            if "=" in argument or len(argument) > len(flag):
                redacted.append("<redacted>")
            else:
                hide_next = has_value
            continue
        if lowered.startswith(("http://", "https://")):
            parsed = urlsplit(argument)
            if parsed.username or parsed.password or parsed.query or parsed.fragment:
                redacted.append("<redacted>")
                continue
        redacted.append(argument)
    return redacted


def _require_exact_outputs(
    actual: Sequence[Path], expected: Sequence[Path], operation: str
) -> None:
    actual_strings = [str(path) for path in actual]
    expected_strings = [str(path) for path in expected]
    if len(actual_strings) != len(expected_strings) or set(actual_strings) != set(
        expected_strings
    ):
        raise DatasetProcessPolicyError(
            f"{operation} must declare its exact committed output set"
        )


def _projection_outputs(cwd: Path, release: str) -> tuple[Path, Path, Path, Path]:
    output = _resolve_path(
        _PROJECTION_ROOT / release / "activity_identity.parquet", cwd
    )
    return (
        output,
        output.with_name(f"{output.stem}.documents.parquet"),
        output.with_suffix(".manifest.json"),
        output.with_suffix(".replay.json"),
    )


def _presence_outputs(cwd: Path) -> tuple[Path, Path, Path]:
    output = _resolve_path(_PRESENCE_OUTPUT, cwd)
    return (
        output,
        output.with_suffix(".manifest.json"),
        output.with_suffix(".replay.json"),
    )


def _a2s_validation_outputs(cwd: Path) -> list[Path]:
    root = _resolve_path(_A2S_VALIDATION_OUTPUT, cwd)
    return [root / name for name in ("rows.parquet", "ligand_features.npy", "target_embeddings.npz", "episodes.parquet", "leakage_report.json", "normalization.json", "manifest.json")]


def _windows_kernel32() -> Any:
    global _WINDOWS_KERNEL32
    if _WINDOWS_KERNEL32 is not None:
        return _WINDOWS_KERNEL32
    if os.name != "nt":  # pragma: no cover - Windows is the project platform.
        raise OSError("Windows process APIs are unavailable on this platform")
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    filetime_pointer = ctypes.POINTER(wintypes.FILETIME)
    kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.GetProcessTimes.argtypes = [
        wintypes.HANDLE,
        filetime_pointer,
        filetime_pointer,
        filetime_pointer,
        filetime_pointer,
    ]
    kernel32.GetProcessTimes.restype = wintypes.BOOL
    kernel32.GetExitCodeProcess.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(wintypes.DWORD),
    ]
    kernel32.GetExitCodeProcess.restype = wintypes.BOOL
    kernel32.QueryFullProcessImageNameW.argtypes = [
        wintypes.HANDLE,
        wintypes.DWORD,
        wintypes.LPWSTR,
        ctypes.POINTER(wintypes.DWORD),
    ]
    kernel32.QueryFullProcessImageNameW.restype = wintypes.BOOL
    kernel32.TerminateProcess.argtypes = [wintypes.HANDLE, wintypes.UINT]
    kernel32.TerminateProcess.restype = wintypes.BOOL
    kernel32.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    kernel32.WaitForSingleObject.restype = wintypes.DWORD
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    kernel32.CloseHandle.restype = wintypes.BOOL
    _WINDOWS_KERNEL32 = kernel32
    return kernel32


def _windows_error(message: str) -> OSError:
    error = ctypes.get_last_error()
    return OSError(error, f"{message}: {ctypes.FormatError(error)}")


def _open_windows_process(pid: int, access: int) -> tuple[Any, Any | None, int]:
    kernel32 = _windows_kernel32()
    ctypes.set_last_error(0)
    handle = kernel32.OpenProcess(access, False, pid)
    return kernel32, handle or None, ctypes.get_last_error()


def _windows_process_details(kernel32: Any, handle: Any, pid: int) -> dict[str, Any]:
    creation = wintypes.FILETIME()
    exit_time = wintypes.FILETIME()
    kernel_time = wintypes.FILETIME()
    user_time = wintypes.FILETIME()
    if not kernel32.GetProcessTimes(
        handle,
        ctypes.byref(creation),
        ctypes.byref(exit_time),
        ctypes.byref(kernel_time),
        ctypes.byref(user_time),
    ):
        raise _windows_error("GetProcessTimes failed")
    image_buffer = ctypes.create_unicode_buffer(32768)
    image_length = wintypes.DWORD(len(image_buffer))
    if not kernel32.QueryFullProcessImageNameW(
        handle, 0, image_buffer, ctypes.byref(image_length)
    ):
        raise _windows_error("QueryFullProcessImageNameW failed")
    exit_code = wintypes.DWORD()
    if not kernel32.GetExitCodeProcess(handle, ctypes.byref(exit_code)):
        raise _windows_error("GetExitCodeProcess failed")
    return {
        "pid": int(pid),
        "creation_time_filetime": (
            int(creation.dwHighDateTime) << 32
        )
        | int(creation.dwLowDateTime),
        "executable_image_path": image_buffer.value,
        "running": exit_code.value == _STILL_ACTIVE,
        "observed_exit_code": int(exit_code.value),
    }


def _normalize_windows_image(path: str) -> str:
    normalized = path.strip().replace("/", "\\")
    if normalized.startswith("\\\\?\\"):
        normalized = normalized[4:]
    return os.path.normpath(normalized).casefold()


def _started_identity_complete(started: Mapping[str, Any]) -> bool:
    return (
        isinstance(started.get("pid"), int)
        and started["pid"] > 0
        and isinstance(started.get("creation_time_filetime"), int)
        and started["creation_time_filetime"] > 0
        and isinstance(started.get("executable_image_path"), str)
        and bool(started["executable_image_path"])
        and isinstance(started.get("command_sha256"), str)
        and re.fullmatch(r"[0-9a-f]{64}", started["command_sha256"]) is not None
    )


def _process_identity_matches(
    expected: Mapping[str, Any], observed: Mapping[str, Any]
) -> bool:
    return (
        int(expected["pid"]) == int(observed["pid"])
        and int(expected["creation_time_filetime"])
        == int(observed["creation_time_filetime"])
        and _normalize_windows_image(str(expected["executable_image_path"]))
        == _normalize_windows_image(str(observed["executable_image_path"]))
    )


def _capture_started_identity(
    process: subprocess.Popen[Any], command: Sequence[str]
) -> dict[str, Any]:
    if os.name == "nt":
        kernel32 = _windows_kernel32()
        details = _windows_process_details(kernel32, process._handle, process.pid)
        return {
            "pid": process.pid,
            "creation_time_filetime": details["creation_time_filetime"],
            "executable_image_path": details["executable_image_path"],
            "command_sha256": _json_hash(list(command)),
        }
    return {  # pragma: no cover - Windows is the project platform.
        "pid": process.pid,
        "creation_time_filetime": None,
        "executable_image_path": str(Path(command[0]).resolve(strict=False)),
        "command_sha256": _json_hash(list(command)),
    }


def _inspect_started_process(
    started: Mapping[str, Any], *, access: int = _PROCESS_QUERY_LIMITED_INFORMATION
) -> tuple[str, dict[str, Any] | None]:
    if not _started_identity_complete(started):
        return "IDENTITY_UNKNOWN", None
    if os.name != "nt":  # pragma: no cover - Windows is the project platform.
        return "IDENTITY_UNKNOWN", None
    pid = int(started["pid"])
    kernel32, handle, error = _open_windows_process(pid, access)
    if handle is None:
        if error == _ERROR_INVALID_PARAMETER:
            return "NOT_RUNNING", None
        return "IDENTITY_UNKNOWN", None
    try:
        observed = _windows_process_details(kernel32, handle, pid)
        if not _process_identity_matches(started, observed):
            return "MISMATCH", observed
        if not observed["running"]:
            return "NOT_RUNNING", observed
        return "MATCH", observed
    except OSError:
        return "IDENTITY_UNKNOWN", None
    finally:
        kernel32.CloseHandle(handle)


def _validate_command_policy(
    command: Sequence[str],
    cwd: Path,
    outputs: Sequence[Path] = (),
    expected_bytes: Mapping[str, int] | None = None,
    expected_hashes: Mapping[str, str] | None = None,
    publish_outputs: Mapping[str, str] | None = None,
) -> None:
    executable = Path(command[0]).resolve(strict=True)
    official_curl = Path(r"C:\Windows\System32\curl.exe").resolve(strict=True)
    if executable == official_curl:
        lowered = [argument.lower() for argument in command[1:]]
        if any(_sensitive_argument(argument) is not None for argument in command[1:]):
            raise DatasetProcessPolicyError("curl credentials/configuration are forbidden")
        urls = [
            argument
            for argument in command[1:]
            if argument.lower().startswith(("http://", "https://"))
        ]
        if len(urls) != 1:
            raise DatasetProcessPolicyError("curl requires one official source URL")
        parsed = urlsplit(urls[0])
        if (
            parsed.scheme != "https"
            or parsed.hostname != "ftp.ebi.ac.uk"
            or not parsed.path.startswith(
                "/pub/databases/chembl/ChEMBLdb/releases/"
            )
            or parsed.username
            or parsed.password
            or parsed.query
            or parsed.fragment
        ):
            raise DatasetProcessPolicyError(
                "curl is limited to official ChEMBL release downloads"
            )
        required_flags = {"--fail", "--location", "--continue-at", "--output"}
        if not required_flags.issubset(lowered):
            raise DatasetProcessPolicyError(
                "curl requires --fail --location --continue-at - --output"
            )
        continue_index = lowered.index("--continue-at")
        output_index = lowered.index("--output")
        if continue_index + 1 >= len(command) or command[continue_index + 2] != "-":
            raise DatasetProcessPolicyError("curl resume value must be '-'")
        if output_index + 1 >= len(command):
            raise DatasetProcessPolicyError("curl output path is missing")
        output = _resolve_path(command[output_index + 2], cwd)
        if list(outputs) != [output] or not output.name.endswith(".partial"):
            raise DatasetProcessPolicyError(
                "curl must write the single declared .partial staging output"
            )
        if set(expected_bytes or {}) != {str(output)} or set(expected_hashes or {}) != {
            str(output)
        }:
            raise DatasetProcessPolicyError(
                "curl requires the official expected byte count and SHA-256"
            )
        publish = publish_outputs or {}
        if set(publish) != {str(output)}:
            raise DatasetProcessPolicyError(
                "curl requires one recorded staging-to-canonical publish mapping"
            )
        destination = Path(publish[str(output)])
        if output.name != f"{destination.name}.partial":
            raise DatasetProcessPolicyError(
                "curl staging name must be the canonical filename plus .partial"
            )
        allowed_values = {"--fail", "--location", "--continue-at", "-", "--output", str(output), urls[0]}
        original_output = command[output_index + 2]
        allowed_values.add(original_output)
        for flag, lower, upper in (
            ("--retry", 1, 10),
            ("--retry-delay", 0, 60),
            ("--speed-limit", 1024, 100000000),
            ("--speed-time", 10, 600),
        ):
            if flag in lowered:
                index = lowered.index(flag)
                if index + 2 >= len(command):
                    raise DatasetProcessPolicyError(f"curl {flag} value is missing")
                value = command[index + 2]
                if not value.isdigit() or not lower <= int(value) <= upper:
                    raise DatasetProcessPolicyError(f"curl {flag} value is outside policy")
                allowed_values.update({flag, value})
        if "--retry-all-errors" in lowered:
            allowed_values.add("--retry-all-errors")
        unexpected = [argument for argument in command[1:] if argument not in allowed_values]
        if unexpected:
            raise DatasetProcessPolicyError(f"unsupported curl arguments: {unexpected}")
        return
    if executable != DRUG_PYTHON.resolve(strict=True):
        raise DatasetProcessPolicyError(
            "Python dataset operations must use D:\\anaconda\\envs\\drug\\python.exe"
        )
    if len(command) < 3 or _resolve_path(command[1], cwd) != (cwd / "main.py"):
        raise DatasetProcessPolicyError("Python must invoke the workspace main.py")
    allowed = {"historical-project", "historical-presence", "h0-inventory", "prepare-a2s-validation"}
    if command[2] not in allowed:
        raise DatasetProcessPolicyError(
            f"command {command[2]!r} is not an allowed metadata operation"
        )
    if not outputs:
        raise DatasetProcessPolicyError(
            "metadata transformations must declare every committed output"
        )
    if publish_outputs:
        raise DatasetProcessPolicyError(
            "Python metadata transformations may not use runner publish mappings"
        )

    arguments = tuple(command[3:])
    operation = command[2]
    if operation == "historical-project":
        primary = (
            len(arguments) == 2
            and arguments[0] == "--release"
            and arguments[1] in _HISTORICAL_RELEASES
        )
        replay = (
            len(arguments) == 3
            and arguments[0] == "--release"
            and arguments[1] in _HISTORICAL_RELEASES
            and arguments[2] == "--replay"
        )
        if not primary and not replay:
            raise DatasetProcessPolicyError(
                "historical-project requires exactly --release RELEASE and optional --replay"
            )
        projection, documents, manifest, certificate = _projection_outputs(
            cwd, arguments[1]
        )
        if replay:
            missing = [
                str(path) for path in (projection, documents, manifest) if not path.is_file()
            ]
            if missing:
                raise DatasetProcessPolicyError(
                    "historical-project replay requires all frozen primary artifacts"
                )
            _require_exact_outputs(outputs, [certificate], operation)
        else:
            _require_exact_outputs(outputs, [projection, documents, manifest], operation)
        return

    if operation == "historical-presence":
        if arguments not in {(), ("--replay",)}:
            raise DatasetProcessPolicyError(
                "historical-presence permits only the primary or --replay invocation"
            )
        output, manifest, certificate = _presence_outputs(cwd)
        if arguments == ("--replay",):
            if not output.is_file() or not manifest.is_file():
                raise DatasetProcessPolicyError(
                    "historical-presence replay requires the frozen primary artifacts"
                )
            _require_exact_outputs(outputs, [certificate], operation)
        else:
            _require_exact_outputs(outputs, [output, manifest], operation)
        return

    if operation == "prepare-a2s-validation":
        if arguments not in {(), ("--output", str(_resolve_path(_A2S_VALIDATION_OUTPUT, cwd)))}:
            raise DatasetProcessPolicyError("prepare-a2s-validation permits only the default output")
        _require_exact_outputs(outputs, _a2s_validation_outputs(cwd), operation)
        return

    inventory_output = _resolve_path(_H0_INVENTORY_OUTPUT, cwd)
    explicit_output = (
        len(arguments) == 2
        and arguments[0] == "--output"
        and _resolve_path(arguments[1], cwd) == inventory_output
    )
    if arguments and not explicit_output:
        raise DatasetProcessPolicyError(
            "h0-inventory permits only the local v3 output and no remote acquisition"
        )
    _require_exact_outputs(outputs, [inventory_output], operation)


def _automatic_bundles(command: Sequence[str], cwd: Path) -> tuple[list[Path], list[Path]]:
    code = [cwd / "main.py"]
    module_by_command = {
        "historical-project": cwd / "research" / "a2s_historical_projection.py",
        "historical-presence": cwd / "research" / "a2s_historical_presence.py",
        "h0-inventory": cwd / "research" / "a2s_h0_metadata_inventory.py",
        "prepare-a2s-validation": cwd / "scripts" / "preprocess.py",
    }
    if len(command) >= 3 and command[2] in module_by_command:
        code.append(module_by_command[command[2]])
        if command[2] in {"historical-presence", "h0-inventory"}:
            code.append(cwd / "research" / "a2s_historical_projection.py")
    config = [cwd / "requirements.txt", cwd / "task.md"]
    conda_history = DRUG_PYTHON.parent / "conda-meta" / "history"
    if conda_history.exists():
        config.append(conda_history)
    return [path for path in code if path.exists()], [path for path in config if path.exists()]


def _new_run(root: Path, operation: str) -> tuple[str, Path]:
    if not _OPERATION.fullmatch(operation):
        raise ValueError("operation must be a 2-80 character lowercase slug")
    root.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    run_id = f"{stamp}_{operation}_{uuid4().hex[:8]}"
    run_dir = root / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_id, run_dir


def require_run_context(argv: Sequence[str], cwd: str | Path = ".") -> None:
    """Fail closed when a protected dataset command bypasses this runner."""
    working_directory = Path(cwd).resolve(strict=True)
    run_dir_value = os.environ.get("A2S_DATASET_RUN_DIR")
    run_id = os.environ.get("A2S_DATASET_RUN_ID")
    if not run_dir_value or not run_id:
        raise SystemExit(
            "dataset command blocked: execute it through `main.py dataset-run execute`"
        )
    run_dir = Path(run_dir_value).resolve(strict=True)
    prepared_path = run_dir / "prepared.json"
    intent_path = run_dir / "intent.json"
    if not prepared_path.is_file() or not intent_path.is_file():
        raise SystemExit("dataset command blocked: run context is incomplete")
    prepared = json.loads(prepared_path.read_text(encoding="utf-8"))
    intent = json.loads(intent_path.read_text(encoding="utf-8"))
    actual = [str(_resolve_path(argv[0], working_directory)), *argv[1:]]
    expected = list(prepared.get("command_argv", [])[1:])
    if expected:
        expected[0] = str(_resolve_path(expected[0], working_directory))
    valid = (
        run_dir.name == run_id
        and run_dir.parent.name == "runs"
        and prepared.get("run_id") == run_id
        and intent.get("run_id") == run_id
        and prepared.get("record_type") == "contemporaneous_execution"
        and prepared.get("cwd") == str(working_directory)
        and actual == expected
    )
    if not valid:
        raise SystemExit("dataset command blocked: run context does not match argv/cwd")


def _content_record(payload: dict[str, Any]) -> dict[str, Any]:
    payload = dict(payload)
    payload["content_sha256"] = _json_hash(payload)
    return payload


def _claim_path_key(path: Path) -> str:
    """Return the case-insensitive canonical key used for output ownership."""
    return os.path.normcase(str(path.resolve(strict=False)))


def _current_process_identity() -> dict[str, Any]:
    """Capture the runner identity used to recover claims after a crash."""
    command = [str(Path(sys.executable).resolve(strict=False)), *sys.argv]
    if os.name == "nt":
        kernel32, handle, error = _open_windows_process(
            os.getpid(), _PROCESS_QUERY_LIMITED_INFORMATION
        )
        if handle is None:
            raise DatasetProcessPolicyError(
                f"cannot identify dataset runner process (Windows error {error})"
            )
        try:
            details = _windows_process_details(kernel32, handle, os.getpid())
        finally:
            kernel32.CloseHandle(handle)
        return {
            "pid": os.getpid(),
            "creation_time_filetime": details["creation_time_filetime"],
            "executable_image_path": details["executable_image_path"],
            "command_sha256": _json_hash(command),
        }
    return {  # pragma: no cover - Windows is the project platform.
        "pid": os.getpid(),
        "creation_time_filetime": None,
        "executable_image_path": command[0],
        "command_sha256": _json_hash(command),
    }


def _claim_outputs(
    *,
    root: Path,
    run_dir: Path,
    run_id: str,
    operation: str,
    outputs: Sequence[Path],
    publish_map: Mapping[str, str],
    runner_identity: Mapping[str, Any],
) -> dict[str, Any]:
    roles: dict[str, set[str]] = {}
    paths: dict[str, Path] = {}

    def add(path: Path, role: str) -> None:
        resolved = path.resolve(strict=False)
        key = _claim_path_key(resolved)
        paths[key] = resolved
        roles.setdefault(key, set()).add(role)

    for path in outputs:
        add(path, "declared_output")
    for source, destination in publish_map.items():
        add(Path(source), "publish_staging_source")
        add(Path(destination), "publish_destination")

    claim_paths = [
        {
            "path": str(paths[key]),
            "normalized_path": key,
            "roles": sorted(roles[key]),
        }
        for key in sorted(paths)
    ]
    ledger = root / "ledger.jsonl"
    with _ledger_lock(root):
        indexed = {entry["run_id"] for entry in _validated_ledger_entries(ledger)}
        for other_dir in sorted((root / "runs").glob("*")):
            if not other_dir.is_dir() or other_dir.name in indexed:
                continue
            other_claim_path = other_dir / "output_claim.json"
            if not other_claim_path.is_file():
                continue
            try:
                other_claim = json.loads(other_claim_path.read_text(encoding="utf-8"))
                unsigned = {
                    key: value
                    for key, value in other_claim.items()
                    if key != "content_sha256"
                }
                if (
                    other_claim.get("content_sha256") != _json_hash(unsigned)
                    or other_claim.get("run_id") != other_dir.name
                ):
                    raise ValueError("invalid output claim content hash or identity")
                other_paths = {
                    str(item["normalized_path"]): str(item["path"])
                    for item in other_claim.get("paths", [])
                }
            except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as caught:
                raise DatasetProcessPolicyError(
                    f"cannot validate unindexed output claim {other_dir.name}: {caught}"
                ) from caught
            overlap = sorted(set(other_paths).intersection(paths))
            if overlap:
                raise DatasetProcessPolicyError(
                    f"output already reserved by {other_dir.name}: {other_paths[overlap[0]]}"
                )

        claim = _content_record(
            {
                "schema": "a2s-dataset-processing-output-claim-v1",
                "run_id": run_id,
                "operation": operation,
                "acquired_utc": _utc_now(),
                "runner_identity": dict(runner_identity),
                "paths": claim_paths,
            }
        )
        _write_json_once(run_dir / "output_claim.json", claim)
        return claim


@contextmanager
def _ledger_lock(root: Path) -> Iterable[None]:
    lock_path = root / ".ledger.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+b") as stream:
        stream.seek(0, os.SEEK_END)
        if stream.tell() == 0:
            stream.write(b"0")
            stream.flush()
        stream.seek(0)
        deadline = time.monotonic() + 30
        if os.name == "nt":
            import msvcrt

            while True:
                try:
                    msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError:
                    if time.monotonic() >= deadline:
                        raise TimeoutError("timed out acquiring processing ledger lock")
                    time.sleep(0.05)
            try:
                yield
            finally:
                stream.seek(0)
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            return
        import fcntl  # pragma: no cover - Windows is the project platform.

        fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _validated_ledger_entries(ledger: Path) -> list[dict[str, Any]]:
    if not ledger.exists():
        return []
    entries: list[dict[str, Any]] = []
    previous: str | None = None
    for index, line in enumerate(ledger.read_bytes().splitlines()):
        if not line.strip():
            continue
        entry = json.loads(line)
        digest = entry.get("entry_sha256")
        unsigned = {key: value for key, value in entry.items() if key != "entry_sha256"}
        if digest != _json_hash(unsigned) or entry.get("previous_entry_sha256") != previous:
            raise ValueError(f"processing ledger chain failure at entry {index}")
        entries.append(entry)
        previous = digest
    return entries


def _seal_and_index(
    root: Path,
    run_dir: Path,
    run_id: str,
    operation: str,
    status: str,
) -> None:
    names = [
        name
        for name in _SEALED_RUN_FILES
        if (run_dir / name).is_file()
    ]
    seal = _content_record(
        {
            "schema": "a2s-dataset-processing-seal-v1",
            "run_id": run_id,
            "files": {name: _file_hash(run_dir / name) for name in names},
        }
    )
    _write_json_once(run_dir / "seal.json", seal)

    _index_sealed_run(root, run_dir, run_id, operation, status)


def _index_sealed_run(
    root: Path,
    run_dir: Path,
    run_id: str,
    operation: str,
    status: str,
) -> None:
    if not (run_dir / "run.json").is_file() or not (run_dir / "seal.json").is_file():
        raise ValueError("run.json and seal.json are required before ledger indexing")

    ledger = root / "ledger.jsonl"
    with _ledger_lock(root):
        entries = _validated_ledger_entries(ledger)
        existing = [entry for entry in entries if entry.get("run_id") == run_id]
        if existing:
            if len(existing) != 1:
                raise ValueError(f"duplicate processing ledger entries for {run_id}")
            entry = existing[0]
            expected = {
                "operation": operation,
                "status": status,
                "run_path": str(run_dir),
                "run_json_sha256": _file_hash(run_dir / "run.json"),
                "seal_sha256": _file_hash(run_dir / "seal.json"),
            }
            if any(entry.get(key) != value for key, value in expected.items()):
                raise ValueError(f"conflicting existing ledger entry for {run_id}")
            return
        previous = entries[-1]["entry_sha256"] if entries else None
        entry = {
            "schema": "a2s-dataset-processing-ledger-entry-v1",
            "indexed_utc": _utc_now(),
            "run_id": run_id,
            "operation": operation,
            "status": status,
            "run_path": str(run_dir),
            "run_json_sha256": _file_hash(run_dir / "run.json"),
            "seal_sha256": _file_hash(run_dir / "seal.json"),
            "previous_entry_sha256": previous,
        }
        entry["entry_sha256"] = _json_hash(entry)
        with ledger.open("ab") as stream:
            stream.write(_canonical_bytes(entry) + b"\n")
            stream.flush()
            os.fsync(stream.fileno())


def execute_run(
    *,
    operation: str,
    command: Sequence[str],
    cwd: str | Path = ".",
    root: str | Path = DEFAULT_ROOT,
    inputs: Sequence[str | Path] = (),
    outputs: Sequence[str | Path] = (),
    code_inputs: Sequence[str | Path] = (),
    config_inputs: Sequence[str | Path] = (),
    expected_output_bytes: Mapping[str | Path, int] | None = None,
    expected_output_sha256: Mapping[str | Path, str] | None = None,
    publish_outputs: Mapping[str | Path, str | Path] | None = None,
    parent_run_id: str | None = None,
    environment_overrides: Mapping[str, str] | None = None,
    _policy_validator: Any = _validate_command_policy,
) -> dict[str, Any]:
    working_directory = Path(cwd).resolve(strict=True)
    audit_root = _resolve_path(root, working_directory)
    run_id, run_dir = _new_run(audit_root, operation)
    events = run_dir / "events.jsonl"
    (run_dir / "stdout.log").touch(exist_ok=False)
    (run_dir / "stderr.log").touch(exist_ok=False)
    intent = {
        "schema": "a2s-dataset-processing-intent-v1",
        "run_id": run_id,
        "record_type": "contemporaneous_execution",
        "operation": operation,
        "created_utc": _utc_now(),
        "cwd": str(working_directory),
        "requested_command_redacted": _redacted_command(command),
        "declared_inputs": [str(path) for path in inputs],
        "declared_outputs": [str(path) for path in outputs],
        "declared_code_inputs": [str(path) for path in code_inputs],
        "declared_config_inputs": [str(path) for path in config_inputs],
        "expected_output_bytes": {
            str(path): value for path, value in (expected_output_bytes or {}).items()
        },
        "expected_output_sha256": {
            str(path): value for path, value in (expected_output_sha256 or {}).items()
        },
        "publish_outputs": {
            str(source): str(destination)
            for source, destination in (publish_outputs or {}).items()
        },
        "parent_run_id": parent_run_id,
    }
    _write_json_once(run_dir / "intent.json", intent)
    _append_event(events, "intent_frozen")
    prepared: dict[str, Any] | None = None
    child_started: str | None = None
    child_ended: str | None = None
    pid: int | None = None
    exit_code: int | None = None
    error: dict[str, Any] | None = None
    output_claim: dict[str, Any] | None = None
    effective_code_inputs = list(map(Path, code_inputs))
    effective_config_inputs = list(map(Path, config_inputs))
    publish_map: dict[str, str] = {}
    start = time.monotonic()
    try:
        environment = _controlled_environment(environment_overrides)
        environment["A2S_DATASET_RUN_DIR"] = str(run_dir)
        environment["A2S_DATASET_RUN_ID"] = run_id
        resolved_command = _resolve_command(command, environment)
        output_paths = [_resolve_path(path, working_directory) for path in outputs]
        expected_bytes = {
            str(_resolve_path(path, working_directory)): int(value)
            for path, value in (expected_output_bytes or {}).items()
        }
        expected_hashes = {
            str(_resolve_path(path, working_directory)): str(value).lower()
            for path, value in (expected_output_sha256 or {}).items()
        }
        publish_map = {
            str(_resolve_path(source, working_directory)): str(
                _resolve_path(destination, working_directory)
            )
            for source, destination in (publish_outputs or {}).items()
        }
        _policy_validator(
            resolved_command,
            working_directory,
            output_paths,
            expected_bytes,
            expected_hashes,
            publish_map,
        )
        runner_identity = _current_process_identity()
        output_claim = _claim_outputs(
            root=audit_root,
            run_dir=run_dir,
            run_id=run_id,
            operation=operation,
            outputs=output_paths,
            publish_map=publish_map,
            runner_identity=runner_identity,
        )
        _append_event(
            events,
            "output_claim_acquired",
            output_claim_sha256=output_claim["content_sha256"],
        )
        automatic_code, automatic_config = _automatic_bundles(
            resolved_command, working_directory
        )
        effective_code_inputs = [*automatic_code, *map(Path, code_inputs)]
        effective_config_inputs = [*automatic_config, *map(Path, config_inputs)]
        _append_event(events, "input_hashing_started")
        before_state = _snapshots(
            [*inputs, *outputs, *publish_map.values()], working_directory
        )
        before_by_path = {item["path"]: item for item in before_state}
        input_paths = [_resolve_path(path, working_directory) for path in inputs]
        input_state = [before_by_path[str(path)] for path in input_paths]
        output_state_before = [before_by_path[str(path)] for path in output_paths]
        publish_targets_before = [
            before_by_path[path] for path in publish_map.values()
        ]
        existing_publish_targets = [
            item["path"]
            for item in publish_targets_before
            if item["kind"] != "missing"
        ]
        if existing_publish_targets:
            raise DatasetProcessPolicyError(
                f"refusing to replace canonical outputs: {existing_publish_targets}"
            )
        input_path_set = set(input_paths)
        preexisting_outputs = [
            item["path"]
            for path, item in zip(output_paths, output_state_before)
            if item["kind"] != "missing" and path not in input_path_set
        ]
        if preexisting_outputs:
            raise DatasetProcessPolicyError(
                f"refusing undeclared overwrite of outputs: {preexisting_outputs}"
            )
        resumed_outputs = [
            item["path"]
            for path, item in zip(output_paths, output_state_before)
            if item["kind"] != "missing" and path in input_path_set
        ]
        if resumed_outputs and not parent_run_id:
            raise DatasetProcessPolicyError(
                "resumed outputs require parent_run_id/resume lineage"
            )
        if parent_run_id and parent_run_id not in _indexed_run_ids(
            audit_root / "ledger.jsonl"
        ):
            raise DatasetProcessPolicyError(
                f"parent_run_id is not present in the processing ledger: {parent_run_id}"
            )
        code_state = _snapshots(
            [Path(__file__), *effective_code_inputs], working_directory
        )
        config_state = _snapshots(effective_config_inputs, working_directory)
        executable_state = _snapshot_file(Path(resolved_command[0]))
        unknown_expectations = set(expected_bytes).union(expected_hashes).difference(
            str(path) for path in output_paths
        )
        if unknown_expectations:
            raise ValueError(
                f"expectations reference undeclared outputs: {sorted(unknown_expectations)}"
            )
        invalid_hashes = [
            value
            for value in expected_hashes.values()
            if not re.fullmatch(r"[0-9a-f]{64}", value)
        ]
        if invalid_hashes:
            raise ValueError("expected SHA-256 values must be 64 lowercase hex digits")
        prepared = _content_record(
            {
                "schema": "a2s-dataset-processing-prepared-v1",
                "run_id": run_id,
                "record_type": "contemporaneous_execution",
                "prepared_utc": _utc_now(),
                "operation": operation,
                "cwd": str(working_directory),
                "command_argv": resolved_command,
                "command_line_windows": subprocess.list2cmdline(resolved_command),
                "environment": environment,
                "environment_lock_sha256": _json_hash(environment),
                "runtime": {
                    "runner_python": sys.version,
                    "platform": platform.platform(),
                    "executable": executable_state,
                },
                "inputs_before": input_state,
                "outputs_before": output_state_before,
                "publish_targets_before": publish_targets_before,
                "publish_outputs": publish_map,
                "parent_run_id": parent_run_id,
                "expected_output_bytes": expected_bytes,
                "expected_output_sha256": expected_hashes,
                "code_bundle": _bundle(code_state),
                "config_bundle": _bundle(config_state),
                "git": _git_state(working_directory),
                "declared_outputs": [
                    str(_resolve_path(path, working_directory)) for path in outputs
                ],
                "output_claim_sha256": output_claim["content_sha256"],
            }
        )
        _write_json_once(run_dir / "prepared.json", prepared)
        _append_event(events, "prepared_frozen", prepared_sha256=prepared["content_sha256"])
        _append_event(events, "child_launch_attempted")
        with (run_dir / "stdout.log").open("ab") as stdout, (run_dir / "stderr.log").open("ab") as stderr:
            process = subprocess.Popen(
                resolved_command,
                cwd=working_directory,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                shell=False,
            )
            pid = process.pid
            child_started = _utc_now()
            try:
                process_identity = _capture_started_identity(process, resolved_command)
                _write_json_once(
                    run_dir / "started.json",
                    {
                        "schema": "a2s-dataset-processing-started-v2",
                        "run_id": run_id,
                        "started_utc": child_started,
                        **process_identity,
                    },
                )
                _append_event(events, "child_started", pid=pid)
                try:
                    exit_code = process.wait()
                except KeyboardInterrupt:
                    process.terminate()
                    try:
                        exit_code = process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        exit_code = process.wait()
                    child_ended = _utc_now()
                    raise
            except BaseException:
                if process.poll() is None:
                    process.terminate()
                    try:
                        exit_code = process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        exit_code = process.wait()
                child_ended = child_ended or _utc_now()
                raise
            finally:
                stdout.flush()
                stderr.flush()
                os.fsync(stdout.fileno())
                os.fsync(stderr.fileno())
        child_ended = child_ended or _utc_now()
        _append_event(events, "child_ended", pid=pid, exit_code=exit_code)
    except BaseException as caught:
        error = {
            "type": type(caught).__name__,
            "message": str(caught),
            "traceback": traceback.format_exc(),
        }
        with (run_dir / "runner_error.log").open(
            "x", encoding="utf-8", newline="\n"
        ) as stream:
            stream.write(error["traceback"])
        _append_event(events, "runner_failure", error_type=error["type"])
    duration = time.monotonic() - start
    try:
        output_state = _snapshots(outputs, working_directory)
        code_state_after = _snapshots(
            [Path(__file__), *effective_code_inputs], working_directory
        )
        config_state_after = _snapshots(effective_config_inputs, working_directory)
    except BaseException as caught:
        output_state = []
        code_state_after = []
        config_state_after = []
        if error is None:
            error = {
                "type": type(caught).__name__,
                "message": str(caught),
                "traceback": traceback.format_exc(),
            }
    outputs_complete = all(item.get("kind") != "missing" for item in output_state)
    expected_bytes = prepared.get("expected_output_bytes", {}) if prepared else {}
    expected_hashes = prepared.get("expected_output_sha256", {}) if prepared else {}
    validation_failures: list[str] = []
    for item in output_state:
        path = item["path"]
        if path in expected_bytes and item.get("bytes") != expected_bytes[path]:
            validation_failures.append(
                f"{path}: bytes {item.get('bytes')} != {expected_bytes[path]}"
            )
        if path in expected_hashes and item.get("sha256") != expected_hashes[path]:
            validation_failures.append(f"{path}: SHA-256 mismatch")
    if not outputs_complete:
        validation_failures.append("one or more declared outputs are missing")
    code_unchanged = prepared is not None and _bundle(code_state_after)[
        "bundle_sha256"
    ] == prepared["code_bundle"]["bundle_sha256"]
    config_unchanged = prepared is not None and _bundle(config_state_after)[
        "bundle_sha256"
    ] == prepared["config_bundle"]["bundle_sha256"]
    if prepared is not None and not code_unchanged:
        validation_failures.append("code bundle changed during execution")
    if prepared is not None and not config_unchanged:
        validation_failures.append("config bundle changed during execution")
    input_only_mutations: list[str] = []
    if prepared is not None:
        output_path_strings = {item["path"] for item in output_state}
        input_only_after = _snapshots(
            [item["path"] for item in prepared["inputs_before"] if item["path"] not in output_path_strings],
            working_directory,
        )
        after_by_path = {item["path"]: item for item in input_only_after}
        input_only_mutations = [
            item["path"]
            for item in prepared["inputs_before"]
            if item["path"] not in output_path_strings
            and _snapshot_identity(item)
            != _snapshot_identity(after_by_path[item["path"]])
        ]
        if input_only_mutations:
            validation_failures.append("input-only artifacts changed during execution")
    publish_status = "not_requested"
    publish_error: dict[str, str] | None = None
    published_outputs: list[dict[str, Any]] = []
    staging_after_publish: list[dict[str, Any]] = []
    if publish_map:
        publish_status = "not_attempted"
        if exit_code == 0 and error is None and not validation_failures:
            try:
                for source, destination in publish_map.items():
                    os.replace(source, destination)
                    _append_event(
                        events,
                        "validated_output_published",
                        source=source,
                        destination=destination,
                    )
                published_outputs = _snapshots(
                    publish_map.values(), working_directory
                )
                staging_after_publish = _snapshots(
                    publish_map.keys(), working_directory
                )
                publish_status = "published"
            except BaseException as caught:
                publish_status = "failed"
                publish_error = {
                    "type": type(caught).__name__,
                    "message": str(caught),
                }
                validation_failures.append("atomic canonical publish failed")
    if prepared is None:
        process_status = "preparation_failed"
    elif error is not None and error["type"] == "KeyboardInterrupt":
        process_status = "interrupted"
    elif pid is None:
        process_status = "spawn_failed"
    elif exit_code == 0:
        process_status = "completed"
    else:
        process_status = "nonzero_exit"
    validation_status = (
        "pass"
        if prepared is not None and not validation_failures
        else "fail"
        if prepared is not None
        else "not_applicable"
    )
    status = (
        "SUCCESS"
        if process_status == "completed" and validation_status == "pass" and error is None
        else "FAILED"
    )
    _append_event(events, "run_finalized", status=status)
    result = _content_record(
        {
            "schema": "a2s-dataset-processing-run-v1",
            "run_id": run_id,
            "record_type": "contemporaneous_execution",
            "operation": operation,
            "status": status,
            "process_status": process_status,
            "validation_status": validation_status,
            "validation_failures": validation_failures,
            "prepared_sha256": prepared["content_sha256"] if prepared else None,
            "child_started_utc": child_started,
            "child_ended_utc": child_ended,
            "duration_seconds": duration,
            "pid": pid,
            "exit_code": exit_code,
            "outputs_complete": outputs_complete,
            "validated_outputs_before_publish": output_state,
            "publish_status": publish_status,
            "publish_error": publish_error,
            "published_outputs_after": published_outputs,
            "staging_outputs_after_publish": staging_after_publish,
            "code_bundle_after": _bundle(code_state_after),
            "config_bundle_after": _bundle(config_state_after),
            "code_unchanged": code_unchanged,
            "config_unchanged": config_unchanged,
            "input_only_mutations": input_only_mutations,
            "stdout": _snapshot_file(run_dir / "stdout.log"),
            "stderr": _snapshot_file(run_dir / "stderr.log"),
            "error": error,
        }
    )
    _write_json_once(run_dir / "run.json", result)
    _seal_and_index(audit_root, run_dir, run_id, operation, status)
    result["run_path"] = str(run_dir)
    return result


def relocate_file_run(
    *,
    operation: str,
    source: str | Path,
    destination: str | Path,
    cwd: str | Path = ".",
    root: str | Path = DEFAULT_ROOT,
    expected_bytes: int | None = None,
) -> dict[str, Any]:
    """Audit and atomically relocate one pre-run legacy artifact without deletion."""
    working_directory = Path(cwd).resolve(strict=True)
    audit_root = _resolve_path(root, working_directory)
    source_path = _resolve_path(source, working_directory)
    destination_path = _resolve_path(destination, working_directory)
    archive_root = (
        working_directory / "dataset" / "public" / "chembl_historical" / "archives"
    ).resolve(strict=True)
    if source_path.parent != archive_root or destination_path.parent != archive_root:
        raise DatasetProcessPolicyError(
            "legacy relocation is restricted to the ChEMBL historical archive directory"
        )
    run_id, run_dir = _new_run(audit_root, operation)
    events = run_dir / "events.jsonl"
    intent = {
        "schema": "a2s-dataset-processing-intent-v1",
        "run_id": run_id,
        "record_type": "contemporaneous_internal_relocation",
        "operation": operation,
        "created_utc": _utc_now(),
        "cwd": str(working_directory),
        "source": str(source_path),
        "destination": str(destination_path),
        "expected_bytes": expected_bytes,
        "policy": "preserve legacy partial bytes; destination must not exist",
    }
    _write_json_once(run_dir / "intent.json", intent)
    _append_event(events, "relocation_intent_frozen")
    source_before = _snapshot_path(source_path)
    destination_before = _snapshot_path(destination_path)
    code_state = _snapshots([Path(__file__), working_directory / "main.py"], working_directory)
    prepared = _content_record(
        {
            "schema": "a2s-dataset-processing-prepared-v1",
            "run_id": run_id,
            "record_type": "contemporaneous_internal_relocation",
            "prepared_utc": _utc_now(),
            "operation": operation,
            "source_before": source_before,
            "destination_before": destination_before,
            "code_bundle": _bundle(code_state),
            "git": _git_state(working_directory),
        }
    )
    _write_json_once(run_dir / "prepared.json", prepared)
    _append_event(events, "relocation_prepared_frozen")
    start = time.monotonic()
    failure: str | None = None
    if source_before["kind"] != "file":
        failure = "source is not a file"
    elif destination_before["kind"] != "missing":
        failure = "destination already exists"
    elif expected_bytes is not None and source_before["bytes"] != expected_bytes:
        failure = (
            f"source bytes {source_before['bytes']} do not match expected {expected_bytes}"
        )
    if failure is None:
        try:
            os.replace(source_path, destination_path)
            _append_event(events, "legacy_artifact_relocated")
        except OSError as caught:
            failure = f"{type(caught).__name__}: {caught}"
    source_after = _snapshot_path(source_path)
    destination_after = _snapshot_path(destination_path)
    code_after = _snapshots([Path(__file__), working_directory / "main.py"], working_directory)
    if failure is None and (
        source_after["kind"] != "missing"
        or destination_after.get("sha256") != source_before.get("sha256")
        or _bundle(code_after)["bundle_sha256"] != prepared["code_bundle"]["bundle_sha256"]
    ):
        failure = "post-relocation hash or code validation failed"
    status = "SUCCESS" if failure is None else "FAILED"
    _append_event(events, "relocation_finalized", status=status)
    result = _content_record(
        {
            "schema": "a2s-dataset-processing-run-v1",
            "run_id": run_id,
            "record_type": "contemporaneous_internal_relocation",
            "operation": operation,
            "status": status,
            "process_status": "completed" if failure is None else "not_completed",
            "validation_status": "pass" if failure is None else "fail",
            "failure": failure,
            "duration_seconds": time.monotonic() - start,
            "source_after": source_after,
            "destination_after": destination_after,
            "code_bundle_after": _bundle(code_after),
            "stdout": {"available": False, "reason": "runner-internal atomic operation"},
            "stderr": {"available": False, "reason": "runner-internal atomic operation"},
        }
    )
    _write_json_once(run_dir / "run.json", result)
    _seal_and_index(audit_root, run_dir, run_id, operation, status)
    result["run_path"] = str(run_dir)
    return result


def terminate_active_run(
    *,
    target_run_id: str,
    reason: str,
    cwd: str | Path = ".",
    root: str | Path = DEFAULT_ROOT,
) -> dict[str, Any]:
    """Terminate a stalled logged child while preserving both audit records."""
    working_directory = Path(cwd).resolve(strict=True)
    audit_root = _resolve_path(root, working_directory)
    target_dir = audit_root / "runs" / target_run_id
    started_path = target_dir / "started.json"
    prepared_path = target_dir / "prepared.json"
    if not started_path.is_file() or not prepared_path.is_file():
        raise DatasetProcessPolicyError("target run has no frozen started/prepared record")
    if (target_dir / "run.json").exists():
        raise DatasetProcessPolicyError("target run is already finalized")
    started = json.loads(started_path.read_text(encoding="utf-8"))
    if not _started_identity_complete(started):
        raise DatasetProcessPolicyError(
            "IDENTITY_UNKNOWN: target started record lacks a complete process identity"
        )
    pid = int(started["pid"])
    if os.name != "nt":  # pragma: no cover - Windows is the project platform.
        raise DatasetProcessPolicyError(
            "IDENTITY_UNKNOWN: verified termination requires Windows process identity"
        )
    access = _PROCESS_QUERY_LIMITED_INFORMATION | _PROCESS_TERMINATE | _SYNCHRONIZE
    kernel32, handle, open_error = _open_windows_process(pid, access)
    if handle is None:
        if open_error == _ERROR_INVALID_PARAMETER:
            raise DatasetProcessPolicyError("target child process is not running")
        raise DatasetProcessPolicyError(
            "IDENTITY_UNKNOWN: target process identity could not be queried"
        )
    try:
        try:
            observed = _windows_process_details(kernel32, handle, pid)
        except OSError as caught:
            raise DatasetProcessPolicyError(
                "IDENTITY_UNKNOWN: target process identity could not be queried"
            ) from caught
        if not _process_identity_matches(started, observed):
            raise DatasetProcessPolicyError(
                "process identity mismatch; PID may have been reused"
            )
        if not observed["running"]:
            raise DatasetProcessPolicyError("target child process is not running")

        run_id, run_dir = _new_run(audit_root, f"terminate-{target_run_id[-8:]}")
        events = run_dir / "events.jsonl"
        intent = {
            "schema": "a2s-dataset-processing-intent-v1",
            "run_id": run_id,
            "record_type": "contemporaneous_control_action",
            "operation": "terminate-active-run",
            "created_utc": _utc_now(),
            "target_run_id": target_run_id,
            "target_pid": pid,
            "reason": reason,
        }
        _write_json_once(run_dir / "intent.json", intent)
        _append_event(events, "termination_intent_frozen")
        prepared = _content_record(
            {
                "schema": "a2s-dataset-processing-prepared-v1",
                "run_id": run_id,
                "record_type": "contemporaneous_control_action",
                "prepared_utc": _utc_now(),
                "target_started": started,
                "target_identity_observed": observed,
                "target_prepared_sha256": _file_hash(prepared_path),
                "runner_sha256": _file_hash(Path(__file__)),
            }
        )
        _write_json_once(run_dir / "prepared.json", prepared)
        _append_event(events, "termination_prepared_frozen")
        termination_error: str | None = None
        if not kernel32.TerminateProcess(handle, 1):
            termination_error = str(_windows_error("TerminateProcess failed"))
            stopped = False
        else:
            wait_result = int(kernel32.WaitForSingleObject(handle, 15_000))
            exit_code = wintypes.DWORD()
            exit_observed = bool(
                kernel32.GetExitCodeProcess(handle, ctypes.byref(exit_code))
            )
            stopped = (
                wait_result == _WAIT_OBJECT_0
                and exit_observed
                and exit_code.value != _STILL_ACTIVE
            )
            if not stopped:
                termination_error = (
                    "verified target did not stop before timeout"
                    if wait_result == _WAIT_TIMEOUT
                    else "verified target termination could not be confirmed"
                )
        status = "SUCCESS" if stopped else "FAILED"
        _append_event(events, "termination_finalized", status=status)
        result = _content_record(
            {
                "schema": "a2s-dataset-processing-run-v1",
                "run_id": run_id,
                "record_type": "contemporaneous_control_action",
                "operation": "terminate-active-run",
                "status": status,
                "process_status": "completed" if stopped else "target_still_running",
                "validation_status": "pass" if stopped else "fail",
                "target_run_id": target_run_id,
                "target_pid": pid,
                "target_identity_observed": observed,
                "reason": reason,
                "termination_error": termination_error,
                "exit_code": None,
                "note": "The target run owns and records its child exit code.",
            }
        )
        _write_json_once(run_dir / "run.json", result)
        _seal_and_index(audit_root, run_dir, run_id, "terminate-active-run", status)
        result["run_path"] = str(run_dir)
        return result
    finally:
        kernel32.CloseHandle(handle)


def record_retrospective(
    *,
    event_file: str | Path,
    event_id: str | None = None,
    cwd: str | Path = ".",
    root: str | Path = DEFAULT_ROOT,
    artifacts: Sequence[str | Path] = (),
) -> dict[str, Any]:
    working_directory = Path(cwd).resolve(strict=True)
    event_path = _resolve_path(event_file, working_directory)
    payload = json.loads(event_path.read_text(encoding="utf-8"))
    if payload.get("schema") == "a2s-dataset-retrospective-events-v1":
        matches = [
            event for event in payload.get("events", []) if event.get("event_id") == event_id
        ]
        if len(matches) != 1:
            raise ValueError(f"retrospective event_id must select one event: {event_id}")
        event = matches[0]
    else:
        event = payload
    operation = str(event.get("operation", ""))
    limitations = event.get("limitations")
    if event.get("schema") != "a2s-dataset-retrospective-event-v1":
        raise ValueError("retrospective event has the wrong schema")
    if not isinstance(limitations, list) or not limitations:
        raise ValueError("retrospective event must state at least one limitation")
    if not isinstance(event.get("evidence"), dict) or not event["evidence"]:
        raise ValueError("retrospective event must contain hashed or cited evidence")
    if not isinstance(event.get("unknown_fields"), list) or not event["unknown_fields"]:
        raise ValueError("retrospective event must enumerate unavailable fields")
    audit_root = _resolve_path(root, working_directory)
    run_id, run_dir = _new_run(audit_root, operation)
    events = run_dir / "events.jsonl"
    intent = {
        "schema": "a2s-dataset-processing-intent-v1",
        "run_id": run_id,
        "record_type": "retrospective_reconstruction",
        "operation": operation,
        "created_utc": _utc_now(),
        "source_event_file": str(event_path),
        "warning": "This record was reconstructed after the event and is not a contemporaneous log.",
    }
    _write_json_once(run_dir / "intent.json", intent)
    _append_event(events, "retrospective_intent_frozen")
    event_artifacts = event.get("observed_artifacts", [])
    if not isinstance(event_artifacts, list):
        raise ValueError("retrospective observed_artifacts must be a list")
    observed = _snapshots([*artifacts, *event_artifacts], working_directory)
    prepared = _content_record(
        {
            "schema": "a2s-dataset-processing-prepared-v1",
            "run_id": run_id,
            "record_type": "retrospective_reconstruction",
            "prepared_utc": _utc_now(),
            "operation": operation,
            "event": event,
            "event_file": _snapshot_file(event_path),
            "artifacts_observed_at_reconstruction": observed,
            "contemporaneous_stdout_available": False,
            "contemporaneous_stderr_available": False,
        }
    )
    _write_json_once(run_dir / "prepared.json", prepared)
    _append_event(events, "retrospective_evidence_frozen")
    result = _content_record(
        {
            "schema": "a2s-dataset-processing-run-v1",
            "run_id": run_id,
            "record_type": "retrospective_reconstruction",
            "operation": operation,
            "status": "RECORDED_RETROSPECTIVE",
            "prepared_sha256": prepared["content_sha256"],
            "historical_event_status": event.get("historical_event_status", "UNKNOWN"),
            "recorded_utc": _utc_now(),
            "limitations": limitations,
            "stdout": {"available": False, "reason": "not contemporaneously captured"},
            "stderr": {"available": False, "reason": "not contemporaneously captured"},
        }
    )
    _append_event(events, "retrospective_record_finalized")
    _write_json_once(run_dir / "run.json", result)
    _seal_and_index(
        audit_root, run_dir, run_id, operation, "RECORDED_RETROSPECTIVE"
    )
    result["run_path"] = str(run_dir)
    return result


def _pid_is_running(pid: int) -> bool:
    if pid <= 0:
        return False
    if os.name == "nt":
        kernel32, handle, _ = _open_windows_process(
            pid, _PROCESS_QUERY_LIMITED_INFORMATION
        )
        if handle is None:
            return False
        try:
            return bool(_windows_process_details(kernel32, handle, pid)["running"])
        except OSError:
            return False
        finally:
            kernel32.CloseHandle(handle)
    try:  # pragma: no cover - Windows is the project platform.
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def _indexed_run_ids(ledger: Path) -> set[str]:
    return {str(entry["run_id"]) for entry in _validated_ledger_entries(ledger)}


def _orphan_artifact_inventory(
    prepared: Mapping[str, Any] | None,
    intent: Mapping[str, Any],
    working_directory: Path,
) -> dict[str, Any]:
    source = prepared or intent
    declared = source.get("declared_outputs", intent.get("declared_outputs", []))
    publish_outputs = source.get(
        "publish_outputs", intent.get("publish_outputs", {})
    )
    if not isinstance(declared, list):
        declared = []
    if not isinstance(publish_outputs, dict):
        publish_outputs = {}
    staging_paths = list(publish_outputs)
    canonical_paths = list(publish_outputs.values())
    partial_paths: list[Path] = []
    for value in declared:
        path = _resolve_path(value, working_directory)
        if not path.name.endswith(".partial"):
            partial_paths.append(path.with_name(f".{path.name}.partial"))

    source_path: str | None = None
    destination_path: str | None = None
    if prepared:
        source_before = prepared.get("source_before")
        destination_before = prepared.get("destination_before")
        if isinstance(source_before, dict):
            source_path = source_before.get("path")
        if isinstance(destination_before, dict):
            destination_path = destination_before.get("path")
    source_path = source_path or intent.get("source")
    destination_path = destination_path or intent.get("destination")

    return {
        "declared_outputs": _snapshots(declared, working_directory),
        "publish_staging_sources": _snapshots(staging_paths, working_directory),
        "publish_canonical_destinations": _snapshots(
            canonical_paths, working_directory
        ),
        "dot_partial_candidates": _snapshots(partial_paths, working_directory),
        "relocation_source": (
            _snapshot_path(_resolve_path(source_path, working_directory))
            if source_path
            else None
        ),
        "relocation_destination": (
            _snapshot_path(_resolve_path(destination_path, working_directory))
            if destination_path
            else None
        ),
    }


def _claim_owner_status(run_dir: Path) -> tuple[str, dict[str, Any] | None]:
    claim_path = run_dir / "output_claim.json"
    if not claim_path.is_file():
        return "NOT_RECORDED", None
    try:
        claim = json.loads(claim_path.read_text(encoding="utf-8"))
        unsigned = {
            key: value for key, value in claim.items() if key != "content_sha256"
        }
        identity = claim.get("runner_identity")
        if (
            claim.get("run_id") != run_dir.name
            or claim.get("content_sha256") != _json_hash(unsigned)
            or not isinstance(identity, dict)
            or not _started_identity_complete(identity)
        ):
            return "IDENTITY_UNKNOWN", None
    except (OSError, TypeError, ValueError, KeyError, json.JSONDecodeError):
        return "IDENTITY_UNKNOWN", None
    return _inspect_started_process(identity)


def reconcile_orphans(
    root: str | Path = DEFAULT_ROOT, cwd: str | Path = "."
) -> dict[str, Any]:
    working_directory = Path(cwd).resolve(strict=True)
    audit_root = _resolve_path(root, working_directory)
    runs_root = audit_root / "runs"
    indexed = _indexed_run_ids(audit_root / "ledger.jsonl")
    recovered: list[str] = []
    still_running: list[str] = []
    identity_unknown: list[str] = []
    for run_dir in sorted(runs_root.glob("*")) if runs_root.exists() else []:
        if run_dir.name in indexed:
            continue
        intent_path = run_dir / "intent.json"
        if not intent_path.exists():
            _write_json_once(
                intent_path,
                {
                    "schema": "a2s-dataset-processing-intent-v1",
                    "run_id": run_dir.name,
                    "record_type": "recovered_empty_run_directory",
                    "operation": "unknown-orphan",
                    "created_utc": None,
                    "recovered_utc": _utc_now(),
                },
            )
        intent = json.loads(intent_path.read_text(encoding="utf-8"))
        existing_run_path = run_dir / "run.json"
        if existing_run_path.exists():
            existing_run = json.loads(existing_run_path.read_text(encoding="utf-8"))
            if not (run_dir / "seal.json").exists():
                _write_json_once(
                    run_dir / "recovery.json",
                    _content_record(
                        {
                            "schema": "a2s-dataset-processing-orphan-recovery-v1",
                            "run_id": run_dir.name,
                            "recovered_utc": _utc_now(),
                            "recovery_point": "run_json_written_before_seal_or_index",
                            "original_run_status_preserved": existing_run.get("status"),
                        }
                    ),
                )
                _append_event(run_dir / "events.jsonl", "completed_run_recovered_for_index")
                _seal_and_index(
                    audit_root,
                    run_dir,
                    run_dir.name,
                    str(existing_run.get("operation", intent.get("operation"))),
                    str(existing_run.get("status", "FAILED")),
                )
            else:
                _index_sealed_run(
                    audit_root,
                    run_dir,
                    run_dir.name,
                    str(existing_run.get("operation", intent.get("operation"))),
                    str(existing_run.get("status", "FAILED")),
                )
            recovered.append(run_dir.name)
            continue
        started_path = run_dir / "started.json"
        pid: int | None = None
        started: dict[str, Any] | None = None
        identity_status = "NOT_RECORDED"
        identity_observed: dict[str, Any] | None = None
        if started_path.exists():
            started = json.loads(started_path.read_text(encoding="utf-8"))
            pid = int(started["pid"])
            identity_status, identity_observed = _inspect_started_process(started)
            if identity_status == "MATCH":
                still_running.append(run_dir.name)
                continue
            if identity_status == "IDENTITY_UNKNOWN":
                identity_unknown.append(run_dir.name)
                continue
        claim_status, _ = _claim_owner_status(run_dir)
        if claim_status == "MATCH":
            still_running.append(run_dir.name)
            continue
        if claim_status == "IDENTITY_UNKNOWN":
            identity_unknown.append(run_dir.name)
            continue
        prepared_path = run_dir / "prepared.json"
        prepared = (
            json.loads(prepared_path.read_text(encoding="utf-8"))
            if prepared_path.exists()
            else None
        )
        observed_artifacts = _orphan_artifact_inventory(
            prepared, intent, working_directory
        )
        recovery = _content_record(
            {
                "schema": "a2s-dataset-processing-orphan-recovery-v1",
                "run_id": run_dir.name,
                "recovered_utc": _utc_now(),
                "pid_observed": pid,
                "pid_running_at_recovery": False,
                "process_identity_status": identity_status,
                "process_identity_observed": identity_observed,
                "exit_code": None,
                "child_ended_utc": None,
                "artifacts_observed_after_recovery": observed_artifacts,
                "publish_success_inferred": False,
                "limitation": (
                    "The original process identity was no longer active at "
                    "reconciliation; its exit code, end time, and publish outcome "
                    "are unknown and were not inferred."
                ),
            }
        )
        _write_json_once(run_dir / "recovery.json", recovery)
        _append_event(run_dir / "events.jsonl", "orphan_reconciled_unknown")
        result = _content_record(
            {
                "schema": "a2s-dataset-processing-run-v1",
                "run_id": run_dir.name,
                "record_type": "recovered_incomplete_execution",
                "operation": intent["operation"],
                "status": "FAILED",
                "process_status": "orphaned_unknown",
                "validation_status": "fail",
                "validation_failures": ["execution ended without a final run record"],
                "prepared_sha256": prepared.get("content_sha256") if prepared else None,
                "child_started_utc": (
                    started.get("started_utc") if started is not None else None
                ),
                "child_ended_utc": None,
                "duration_seconds": None,
                "pid": pid,
                "exit_code": None,
                "process_identity_status": identity_status,
                "process_identity_observed": identity_observed,
                "outputs_after_recovery": observed_artifacts["declared_outputs"],
                "artifacts_observed_after_recovery": observed_artifacts,
                "publish_success_inferred": False,
                "recovery_sha256": recovery["content_sha256"],
            }
        )
        _write_json_once(run_dir / "run.json", result)
        _seal_and_index(
            audit_root, run_dir, run_dir.name, intent["operation"], "FAILED"
        )
        recovered.append(run_dir.name)
    status = (
        "UNRESOLVED_OPERATIONS_PRESENT"
        if still_running and identity_unknown
        else "RUNNING_OPERATIONS_PRESENT"
        if still_running
        else "IDENTITY_UNKNOWN_OPERATIONS_PRESENT"
        if identity_unknown
        else "PASS"
    )
    return {
        "schema": "a2s-dataset-processing-reconciliation-v1",
        "status": status,
        "recovered": recovered,
        "still_running": still_running,
        "identity_unknown": identity_unknown,
    }


def verify_ledger(root: str | Path = DEFAULT_ROOT, cwd: str | Path = ".") -> dict[str, Any]:
    working_directory = Path(cwd).resolve(strict=True)
    audit_root = _resolve_path(root, working_directory)
    ledger = audit_root / "ledger.jsonl"
    entries = _validated_ledger_entries(ledger)
    runs_root = audit_root / "runs"
    if runs_root.exists():
        resolved_runs_root = runs_root.resolve(strict=True)
        if resolved_runs_root.parent != audit_root:
            raise ValueError("processing runs directory escapes the audit root")
        runs_root = resolved_runs_root
    previous: str | None = None
    checked = 0
    for entry in entries:
        digest = entry["entry_sha256"]
        run_id = entry.get("run_id")
        if not isinstance(run_id, str) or Path(run_id).name != run_id:
            raise ValueError(f"invalid processing run id: {run_id!r}")
        expected_run_dir = runs_root / run_id
        run_dir = Path(entry["run_path"]).resolve(strict=True)
        if run_dir != expected_run_dir or run_dir.parent != runs_root:
            raise ValueError(f"run_path escapes processing history for {run_id}")

        def sealed_path(name: str) -> Path:
            if name not in _SEALED_RUN_FILES:
                raise ValueError(f"invalid sealed file name for {run_id}: {name!r}")
            path = (run_dir / name).resolve(strict=True)
            if path.parent != run_dir:
                raise ValueError(f"sealed file escapes run directory for {run_id}: {name!r}")
            return path

        run_path = sealed_path("run.json")
        seal_path = (run_dir / "seal.json").resolve(strict=True)
        if seal_path.parent != run_dir:
            raise ValueError(f"seal.json escapes run directory for {run_id}")
        if _file_hash(run_path) != entry["run_json_sha256"]:
            raise ValueError(f"run.json hash mismatch for {entry['run_id']}")
        if _file_hash(seal_path) != entry["seal_sha256"]:
            raise ValueError(f"seal hash mismatch for {entry['run_id']}")
        run = json.loads(run_path.read_text(encoding="utf-8"))
        seal = json.loads(seal_path.read_text(encoding="utf-8"))
        if run.get("run_id") != run_id or seal.get("run_id") != run_id:
            raise ValueError(f"sealed run identity mismatch for {run_id}")
        if run.get("operation") != entry.get("operation"):
            raise ValueError(f"sealed run operation mismatch for {run_id}")
        if run.get("status") != entry.get("status"):
            raise ValueError(f"sealed run status mismatch for {run_id}")
        seal_digest = seal.get("content_sha256")
        seal_unsigned = {key: value for key, value in seal.items() if key != "content_sha256"}
        if seal_digest != _json_hash(seal_unsigned):
            raise ValueError(f"seal content hash mismatch for {entry['run_id']}")
        for name, expected in seal["files"].items():
            if _file_hash(sealed_path(name)) != expected:
                raise ValueError(f"sealed file mismatch for {entry['run_id']}/{name}")
        previous = digest
        checked += 1
    indexed = _indexed_run_ids(ledger)
    unindexed_paths = (
        [
            path
            for path in sorted(runs_root.glob("*"))
            if path.is_dir() and path.name not in indexed
        ]
        if runs_root.exists()
        else []
    )
    running_unindexed: list[str] = []
    unresolved_unindexed: list[str] = []
    for path in unindexed_paths:
        started_path = path / "started.json"
        if (path / "run.json").exists() or (path / "seal.json").exists():
            unresolved_unindexed.append(path.name)
            continue
        if not (path / "intent.json").is_file():
            unresolved_unindexed.append(path.name)
            continue
        if not started_path.is_file() or not (path / "prepared.json").is_file():
            claim_status, _ = _claim_owner_status(path)
            if claim_status == "MATCH":
                running_unindexed.append(path.name)
            else:
                unresolved_unindexed.append(path.name)
            continue
        try:
            intent = json.loads((path / "intent.json").read_text(encoding="utf-8"))
            prepared = json.loads(
                (path / "prepared.json").read_text(encoding="utf-8")
            )
            started = json.loads(started_path.read_text(encoding="utf-8"))
            command = prepared.get("command_argv")
            valid_active_context = (
                intent.get("run_id") == path.name
                and prepared.get("run_id") == path.name
                and started.get("run_id") == path.name
                and intent.get("record_type") == "contemporaneous_execution"
                and prepared.get("record_type") == "contemporaneous_execution"
                and isinstance(command, list)
                and started.get("command_sha256") == _json_hash(command)
            )
            if not valid_active_context:
                unresolved_unindexed.append(path.name)
                continue
            identity_status, _ = _inspect_started_process(started)
        except (OSError, ValueError, KeyError, json.JSONDecodeError):
            identity_status = "IDENTITY_UNKNOWN"
        if identity_status == "MATCH":
            running_unindexed.append(path.name)
            continue
        if identity_status == "IDENTITY_UNKNOWN":
            unresolved_unindexed.append(path.name)
            continue
        claim_status, _ = _claim_owner_status(path)
        if claim_status == "MATCH":
            running_unindexed.append(path.name)
        elif claim_status == "IDENTITY_UNKNOWN":
            unresolved_unindexed.append(path.name)
        else:
            unresolved_unindexed.append(path.name)
    if unresolved_unindexed:
        raise ValueError(
            "unindexed processing run directories require reconciliation: "
            f"{unresolved_unindexed}"
        )
    status = "RUNNING_OPERATIONS_PRESENT" if running_unindexed else "PASS"
    return {
        "schema": "a2s-dataset-processing-ledger-verification-v1",
        "status": status,
        "entries": checked,
        "head_entry_sha256": previous,
        "running_unindexed_runs": running_unindexed,
    }


def _environment_values(values: Sequence[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"environment override must be KEY=VALUE: {value}")
        key, setting = value.split("=", 1)
        result[key] = setting
    return result


def _path_values(values: Sequence[str], converter: Any) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"path expectation must be PATH=VALUE: {value}")
        path, expectation = value.rsplit("=", 1)
        result[path] = converter(expectation)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    execute = subparsers.add_parser("execute")
    execute.add_argument("--operation", required=True)
    execute.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    execute.add_argument("--cwd", type=Path, default=Path("."))
    execute.add_argument("--input", action="append", default=[])
    execute.add_argument("--output", action="append", default=[])
    execute.add_argument("--code-input", action="append", default=[])
    execute.add_argument("--config-input", action="append", default=[])
    execute.add_argument("--expected-bytes", action="append", default=[])
    execute.add_argument("--expected-sha256", action="append", default=[])
    execute.add_argument("--publish", action="append", default=[])
    execute.add_argument("--parent-run-id")
    execute.add_argument("--env", action="append", default=[])
    execute.add_argument("command", nargs=argparse.REMAINDER)
    record = subparsers.add_parser("record-retrospective")
    record.add_argument("--event-file", type=Path, required=True)
    record.add_argument("--event-id")
    record.add_argument("--artifact", action="append", default=[])
    record.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    record.add_argument("--cwd", type=Path, default=Path("."))
    record_batch = subparsers.add_parser("record-retrospective-batch")
    record_batch.add_argument("--event-file", type=Path, required=True)
    record_batch.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    record_batch.add_argument("--cwd", type=Path, default=Path("."))
    verify = subparsers.add_parser("verify")
    verify.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    verify.add_argument("--cwd", type=Path, default=Path("."))
    reconcile = subparsers.add_parser("reconcile")
    reconcile.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    reconcile.add_argument("--cwd", type=Path, default=Path("."))
    relocate = subparsers.add_parser("relocate-existing")
    relocate.add_argument("--operation", required=True)
    relocate.add_argument("--source", type=Path, required=True)
    relocate.add_argument("--destination", type=Path, required=True)
    relocate.add_argument("--expected-bytes", type=int)
    relocate.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    relocate.add_argument("--cwd", type=Path, default=Path("."))
    terminate = subparsers.add_parser("terminate-active")
    terminate.add_argument("--target-run-id", required=True)
    terminate.add_argument("--reason", required=True)
    terminate.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    terminate.add_argument("--cwd", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.action == "execute":
        command = args.command[1:] if args.command[:1] == ["--"] else args.command
        result = execute_run(
            operation=args.operation,
            command=command,
            cwd=args.cwd,
            root=args.root,
            inputs=args.input,
            outputs=args.output,
            code_inputs=args.code_input,
            config_inputs=args.config_input,
            expected_output_bytes=_path_values(args.expected_bytes, int),
            expected_output_sha256=_path_values(args.expected_sha256, str),
            publish_outputs=_path_values(args.publish, str),
            parent_run_id=args.parent_run_id,
            environment_overrides=_environment_values(args.env),
        )
        print(json.dumps({"status": result["status"], "run_path": result["run_path"]}, indent=2))
        raise SystemExit(0 if result["status"] == "SUCCESS" else 1)
    if args.action == "record-retrospective":
        result = record_retrospective(
            event_file=args.event_file,
            event_id=args.event_id,
            artifacts=args.artifact,
            cwd=args.cwd,
            root=args.root,
        )
        print(json.dumps({"status": result["status"], "run_path": result["run_path"]}, indent=2))
        return
    if args.action == "record-retrospective-batch":
        payload = json.loads(args.event_file.read_text(encoding="utf-8"))
        event_ids = [event["event_id"] for event in payload.get("events", [])]
        results = [
            record_retrospective(
                event_file=args.event_file,
                event_id=event_id,
                cwd=args.cwd,
                root=args.root,
            )
            for event_id in event_ids
        ]
        print(json.dumps({"status": "RECORDED_RETROSPECTIVE", "runs": len(results)}, indent=2))
        return
    if args.action == "reconcile":
        print(json.dumps(reconcile_orphans(args.root, args.cwd), indent=2, sort_keys=True))
        return
    if args.action == "relocate-existing":
        result = relocate_file_run(
            operation=args.operation,
            source=args.source,
            destination=args.destination,
            expected_bytes=args.expected_bytes,
            root=args.root,
            cwd=args.cwd,
        )
        print(json.dumps({"status": result["status"], "run_path": result["run_path"]}, indent=2))
        raise SystemExit(0 if result["status"] == "SUCCESS" else 1)
    if args.action == "terminate-active":
        result = terminate_active_run(
            target_run_id=args.target_run_id,
            reason=args.reason,
            root=args.root,
            cwd=args.cwd,
        )
        print(json.dumps({"status": result["status"], "run_path": result["run_path"]}, indent=2))
        raise SystemExit(0 if result["status"] == "SUCCESS" else 1)
    print(json.dumps(verify_ledger(args.root, args.cwd), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
