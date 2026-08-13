"""Unified command-line entry point for retained MetaSieve capabilities."""
from __future__ import annotations

from dataclasses import dataclass
import importlib
import sys
from typing import Sequence


@dataclass(frozen=True)
class Command:
    module: str
    injected_args: tuple[str, ...] = ()
    cuda_required: bool = False
    description: str = ""


COMMANDS: dict[tuple[str, ...], Command] = {
    ("status",): Command(
        "scripts.project_status", description="show governed project status"),
    ("archive", "status"): Command(
        "scripts.project_status", ("--archive",),
        description="show the read-only research archive policy"),
    ("verify", "tests"): Command(
        "scripts.verify_project", description="run the maintained test suite"),
    ("qpsmp", "train"): Command(
        "scripts.train_qpsmp",
        cuda_required=True,
        description="train the active BPSF meta-learner"),
    ("qpsmp", "evaluate"): Command(
        "scripts.evaluate_qpsmp",
        cuda_required=True,
        description="run governed nested-k BPSF evaluation"),
    ("geometry", "pretrain"): Command(
        "scripts.pretrain_pair_geometry",
        cuda_required=True,
        description="pretrain the source-only BPSF geometry head"),
    ("geometry", "evaluate"): Command(
        "scripts.evaluate_pair_geometry",
        cuda_required=True,
        description="evaluate the source-only BPSF geometry head"),
    ("data", "prepare"): Command(
        "scripts.preprocess_dataset",
        description="compile and audit a governed canonical dataset"),
    ("data", "verify"): Command(
        "scripts.verify_dataset",
        description="read-only audit of a canonical dataset"),
}


def _render_help(prefix: tuple[str, ...] = ()) -> str:
    rows = []
    for path, command in sorted(COMMANDS.items()):
        if path[:len(prefix)] != prefix:
            continue
        suffix = path[len(prefix):]
        if not suffix:
            continue
        rows.append((" ".join(suffix), command.description))
    heading = "MetaSieve retained command surface"
    usage_prefix = " ".join(("python main.py", *prefix))
    lines = [heading, "", f"usage: {usage_prefix} <command> [options]", "", "commands:"]
    width = max((len(name) for name, _ in rows), default=0)
    lines.extend(f"  {name:<{width}}  {description}" for name, description in rows)
    lines.extend((
        "",
        "Run a leaf command with --help for its module-specific arguments.",
        "Model execution is development-only and requires CUDA device cuda:0.",
    ))
    return "\n".join(lines)


def _match(argv: Sequence[str]) -> tuple[tuple[str, ...], Command] | None:
    matches = [
        (path, command) for path, command in COMMANDS.items()
        if tuple(argv[:len(path)]) == path
    ]
    return max(matches, key=lambda item: len(item[0]), default=None)


def _cuda_args(args: list[str]) -> list[str]:
    if "-h" in args or "--help" in args:
        return args
    values = []
    for index, value in enumerate(args):
        if value == "--device":
            if index + 1 >= len(args):
                return args  # Let the leaf argparse report the missing value.
            values.append(args[index + 1])
        elif value.startswith("--device="):
            values.append(value.split("=", 1)[1])
    if len(values) > 1:
        raise ValueError("--device may be specified only once")
    if values and values[0] != "cuda:0":
        raise ValueError(
            f"MetaSieve model commands require --device cuda:0; received {values[0]!r}")
    return args if values else [*args, "--device", "cuda:0"]


def _invoke(command: Command, args: list[str]) -> int:
    forwarded = [*command.injected_args, *args]
    if command.cuda_required:
        forwarded = _cuda_args(forwarded)
    module = importlib.import_module(command.module)
    prior = sys.argv
    sys.argv = [command.module, *forwarded]
    try:
        result = module.main()
    finally:
        sys.argv = prior
    return 0 if result is None else int(result)


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help"}:
        print(_render_help())
        return 0
    if len(args) >= 2 and args[-1] in {"-h", "--help"}:
        prefix = tuple(args[:-1])
        if prefix not in COMMANDS and any(path[:len(prefix)] == prefix for path in COMMANDS):
            print(_render_help(prefix))
            return 0
    matched = _match(args)
    if matched is None:
        print(f"unknown command: {' '.join(args)}", file=sys.stderr)
        print(_render_help(), file=sys.stderr)
        return 2
    path, command = matched
    try:
        return _invoke(command, args[len(path):])
    except (FileExistsError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
