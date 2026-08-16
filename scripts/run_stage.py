"""Robust stage runner (engineering contract 2026-08-16).

Reads a stage spec JSON and executes it with the required discipline:

1. runs the one-step end-to-end smoke first and refuses to proceed on
   failure;
2. runs every arm sequentially, checking each exit code — the first
   failure stops the stage and the comparison script is NOT run;
3. rejects an arm whose output directory is missing or has no RESULT.json
   after a nominally successful exit;
4. records every actual command line in `commands.jsonl`.

Stage spec schema:

{
  "name": "...",
  "smoke": {"output_root": "...", "trainers": ["scripts.train_reltransport"]},
  "runs": [
    {"name": "arm_seed", "trainer": "scripts.train_reltransport",
     "args": ["--no-gate", ...], "output": "..."},
    ...
  ],
  "compare": {"trainer": "scripts.stageR6_compare_arms",
              "args": ["--arm", "A0=...", ...]}
}
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable


def execute(command: list[str], log_path: Path) -> int:
    completed = subprocess.run(
        command, cwd=ROOT, text=True, timeout=None,
        env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(
            {"command": " ".join(command), "exit_code": completed.returncode})
            + "\n")
    return completed.returncode


def _completed(output: Path) -> bool:
    """A run counts as completed only when its RESULT.json exists."""
    return output.is_dir() and (output / "RESULT.json").is_file()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, required=True)
    args = parser.parse_args()
    spec = json.loads(args.spec.read_text(encoding="utf-8"))
    log_path = args.spec.with_suffix(".commands.jsonl")
    if log_path.exists():
        log_path.unlink()

    # 1. Smoke first.
    smoke = spec.get("smoke")
    if smoke is None:
        raise ValueError("stage spec requires a smoke block")
    smoke_output = Path(smoke["output_root"])
    if smoke_output.exists():
        shutil.rmtree(smoke_output)
    code = execute(
        [PYTHON, "-m", "scripts.stage_smoke",
         "--output-root", str(smoke_output),
         "--trainers", *smoke.get("trainers",
                                  ["scripts.train_qpsmp",
                                   "scripts.train_level_shape",
                                   "scripts.train_reltransport",
                                   "scripts.train_grammar_shape"])],
        log_path)
    if code != 0:
        print("SMOKE FAILED — stage aborted before any training run.",
              file=sys.stderr)
        sys.exit(code)

    # 2. Runs, sequentially, first failure stops the stage.
    for run_spec in spec.get("runs", []):
        if "output" in run_spec:
            output = Path(run_spec["output"])
        elif "--output" in run_spec.get("args", []):
            args_list = run_spec["args"]
            output = Path(args_list[args_list.index("--output") + 1])
        else:
            raise ValueError(
                f"run {run_spec.get('name')} declares no output directory")
        if _completed(output):
            print(f"skipping completed run: {run_spec['name']}", flush=True)
            continue
        if output.is_dir():
            # A partial directory from an aborted attempt is not a result.
            shutil.rmtree(output)
            print(f"removing partial output and rerunning: {run_spec['name']}",
                  flush=True)
        command = [PYTHON, "-m", run_spec["trainer"], *run_spec["args"]]
        print(f"running {run_spec['name']}: {run_spec['trainer']}", flush=True)
        code = execute(command, log_path)
        if code != 0:
            print(f"RUN FAILED ({code}): {run_spec['name']} — stage stopped; "
                  "comparison not executed.", file=sys.stderr)
            sys.exit(code)
        if not output.is_dir() or not (output / "RESULT.json").is_file():
            print(f"RUN PRODUCED NO RESULT: {run_spec['name']} ({output}) — "
                  "stage stopped; comparison not executed.", file=sys.stderr)
            sys.exit(1)

    # 3. Comparison, only after every run succeeded.
    compare = spec.get("compare")
    if compare is not None:
        command = [PYTHON, "-m", compare["trainer"], *compare["args"]]
        print(f"running comparison: {compare['trainer']}", flush=True)
        code = execute(command, log_path)
        if code != 0:
            print(f"COMPARISON FAILED ({code}).", file=sys.stderr)
            sys.exit(code)
    print(f"stage {spec.get('name')} complete; commands recorded in {log_path}")


if __name__ == "__main__":
    main()
