"""End-to-end one-step smoke for every training entry point.

Runs each trainer for 2 steps on CPU and verifies the full pipeline surface:
forward API, loss, backward, optimizer step, checkpoint save, RESULT.json
schema, PREDICTIONS jsonl, progress.jsonl, donors block, checkpoint sha256
and the meta_test seal fields. A stage runner must pass this smoke before
launching real training (engineering contract 2026-08-16).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable
SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"


def run(command: list[str], output: Path) -> dict:
    completed = subprocess.run(
        command, cwd=ROOT, capture_output=True, text=True,
        env={**__import__("os").environ, "PYTHONIOENCODING": "utf-8"},
        timeout=3600)
    if completed.returncode != 0:
        raise RuntimeError(
            f"smoke command failed ({completed.returncode}): {' '.join(command)}\n"
            f"--- stdout tail ---\n{completed.stdout[-2000:]}\n"
            f"--- stderr tail ---\n{completed.stderr[-2000:]}")
    return {"command": " ".join(command), "output": str(output)}


def check_artifacts(output: Path, expect_meta_val: bool = True) -> dict:
    problems = []
    for name in ("checkpoint.pt", "RESULT.json", "PREDICTIONS_meta_val.jsonl",
                 "progress.jsonl"):
        path = output / name
        if not path.is_file() or path.stat().st_size == 0:
            problems.append(f"missing or empty: {name}")
    payload = json.loads((output / "RESULT.json").read_text(encoding="utf-8"))
    for field in ("config", "split_assignment_sha256", "checkpoint_sha256",
                  "donors", "meta_test"):
        if field not in payload:
            problems.append(f"RESULT.json missing field: {field}")
    if payload.get("meta_test", {}).get("evaluated") is not False:
        problems.append("meta_test.evaluated is not False")
    if payload.get("meta_test", {}).get("included") is not False:
        problems.append("meta_test.included is not False")
    if expect_meta_val:
        try:
            rows = [json.loads(line) for line in
                    (output / "PREDICTIONS_meta_val.jsonl").read_text(
                        encoding="utf-8").splitlines() if line.strip()]
        except json.JSONDecodeError:
            rows = []
        if not rows or "target" not in rows[0] or "k" not in rows[0]:
            problems.append("PREDICTIONS_meta_val.jsonl rows malformed")
    return {"ok": not problems, "problems": problems}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--trainers", nargs="+",
                        default=("scripts.train_qpsmp", "scripts.train_level_shape",
                                 "scripts.train_reltransport",
                                 "scripts.train_grammar_shape"))
    parser.add_argument("--device", default=None,
                        help="smoke device; defaults to the training device "
                             "(cuda when available) so the smoke mirrors the "
                             "real run's execution path")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    device = args.device or (
        "cuda" if __import__("torch").cuda.is_available() else "cpu")
    results = {}
    for trainer in args.trainers:
        name = trainer.split(".")[-1]
        variants: list[tuple[str, list[str]]] = []
        if name == "train_qpsmp":
            variants.append((name, ["--arch", "similarity_only"]))
        elif name == "train_reltransport":
            variants.extend([
                (f"{name}_nogate", ["--no-gate"]),
                (name, []),
                (f"{name}_ordinary", ["--ordinary"]),
            ])
        elif name == "train_grammar_shape":
            variants.extend([
                (name, []),
                (f"{name}_nocliff", ["--cliff-pair-weight", "1.0"]),
            ])
        else:
            variants.append((name, []))
        for tag, extra in variants:
            output = args.output_root / f"smoke_{tag}"
            if output.exists():
                import shutil
                shutil.rmtree(output)
            command = [PYTHON, "-m", trainer,
                       "--split-directory", str(SPLIT),
                       "--output", str(output),
                       "--steps", "2", "--episodes-per-step", "1",
                       "--val-interval", "1", "--device", device,
                       *extra]
            record = run(command, args.output_root)
            check = check_artifacts(output,
                                    expect_meta_val=(name != "train_qpsmp"))
            record.update(check)
            results[tag] = record
            print(f"smoke {tag}: {'OK' if check['ok'] else check['problems']}",
                  flush=True)
            if not check["ok"]:
                sys.exit(1)
    summary = args.output_root / "SMOKE_SUMMARY.json"
    summary.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n",
                       encoding="utf-8")
    print("all smokes passed; summary written to", summary)


if __name__ == "__main__":
    main()
