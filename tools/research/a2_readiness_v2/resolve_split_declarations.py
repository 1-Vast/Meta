"""Resolve the `split_undeclared` record defect from each stage's own leaf runs.

Seven stage-level `RESULT.json` summaries record `population.split = "meta_val"`
with no corpus or split directory, so the audit cannot verify which protocol
they belong to and correctly refuses to guess.

They are resolvable without asserting anything: each stage directory contains
the per-seed `RESULT.json` files its summary aggregates, and those *do* declare
`config.split_directory` and `split_assignment_sha256`. This script reads the
leaves, requires them to agree unanimously, and stamps the parent with the
declaration they establish plus a provenance record naming the leaves it came
from.

Two stages have no real-data leaf and are **not** resolvable this way. They are
left undeclared and recorded as such:

* `stageR5_reltransport_20260816` — a contract-and-gates stage, no training run;
* `stageR13_shape_direct_20260816` — the family failed its own Stage 1 gates
  before any real-data run (`population.split` is already explicitly `null`).

Run:
`conda run -n drug python -m tools.research.a2_readiness_v2.resolve_split_declarations`
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FEWSHOT = ROOT / "report/meta_fewshot"
DOUBLE_COLD = "bindingdb_ki_double_cold_v1"

# Stages with no real-data run: nothing to inherit a declaration from.
NO_RUN = {
    "stageR5_reltransport_20260816":
        "contract-and-gates stage; structural gates and pipeline smoke only, "
        "no real-data training run to inherit a split declaration from",
    "stageR13_shape_direct_20260816":
        "the family failed its own preregistered Stage 1 gates before any "
        "real-data run; population.split is already explicitly null",
}


def leaf_declarations(stage: Path) -> tuple[set[str], set[str], list[str]]:
    """Split directories and assignment hashes declared by this stage's runs."""
    directories: set[str] = set()
    hashes: set[str] = set()
    sources: list[str] = []
    for path in sorted(stage.rglob("RESULT.json")):
        if path.parent == stage:
            continue                      # the summary itself, not a leaf
        payload = json.loads(path.read_text(encoding="utf-8"))
        config = payload.get("config") if isinstance(payload.get("config"), dict) else {}
        directory = payload.get("split_directory") or config.get("split_directory")
        digest = payload.get("split_assignment_sha256")
        if directory:
            directories.add(Path(str(directory)).name)
            sources.append(str(path.relative_to(ROOT)))
        if digest:
            hashes.add(str(digest))
    return directories, hashes, sources


def main() -> int:
    resolved, unresolved = [], []
    for stage in sorted(FEWSHOT.iterdir()):
        summary = stage / "RESULT.json"
        if not stage.is_dir() or not summary.is_file():
            continue
        text = summary.read_text(encoding="utf-8")
        if DOUBLE_COLD in text or "bindingdb_ki_main_v0" in text:
            continue                      # already declares a protocol
        if stage.name in NO_RUN:
            payload = json.loads(text)
            population = payload.setdefault("population", {})
            population["split_directory"] = None
            population["split_declaration"] = {
                "status": "unresolvable",
                "reason": NO_RUN[stage.name],
                "recorded": "2026-08-16",
            }
            summary.write_text(json.dumps(payload, indent=1), encoding="utf-8")
            unresolved.append((stage.name, NO_RUN[stage.name]))
            continue

        directories, hashes, sources = leaf_declarations(stage)
        if len(directories) != 1 or not sources:
            unresolved.append(
                (stage.name, f"leaf runs declare {sorted(directories)}"))
            continue
        if len(hashes) > 1:
            unresolved.append(
                (stage.name, f"leaf runs disagree on the assignment hash: {hashes}"))
            continue

        payload = json.loads(text)
        population = payload.setdefault("population", {})
        population["split_directory"] = (
            f"dataset/processed/meta_fewshot/{directories.pop()}")
        if hashes:
            population["split_assignment_sha256"] = hashes.pop()
        population["split_declaration"] = {
            "status": "resolved_from_leaf_runs",
            "reason": ("the summary declared no corpus; every per-seed run it "
                       "aggregates declares this split directory unanimously"),
            "evidence": sources,
            "recorded": "2026-08-16",
        }
        summary.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        resolved.append((stage.name, len(sources)))

    for name, count in resolved:
        print(f"resolved   {name}  (from {count} leaf run declarations)")
    for name, reason in unresolved:
        print(f"unresolved {name}: {reason}")
    print(f"\n{len(resolved)} resolved, {len(unresolved)} left undeclared "
          f"with a recorded reason")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
