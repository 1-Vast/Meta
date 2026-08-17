"""Replace false `meta_test` seal claims in recorded artifacts. Numbers untouched.

Scope, established by producer, not by string match. The literal claim
`"physical: QPSMPData include_meta_test=False"` appears in 37 recorded JSON
artifacts. In 30 of them it is **true**: `train_level_shape.py`,
`train_reltransport.py` and `train_grammar_shape.py` hardcode
`include_meta_test=False`, so the artifact matches the code that produced it.

It is **false** in seven, all downstream of `scripts/r14_dispersion_audit.py`,
which constructed `QPSMPData` without the flag while writing the claim:

* 5 artifacts written directly by that script (`MetaSieve.R14DispersionAudit.v1`);
* 2 stage-level `RESULT.json` summaries that copied the string from them.

This script rewrites only the `meta_test` block of those seven, and adds a
`seal_correction` block recording what was actually true. Every metric, row and
hash in the files is preserved byte-for-byte; the correction is verified by
`SEAL_REPAIR_REPRODUCTION.json`, a re-run of the same audit under the repaired
seal whose 105 numeric A0 fields are bit-identical to the recorded ones.

Run: `conda run -n drug python -m tools.research.a2_readiness_v2.correct_seal_claims`
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
REPORT = ROOT / "report/meta_fewshot"

FALSE_CLAIM = "physical: QPSMPData include_meta_test=False"

TARGETS = (
    "stageR14_diagnostics_20260816/DISPERSION_meta_val.json",
    "stageR14_diagnostics_20260816/DISPERSION_attribution_meta_val.json",
    "stageR14_diagnostics_20260816/DISPERSION_rows.json",
    "stageR14_diagnostics_20260816/RESULT.json",
    "stageR14_screening_20260816/DISPERSION_R14_3seed.json",
    "stageR14_screening_20260816/DISPERSION_R14_meta_val.json",
    "stageR14_screening_20260816/RESULT.json",
)

CORRECTED_SEAL = (
    "population sealed, process not sealed: scripts/r14_dispersion_audit.py "
    "constructed QPSMPData without include_meta_test, so meta_test cells were "
    "parsed and indexed in this process. No meta_test value entered any "
    "recorded computation — every bank is fixed_nested_episode_banks('meta_val', "
    "...), donors are meta_val, the label scale is meta_train-only. Verified by "
    "re-running the identical audit under the repaired fail-closed seal: 105/105 "
    "numeric A0 fields bit-identical."
)

CORRECTION = {
    "date": "2026-08-16",
    "incident": "tools/research/a2_readiness_v2/GOVERNANCE_INCIDENT.md",
    "original_claim": FALSE_CLAIM,
    "why_it_was_false": (
        "the producing script never passed include_meta_test, and the class "
        "default at the time admitted the sealed split"),
    "numerical_impact": "none; no field in this artifact was modified",
    "verification": (
        "tools/research/a2_readiness_v2/SEAL_REPAIR_REPRODUCTION.json"),
    "repair": (
        "QPSMPData now defaults to include_meta_test=False, requires a written "
        "meta_test_authorization to open, and emits seal_record() so an "
        "artifact cannot contradict its constructor; contracts in "
        "tools/tests/test_meta_test_seal_contract.py"),
}


def corrected_block(existing: dict) -> dict:
    return {
        "evaluated": False,
        "included": True,
        "authorization": None,
        "seal": CORRECTED_SEAL,
        "seal_correction": CORRECTION,
    }


RELABELLED_SEAL = (
    "logical exclusion after parsing: QPSMPData(include_meta_test=False). "
    "Sealed rows are parsed from the all-label corpus and discarded during "
    "construction; they are unreachable from cells, tasks, components and "
    "materialize afterwards. This is NOT a physical label seal.")

RELABEL_NOTE = {
    "date": "2026-08-16",
    "was": FALSE_CLAIM,
    "why": ("the producing script did pass include_meta_test=False, so the "
            "substance of the claim was true, but 'physical' overstated it: "
            "cells.jsonl.gz is a single all-label artifact and every meta_test "
            "label is decompressed and parsed on every construction"),
    "numerical_impact": "none; no field in this artifact was modified",
    "specification_for_a_real_seal": (
        "tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md"),
}


def relabel_true_claims() -> int:
    """Correct the *wording* of artifacts whose seal claim was substantively true.

    Thirty artifacts written by `train_level_shape.py`, `train_reltransport.py`
    and `train_grammar_shape.py` carry the same literal string. Those scripts
    hardcode `include_meta_test=False`, so unlike the seven in `TARGETS` their
    claim matched their code. What was wrong is only the word "physical",
    which asserts an isolation property the lineage has never had.

    This pass rewrites the `seal` string and adds a `seal_relabel` note. It
    does not touch `included`, `evaluated`, or any other field, and it asserts
    non-seal content is unchanged before writing.
    """
    changed = 0
    for path in sorted(REPORT.rglob("RESULT.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        block = payload.get("meta_test")
        # Three literal variants of the same overstatement were emitted by
        # different producers. All are substantively true and all use the wrong
        # word; an already-relabelled block is skipped.
        seal = block.get("seal") if isinstance(block, dict) else None
        if not isinstance(seal, str) or not seal.startswith("physical:"):
            continue
        if block.get("included") is not False:
            continue          # an open run is not this pass's business
        before = {k: v for k, v in payload.items() if k != "meta_test"}
        block["seal"] = RELABELLED_SEAL + (
            "" if "never reached" not in seal
            else " (never reached; no training run)")
        block["seal_relabel"] = dict(RELABEL_NOTE, was=seal)
        payload["meta_test"] = block
        assert {k: v for k, v in payload.items()
                if k != "meta_test"} == before, f"non-seal content changed in {path}"
        path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        changed += 1
    print(f"{changed} artifact(s) relabelled (claim was true, wording was not)")
    return changed


def main() -> int:
    changed = 0
    for name in TARGETS:
        path = REPORT / name
        payload = json.loads(path.read_text(encoding="utf-8"))
        block = payload.get("meta_test")
        if not isinstance(block, dict) or block.get("seal") != FALSE_CLAIM:
            print(f"skip (already corrected or unexpected shape): {name}")
            continue
        before = {k: v for k, v in payload.items() if k != "meta_test"}
        payload["meta_test"] = corrected_block(block)
        after = {k: v for k, v in payload.items() if k != "meta_test"}
        assert before == after, f"non-seal content changed in {name}"
        path.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        changed += 1
        print(f"corrected {name}")
    print(f"\n{changed} artifact(s) corrected; no numeric field modified")
    relabel_true_claims()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
