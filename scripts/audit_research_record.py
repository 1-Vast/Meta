"""Verify the narrative research record against the leaf artifacts.

Every headline number in `report/BOUNDARY_20260816.md`,
`report/CURRENT_MODEL_EVIDENCE.md` and `report/EVIDENCE_LEDGER.md` is supposed
to have a `RESULT.json` / `COMPARE_*.json` authority under `report/`. This
command recomputes those numbers from the artifacts and reports every
disagreement, so a stale narrative cannot silently outlive its evidence.

It also checks the two structural invariants that protect the confirmation
population:

* every recorded run declares the double-cold `meta_test` unevaluated;
* every recorded `checkpoint_sha256` matches the retained checkpoint bytes
  (when the checkpoint payload is present locally; payloads are gitignored).

Usage::

    python -m scripts.audit_research_record
    python -m scripts.audit_research_record --json report/RECORD_AUDIT.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

# Both invocations must work and must fail for the same reasons:
#   python scripts/audit_research_record.py
#   python -m scripts.audit_research_record
# Direct invocation puts `scripts/` on `sys.path`, not the repository root, so
# `check_strict_loading`'s lazy `from scripts.stageR6_compare_arms import ...`
# raised ModuleNotFoundError and the process exited 1 — the same exit code the
# open governance incident produces. A green-looking failure and a real finding
# were indistinguishable from the shell.
if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"
FEWSHOT = REPORT / "meta_fewshot"

# The double-cold protocol. Absolute numbers are only comparable within it.
DOUBLE_COLD = "bindingdb_ki_double_cold_v1"

# The older, already-consumed protein-only-cold corpus. An artifact must *name*
# this to be called older-protocol; silence is not evidence of it.
MAIN_V0 = "bindingdb_ki_main_v0"

# Arm summaries live in these comparison authorities. Each entry maps a stage
# label to the COMPARE file that carries the three-seed arm means.
COMPARE_SOURCES = {
    "R3R4": FEWSHOT / "stageR3R4_level_shape_20260815" / "COMPARE_v2_3seed.json",
    "R7": FEWSHOT / "stageR7_reltransport_3seed_20260816" / "COMPARE_R7_meta_val.json",
    "R8": FEWSHOT / "stageR8_stronger_shape_20260816" / "COMPARE_R8_meta_val.json",
    "R9": FEWSHOT / "stageR9_cliffweight_20260816" / "COMPARE_R9_meta_val.json",
    "R10": FEWSHOT / "stageR10_variance_20260816" / "COMPARE_R10_meta_val.json",
    "R11": FEWSHOT / "stageR11_grammar_shape_20260816" / "COMPARE_R11_meta_val.json",
    "R12": FEWSHOT / "stageR12_margin_20260816" / "COMPARE_R12_meta_val.json",
}

# Canonical arm identity: the same arm is re-reported by several stages against
# the frozen incumbent, so the record must name one owning stage per arm.
ARM_OWNER = {
    "A0": "R3R4", "B1_R3R4": "R3R4", "B2": "R3R4", "B3": "R3R4",
    "A1": "R7", "A2": "R7", "A3": "R7",
    "B1": "R8",
    "C1": "R9", "C2": "R9",
    "D1": "R10", "G1": "R11", "D2": "R12",
}

METRICS = {
    "mse": "full_mse_pk",
    "ci": "full_ci",
    "spearman": "full_spearman",
    "calibration": "calibration_pk",
    "shape": "shape_pk",
    "cliff_sign": "full_cliff_sign_accuracy",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_arms() -> dict[str, dict]:
    """Return {arm_label: metrics} using one owning stage per arm."""
    arms: dict[str, dict] = {}
    for stage, path in COMPARE_SOURCES.items():
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        for arm, blocks in payload.get("summary", {}).items():
            if not isinstance(blocks, dict):
                continue
            label = arm if ARM_OWNER.get(arm) == stage else f"{arm}@{stage}"
            if ARM_OWNER.get(arm) not in (None, stage):
                continue        # a re-report of another stage's arm
            row = {"stage": stage, "source": str(path.relative_to(ROOT))}
            for k in ("0", "1", "2", "3", "5"):
                block = blocks.get(k)
                if not isinstance(block, dict):
                    continue
                row[k] = {name: block.get(key) for name, key in METRICS.items()}
            arms[label] = row
    return arms


def pareto(arms: dict[str, dict]) -> list[str]:
    """Non-dominated arms on (k=0 MSE minimise, k=0 CI maximise).

    A frontier is a set of whole configurations. Reporting the best MSE and the
    best CI as one interval would describe a model that does not exist.
    """
    points = {
        label: (row["0"]["mse"], row["0"]["ci"])
        for label, row in arms.items()
        if "0" in row
        and isinstance(row["0"].get("mse"), (int, float))
        and isinstance(row["0"].get("ci"), (int, float))
    }
    front = []
    for label, (mse, ci) in points.items():
        dominated = any(
            (m <= mse and c >= ci) and (m < mse or c > ci)
            for other, (m, c) in points.items() if other != label
        )
        if not dominated:
            front.append(label)
    return sorted(front, key=lambda label: points[label][0])


def check_seals() -> dict:
    """Classify every artifact's relationship to the double-cold `meta_test`.

    Three states, and only the third is a violation:

    ``sealed_explicit``
        a double-cold run carrying ``meta_test.evaluated == False``. This is
        **logical exclusion after parsing**, not a physical label seal: the
        corpus is a single all-label artifact and every sealed row is
        decompressed and parsed before being discarded. See
        `QPSMPData.seal_record()["isolation"]`.
    ``sealed_implicit``
        a double-cold run from before R5 that records only ``meta_val`` and
        never mentions ``meta_test`` at all.
    ``sealed_quarantined``
        a pre-R5 double-cold run whose trainer *did* compute `meta_test`
        metrics before the seal was authorised, with those numbers moved into
        a `SEALED_meta_test_DO_NOT_OPEN.json` sidecar and the `test` field
        replaced by a pointer. The values exist on disk. They are never read
        and were used for no decision, but the honest claim is
        "metric-unconsumed", not "never computed" and not "never read".
    ``process_unsealed``
        a double-cold run whose *population* was sealed — no `meta_test` value
        entered any recorded metric — but whose *process* was not: the
        producing script omitted `include_meta_test`, so the sealed cells were
        parsed **and indexed** in memory, and an episode could in principle have
        been drawn from them. Corrected on 2026-08-16 and carrying a
        `meta_test.seal_correction` block. This is tracked separately because
        folding it into `sealed_explicit` would overstate the seal, and folding
        it into `violations` would overstate the damage: the numbers are
        verified unchanged (`SEAL_REPAIR_REPRODUCTION.json`).

        **A non-empty `process_unsealed` list is a standing open incident.**
        `violations == 0` does not clear it and must never be reported as
        though it did; `main()` prints it as an explicit banner and returns a
        non-zero status while any entry remains.
    ``older_protocol``
        an artifact that **names** the pre-double-cold `bindingdb_ki_main_v0`
        population. Stages 4/6/7 legitimately report a `meta_test` there; it is
        a *different, consumed* population and must never be conflated with the
        sealed double-cold confirmation split.
    ``split_undeclared``
        an artifact that names neither. Until 2026-08-16 these were swept into
        `older_protocol`, which asserted a consumed population the file never
        claimed — every one of them is in fact a double-cold stage summary that
        records `population.split = "meta_val"` and no directory. They are
        reported as their own state rather than assigned a protocol the
        evidence does not support.
    """
    explicit, implicit, quarantined, older, violations = [], [], [], [], []
    process_unsealed, undeclared = [], []
    for path in sorted(FEWSHOT.rglob("RESULT.json")):
        name = str(path.relative_to(ROOT))
        raw = path.read_text(encoding="utf-8")
        try:
            payload = json.loads(raw)
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            violations.append({"artifact": name, "problem": f"unreadable: {error}"})
            continue
        # Artifacts record the split in at least four places: top level
        # (train_reltransport), inside `config` (train_qpsmp), inside
        # `population` (the hand-written stage summaries), or not at all.
        # Scanning the whole document is the only detection that cannot
        # silently skip an artifact's seal check by missing a nesting level.
        if DOUBLE_COLD not in raw:
            (older if MAIN_V0 in raw else undeclared).append(name)
            continue
        seal = payload.get("meta_test")
        sidecars = sorted(path.parent.glob("SEALED_meta_test*"))
        if isinstance(seal, dict) and seal.get("seal_correction"):
            process_unsealed.append({
                "artifact": name,
                "incident": seal["seal_correction"].get("incident"),
                "numerical_impact": seal["seal_correction"].get("numerical_impact")})
        elif isinstance(seal, dict) and seal.get("evaluated") is False:
            if seal.get("included"):
                violations.append({
                    "artifact": name,
                    "problem": "meta_test included without a seal_correction "
                               f"or authorization: {seal!r}"})
            else:
                explicit.append(name)
        elif sidecars and str(payload.get("test", "")).startswith("SEALED"):
            quarantined.append({"artifact": name,
                                "sidecar": str(sidecars[0].relative_to(ROOT))})
        elif seal is None and "meta_test" not in raw:
            implicit.append(name)
        elif isinstance(seal, dict) and seal.get("opened") is False and (
                seal.get("evaluated") is False):
            # The hand-written stage summaries use `opened`/`evaluated` rather
            # than the trainer's `included`/`evaluated` vocabulary.
            explicit.append(name)
        else:
            violations.append({"artifact": name,
                               "problem": f"double-cold meta_test not sealed: {seal!r}"})
    return {"sealed_explicit": explicit, "sealed_implicit": implicit,
            "sealed_quarantined": quarantined,
            "process_unsealed": process_unsealed, "older_protocol": older,
            "split_undeclared": undeclared, "violations": violations}


def check_checkpoint_hashes() -> dict:
    """Recompute every recorded checkpoint sha256 whose payload is present."""
    checked, mismatched, absent = 0, [], 0
    for path in sorted(FEWSHOT.rglob("RESULT.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        recorded = payload.get("checkpoint_sha256")
        if not isinstance(recorded, str):
            continue
        candidates = [path.parent / name
                      for name in ("checkpoint.pt", "best.pt", "model.pt")]
        checkpoint = next((c for c in candidates if c.exists()), None)
        if checkpoint is None:
            absent += 1
            continue
        checked += 1
        actual = sha256_file(checkpoint)
        if actual != recorded:
            mismatched.append({"artifact": str(path.relative_to(ROOT)),
                               "recorded": recorded, "actual": actual})
    return {"verified": checked, "payload_absent": absent, "mismatched": mismatched}


def check_strict_loading(stages: list[str] | None = None) -> dict:
    """Reload every retained checkpoint through the production loader.

    `load_state_dict` defaults to `strict=True`, so a renamed or dropped
    parameter surfaces here and nowhere else once training has finished.

    A failure is **not automatically a defect**. Checkpoints written before a
    deliberate architecture change belong to the superseded architecture; they
    are retained because their `RESULT.json` metrics are the evidence, not
    because the bytes are re-runnable. Those are reported as `superseded`, and
    the caller must be able to name the change. Anything else is a real break.
    """
    import types

    import torch

    from scripts.stageR6_compare_arms import load_arm

    manifest = (ROOT / "dataset/processed/meta_fewshot"
                / "bindingdb_ki_main_v0_protein_bank" / "manifest.json")
    if not manifest.exists():
        return {"skipped": "protein bank manifest absent (dataset is local-only)"}
    hidden = int(json.loads(manifest.read_text(encoding="utf-8"))["hidden_dim"])
    stub = types.SimpleNamespace(
        protein_bank=types.SimpleNamespace(manifest={"hidden_dim": hidden}))

    # Retained pre-fix arms of the R3R4 ladder. Their `TypedLigandChannels`
    # pooled each ligand to five vectors before any protein contact; the fix
    # replaced it with 16 query-slot tokens and changed the anchor shape, so
    # these state dicts cannot enter the current LevelShapeModel by design.
    # They are evidence for the identifiability and capacity defects the R3R4
    # report documents, not comparators.
    SUPERSEDED = ("stageR3R4_level_shape_20260815/A1_shared",
                  "stageR3R4_level_shape_20260815/A2_routed",
                  "stageR3R4_level_shape_20260815/A3_full")

    stages = stages or ["stageR3R4_level_shape_20260815",
                        "stageR7_reltransport_3seed_20260816",
                        "stageR8_stronger_shape_20260816",
                        "stageR9_cliffweight_20260816",
                        "stageR10_variance_20260816",
                        "stageR11_grammar_shape_20260816",
                        "stageR12_margin_20260816"]

    loaded, superseded, broken = [], [], []
    for stage in stages:
        for result in sorted((FEWSHOT / stage).glob("*/RESULT.json")):
            checkpoint = result.parent / "checkpoint.pt"
            if not checkpoint.exists():
                continue
            label = f"{stage}/{result.parent.name}"
            try:
                model, kind, _ = load_arm(checkpoint, stub, "cpu")
                loaded.append({"artifact": label, "arch": kind,
                               "parameters": sum(p.numel() for p in model.parameters())})
            except Exception as error:                             # noqa: BLE001
                row = {"artifact": label, "error": repr(error)[:200]}
                (superseded if label.startswith(SUPERSEDED) else broken).append(row)
    return {"loaded": loaded, "superseded_architecture": superseded,
            "broken": broken}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=None,
                        help="also write the audit as JSON to this path")
    parser.add_argument("--skip-loading", action="store_true",
                        help="skip the strict checkpoint reload (needs torch "
                             "and the local protein-bank manifest)")
    arguments = parser.parse_args()

    arms = collect_arms()
    frontier = pareto(arms)
    seals = check_seals()
    hashes = check_checkpoint_hashes()
    loading = {"skipped": "--skip-loading"} if arguments.skip_loading \
        else check_strict_loading()

    print(f"{'arm':<8}{'stage':<7}{'k0 MSE':>9}{'k0 CI':>8}{'rho':>8}"
          f"{'calib':>8}{'shape':>8}{'cliff0':>8}{'cliff5':>8}")
    for label, row in sorted(arms.items()):
        zero, five = row.get("0", {}), row.get("5", {})

        def fmt(value):
            return f"{value:.4f}" if isinstance(value, (int, float)) else "-"

        print(f"{label:<8}{row['stage']:<7}{fmt(zero.get('mse')):>9}"
              f"{fmt(zero.get('ci')):>8}{fmt(zero.get('spearman')):>8}"
              f"{fmt(zero.get('calibration')):>8}{fmt(zero.get('shape')):>8}"
              f"{fmt(zero.get('cliff_sign')):>8}{fmt(five.get('cliff_sign')):>8}")

    print("\nk=0 Pareto frontier (MSE down, CI up) — whole configurations:")
    for label in frontier:
        zero = arms[label]["0"]
        print(f"  {label:<6} MSE {zero['mse']:.4f}  CI {zero['ci']:.4f}"
              f"  ({arms[label]['source']})")

    print(f"\ndouble-cold meta_test seal: {len(seals['sealed_explicit'])} explicit, "
          f"{len(seals['sealed_implicit'])} implicit (pre-R5), "
          f"{len(seals['process_unsealed'])} population-sealed but "
          f"process-unsealed (2026-08-16 incident), "
          f"{len(seals['split_undeclared'])} split-undeclared, "
          f"{len(seals['sealed_quarantined'])} quarantined "
          f"(computed pre-authorisation, moved to a sidecar, never read), "
          f"{len(seals['older_protocol'])} older-protocol artifacts "
          f"(bindingdb_ki_main_v0 — a different, consumed population), "
          f"{len(seals['violations'])} violation(s)")
    for violation in seals["violations"]:
        print(f"  {violation['artifact']}: {violation['problem']}")

    # `violations == 0` is not a clean bill of health while an incident stands.
    # Print the open states as their own banner so a reader who skims the
    # violation count cannot conclude the seal was never breached.
    print("\nseal isolation level for the artifacts audited here: LOGICAL "
          "EXCLUSION AFTER PARSING — they were produced against "
          "bindingdb_ki_double_cold_v1, a single all-label corpus, so every "
          "meta_test label is decompressed and parsed on every construction. "
          "This is not a physical label seal. A physically isolated surface "
          "now exists (scripts/build_governed_split_views.py -> "
          "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views, "
          "spec: tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md), but "
          "the trainers still load the all-label corpus, so no recorded "
          "artifact yet claims physical isolation.")
    if seals["process_unsealed"]:
        print(f"\n*** OPEN INCIDENT: {len(seals['process_unsealed'])} artifact(s) "
              f"were produced by a process that parsed AND INDEXED the sealed "
              f"split. Population sealed, process unsealed. `violations=0` "
              f"does not clear this. ***")
        for row in seals["process_unsealed"]:
            print(f"  {row['artifact']}  ->  {row['incident']}")
    if seals["split_undeclared"]:
        print(f"\nrecord defect: {len(seals['split_undeclared'])} artifact(s) "
              f"declare no corpus, so their protocol cannot be verified from "
              f"the file:")
        for name in seals["split_undeclared"]:
            print(f"  {name}")

    print(f"checkpoint sha256: {hashes['verified']} verified, "
          f"{len(hashes['mismatched'])} mismatched, "
          f"{hashes['payload_absent']} payload absent (gitignored)")
    for mismatch in hashes["mismatched"]:
        print(f"  {mismatch['artifact']}: recorded {mismatch['recorded'][:12]}… "
              f"actual {mismatch['actual'][:12]}…")

    if "skipped" in loading:
        print(f"strict checkpoint reload: skipped ({loading['skipped']})")
    else:
        print(f"strict checkpoint reload: {len(loading['loaded'])} loaded, "
              f"{len(loading['superseded_architecture'])} superseded-architecture "
              f"(expected), {len(loading['broken'])} broken")
        for row in loading["broken"]:
            print(f"  BROKEN {row['artifact']}: {row['error']}")

    if arguments.json is not None:
        arguments.json.write_text(json.dumps({
            "schema": "MetaSieve.RecordAudit.v2",
            "arms": arms,
            "k0_pareto_frontier": frontier,
            "meta_test_seal": seals,
            "seal_isolation": {
                "level": "logical_exclusion_after_parsing",
                "physically_isolated": False,
                "why_not": ("the audited artifacts were produced against "
                            "bindingdb_ki_double_cold_v1, a single all-label "
                            "cells.jsonl.gz"),
                "specification": ("tools/research/a2_readiness_v2/"
                                  "SPLIT_ISOLATION_SPEC.md"),
                "specification_status": "implemented, trainers not migrated",
                "isolated_surface": ("dataset/processed/meta_fewshot/"
                                     "bindingdb_ki_double_cold_v1_views"),
                "isolated_surface_builder": (
                    "scripts/build_governed_split_views.py"),
            },
            "open_incidents": {
                "process_unsealed": seals["process_unsealed"],
                "note": ("population sealed, process unsealed; "
                         "violations==0 does not clear this"),
            },
            "record_defects": {"split_undeclared": seals["split_undeclared"]},
            "checkpoint_hashes": hashes,
            "strict_loading": loading,
        }, indent=1), encoding="utf-8")
        print(f"\nwrote {arguments.json}")

    # Exit non-zero while the incident stands. A green audit must mean the
    # record is clean, not merely that no *new* violation was found.
    return 1 if (seals["violations"] or hashes["mismatched"]
                 or loading.get("broken") or seals["process_unsealed"]) else 0


if __name__ == "__main__":
    raise SystemExit(main())
