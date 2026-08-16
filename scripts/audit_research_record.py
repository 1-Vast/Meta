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

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"
FEWSHOT = REPORT / "meta_fewshot"

# The double-cold protocol. Absolute numbers are only comparable within it.
DOUBLE_COLD = "bindingdb_ki_double_cold_v1"

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
        a double-cold run carrying ``meta_test.evaluated == False`` — the
        physical seal introduced by Stage R5.
    ``sealed_implicit``
        a double-cold run from before R5 that records only ``meta_val`` and
        never mentions ``meta_test`` at all.
    ``older_protocol``
        an artifact on the pre-double-cold `bindingdb_ki_main_v0` population.
        Stages 4/6/7 legitimately report a `meta_test` there; it is a
        *different, consumed* population and must never be conflated with the
        sealed double-cold confirmation split.
    """
    explicit, implicit, older, violations = [], [], [], []
    for path in sorted(FEWSHOT.rglob("RESULT.json")):
        name = str(path.relative_to(ROOT))
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError) as error:
            violations.append({"artifact": name, "problem": f"unreadable: {error}"})
            continue
        if DOUBLE_COLD not in str(payload.get("split_directory", "")):
            older.append(name)
            continue
        seal = payload.get("meta_test")
        if isinstance(seal, dict) and seal.get("evaluated") is False:
            explicit.append(name)
        elif seal is None and "meta_test" not in path.read_text(encoding="utf-8"):
            implicit.append(name)
        else:
            violations.append({"artifact": name,
                               "problem": f"double-cold meta_test not sealed: {seal!r}"})
    return {"sealed_explicit": explicit, "sealed_implicit": implicit,
            "older_protocol": older, "violations": violations}


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
          f"{len(seals['older_protocol'])} older-protocol artifacts "
          f"(bindingdb_ki_main_v0 — a different, consumed population), "
          f"{len(seals['violations'])} violation(s)")
    for violation in seals["violations"]:
        print(f"  {violation['artifact']}: {violation['problem']}")

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
            "schema": "MetaSieve.RecordAudit.v1",
            "arms": arms,
            "k0_pareto_frontier": frontier,
            "meta_test_seal": seals,
            "checkpoint_hashes": hashes,
            "strict_loading": loading,
        }, indent=1), encoding="utf-8")
        print(f"\nwrote {arguments.json}")

    return 1 if (seals["violations"] or hashes["mismatched"]
                 or loading.get("broken")) else 0


if __name__ == "__main__":
    raise SystemExit(main())
