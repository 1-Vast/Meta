"""Keep the narrative research record consistent with its leaf artifacts.

These are cheap, no-training checks. They exist because the R13 gate-count
inconsistency (`15/16 passed` alongside two recorded xfails) and the mixed
k=0 "frontier" both survived several consolidation passes: a number that no
test recomputes will drift away from its evidence.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.audit_research_record import (
    DOUBLE_COLD, MAIN_V0, check_seals, check_strict_loading, collect_arms,
    pareto,
)

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "report"
BOUNDARY = REPORT / "BOUNDARY_20260816.md"


def test_double_cold_meta_test_is_never_reported_as_evaluated():
    """The confirmation split opens once, and it has not opened yet."""
    seals = check_seals()
    assert seals["violations"] == [], seals["violations"]
    assert seals["sealed_explicit"], "no explicitly sealed double-cold artifact found"


def test_an_artifact_is_only_called_older_protocol_if_it_says_so():
    """The 2026-08-16 misclassification: silence was read as a consumed split.

    Every artifact that declares neither corpus was being reported as
    `bindingdb_ki_main_v0` older-protocol evidence. All of them were in fact
    double-cold stage summaries, so the audit was asserting a protocol the
    files never claimed. Undeclared is now its own state.
    """
    seals = check_seals()
    for name in seals["older_protocol"]:
        assert MAIN_V0 in (ROOT / name).read_text(encoding="utf-8"), (
            f"{name} is called older-protocol but never names {MAIN_V0}")
    for name in seals["split_undeclared"]:
        text = (ROOT / name).read_text(encoding="utf-8")
        assert MAIN_V0 not in text and DOUBLE_COLD not in text, (
            f"{name} declares a protocol and should not be undeclared")


def test_both_audit_invocations_reach_the_audit_and_fail_for_the_right_reason():
    """A nonzero exit must mean a finding, not a broken import.

    `python scripts/audit_research_record.py` puts `scripts/` on `sys.path`,
    not the repository root, so the lazy
    `from scripts.stageR6_compare_arms import load_arm` inside
    `check_strict_loading` raised `ModuleNotFoundError` and the process exited
    1 — the same code the open governance incident produces. From the shell the
    two were indistinguishable.

    Both forms must now run to completion, print the incident banner, and exit
    1 because of it.
    """
    import subprocess
    import sys as _sys
    for command in (["-m", "scripts.audit_research_record"],
                    ["scripts/audit_research_record.py"]):
        # Decode explicitly as UTF-8: the audit prints em dashes, and on a
        # CJK-locale Windows console `text=True` would decode the pipe as
        # cp936 and drop the captured output entirely, failing this test for a
        # reason that has nothing to do with the record.
        finished = subprocess.run(
            [_sys.executable, *command], cwd=str(ROOT),
            capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=900)
        assert "OPEN INCIDENT" in finished.stdout, (
            f"{command} did not reach the seal report:\n"
            f"{finished.stderr[-2000:]}")
        assert "Traceback" not in finished.stderr, finished.stderr[-2000:]
        assert "strict checkpoint reload:" in finished.stdout, (
            f"{command} skipped the loading check:\n{finished.stderr[-2000:]}")
        assert finished.returncode == 1, (
            f"{command} exited {finished.returncode}; expected 1 from the "
            "open governance incident")


def test_no_recorded_artifact_asserts_a_physical_label_seal():
    """The terminology guard, at repository scale.

    The lineage has never had physical label isolation: `cells.jsonl.gz` is a
    single all-label artifact, so every sealed label is parsed on every load.
    An artifact whose `meta_test.seal` begins "physical:" is making the same
    class of overstatement as the 2026-08-16 incident, even when the flag it
    names was genuinely passed.

    Historical strings are preserved inside `seal_relabel.was` as provenance;
    only the live `seal` field is checked.
    """
    offenders = []
    for path in (ROOT / "report").rglob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        block = payload.get("meta_test")
        if isinstance(block, dict) and str(block.get("seal", "")).startswith(
                "physical:"):
            offenders.append(str(path.relative_to(ROOT)))
    assert offenders == [], (
        "artifacts claim a physical label seal the lineage does not have: "
        f"{offenders}")


def test_the_split_isolation_specification_describes_the_implementation():
    """The honest state: built, and not retroactive.

    Until 2026-08-18 this test asserted the opposite — that the specification
    was a proposal — and required it to be rewritten once an isolated artifact
    appeared. It has appeared: `scripts/build_governed_split_views.py` produced
    the development surface and `QPSMPData(split_view=…)` mounts it. What the
    test now guards is the half that must never drift: the implementation is
    the surface a *future* run may mount, and no recorded artifact was produced
    on it.
    """
    spec = ROOT / "tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md"
    assert spec.exists()
    text = spec.read_text(encoding="utf-8")
    assert "IMPLEMENTED" in text, "the specification still reads as a proposal"
    assert "not retroactive" in text.lower(), (
        "the specification must state that recorded artifacts were produced "
        "on the weaker all-label surface")
    builder = ROOT / "scripts/build_governed_split_views.py"
    assert builder.is_file(), "the specification claims a builder that is absent"
    for contract in ("tools/tests/test_governed_split_views.py",
                     "tools/tests/test_physical_meta_test_seal.py"):
        assert (ROOT / contract).is_file(), contract
        assert contract in text, f"{contract} is not cited by the specification"


def test_process_unsealed_artifacts_carry_their_correction():
    """A corrected artifact must keep its incident pointer, not just a new string.

    These are the R14 artifacts whose producer omitted `include_meta_test`:
    the *population* was sealed (no meta_test value entered any metric, verified
    by bit-identical re-run) but the *process* was not. Downgrading them to
    plain `sealed_explicit` would overstate the seal.
    """
    seals = check_seals()
    for row in seals["process_unsealed"]:
        assert row["incident"], row
        assert (ROOT / row["incident"]).exists(), (
            f"{row['artifact']} points at a missing incident report")
        assert row["numerical_impact"] == "none; no field in this artifact was modified"


def test_k0_frontier_in_the_boundary_document_is_the_real_pareto_set():
    """A frontier is a set of whole configurations, not a mix of metrics."""
    arms = collect_arms()
    front = pareto(arms)
    assert front == ["B3", "C2", "A0"], front

    text = BOUNDARY.read_text(encoding="utf-8")
    for label in front:
        mse = arms[label]["0"]["mse"]
        ci = arms[label]["0"]["ci"]
        assert f"{mse:.3f}" in text, f"{label} MSE {mse:.3f} absent from BOUNDARY"
        assert f"{ci:.3f}" in text, f"{label} CI {ci:.3f} absent from BOUNDARY"

    # Every arm on the frontier must be beaten on the other metric by another
    # frontier member; that is what makes the trade real rather than rhetorical.
    best_mse = min(arms[label]["0"]["mse"] for label in front)
    best_ci = max(arms[label]["0"]["ci"] for label in front)
    assert not any(arms[label]["0"]["mse"] == best_mse
                   and arms[label]["0"]["ci"] == best_ci for label in front), (
        "one arm dominates both metrics; the frontier language must be rewritten")


def test_best_cliff_record_is_scoped_to_meta_val_development():
    """0.782 is a development record on a Pareto-dominated arm; say so."""
    arms = collect_arms()
    best = max(arms.items(), key=lambda kv: kv[1].get("5", {}).get("cliff_sign") or 0)
    label, row = best
    assert label == "C1" and abs(row["5"]["cliff_sign"] - 0.782) < 5e-4, best

    assert label not in pareto(arms), (
        "C1 is no longer Pareto-dominated; the BOUNDARY caveat needs rewriting")

    text = BOUNDARY.read_text(encoding="utf-8")
    heading = "## The activity-cliff record, stated at its true scope"
    assert heading in text, "BOUNDARY has no scoped activity-cliff section"
    section = text[text.index(heading):]
    section = section[:section.index("\n## ", 1)]
    assert "0.782" in section, "the cliff record section omits the number"
    assert "meta_val" in section, "the cliff record is not scoped to meta_val"
    assert "development" in section, "the cliff record is not graded development"
    assert "dominated" in section, "the cliff record omits its dominated status"


def test_every_comparator_checkpoint_still_reloads_strictly():
    """A comparator you cannot reload is not a comparator.

    Six R3R4 pre-fix arms are expected to fail: their `TypedLigandChannels`
    pooled each ligand to five vectors before any protein contact, and the
    documented fix replaced it with 16 query-slot tokens. They are retained as
    evidence for the identifiability and capacity defects, not as comparators.
    Any *other* failure means an architecture change silently orphaned a
    checkpoint the record still cites.
    """
    loading = check_strict_loading()
    if "skipped" in loading:
        pytest.skip(loading["skipped"])
    assert loading["broken"] == [], loading["broken"]
    assert loading["loaded"], "no checkpoint was reloaded at all"

    superseded = {row["artifact"] for row in loading["superseded_architecture"]}
    assert all("stageR3R4" in name for name in superseded), superseded

    # The three frontier arms must be among the checkpoints that do reload.
    reloadable = {row["artifact"] for row in loading["loaded"]}
    for arm in ("B3_full", "C2_cliffw2", "A0_incumbent"):
        assert any(arm in name for name in reloadable), (
            f"{arm} is on the k=0 Pareto frontier but no checkpoint reloads")


@pytest.mark.parametrize("path,expected", [
    ("docs/PROJECT_FILE_ORGANIZATION.md", "stageR0"),
    ("report/BOUNDARY_20260816.md", "no fitting, selection or reported metric"),
])
def test_documents_carry_the_current_cycle_scope(path, expected):
    assert expected in (ROOT / path).read_text(encoding="utf-8")


CLOSING_AUTHORITIES = [
    "report/POST_COMPLETION_REVIEW_20260818.md",
    "report/BOUNDARY_20260817_NIGHT.md",
    "report/FINAL_STATE_20260818.md",
    "report/COMPLETION_STATEMENT_20260818.md",
    "report/CURRENT_MODEL_EVIDENCE.md",
    "report/EVIDENCE_LEDGER.md",
    "task.md",
    "history.md",
]


@pytest.mark.parametrize("document", [
    "report/BOUNDARY_20260816.md",
    "report/CURRENT_MODEL_EVIDENCE.md",
    "report/EVIDENCE_LEDGER.md",
    "task.md",
    "history.md",
])
def test_authorities_do_not_call_meta_test_unopened_or_physically_sealed(document):
    """Two claims the evidence does not support, in the authorities that matter.

    `meta_test` was parsed and indexed by six analysis scripts. The supportable
    statement is that no sealed label entered any fitting, selection or
    reported metric, and that a process-isolation incident remains open.
    "Never opened" and "physically sealed" both overstate it.

    Prose *about* the prohibition is allowed; the banned phrasing is only a
    violation when it is used as an assertion, so a line that also contains a
    negation or the word "not" is exempt.
    """
    banned = ("never been opened", "never opened", "unopened",
              "physically sealed", "physically excluded")
    offending = []
    for number, line in enumerate(
            (ROOT / document).read_text(encoding="utf-8").splitlines(), 1):
        lowered = line.lower()
        if any(phrase in lowered for phrase in banned) and not any(
                marker in lowered for marker in
                ("not ", "never write", "do not", "n't", "overstate",
                 "instead of", "rather than")):
            offending.append(f"{document}:{number}: {line.strip()[:100]}")
    assert offending == [], offending


# --- the post-completion narrowing, pinned so it cannot drift back --------

@pytest.mark.parametrize("document", CLOSING_AUTHORITIES)
def test_no_authority_turns_a_measured_failure_into_a_bound(document):
    """Four overstatements the post-completion review retired.

    Each is a claim the experiments cannot support: a measured probe result
    stated as a ceiling, an unreached target stated as impossible, a
    trade-off seen four times stated as a property of every architecture, and
    a proxy experiment stated as a falsification of the method it stood in
    for. The measurements behind all four are unchanged and still recorded;
    only the quantifier was wrong.

    As with the seal guard above, prose *about* the overstatement is allowed:
    a passage that negates, quotes or corrects the phrase is exempt. The scan
    is per paragraph rather than per line, because these documents are
    hard-wrapped and a negation routinely lands on the line above its claim.
    """
    banned = (
        "at most 26%", "at most ~26%", "<=26%", "at most 26 %",
        "is not achievable", "is not reachable", "remains unreachable",
        "mathematically impossible",
        "fundamental to single-stage", "fundamental to all",
    )
    exempt = ("not ", "never write", "do not", "n't", "overstate",
              "instead of", "rather than", "originally read", "corrected",
              "no longer", "superseded", "was wrong", "retired", "banned",
              "replace", "must become", "overstatement")
    text = (ROOT / document).read_text(encoding="utf-8")
    offending = []
    for paragraph in text.split("\n\n"):
        lowered = paragraph.lower()
        if any(marker in lowered for marker in exempt):
            continue
        for phrase in banned:
            if phrase in lowered:
                offending.append(
                    f"{document}: {phrase!r} in: "
                    f"{' '.join(paragraph.split())[:140]}")
    assert offending == [], offending


@pytest.mark.parametrize("document", [
    "report/BOUNDARY_20260817_NIGHT.md",
    "report/FINAL_STATE_20260818.md",
    "report/COMPLETION_STATEMENT_20260818.md",
    "task.md",
])
def test_closing_authorities_record_the_oracle_level_floor(document):
    """The decomposition makes MSE <= 1.00 arithmetically possible; say so.

    Omitting it is how "no tested candidate reached 1.00" silently becomes
    "1.00 cannot be reached". The centered term at k=0 is 0.8648, so an oracle
    level predictor lands near 0.865 — below the target.
    """
    text = (ROOT / document).read_text(encoding="utf-8")
    assert "0.865" in text, (
        f"{document} omits the oracle-level k=0 floor (~0.865)")
    assert "arithmetically possible" in text, (
        f"{document} does not state that the target is arithmetically "
        "possible")


def test_the_closure_map_separates_proxies_from_direct_implementations():
    """A named method is only falsified if its defining operator was built."""
    text = (ROOT / "tools/research/method_ladder/CLOSURE_MAP.md").read_text(
        encoding="utf-8")
    assert "proxy negative; direct method not instantiated" in text
    for method in ("OGM", "Gradient Blending", "Disentangled Gradient "
                   "Learning", "Set Transformer", "DrugBAN", "FS-CAP"):
        assert method in text, method
    # No *table row* may claim measurement closure. The phrase may still
    # appear in the prose that explains why the earlier version was wrong.
    rows = [line for line in text.splitlines()
            if line.startswith("|") and "closed by measurement" in line.lower()]
    assert rows == [], (
        "the closure map still claims measurement closure for families whose "
        f"defining operator was never implemented: {rows}")
