"""Stage W0b audit tests."""
from __future__ import annotations
import ast, hashlib, json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

STAGE = Path(__file__).resolve().parents[1]
PREREG_SHA = "ff23c408d20cc79b1bd5fcd20854a0443280d32d6fc3dbb8abf0733a9a70631f"


def test_preregistration_frozen():
    actual = hashlib.sha256((STAGE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    assert actual == PREREG_SHA


def test_audit_records_no_go_and_censoring():
    d = json.loads((STAGE / "W0B_AUDIT.json").read_text(encoding="utf-8"))
    assert d["preregistration_sha256"] == PREREG_SHA
    assert d["go_no_go"]["w1_biological_interpretation"] == "NO-GO"
    assert d["positive_control_audit"]["runnable"] is False
    assert d["censoring"]["davis"]["censored_fraction"] > 0.5
    assert d["censoring"]["metz"]["censored_fraction"] > 0.5
    assert d["censoring"]["klaeger"]["censored_fraction"] > 0.9
    assert d["stage_w_audit"]["w1_training_status"].startswith("PAUSED")


def test_decision_artifact():
    d = json.loads((STAGE / "W0B_DECISION.json").read_text(encoding="utf-8"))
    assert d["decision"].startswith("NO-GO")


def test_no_hash_or_forbidden_split_name_in_source():
    forbidden = "meta" + "_val"
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "hash":
                raise AssertionError(f"{path.name}:{node.lineno} hash()")
            if isinstance(node, ast.Constant) and node.value == forbidden:
                raise AssertionError(f"{path.name}:{node.lineno}")

def test_final_decision_no_go():
    import json
    d=json.loads((STAGE/'W0B_W0P_FINAL_DECISION.json').read_text(encoding='utf-8'))
    assert d['decision'].startswith('NO-GO')
    assert d['evidence']['w0p_result']['verdict'].startswith('FAIL')

def test_censored_recensus_updates_support_map():
    d=json.loads((STAGE/'W0B_CENSORED_RECENSUS.json').read_text(encoding='utf-8'))
    davis=d['datasets']['davis']['layers']['all_pairs']
    assert davis['cross_component_D_rows'] > 500
    assert d['datasets']['davis']['components'] >= 10
    metz=d['datasets']['metz']['layers']['all_pairs']
    assert metz['cross_component_D_rows'] > 500
    assert d['datasets']['metz']['components'] >= 10

def test_crossplatform_is_transfer_gate_only():
    d=json.loads((STAGE/'W0B_CROSSPLATFORM.json').read_text(encoding='utf-8'))
    mk=d['platforms']['metz_vs_klaeger']
    assert mk['spearman'] > 0.4
    md=d['platforms']['metz_vs_davis']
    assert md['spearman'] < 0.0
    assert md['interpretation'] == 'cross-platform transfer gate only'
