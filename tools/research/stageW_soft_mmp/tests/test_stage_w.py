"""Stage W W0 tests: soft-family construction and the frozen admission result."""
from __future__ import annotations
import ast, hashlib, json
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from tools.research.stageW_soft_mmp.soft_mmp import category, murcko_core

STAGE = Path(__file__).resolve().parents[1]
PREREG_SHA = "ae96762e319521f30aa09eb1a79fb8bb0e3ea324b21d4b40868aa6826a45dc71"


def test_preregistration_is_frozen():
    actual = hashlib.sha256((STAGE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    assert actual == PREREG_SHA


def test_w0_gates_davis_fail_kiba_pass():
    d = json.loads((STAGE / "W0_SOFT_MMP_CENSUS.json").read_text(encoding="utf-8"))
    assert d["preregistration_sha256"] == PREREG_SHA
    assert d["datasets"]["davis"]["all_pass"] is False
    assert d["datasets"]["kiba"]["all_pass"] is True
    davis = d["datasets"]["davis"]["census"]
    assert davis["rich_families"] == 7
    kiba = d["datasets"]["kiba"]["census"]
    assert kiba["rich_families"] >= 20
    assert d["datasets"]["kiba"]["same_core_residual"]["median"] <= 1.0


def test_soft_family_helpers_are_structure_only():
    assert murcko_core("Clc1ccc(C[*:1])cc1").startswith("c1cc")
    cat = category("C[*:1]")
    assert cat == (1, 0, 0, 0, 0, 1)
    charged = category("[O-][*:1]")
    assert charged[-1] == 0


def test_stage_w_source_has_no_hash_or_forbidden_split_name():
    forbidden = "meta" + "_val"
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "hash":
                raise AssertionError(f"{path.name}:{node.lineno} hash()")
            if isinstance(node, ast.Constant) and node.value == forbidden:
                raise AssertionError(f"{path.name}:{node.lineno} forbidden split")


def test_w1_preregistration_and_split_admission_pass():
    prereg_sha = "038f4d97f74841023c48a2e9b3bab5592a0bad2bb9fa54a464d5290641549082"
    actual = hashlib.sha256(
        (STAGE / "W1_PREREGISTRATION.md").read_bytes()).hexdigest()
    assert actual == prereg_sha
    d = json.loads((STAGE / "W1_DATA.json").read_text(encoding="utf-8"))
    assert d["w1_preregistration_sha256"] == prereg_sha
    assert d["all_pass"] is True
    assert d["counts"]["heldout_repeated_rows"] >= 500
    assert d["counts"]["heldout_repeated_components"] >= 10
    assert d["counts"]["repeated_families"] >= 50
