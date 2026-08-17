"""Stage V structural and leakage verification.

Stage V stopped at its V1 gate, so no neural arm exists to test. What must still
be verified is everything the census and the estimand rest on — above all the
repair itself: that the shared core is present in the transformation key *and*
in the model-facing descriptor, and that Stage T's key genuinely lacked it.
"""
from __future__ import annotations

import ast
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:  # pragma: no cover
    sys.path.insert(0, str(REPO))

from tools.research.stageT_mmp import mmp as stage_t_mmp  # noqa: E402
from tools.research.stageU_mmp_interaction import mmp as chem  # noqa: E402
from tools.research.stageV_core_mmp import core_mmp  # noqa: E402

STAGE = Path(__file__).resolve().parents[1]
SLOW = os.environ.get("RUN_SLOW") == "1"


def _one(left_smiles: str, right_smiles: str, module=chem):
    for left in module.fragment(left_smiles):
        for right in module.fragment(right_smiles):
            built = module.transformation(left, right)
            if built is not None:
                return built
    raise AssertionError("no matched pair")


# -- the repair itself ------------------------------------------------------


def test_stage_t_key_omitted_the_core_regression_pin():
    """Pins the forensic claim: Stage T's key is core-blind.

    If this ever fails, `CORRECTION_20260817_CORE_KEY.md` is describing code
    that no longer exists and must be revisited rather than silently trusted.
    """
    item, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1", module=stage_t_mmp)
    assert item.core not in item.exact_key
    assert item.exact_key == f"{item.context}|{item.r_a}>>{item.r_b}"


def test_core_inclusive_key_separates_the_same_edit_on_different_cores():
    """The central repair. Same nominal Cl->Br edit, two different scaffolds."""
    first, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    second, _ = _one("CCCCc1ccc(Cl)cc1", "CCCCc1ccc(Br)cc1")
    assert first.r_a == second.r_a and first.r_b == second.r_b
    assert first.context == second.context
    assert first.core != second.core
    # Stage T would have merged these; Stage V must not.
    stage_t_first, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1",
                            module=stage_t_mmp)
    stage_t_second, _ = _one("CCCCc1ccc(Cl)cc1", "CCCCc1ccc(Br)cc1",
                             module=stage_t_mmp)
    assert stage_t_first.exact_key == stage_t_second.exact_key
    assert first.exact_key != second.exact_key


def test_descriptor_consumes_the_core():
    first, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    second, _ = _one("CCCCc1ccc(Cl)cc1", "CCCCc1ccc(Br)cc1")
    assert chem.descriptor(first) != chem.descriptor(second)
    assert chem.edit_features(first) != chem.edit_features(second)
    assert len(chem.edit_features(first)) == chem.EDIT_WIDTH


def test_stage_t_descriptor_did_not_consume_the_core():
    first, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1", module=stage_t_mmp)
    second, _ = _one("CCCCc1ccc(Cl)cc1", "CCCCc1ccc(Br)cc1", module=stage_t_mmp)
    assert stage_t_mmp.descriptor(first) == stage_t_mmp.descriptor(second)


# -- chemistry contracts ----------------------------------------------------


def test_key_is_deterministic_and_stereo_sensitive():
    first, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    again, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    assert first.exact_key == again.exact_key
    stereo, _ = _one("C[C@H](N)c1ccccc1", "C[C@@H](N)c1ccccc1")
    assert stereo.stereo_edit is True
    assert stereo.exact_key != stereo.inverse.exact_key


def test_canonical_direction_and_inverse_negation():
    item, flipped = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    reverse, reverse_flipped = _one("CCc1ccc(Br)cc1", "CCc1ccc(Cl)cc1")
    assert item.exact_key == reverse.exact_key
    assert flipped != reverse_flipped
    assert item.r_a <= item.r_b
    assert item.inverse.inverse == item
    assert item.inverse.charge_change == -item.charge_change


def test_charge_change_is_recorded():
    item, _ = _one("CCc1ccc(CC(=O)O)cc1", "CCc1ccc(CC(=O)[O-])cc1")
    assert item.charge_change != 0


def test_attachment_environment_includes_hybridization():
    item, _ = _one("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    assert len(item.context) == 6
    assert isinstance(item.context[5], str)


# -- double-difference algebra ---------------------------------------------


def _effect(key, target, component, value):
    return core_mmp.TargetEffect(key, f"coarse::{key}", target, component,
                                 value, 1, 1, False, (0.0,))


def test_double_difference_cancels_level_and_generic_effect():
    mu = 1.75
    specific = {"t1": +0.4, "t2": -0.1}
    effects = [_effect("tau", t, f"c_{t}", mu + specific[t]) for t in ("t1", "t2")]
    rows = core_mmp.double_differences(effects)
    assert len(rows) == 1
    assert abs(rows[0].value - (specific["t1"] - specific["t2"])) < 1e-12


def test_double_differences_never_cross_keys():
    effects = [_effect("tauA", "t1", "c1", 1.0), _effect("tauB", "t2", "c2", 2.0)]
    assert core_mmp.double_differences(effects) == []


def test_double_difference_target_order_is_canonical():
    effects = [_effect("tau", "t2", "c2", 1.0), _effect("tau", "t1", "c1", 0.25)]
    rows = core_mmp.double_differences(effects)
    assert rows[0].target_left == "t1" and rows[0].target_right == "t2"
    assert abs(rows[0].value - (0.25 - 1.0)) < 1e-12


# -- leakage and hygiene ----------------------------------------------------


def _modules():
    for path in sorted(STAGE.rglob("*.py")):
        yield path, ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_no_python_hash_anywhere():
    for path, tree in _modules():
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "hash"):
                raise AssertionError(f"{path.name}:{node.lineno} calls hash()")


def test_every_rng_is_seeded_deterministically():
    offenders = []
    for path, tree in _modules():
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call)
                    and getattr(node.func, "attr", None) == "default_rng"):
                continue
            if not node.args:
                offenders.append(f"{path.name}:{node.lineno} unseeded")
                continue
            argument = node.args[0]
            ok = (isinstance(argument, ast.Call)
                  and getattr(argument.func, "id", None) == "stable_seed")
            ok = ok or isinstance(argument, (ast.Constant, ast.Name))
            if not ok:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, offenders


def test_stage_never_names_the_development_validation_split():
    forbidden = "meta" + "_val"
    offenders = []
    for path, tree in _modules():
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == forbidden:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, offenders


def test_mmp_construction_never_reads_a_label():
    path = Path(chem.__file__)
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and node.value in {"pK", "delta_y"}:
            raise AssertionError(f"{path.name}:{node.lineno} reads a label")
        if isinstance(node, ast.Attribute) and node.attr in {"pK", "delta_y"}:
            raise AssertionError(f"{path.name}:{node.lineno} reads a label")


# -- governed-data contracts ------------------------------------------------


@pytest.fixture(scope="module")
def governed():
    return core_mmp.load_governed()


def test_physical_meta_test_seal(governed):
    data, seal = governed
    assert seal["included"] is False
    assert seal["evaluated"] is False
    assert seal["isolation"]["physically_isolated"] is True
    assert "meta_test" not in data.tasks
    assert "meta_test" not in data.components


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to build the observation bank")
def test_observations_are_within_target_same_core_and_same_split(governed):
    from scripts.internal_validation import partition_components

    data, _seal = governed
    _fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    targets = sorted(t for t, c in component_of.items() if c in set(internal))
    rows = core_mmp.build_observations(data, targets)["observations"]
    assert len(rows) > 500
    for row in rows:
        assert data.cells[row.cell_a]["target_id"] == row.target
        assert data.cells[row.cell_b]["target_id"] == row.target
        assert data.cells[row.cell_a]["split"] == "meta_train"
        assert data.cells[row.cell_b]["split"] == "meta_train"
        expected = data.cells[row.cell_b]["pK"] - data.cells[row.cell_a]["pK"]
        assert abs(row.delta_y - expected) < 1e-12
        shared = bool(set(data.cells[row.cell_a]["panel_ids"])
                      & set(data.cells[row.cell_b]["panel_ids"]))
        assert shared == row.same_panel


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to check split isolation")
def test_no_double_difference_crosses_the_fit_internal_boundary(governed):
    from scripts.internal_validation import partition_components

    data, _seal = governed
    fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    fit_targets = sorted(t for t, c in component_of.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component_of.items()
                              if c in set(internal))
    rows = core_mmp.build_observations(
        data, fit_targets + internal_targets)["observations"]
    same_panel = [o for o in rows if o.same_panel]
    fit_set = set(fit_targets)
    for population in (set(fit_targets), set(internal_targets)):
        effects = core_mmp.target_effects(
            [o for o in same_panel if o.target in population])
        for row in core_mmp.double_differences(effects):
            assert (row.target_left in population
                    and row.target_right in population)
    # And the two populations share no protein component.
    assert not (set(fit) & set(internal))
    assert fit_set.isdisjoint(internal_targets)


DIGEST_SCRIPT = r"""
import hashlib, json, sys
sys.path.insert(0, r"{root}")
from scripts.internal_validation import partition_components
from tools.research.stageV_core_mmp import core_mmp
data, _ = core_mmp.load_governed()
_fit, internal = partition_components(data)
component = {{c["target_id"]: c["protein_group_40"] for c in data.cells}}
targets = sorted(t for t, c in component.items() if c in set(internal))
rows = core_mmp.build_observations(data, targets)["observations"]
payload = [(r.target, r.exact_key, r.cell_a, r.cell_b, round(r.delta_y, 12),
            r.same_panel, r.stratum) for r in rows]
print(hashlib.sha256(json.dumps(payload).encode()).hexdigest())
"""


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to run the subprocess test")
def test_bank_is_stable_across_pythonhashseed(tmp_path):
    script = tmp_path / "digest.py"
    script.write_text(DIGEST_SCRIPT.format(root=str(REPO)), encoding="utf-8")
    digests = []
    for salt in ("0", "1", "12345"):
        environment = dict(os.environ, PYTHONHASHSEED=salt)
        result = subprocess.run([sys.executable, str(script)], check=True,
                                capture_output=True, text=True,
                                env=environment, cwd=str(REPO))
        digests.append(result.stdout.strip())
    assert len(set(digests)) == 1, digests
