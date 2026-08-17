"""Stage U structural verification.

Because U0's frozen admission gate failed, no U1/U2 training was authorized and
the model-contract tests (antisymmetry/identity/cycle, gradient coverage) are
recorded as not applicable rather than pretended. Every test that guards the
measured negative decision is here: the MMP decomposition, canonical direction,
provenance seal, deterministic banks, PYTHONHASHSEED independence, and no label
path into construction or splits.
"""
from __future__ import annotations

import ast
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import torch

REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:  # pragma: no cover
    sys.path.insert(0, str(REPO))

from tools.research.stageU_mmp_interaction import mmp as mmp_module  # noqa: E402
from tools.research.stageU_mmp_interaction.mmp import (  # noqa: E402
    descriptor, edit_features, fragment, strip_stereochemistry,
    transformation,
)
from tools.research.stageU_mmp_interaction.observation_cache import (  # noqa: E402
    load_observations,
)

STAGE = Path(__file__).resolve().parents[1]
SLOW = os.environ.get("RUN_SLOW") == "1"

BENZENE_CL = "CCc1ccc(Cl)cc1"
BENZENE_BR = "CCc1ccc(Br)cc1"


# -- 1. deterministic, non-empty MMP decomposition --------------------------

def test_fragmentation_is_non_empty_positive_control():
    pieces = fragment(BENZENE_CL)
    assert len(pieces) >= 3, pieces
    assert all("[*:1]" in piece.core and "[*:1]" in piece.r_group
               for piece in pieces)


def test_fragmentation_is_deterministic_and_order_stable():
    first = fragment(BENZENE_CL)
    fragment.cache_clear()
    second = fragment(BENZENE_CL)
    assert first == second


def test_attachment_context_returns_plain_data_not_rdkit_objects():
    context = mmp_module._attachment_context("Clc1ccc(C[*:1])cc1")
    assert isinstance(context, tuple)
    for value in context:
        assert isinstance(value, (str, bool, int)), type(value)


def test_core_is_the_larger_fragment():
    for piece in fragment(BENZENE_CL):
        assert (mmp_module._heavy_atoms(piece.core)
                >= mmp_module._heavy_atoms(piece.r_group))


# -- 2. inverse transformation and sign consistency -------------------------

def _one_transformation(left_smiles: str, right_smiles: str):
    for left in fragment(left_smiles):
        for right in fragment(right_smiles):
            built = transformation(left, right)
            if built is not None:
                return built
    raise AssertionError("no matched pair")


def test_canonical_direction_is_a_function_of_structure_only():
    item, flipped = _one_transformation(BENZENE_CL, BENZENE_BR)
    reverse, reverse_flipped = _one_transformation(BENZENE_BR, BENZENE_CL)
    assert item.exact_key == reverse.exact_key
    assert flipped != reverse_flipped
    assert item.r_a <= item.r_b


def test_inverse_transformation_negates_the_key_direction():
    item, _ = _one_transformation(BENZENE_CL, BENZENE_BR)
    assert item.inverse.r_a == item.r_b
    assert item.inverse.r_b == item.r_a
    assert item.inverse.charge_change == -item.charge_change
    assert item.inverse.inverse == item
    assert item.inverse.exact_key != item.exact_key


def test_exact_key_includes_the_shared_core():
    aromatic, _ = _one_transformation(BENZENE_CL, BENZENE_BR)
    aliphatic, _ = _one_transformation("CCCCCl", "CCCCBr")
    assert aromatic.core != aliphatic.core
    assert aromatic.context != aliphatic.context
    assert aromatic.exact_key != aliphatic.exact_key


def test_attachment_stereo_and_charge_are_preserved():
    item, _ = _one_transformation("C[C@H](N)c1ccccc1", "C[C@@H](N)c1ccccc1")
    assert item.stereo_edit is True
    assert strip_stereochemistry(item.r_a) == strip_stereochemistry(item.r_b)
    assert item.r_a != item.r_b

    charged, _ = _one_transformation("CCc1ccc(CC(=O)O)cc1",
                                     "CCc1ccc(CC(=O)[O-])cc1")
    assert charged.charge_change != 0


def test_descriptor_and_edit_features_are_finite_and_fixed_width():
    item, _ = _one_transformation(BENZENE_CL, BENZENE_BR)
    values = descriptor(item)
    features = edit_features(item)
    assert len(values) == mmp_module.DESCRIPTOR_WIDTH
    assert len(features) == mmp_module.EDIT_WIDTH
    assert all(np.isfinite(values))
    assert all(np.isfinite(features))


# -- 3. no label path into construction ------------------------------------

LABEL_BLIND_MODULES = ("mmp.py",)


def test_mmp_construction_never_reads_a_label():
    for name in LABEL_BLIND_MODULES:
        path = STAGE / name
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value in {"pK", "delta_y"}:
                raise AssertionError(f"{path.name}:{node.lineno} reads a label")
            if isinstance(node, ast.Attribute) and node.attr in {"pK", "delta_y"}:
                raise AssertionError(f"{path.name}:{node.lineno} reads a label")


def test_coverage_function_never_reads_a_query_label():
    path = STAGE / "u0_census.py"
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    found = False
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef)
                and node.name in {"_pairs_form_mmp", "_coarse_relation"}):
            continue
        found = True
        for child in ast.walk(node):
            if isinstance(child, ast.Attribute) and child.attr in {"pK", "delta_y"}:
                raise AssertionError(f"{path.name}:{child.lineno} reads a label")
            if isinstance(child, ast.Constant) and child.value in {"pK", "delta_y"}:
                raise AssertionError(f"{path.name}:{child.lineno} reads a label")
    assert found


def test_no_python_hash_anywhere_in_the_stage():
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "hash"):
                raise AssertionError(f"{path.name}:{node.lineno} calls hash()")


def test_every_default_rng_is_seeded_by_stable_seed_or_a_constant():
    offenders = []
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
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
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == forbidden:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, offenders


# -- 4. frozen artifacts ----------------------------------------------------

def test_preregistration_is_frozen_and_matches_recorded_sha256():
    import hashlib
    actual = hashlib.sha256(
        (STAGE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    assert actual == ("fdc0a830aa92882d07b9aea50f22a4c72fc6d93f92c55a3be6bc15cd6a645c11")
    census = json.loads((STAGE / "U0_CENSUS.json").read_text(encoding="utf-8"))
    assert census["preregistration_sha256"] == actual
    assert census["admission"]["all_pass"] is False


def test_u0_reliability_never_claims_an_mse_floor():
    report = json.loads((STAGE / "U0_RELIABILITY.json").read_text(
        encoding="utf-8"))
    assert "not_a_claim" in report
    assert "NOT" in report["not_a_claim"]


# -- 5. governed-data contracts ---------------------------------------------

@pytest.fixture(scope="module")
def governed():
    from tools.research.stageU_mmp_interaction.observations import load_governed
    return load_governed()


def test_physical_meta_test_seal(governed):
    data, seal = governed
    assert seal["included"] is False
    assert seal["evaluated"] is False
    assert seal["isolation"]["physically_isolated"] is True
    assert "meta_test" not in data.tasks
    assert "meta_test" not in data.components


def test_provenance_cache_contains_only_meta_train_rows(governed):
    from tools.research.stageU_mmp_interaction.provenance import (
        load_cache, train_allow_list,
    )
    allow = train_allow_list()
    rows = load_cache()
    assert rows, "run u0_reliability first"
    assert {row["source_row_id"] for row in rows} <= set(allow.row_ids)
    assert all(row["endpoint"] == "Ki" for row in rows)


def test_observation_cache_rows_are_within_target_and_same_core(governed):
    data, _seal = governed
    rows = load_observations()
    assert len(rows) > 1000
    for row in rows:
        assert data.cells[row.cell_a]["target_id"] == row.target
        assert data.cells[row.cell_b]["target_id"] == row.target
        assert data.cells[row.cell_a]["split"] == "meta_train"
        assert row.ligand_a != row.ligand_b
        expected = data.cells[row.cell_b]["pK"] - data.cells[row.cell_a]["pK"]
        assert abs(row.delta_y - expected) < 1e-12
        shared = bool(set(data.cells[row.cell_a]["panel_ids"])
                      & set(data.cells[row.cell_b]["panel_ids"]))
        assert shared == row.same_panel


def test_no_cross_population_observation_in_cached_bank(governed):
    from scripts.internal_validation import partition_components
    data, _seal = governed
    fit, internal = partition_components(data)
    fit_set, internal_set = set(fit), set(internal)
    for row in load_observations():
        assert (row.component in fit_set) != (row.component in internal_set)


BANK_DIGEST_SCRIPT = r"""
import hashlib, json, sys
sys.path.insert(0, r"{root}")
from tools.research.stageU_mmp_interaction.observation_cache import load_observations
rows = load_observations()
payload = [(r.target, r.exact_key, r.cell_a, r.cell_b, round(r.delta_y, 12),
            r.same_panel, r.stratum) for r in rows]
print(hashlib.sha256(json.dumps(payload).encode()).hexdigest())
"""


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 for the subprocess test")
def test_observation_bank_is_stable_across_pythonhashseed(tmp_path):
    script = tmp_path / "digest.py"
    script.write_text(BANK_DIGEST_SCRIPT.format(root=str(REPO)), encoding="utf-8")
    digests = []
    for salt in ("0", "1", "12345"):
        environment = dict(os.environ, PYTHONHASHSEED=salt)
        result = subprocess.run([sys.executable, str(script)], check=True,
                                capture_output=True, text=True,
                                env=environment, cwd=str(REPO))
        digests.append(result.stdout.strip())
    assert len(set(digests)) == 1, digests


# -- 6. the U2 model contract was gated and is therefore not instantiated ----

def test_no_neural_model_was_trained_before_the_u0_gate():
    runs = STAGE / "runs"
    trained = [p for p in runs.rglob("RUN.json")] if runs.exists() else []
    assert trained == [], "U0 failed: no arm may be trained"
    assert not (STAGE / "U1_VARIANCE.json").exists()
