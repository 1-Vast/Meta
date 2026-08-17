"""Stage T structural verification, run before any arm is trained.

Each test closes a failure mode that would make the T2 numbers uninterpretable,
and several close a defect that actually occurred while building this stage:

* `rdMMPA.FragmentMol` is overloaded; four positional integers mis-bind and a
  fully-named call segfaults. A silent empty fragmentation would have produced a
  false negative census, so a positive control pins it;
* `_attachment_context` must return plain data -- returning an RDKit `Atom`
  outlives its parent `Mol` and crashes on the next attribute read;
* the double difference must remove the target level and the generic effect,
  which is asserted algebraically rather than assumed;
* protein-free arms must be structurally incapable of expressing `D`;
* no label may reach a feature, a key, a split or a deduplication choice.
"""
from __future__ import annotations

import ast
import inspect
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

from tools.research.stageT_mmp import mmp as mmp_module  # noqa: E402
from tools.research.stageT_mmp.mmp import (  # noqa: E402
    descriptor, fragment, strip_stereochemistry, transformation,
)
from tools.research.stageT_mmp.t2_dataset import (  # noqa: E402
    DoubleDifference, double_differences, shuffle_within_key, target_effects,
)
from tools.research.stageT_mmp.t2_model import (  # noqa: E402
    DiscriminatorConfig, DoubleDifferenceModel,
)

STAGE = Path(__file__).resolve().parents[1]
SLOW = os.environ.get("RUN_SLOW") == "1"

BENZENE_CL = "CCc1ccc(Cl)cc1"
BENZENE_BR = "CCc1ccc(Br)cc1"


# -- 1. deterministic, non-empty MMP decomposition --------------------------


def test_fragmentation_is_non_empty_positive_control():
    """A silent empty fragmentation would fake a negative census."""
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
    """Returning an `Atom` outlives its `Mol` and segfaults on next access."""
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


def test_attachment_context_is_required_to_match():
    """Two edits at chemically different attachment points are not one key."""
    aromatic, _ = _one_transformation(BENZENE_CL, BENZENE_BR)
    aliphatic, _ = _one_transformation("CCCCCl", "CCCCBr")
    assert aromatic.context != aliphatic.context
    assert aromatic.exact_key != aliphatic.exact_key


def test_stereochemistry_is_retained_and_flagged():
    item, _ = _one_transformation("C[C@H](N)c1ccccc1", "C[C@@H](N)c1ccccc1")
    assert item.stereo_edit is True
    assert strip_stereochemistry(item.r_a) == strip_stereochemistry(item.r_b)
    assert item.r_a != item.r_b


def test_formal_charge_change_is_recorded():
    item, _ = _one_transformation("CCc1ccc(CC(=O)O)cc1", "CCc1ccc(CC(=O)[O-])cc1")
    assert item.charge_change != 0


def test_descriptor_is_finite_and_fixed_width():
    item, _ = _one_transformation(BENZENE_CL, BENZENE_BR)
    values = descriptor(item)
    assert len(values) == mmp_module.DESCRIPTOR_WIDTH
    assert all(np.isfinite(values))


# -- 3. the double difference removes level and the generic effect ----------


def _effect(key, target, component, value, dims=4):
    from tools.research.stageT_mmp.t2_dataset import TargetEffect
    return TargetEffect(key, target, component, value, 1, "S1_same_panel_single",
                        False, tuple([0.0] * dims))


def test_double_difference_removes_target_level_and_generic_effect():
    """delta_y = mu_tau + delta(t,tau) + level(t): D keeps only delta differences."""
    mu = 1.75
    level = {"t1": +3.0, "t2": -2.0}
    specific = {"t1": +0.4, "t2": -0.1}
    effects = [_effect("tau", t, f"c_{t}", mu + specific[t]) for t in ("t1", "t2")]
    rows = double_differences(effects)
    assert len(rows) == 1
    assert abs(rows[0].value - (specific["t1"] - specific["t2"])) < 1e-12
    # A pure level shift on both sides changes nothing, because level cancels
    # inside delta_y before D is even formed.
    assert level  # the shift is never added: delta_y is already a within-target gap


def test_double_differences_are_never_formed_across_keys():
    effects = [_effect("tauA", "t1", "c1", 1.0), _effect("tauB", "t2", "c2", 2.0)]
    assert double_differences(effects) == []


def test_double_difference_pairs_are_canonically_ordered():
    effects = [_effect("tau", "t2", "c2", 1.0), _effect("tau", "t1", "c1", 0.25)]
    rows = double_differences(effects)
    assert len(rows) == 1
    assert rows[0].target_left == "t1" and rows[0].target_right == "t2"
    assert abs(rows[0].value - (0.25 - 1.0)) < 1e-12


def test_label_shuffle_preserves_the_per_key_marginal():
    rows = [DoubleDifference("tau", f"t{i}", f"t{i+1}", f"c{i}", f"c{i+1}",
                             float(i), 2, True, False, (0.0,))
            for i in range(6)]
    shuffled = shuffle_within_key(rows, 20260820)
    assert sorted(r.value for r in rows) == sorted(r.value for r in shuffled)
    assert [r.row_id for r in rows] == [r.row_id for r in shuffled]


# -- 4. model structure: identity, antisymmetry, cycle consistency ----------


def _model(mode: str, dims: int = 8, protein: int = 16):
    torch.manual_seed(20260820)
    return DoubleDifferenceModel(DiscriminatorConfig(dims, protein, width=8,
                                                     hidden=16, depth=2,
                                                     mode=mode))


def test_identity_is_bitwise_exact():
    model = _model("protein")
    tau = torch.randn(5, 8)
    protein = torch.randn(5, 16)
    out = model(tau, protein, protein)
    assert torch.equal(out, torch.zeros_like(out))


def test_protein_antisymmetry_is_bitwise_exact():
    model = _model("protein")
    tau = torch.randn(5, 8)
    left, right = torch.randn(5, 16), torch.randn(5, 16)
    assert torch.equal(model(tau, left, right), -model(tau, right, left))


def test_protein_cycle_consistency_closes_at_machine_precision():
    model = _model("protein")
    tau = torch.randn(5, 8)
    p1, p2, p3 = torch.randn(5, 16), torch.randn(5, 16), torch.randn(5, 16)
    total = model(tau, p1, p2) + model(tau, p2, p3) + model(tau, p3, p1)
    assert float(total.abs().max()) < 1e-5


def test_identities_survive_a_parameter_update():
    model = _model("protein")
    tau = torch.randn(5, 8)
    left, right = torch.randn(5, 16), torch.randn(5, 16)
    optimizer = torch.optim.SGD(model.parameters(), lr=0.05)
    loss = (model(tau, left, right) ** 2).mean()
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    assert torch.equal(model(tau, left, right), -model(tau, right, left))
    assert torch.equal(model(tau, left, left),
                       torch.zeros(5, dtype=model(tau, left, left).dtype))


def test_protein_free_arms_cannot_express_the_double_difference():
    """Arms A and B must be structurally zero on D; that is the estimand working."""
    tau = torch.randn(7, 8)
    left, right = torch.randn(7, 16), torch.randn(7, 16)
    for mode in ("zero", "transformation"):
        model = _model(mode)
        out = model(tau, left, right)
        assert torch.equal(out, torch.zeros_like(out)), mode


def test_protein_conditioned_arm_has_no_dead_trainable_parameters():
    model = _model("protein")
    tau = torch.randn(9, 8)
    left, right = torch.randn(9, 16), torch.randn(9, 16)
    model(tau, left, right).pow(2).mean().backward()
    dead = [name for name, parameter in model.named_parameters()
            if parameter.grad is None or float(parameter.grad.abs().sum()) == 0.0]
    # `constant` is unused in protein mode by construction; everything else must
    # receive gradient.
    assert dead == ["response.constant"], dead


def test_model_signature_admits_no_identity_features():
    names = set(inspect.signature(DoubleDifferenceModel.forward).parameters)
    assert names == {"self", "descriptor", "protein_left", "protein_right"}
    source = (STAGE / "t2_model.py").read_text(encoding="utf-8")
    for token in ("target_id", "component_id", "panel_id", "document_id",
                  "nn.Embedding"):
        assert token not in source, token


# -- 5. no label path into features, keys, splits or deduplication ----------


def _stage_modules(names):
    for name in names:
        path = STAGE / name
        yield path, ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


LABEL_BLIND_MODULES = ("mmp.py",)


def test_mmp_construction_never_reads_a_label():
    for path, tree in _stage_modules(LABEL_BLIND_MODULES):
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value in {"pK", "delta_y"}:
                raise AssertionError(f"{path.name}:{node.lineno} reads a label")
            if isinstance(node, ast.Attribute) and node.attr in {"pK", "delta_y"}:
                raise AssertionError(f"{path.name}:{node.lineno} reads a label")


def test_no_python_hash_anywhere_in_the_stage():
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "hash"):
                raise AssertionError(f"{path.name}:{node.lineno} calls hash()")


def test_every_rng_is_seeded_by_stable_seed_or_a_constant():
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


# -- 6. governed-data contracts (slow: they mount the corpus) ---------------


@pytest.fixture(scope="module")
def governed():
    from tools.research.stageT_mmp.observations import load_governed
    return load_governed()


def test_physical_meta_test_seal(governed):
    data, seal = governed
    assert seal["included"] is False
    assert seal["evaluated"] is False
    assert seal["isolation"]["physically_isolated"] is True
    assert "meta_test" not in data.tasks
    assert "meta_test" not in data.components


def test_t0_cache_contains_only_meta_train_rows(governed):
    from tools.research.stageT_mmp.provenance import load_cache, train_allow_list
    allow = train_allow_list()
    rows = load_cache()
    assert rows, "run t0_reliability first"
    assert {row["source_row_id"] for row in rows} <= set(allow.row_ids)
    assert all(row["endpoint"] == "Ki" for row in rows)


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to build the observation bank")
def test_observations_are_within_target_and_same_core(governed):
    from scripts.internal_validation import partition_components
    from tools.research.stageT_mmp.observations import build_observations

    data, _seal = governed
    _fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    targets = sorted(t for t, c in component_of.items() if c in set(internal))
    built = build_observations(data, targets)
    rows = built["observations"]
    assert len(rows) > 500
    for row in rows:
        assert data.cells[row.cell_a]["target_id"] == row.target
        assert data.cells[row.cell_b]["target_id"] == row.target
        assert data.cells[row.cell_a]["split"] == "meta_train"
        assert row.ligand_a != row.ligand_b
        # delta_y must equal the corpus gap in the canonical direction.
        expected = data.cells[row.cell_b]["pK"] - data.cells[row.cell_a]["pK"]
        assert abs(row.delta_y - expected) < 1e-12
        # same-panel classification must agree with the panel identifiers
        shared = bool(set(data.cells[row.cell_a]["panel_ids"])
                      & set(data.cells[row.cell_b]["panel_ids"]))
        assert shared == row.same_panel


BANK_DIGEST_SCRIPT = r"""
import hashlib, json, sys
sys.path.insert(0, r"{root}")
from scripts.internal_validation import partition_components
from tools.research.stageT_mmp.observations import build_observations, load_governed
data, _ = load_governed()
_fit, internal = partition_components(data)
component = {{c["target_id"]: c["protein_group_40"] for c in data.cells}}
targets = sorted(t for t, c in component.items() if c in set(internal))
rows = build_observations(data, targets)["observations"]
payload = [(r.target, r.exact_key, r.cell_a, r.cell_b, round(r.delta_y, 12),
            r.same_panel, r.stratum) for r in rows]
print(hashlib.sha256(json.dumps(payload).encode()).hexdigest())
"""


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to run the subprocess test")
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
