"""Phase 1: structural falsification of the SAR field, before any training.

These tests exist to make the Phase 2 numbers interpretable.  Each one closes a
failure mode that has actually occurred in this repository or that would make
the hypothesis untestable:

* the three algebraic identities are the whole claim of the construction;
* a `hash()`-derived seed silently changes the bank per process (Stage R,
  Stage A);
* a cross-target pair would smuggle between-target level into a within-target
  quantity;
* a protein swap that also perturbs the ligand path is not a counterfactual;
* a query label reaching a model input is the leak this protocol forbids;
* a dead gradient into the protein path would make a null result vacuous;
* a synthetic positive task shows the machine can find protein-dependent
  transformation effects when they exist, and a synthetic ligand-only task
  shows it does not invent them when they do not.
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

from scripts.internal_validation import partition_components  # noqa: E402
from tools.research.stageS_sar_field import features as features_module  # noqa: E402
from tools.research.stageS_sar_field import pairs as pairs_module  # noqa: E402
from tools.research.stageS_sar_field.features import (  # noqa: E402
    LigandFeatureStore, ProteinFeatureStore, cross_component_permutation,
    hard_wrong_protein_map, within_target_label_shuffle,
)
from tools.research.stageS_sar_field.field import (  # noqa: E402
    FieldConfig, SARField, build_field,
)
from tools.research.stageS_sar_field.pairs import (  # noqa: E402
    build_target_pairs, component_of_target, expand_orientations, load_data,
    target_balanced_bank,
)

STAGE = Path(__file__).resolve().parents[1]
SLOW = os.environ.get("RUN_SLOW") == "1"


# -- fixtures ---------------------------------------------------------------


@pytest.fixture(scope="module")
def data():
    return load_data()


@pytest.fixture(scope="module")
def small_field():
    torch.manual_seed(20260819)
    return build_field(FieldConfig(hidden=16, graph_layers=2, coordinate=8,
                                   response=8, potential_width=16,
                                   potential_depth=2, dtype=torch.float64))


def synthetic_ligands(count: int, atoms: int = 9, seed: int = 7):
    generator = torch.Generator().manual_seed(seed)
    atom_features = torch.rand(count, atoms, 32, generator=generator,
                               dtype=torch.float64)
    bonds = torch.zeros(count, atoms, atoms, 12, dtype=torch.float64)
    for molecule in range(count):
        for left in range(atoms - 1):
            value = torch.rand(12, generator=generator, dtype=torch.float64)
            bonds[molecule, left, left + 1] = value
            bonds[molecule, left + 1, left] = value
    mask = torch.ones(count, atoms, dtype=torch.float64)
    fingerprint = (torch.rand(count, 1024, generator=generator,
                              dtype=torch.float64) > 0.97).to(torch.float64)
    return atom_features, bonds, mask, fingerprint


def synthetic_proteins(count: int, dim: int = 640, slots: int = 8, seed: int = 11):
    generator = torch.Generator().manual_seed(seed)
    pooled = torch.randn(count, dim, generator=generator, dtype=torch.float64)
    residues = torch.randn(count, slots, dim, generator=generator,
                           dtype=torch.float64)
    mask = torch.ones(count, slots, dtype=torch.float64)
    return pooled, residues, mask


# -- 1. the three algebraic identities --------------------------------------


def test_antisymmetry_is_bitwise_exact(small_field):
    atoms, bonds, mask, prints = synthetic_ligands(6)
    pooled, residues, protein_mask = synthetic_proteins(6)
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    forward = small_field(phi[:3], phi[3:], response[:3])
    backward = small_field(phi[3:], phi[:3], response[:3])
    assert torch.equal(forward, -backward)


def test_identity_is_bitwise_exact(small_field):
    atoms, bonds, mask, prints = synthetic_ligands(4)
    pooled, residues, protein_mask = synthetic_proteins(4)
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    same = small_field(phi, phi, response)
    assert torch.equal(same, torch.zeros_like(same))


def test_cycle_consistency_closes_at_machine_precision(small_field):
    atoms, bonds, mask, prints = synthetic_ligands(9)
    pooled, residues, protein_mask = synthetic_proteins(3)
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    a, b, c = phi[0:3], phi[3:6], phi[6:9]
    total = (small_field(a, b, response) + small_field(b, c, response)
             + small_field(c, a, response))
    # Exact in exact arithmetic (the sum telescopes); in IEEE-754 it closes to
    # the rounding residual, which is what is asserted and reported.
    assert float(total.abs().max()) < 1e-12


def test_identities_hold_after_a_parameter_update(small_field):
    """The identities are properties of the construction, not of the fit."""
    atoms, bonds, mask, prints = synthetic_ligands(6)
    pooled, residues, protein_mask = synthetic_proteins(3)
    optimizer = torch.optim.SGD(small_field.parameters(), lr=0.1)
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    loss = (small_field(phi[:3], phi[3:], response) ** 2).mean()
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    assert torch.equal(small_field(phi[:3], phi[3:], response),
                       -small_field(phi[3:], phi[:3], response))
    assert torch.equal(small_field(phi, phi, torch.cat([response] * 2)),
                       torch.zeros(6, dtype=torch.float64))


def test_no_even_quadratic_term_in_the_potential(small_field):
    """The rejected `e^T H e` term must not be present in any form.

    A term even under `e -> -e` cannot appear in a signed prediction.  The
    construction excludes it structurally -- the potential is a function of one
    ligand coordinate and never sees a difference at all -- and the behavioural
    consequence is asserted directly: negating the direction negates the whole
    prediction, which no even term can survive.
    """
    signature = inspect.signature(SARField.forward)
    assert list(signature.parameters) == ["self", "phi_a", "phi_b", "response"]
    body = inspect.getsource(SARField.forward)
    assert "potential_module(phi_b, response)" in body
    assert "potential_module(phi_a, response)" in body
    # The potential's own signature takes a coordinate, never a pair.
    assert list(inspect.signature(SARField.potential).parameters) == [
        "self", "coordinate", "response"]
    # Behavioural check: a nonzero even component would break exact negation.
    atoms, bonds, mask, prints = synthetic_ligands(6)
    pooled, residues, protein_mask = synthetic_proteins(3)
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    forward = small_field(phi[:3], phi[3:], response)
    assert torch.equal(forward + small_field(phi[3:], phi[:3], response),
                       torch.zeros(3, dtype=torch.float64))
    assert float(forward.abs().max()) > 0.0


# -- 2. support/query and pair-order isolation ------------------------------


def test_there_is_no_support_pathway():
    """No label, support set or episode object can reach a model input."""
    for function in (SARField.forward, SARField.phi, SARField.alpha,
                     SARField.potential):
        names = set(inspect.signature(function).parameters)
        assert not names & {"support", "support_y", "labels", "y", "delta_y",
                            "query_y", "episode"}


def test_pair_order_and_batch_composition_do_not_change_predictions(small_field):
    atoms, bonds, mask, prints = synthetic_ligands(8)
    pooled, residues, protein_mask = synthetic_proteins(4)
    phi = small_field.phi(atoms, bonds, mask, prints)
    response = small_field.alpha(pooled, residues, protein_mask)
    full = small_field(phi[:4], phi[4:], response)
    permutation = torch.tensor([2, 0, 3, 1])
    permuted = small_field(phi[:4][permutation], phi[4:][permutation],
                           response[permutation])
    assert torch.allclose(full[permutation], permuted, atol=0, rtol=0)
    single = small_field(phi[1:2], phi[5:6], response[1:2])
    assert torch.allclose(full[1:2], single, atol=1e-12)


def test_graph_padding_does_not_change_a_ligand_coordinate(small_field):
    """A prediction must not depend on which other ligands share its batch."""
    atoms, bonds, mask, prints = synthetic_ligands(3, atoms=9)
    wide_atoms = torch.nn.functional.pad(atoms, (0, 0, 0, 7))
    wide_bonds = torch.nn.functional.pad(bonds, (0, 0, 0, 7, 0, 7))
    wide_mask = torch.nn.functional.pad(mask, (0, 7))
    narrow = small_field.phi(atoms, bonds, mask, prints)
    wide = small_field.phi(wide_atoms, wide_bonds, wide_mask, prints)
    assert float((narrow - wide).abs().max()) < 1e-10


# -- 3. no Python hash() seeds ----------------------------------------------


def _stage_modules():
    for path in sorted(STAGE.rglob("*.py")):
        yield path, ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def test_stage_source_never_calls_python_hash():
    """`hash()` is salted per process by PYTHONHASHSEED.

    Parsed rather than grepped, so a docstring that *mentions* the defect does
    not count as committing it.
    """
    offenders = []
    for path, tree in _stage_modules():
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "hash"):
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, offenders


def test_every_string_derived_seed_goes_through_stable_seed():
    """`default_rng` must be seeded by a frozen int or by `stable_seed(...)`."""
    offenders = []
    for path, tree in _stage_modules():
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


# -- 4. stable pair banks across processes ----------------------------------


BANK_DIGEST_SCRIPT = r"""
import hashlib, json, sys
sys.path.insert(0, r"{root}")
from scripts.internal_validation import partition_components
from tools.research.stageS_sar_field.pairs import (
    build_target_pairs, component_of_target, load_data, target_balanced_bank)
data = load_data()
fit, internal = partition_components(data)
component = component_of_target(data)
targets = sorted(t for t, c in component.items() if c in set(internal))
by_target = build_target_pairs(data, "meta_train", targets)
bank = target_balanced_bank(by_target, 20260819, 8)
payload = [(s.target, s.a, s.b, round(s.delta_y, 12), round(s.tanimoto, 12),
            s.same_panel, s.stratum) for s in bank]
print(hashlib.sha256(json.dumps(payload).encode()).hexdigest())
"""


@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to run the subprocess test")
def test_pair_bank_is_identical_across_processes(tmp_path):
    script = tmp_path / "digest.py"
    script.write_text(BANK_DIGEST_SCRIPT.format(root=str(REPO)),
                      encoding="utf-8")
    digests = []
    for salt in ("0", "1", "12345"):
        environment = dict(os.environ, PYTHONHASHSEED=salt)
        result = subprocess.run([sys.executable, str(script)], check=True,
                                capture_output=True, text=True,
                                env=environment, cwd=str(REPO))
        digests.append(result.stdout.strip())
    assert len(set(digests)) == 1, digests


# -- 5. no cross-target pair construction -----------------------------------


def test_every_constructed_pair_is_within_one_target(data):
    fit, internal = partition_components(data)
    component = component_of_target(data)
    targets = sorted(t for t, c in component.items() if c in set(internal))
    by_target = build_target_pairs(data, "meta_train", targets)
    checked = 0
    for target, specs in by_target.items():
        for spec in specs:
            assert data.cells[spec.a]["target_id"] == target
            assert data.cells[spec.b]["target_id"] == target
            assert data.cells[spec.a]["ligand_id"] != data.cells[spec.b]["ligand_id"]
            assert data.cells[spec.a]["split"] == "meta_train"
            checked += 1
    assert checked > 1000


def test_delta_y_matches_the_corpus_labels(data):
    fit, internal = partition_components(data)
    component = component_of_target(data)
    targets = sorted(t for t, c in component.items() if c in set(internal))
    by_target = build_target_pairs(data, "meta_train", targets)
    for specs in by_target.values():
        for spec in specs[:50]:
            expected = data.cells[spec.b]["pK"] - data.cells[spec.a]["pK"]
            assert abs(spec.delta_y - expected) < 1e-12


def test_reversed_orientation_negates_delta_y(data):
    fit, internal = partition_components(data)
    component = component_of_target(data)
    targets = sorted(t for t, c in component.items() if c in set(internal))[:3]
    by_target = build_target_pairs(data, "meta_train", targets)
    specs = [s for group in by_target.values() for s in group][:200]
    balanced = expand_orientations(specs)
    assert len(balanced) == 2 * len(specs)
    assert abs(sum(s.delta_y for s in balanced)) < 1e-9


# -- 6. protein replacement changes only the protein input ------------------


def test_protein_swap_leaves_the_ligand_coordinate_bitwise_identical(small_field):
    atoms, bonds, mask, prints = synthetic_ligands(4)
    pooled, residues, protein_mask = synthetic_proteins(2)
    phi_once = small_field.phi(atoms, bonds, mask, prints)
    phi_twice = small_field.phi(atoms, bonds, mask, prints)
    assert torch.equal(phi_once, phi_twice)
    left = small_field.alpha(pooled[:1], residues[:1], protein_mask[:1])
    right = small_field.alpha(pooled[1:], residues[1:], protein_mask[1:])
    assert not torch.equal(left, right)
    # `phi` has no protein argument at all, so a swap cannot reach it.
    assert "pooled" not in inspect.signature(SARField.phi).parameters
    assert "residues" not in inspect.signature(SARField.phi).parameters


def test_hard_wrong_protein_rule_is_admissible(data):
    proteins = ProteinFeatureStore(data)
    fit, internal = partition_components(data)
    component = component_of_target(data)
    fit_targets = sorted(t for t, c in component.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component.items()
                              if c in set(internal))
    donors = hard_wrong_protein_map(data, proteins, internal_targets, fit_targets)
    documents = pairs_module.target_panel_documents(data)
    assert set(donors) == set(internal_targets)
    for recipient, donor in donors.items():
        assert component[donor] != component[recipient]
        assert donor in set(fit_targets)
        assert not (documents[donor] & documents[recipient])


def test_protein_permutation_always_crosses_components(data):
    component = component_of_target(data)
    targets = sorted(data.tasks["meta_train"])
    mapping = cross_component_permutation(data, targets, 20260819)
    assert set(mapping) == set(targets)
    assert sorted(mapping.values()) == sorted(targets)
    for key, value in mapping.items():
        assert component[key] != component[value]


def test_label_shuffle_preserves_the_target_level(data):
    targets = sorted(data.tasks["meta_train"])[:20]
    shuffled = within_target_label_shuffle(data, targets, 20260819)
    by_target: dict[str, list[int]] = {}
    for index, cell in enumerate(data.cells):
        if cell["target_id"] in set(targets):
            by_target.setdefault(cell["target_id"], []).append(index)
    for target, indices in by_target.items():
        original = sorted(round(float(data.cells[i]["pK"]), 9) for i in indices)
        permuted = sorted(round(shuffled[i], 9) for i in indices)
        assert original == permuted


# -- 7. no query labels enter model inputs ----------------------------------


def test_predictions_are_invariant_to_the_label_field(data, small_field):
    """Changing `delta_y` must not change any prediction.

    The label is only ever read to build the loss, but this asserts it on the
    real batching path rather than trusting the reading of the code.
    """
    from dataclasses import replace

    from tools.research.stageS_sar_field.train import Trainer, ARMS, build_banks

    proteins = ProteinFeatureStore(data)
    ligands = LigandFeatureStore(data)
    banks = build_banks(data, proteins, 20260819)
    trainer = Trainer(data, ligands, proteins, banks, ARMS["B_protein"],
                      torch.device("cpu"), 20260819)
    specs = banks.internal_same_panel[:24]
    with torch.no_grad():
        before = trainer.predict(specs, [s.target for s in specs]).clone()
    corrupted = [replace(spec, delta_y=spec.delta_y + 12.5) for spec in specs]
    with torch.no_grad():
        after = trainer.predict(corrupted, [s.target for s in corrupted])
    assert torch.equal(before, after)


# -- 8. gradients reach both paths ------------------------------------------


def test_gradient_reaches_the_ligand_and_protein_paths():
    torch.manual_seed(3)
    field = build_field(FieldConfig(hidden=16, graph_layers=2, coordinate=8,
                                    response=8, potential_width=16,
                                    potential_depth=2, dtype=torch.float64))
    atoms, bonds, mask, prints = synthetic_ligands(6)
    pooled, residues, protein_mask = synthetic_proteins(3)
    phi = field.phi(atoms, bonds, mask, prints)
    response = field.alpha(pooled, residues, protein_mask)
    prediction = field(phi[:3], phi[3:], response)
    prediction.pow(2).mean().backward()
    ligand = sum(float(p.grad.abs().sum()) for p in field.coordinate.parameters()
                 if p.grad is not None)
    protein = sum(float(p.grad.abs().sum())
                  for p in field.response_module.parameters() if p.grad is not None)
    potential = sum(float(p.grad.abs().sum())
                    for p in field.potential_module.parameters()
                    if p.grad is not None)
    assert ligand > 0.0, "no gradient into the ligand transformation path"
    assert protein > 0.0, "no gradient into the protein-response path"
    assert potential > 0.0
    assert all(p.grad is not None for p in field.parameters())


# -- 9/10. synthetic recovery and synthetic restraint ------------------------


def _fit_synthetic(delta_builder, steps: int = 400, seed: int = 5):
    """Train the real field on a synthetic task and score a protein swap.

    Returns (correct MSE, wrong-protein MSE, label variance, correlation between
    the correct and wrong-protein predictions).
    """
    torch.manual_seed(seed)
    ligand_count, protein_count = 24, 6
    atoms, bonds, mask, prints = synthetic_ligands(ligand_count, seed=seed)
    pooled, residues, protein_mask = synthetic_proteins(protein_count, slots=8,
                                                        seed=seed + 1)
    generator = torch.Generator().manual_seed(seed + 2)
    # A fixed, unlearnable ligand descriptor and a fixed protein response.
    descriptor = atoms.mean(dim=1) @ torch.randn(32, 4, generator=generator,
                                                 dtype=torch.float64)
    weights = torch.randn(protein_count, 4, generator=generator,
                          dtype=torch.float64)
    left, right, who = [], [], []
    for protein in range(protein_count):
        for a in range(ligand_count):
            for b in range(ligand_count):
                if a == b:
                    continue
                left.append(a)
                right.append(b)
                who.append(protein)
    left = torch.tensor(left)
    right = torch.tensor(right)
    who = torch.tensor(who)
    truth = delta_builder(descriptor, weights, left, right, who)

    field = build_field(FieldConfig(hidden=24, graph_layers=2, coordinate=8,
                                    response=8, potential_width=32,
                                    potential_depth=2, dtype=torch.float64))
    optimizer = torch.optim.AdamW(field.parameters(), lr=3e-3)
    for _ in range(steps):
        phi = field.phi(atoms, bonds, mask, prints)
        response = field.alpha(pooled, residues, protein_mask)
        prediction = field(phi[left], phi[right], response[who])
        loss = torch.nn.functional.mse_loss(prediction, truth)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    with torch.no_grad():
        phi = field.phi(atoms, bonds, mask, prints)
        response = field.alpha(pooled, residues, protein_mask)
        correct = field(phi[left], phi[right], response[who])
        swapped = (who + 1) % protein_count
        wrong = field(phi[left], phi[right], response[swapped])
    correlation = float(np.corrcoef(correct.numpy(), wrong.numpy())[0, 1])
    return (float(((correct - truth) ** 2).mean()),
            float(((wrong - truth) ** 2).mean()),
            float(truth.var()), correlation)


def test_synthetic_positive_task_recovers_protein_dependent_effects():
    """delta_y = <phi_b - phi_a, w(P)>: the effect is real and protein-specific."""
    def builder(descriptor, weights, left, right, who):
        return ((descriptor[right] - descriptor[left]) * weights[who]).sum(-1)

    correct, wrong, variance, correlation = _fit_synthetic(builder)
    assert correct < 0.20 * variance, (correct, variance)
    assert wrong > 3.0 * correct, (correct, wrong)
    assert correlation < 0.9, correlation


def test_synthetic_ligand_only_task_does_not_invent_protein_dependence():
    """delta_y = <phi_b - phi_a, w0>: no protein effect exists to be found."""
    def builder(descriptor, weights, left, right, who):
        shared = weights[0]
        return ((descriptor[right] - descriptor[left]) * shared).sum(-1)

    correct, wrong, variance, correlation = _fit_synthetic(builder)
    assert correct < 0.20 * variance, (correct, variance)
    assert (wrong - correct) < 0.30 * variance, (correct, wrong, variance)
    assert correlation > 0.9, correlation


# -- bookkeeping ------------------------------------------------------------


def test_stage_never_mounts_meta_test(data):
    seal = data.seal_record()
    assert seal["included"] is False
    assert seal["evaluated"] is False
    assert seal["isolation"]["physically_isolated"] is True
    assert "meta_test" not in data.tasks
    assert "meta_test" not in data.components


def test_stage_source_never_names_the_development_validation_split():
    """No module in this stage may reach for `meta_val` at all.

    Parsed, so the prose in a docstring that explains *why* it is off limits is
    not mistaken for a read.
    """
    # Assembled at run time so this file does not itself contain the literal
    # it forbids -- otherwise the test fails on its own source.
    forbidden = "meta" + "_val"
    offenders = []
    for path, tree in _stage_modules():
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == forbidden:
                offenders.append(f"{path.name}:{node.lineno}")
    assert not offenders, offenders
