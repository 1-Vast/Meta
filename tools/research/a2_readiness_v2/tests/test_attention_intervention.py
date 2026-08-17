"""The causal audit re-implements `ContactGrammar.forward` so that the two
protein channels can be driven separately. That re-implementation is only
trustworthy if it is *bit-equivalent* to the module when both channels get the
same input — otherwise every intervention result measures my transcription
rather than the model.

These tests establish that equivalence, then check the intervention itself:
driving both channels from the same donor must equal a plain donor forward, and
driving neither must equal a plain correct forward.
"""
from __future__ import annotations

from pathlib import Path
import sys

import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM       # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel           # noqa: E402
from tools.research.a2_readiness_v2.attention_causal_audit import (   # noqa: E402
    encode_residues, grammar_split, differential, level_change,
    relative_change,
)

PROTEIN_DIM, SLOTS, ATOMS, QUERIES = 32, 12, 7, 6


def build(seed: int = 0) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12, ligand_layers=2,
        pair_dim=24, pair_latents=6, pair_heads=2,
        use_learned_key=False).double().eval()


def protein_parts(seed: int) -> list:
    generator = torch.Generator().manual_seed(seed)
    return [
        torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        torch.ones(1, SLOTS, dtype=torch.float64),
        torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64),
    ]


def ligands(count: int, seed: int):
    generator = torch.Generator().manual_seed(seed)
    atoms = torch.rand(count, ATOMS, ATOM_FEAT_DIM, generator=generator,
                       dtype=torch.float64)
    bonds = torch.rand(count, ATOMS, ATOMS, BOND_FEAT_DIM, generator=generator,
                       dtype=torch.float64)
    bonds = bonds * (bonds > 0.7)
    bonds = 0.5 * (bonds + bonds.transpose(1, 2))
    mask = torch.ones(count, ATOMS, dtype=torch.float64)
    return atoms, bonds, mask


def stages(model, parts, key_parts=None, value_parts=None):
    atoms, bonds, mask = ligands(QUERIES, seed=5)
    _, atom_states = model.ligand_encoder(atoms, bonds, mask)
    key_residues, _, residue_mask = encode_residues(model, key_parts or parts)
    value_residues, _, _ = encode_residues(model, value_parts or parts)
    return grammar_split(
        model, atom_states, mask, model.atom_chemistry(atoms),
        key_residues.expand(QUERIES, -1, -1),
        value_residues.expand(QUERIES, -1, -1),
        residue_mask.expand(QUERIES, -1))


def test_the_reimplementation_matches_the_module_exactly():
    """Both channels from one protein must reproduce `ContactGrammar.forward`."""
    model = build(1)
    parts = protein_parts(11)
    atoms, bonds, mask = ligands(QUERIES, seed=5)
    _, atom_states = model.ligand_encoder(atoms, bonds, mask)
    residues, _, residue_mask = encode_residues(model, parts)

    reference = model.grammar(
        atom_states, mask, model.atom_chemistry(atoms),
        residues.expand(QUERIES, -1, -1), residue_mask.expand(QUERIES, -1))
    mine = stages(model, parts)

    for index, name in enumerate(("occupancy", "mean_state", "max_state")):
        assert torch.equal(reference[index], mine[name]), name


def test_driving_both_channels_from_the_donor_equals_a_plain_donor_forward():
    model = build(2)
    donor = protein_parts(22)
    swapped = stages(model, protein_parts(11), key_parts=donor, value_parts=donor)
    plain = stages(model, donor)
    for name in ("weight", "context", "occupancy", "mean_state", "max_state"):
        assert torch.equal(swapped[name], plain[name]), name


def test_the_two_channels_are_genuinely_separable():
    """Routing-only and content-only must differ from each other and from both."""
    model = build(3)
    correct, donor = protein_parts(11), protein_parts(22)
    both = stages(model, correct, key_parts=donor, value_parts=donor)
    routing = stages(model, correct, key_parts=donor, value_parts=correct)
    content = stages(model, correct, key_parts=correct, value_parts=donor)

    # Routing changes the attention weights; content does not touch them.
    assert not torch.allclose(routing["weight"], stages(model, correct)["weight"])
    assert torch.equal(content["weight"], stages(model, correct)["weight"])
    assert not torch.allclose(routing["context"], content["context"])
    assert not torch.allclose(both["mean_state"], routing["mean_state"])


def test_content_only_leaves_the_attention_weights_bit_identical():
    """`residue_value` must not reach the softmax; otherwise the split is a lie."""
    model = build(4)
    correct, donor = protein_parts(31), protein_parts(41)
    assert torch.equal(stages(model, correct)["weight"],
                       stages(model, correct, key_parts=correct,
                              value_parts=donor)["weight"])


def test_a_null_intervention_registers_as_zero_change():
    model = build(5)
    correct = protein_parts(51)
    left = stages(model, correct)["mean_state"]
    assert relative_change(left, left) == pytest.approx(0.0, abs=1e-12)
    assert level_change(left, left) == pytest.approx(0.0, abs=1e-12)


# `differential` casts to float32 to match the dtype the audit actually runs
# in, so the numerical floor of every relative-change measurement is float32
# epsilon, not float64. Measured below at ~2e-7 relative — four to six orders
# of magnitude below the 1e-2..1e-1 effects the audit reports.
FLOAT32_FLOOR = 1e-6


def test_the_differential_removes_exactly_the_level():
    generator = torch.Generator().manual_seed(6)
    values = torch.randn(QUERIES, 8, generator=generator, dtype=torch.float64)
    assert abs(differential(values).mean(0)).max() < FLOAT32_FLOOR
    shifted = values + 3.5
    assert abs(differential(values) - differential(shifted)).max() < FLOAT32_FLOOR


def test_relative_change_detects_a_planted_differential_perturbation():
    """The instrument must register a change it is given, and only that change."""
    generator = torch.Generator().manual_seed(7)
    values = torch.randn(QUERIES, 8, generator=generator, dtype=torch.float64)
    level_shifted = values + 2.0            # pure level: differential unchanged
    rotated = values.flip(0)                # pure differential change
    floor = relative_change(values, level_shifted)
    assert floor == pytest.approx(0.0, abs=FLOAT32_FLOOR)
    assert relative_change(values, rotated) > 0.5
    assert level_change(values, level_shifted) > 0.1
    # The reported effects must clear the instrument floor by orders of
    # magnitude, not by a margin.
    assert relative_change(values, rotated) > 1e4 * max(floor, 1e-12)
