"""Stage 2 structural gates for sequence-derived locality-aware refinement."""
import pytest
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.locality_grammar import (
    LocalityGrammarModel, LocalizedContactGrammar, SlotLocalityRefiner, band_mask,
)
from model.similarity_grammar import SimilarityGrammarModel


PROTEIN_DIM = 16
RESIDUES = 12
BITS = 32


def _kwargs(**extra):
    base = dict(protein_dim=PROTEIN_DIM, hidden_dim=8, task_dim=4,
                ligand_layers=1, pair_dim=8, pair_latents=4, pair_heads=2)
    base.update(extra)
    return base


def _episode(generator, batch=2, support=3, query=4, atoms=5, residues=RESIDUES):
    def graph(count):
        raw = torch.rand(batch, count, atoms, ATOM_FEAT_DIM, generator=generator)
        bonds = torch.rand(batch, count, atoms, atoms, BOND_FEAT_DIM,
                           generator=generator)
        return raw, bonds * (bonds > 0.7), torch.ones(batch, count, atoms)
    s, q = graph(support), graph(query)
    return {
        "protein_pooled": torch.randn(batch, PROTEIN_DIM, generator=generator),
        "protein_tokens": torch.randn(batch, residues, PROTEIN_DIM,
                                      generator=generator),
        "protein_mask": torch.ones(batch, residues),
        "protein_chemistry": torch.rand(batch, residues, 4, generator=generator),
        "support_atoms": s[0], "support_bonds": s[1], "support_mask": s[2],
        "support_y": torch.randn(batch, support, generator=generator),
        "query_atoms": q[0], "query_bonds": q[1], "query_mask": q[2],
        "support_fingerprint": (torch.rand(batch, support, BITS,
                                           generator=generator) > 0.6).float(),
        "query_fingerprint": (torch.rand(batch, query, BITS,
                                         generator=generator) > 0.6).float(),
    }


def _forward(model, episode, **kwargs):
    return model(
        episode["protein_pooled"], episode["protein_tokens"],
        episode["protein_mask"], episode["support_atoms"],
        episode["support_bonds"], episode["support_mask"], episode["support_y"],
        episode["query_atoms"], episode["query_bonds"], episode["query_mask"],
        protein_chemistry=episode["protein_chemistry"],
        support_fingerprint=episode["support_fingerprint"],
        query_fingerprint=episode["query_fingerprint"], **kwargs)


def test_band_mask_is_a_sequence_window():
    mask = band_mask(6, 2)
    assert bool(mask[0, 0]) and bool(mask[0, 2]) and not bool(mask[0, 3])
    assert torch.equal(mask, mask.T)


def test_zero_gate_reproduces_the_accepted_baseline_exactly():
    """At initialisation the locality model must equal the similarity model."""
    torch.manual_seed(5501)
    baseline = SimilarityGrammarModel(**_kwargs(), use_learned_key=False)
    torch.manual_seed(5501)
    locality = LocalityGrammarModel(**_kwargs(), use_learned_key=False)
    locality.load_state_dict(baseline.state_dict(), strict=False)
    assert float(locality.slot_locality.gate) == 0.0
    episode = _episode(torch.Generator().manual_seed(5601))
    with torch.no_grad():
        left = _forward(baseline, episode).prediction
        right = _forward(locality, episode).prediction
    assert torch.allclose(left, right, atol=1e-6, rtol=1e-5)


def test_refiner_is_a_no_op_at_zero_gate_and_active_when_opened():
    torch.manual_seed(5502)
    refiner = SlotLocalityRefiner(8, heads=2, band=3)
    slots = torch.randn(2, RESIDUES, 8)
    mask = torch.ones(2, RESIDUES)
    assert torch.allclose(refiner(slots, mask), slots, atol=1e-7)
    with torch.no_grad():
        refiner.gate.fill_(1.0)
    assert not torch.allclose(refiner(slots, mask), slots, atol=1e-4)


def test_padded_slots_stay_zero_and_do_not_leak():
    torch.manual_seed(5503)
    refiner = SlotLocalityRefiner(8, heads=2, band=3)
    with torch.no_grad():
        refiner.gate.fill_(1.0)
    slots = torch.randn(2, RESIDUES, 8)
    mask = torch.ones(2, RESIDUES)
    mask[:, RESIDUES // 2:] = 0.0
    out = refiner(slots, mask)
    assert torch.count_nonzero(out[:, RESIDUES // 2:]) == 0
    # Changing padded slot content must not change live slot outputs.
    other = slots.clone()
    other[:, RESIDUES // 2:] = torch.randn_like(other[:, RESIDUES // 2:])
    assert torch.allclose(out[:, :RESIDUES // 2],
                          refiner(other, mask)[:, :RESIDUES // 2], atol=1e-6)


def test_localizer_never_selects_padding_and_keeps_the_atom_union():
    torch.manual_seed(5504)
    grammar = LocalizedContactGrammar(8, 4, 2, 8, slot_topk=3)
    pairs, atoms, residues = 2, 5, RESIDUES
    atom_states = torch.randn(pairs, atoms, 8)
    atom_mask = torch.ones(pairs, atoms)
    chemistry = torch.rand(pairs, atoms, 4)
    slots = torch.randn(pairs, residues, 8)
    slot_mask = torch.ones(pairs, residues)
    slot_mask[:, -4:] = 0.0
    occupancy, mean_state, max_state = grammar(
        atom_states, atom_mask, chemistry, slots, slot_mask)
    assert torch.isfinite(occupancy).all()
    assert occupancy.shape == (pairs, 4)
    # Perturbing masked-out slots must not change the output.
    other = slots.clone()
    other[:, -4:] = torch.randn_like(other[:, -4:])
    again = grammar(atom_states, atom_mask, chemistry, other, slot_mask)[0]
    assert torch.allclose(occupancy, again, atol=1e-6)


def test_topk_at_or_above_slot_count_is_the_dense_baseline():
    torch.manual_seed(5505)
    dense = LocalizedContactGrammar(8, 4, 2, 8, slot_topk=999)
    pairs, atoms, residues = 2, 4, RESIDUES
    args = (torch.randn(pairs, atoms, 8), torch.ones(pairs, atoms),
            torch.rand(pairs, atoms, 4), torch.randn(pairs, residues, 8),
            torch.ones(pairs, residues))
    from model.interaction_grammar import ContactGrammar
    torch.manual_seed(5505)
    plain = ContactGrammar(8, 4, 2, 8)
    plain.load_state_dict(dense.state_dict())
    assert torch.allclose(dense(*args)[0], plain(*args)[0], atol=1e-6)


@pytest.mark.parametrize("residues,atoms", [(12, 5), (17, 3), (33, 9)])
def test_variable_protein_and_ligand_lengths(residues, atoms):
    torch.manual_seed(5506)
    model = LocalityGrammarModel(**_kwargs(), use_learned_key=False,
                                 locality_band=4, slot_topk=6)
    episode = _episode(torch.Generator().manual_seed(5606),
                       atoms=atoms, residues=residues)
    out = _forward(model, episode)
    assert out.prediction.shape == (2, 4)
    assert torch.isfinite(out.prediction).all()


def test_zero_support_is_exactly_the_zero_shot_endpoint():
    torch.manual_seed(5507)
    model = LocalityGrammarModel(**_kwargs(), use_learned_key=False)
    with torch.no_grad():
        model.slot_locality.gate.fill_(0.7)
    episode = _episode(torch.Generator().manual_seed(5607))
    empty = dict(episode)
    for key in ("support_atoms", "support_bonds", "support_mask", "support_y",
                "support_fingerprint"):
        empty[key] = episode[key][:, :0]
    out = _forward(model, empty)
    assert torch.equal(out.prediction, out.zero_shot)


def test_support_permutation_invariance_and_query_equivariance():
    torch.manual_seed(5508)
    model = LocalityGrammarModel(**_kwargs(), use_learned_key=False)
    with torch.no_grad():
        model.slot_locality.gate.fill_(0.5)
    episode = _episode(torch.Generator().manual_seed(5608), support=4, query=4)
    base = _forward(model, episode)
    order = torch.tensor([2, 0, 3, 1])
    permuted = dict(episode)
    for key in ("support_atoms", "support_bonds", "support_mask", "support_y",
                "support_fingerprint"):
        permuted[key] = episode[key].index_select(1, order)
    assert torch.allclose(base.prediction, _forward(model, permuted).prediction,
                          atol=1e-5, rtol=1e-5)
    shifted = dict(episode)
    for key in ("query_atoms", "query_bonds", "query_mask", "query_fingerprint"):
        shifted[key] = episode[key].index_select(1, order)
    assert torch.allclose(base.prediction.index_select(1, order),
                          _forward(model, shifted).prediction,
                          atol=1e-5, rtol=1e-5)


def test_protein_shuffle_changes_the_zero_shot_prediction():
    torch.manual_seed(5509)
    model = LocalityGrammarModel(**_kwargs(), use_learned_key=False)
    with torch.no_grad():
        model.slot_locality.gate.fill_(1.0)
    episode = _episode(torch.Generator().manual_seed(5609))
    base = _forward(model, episode, adapt=False).prediction
    shuffled = dict(episode)
    order = torch.randperm(RESIDUES)
    shuffled["protein_tokens"] = episode["protein_tokens"].index_select(1, order)
    shuffled["protein_chemistry"] = episode["protein_chemistry"].index_select(1, order)
    moved = (base - _forward(model, shuffled, adapt=False).prediction).abs()
    assert float(moved.max()) > 1e-5, "slot order carries no information"


def test_every_trainable_tensor_receives_gradient():
    torch.manual_seed(5510)
    model = LocalityGrammarModel(**_kwargs(), use_learned_key=False)
    with torch.no_grad():
        model.slot_locality.gate.fill_(0.3)
    episode = _episode(torch.Generator().manual_seed(5610), support=3)
    _forward(model, episode).prediction.square().mean().backward()
    dead = [name for name, parameter in model.named_parameters()
            if parameter.requires_grad and (
                parameter.grad is None or float(parameter.grad.abs().sum()) == 0.0)]
    assert dead == [], f"dead trainable branches: {dead}"


def test_query_labels_are_not_an_input():
    torch.manual_seed(5511)
    model = LocalityGrammarModel(**_kwargs(), use_learned_key=False)
    episode = _episode(torch.Generator().manual_seed(5611))
    with pytest.raises(TypeError):
        _forward(model, episode, query_y=torch.zeros(1))
