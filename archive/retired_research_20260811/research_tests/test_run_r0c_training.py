import numpy as np
import torch

from research.correspondence_router.run_r0c_training import (
    System,
    _batches,
    _new_model,
    _pack,
    _prepare_n1,
)


def system(name: str, atoms: int = 3, residues: int = 4) -> System:
    rng = np.random.default_rng(abs(hash(name)) % (2**32))
    prior = rng.random((atoms, 128, 5), dtype=np.float32)
    prior /= prior.sum(axis=-1, keepdims=True)
    return System(
        entry=name,
        component=name,
        atom_state=rng.normal(size=(atoms, 128)).astype(np.float32),
        residue_state=rng.normal(size=(residues, 128)).astype(np.float32),
        distance_prior_slot=prior,
        atom_class=np.arange(atoms, dtype=np.int8) % 8,
        residue_class=np.arange(residues, dtype=np.int8) % 6,
        slots=np.asarray([0, 0, 1, 1][:residues], dtype=np.int64),
        labels=np.zeros((atoms, residues), dtype=np.int64),
    )


def test_n1_pack_removes_within_slot_state_and_chemistry_identity():
    value = system("a")
    packed = _pack([value], "n1", "cpu")
    assert torch.equal(packed["residue_states"][0, 0], packed["residue_states"][0, 1])
    assert torch.equal(packed["residue_states"][0, 2], packed["residue_states"][0, 3])
    assert torch.equal(packed["compatibility"][0, :, 0],
                       packed["compatibility"][0, :, 1])
    assert torch.equal(packed["compatibility"][0, :, 2],
                       packed["compatibility"][0, :, 3])
    cached_state = value.n1_residue_state
    cached_compatibility = value.n1_compatibility
    _prepare_n1(value)
    assert value.n1_residue_state is cached_state
    assert value.n1_compatibility is cached_compatibility


def test_arms_share_initial_parameters_and_zero_distance_head():
    full = _new_model(1701, "full", "cpu")
    additive = _new_model(1701, "n2", "cpu")
    assert full.interaction_mode == "bilinear"
    assert additive.interaction_mode == "additive"
    for name, value in full.state_dict().items():
        assert torch.equal(value, additive.state_dict()[name])
    assert torch.count_nonzero(full.distance_residual.weight) == 0


def test_pair_budget_batches_every_system_once():
    systems = [system(str(index), atoms=2 + index % 2) for index in range(10)]
    batches = _batches(systems, np.arange(len(systems)))
    assert [value.entry for batch in batches for value in batch] == [
        value.entry for value in systems]


def test_pack_derives_bounded_contact_from_rounded_distance_prior():
    value = system("rounded-prior", atoms=1, residues=1)
    value.distance_prior_slot[0, 0] = np.asarray(
        [0.6, 0.4001, 0.0, 0.0, 0.0], dtype=np.float32)
    packed = _pack([value], "full", "cpu")
    assert packed["contact_prior"][0, 0, 0] == 1.0
    prediction = _new_model(1701, "full", "cpu")(
        **{key: item for key, item in packed.items() if key != "labels"})
    assert torch.isfinite(prediction.distance_prob).all()
