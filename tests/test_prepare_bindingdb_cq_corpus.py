from research.crossed_interaction.prepare_bindingdb_cq_corpus import (
    UnionFind,
    canonical_smiles,
    sequence_identity,
)


def test_union_find_is_deterministic():
    union = UnionFind(["c", "a", "b"])
    union.union("c", "b")
    union.union("b", "a")
    assert {union.find(value) for value in ("a", "b", "c")} == {"a"}


def test_sequence_identity_contract():
    assert sequence_identity("ACDEFG", "ACDEFG") == 1.0
    assert sequence_identity("AAAAAA", "CCCCCC") == 0.0


def test_smiles_identity_is_canonicalized():
    assert canonical_smiles("C(C)O") == canonical_smiles("OCC")
