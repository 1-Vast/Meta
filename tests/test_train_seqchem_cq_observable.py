import numpy as np

from research.crossed_interaction.train_seqchem_cq_observable import (
    donor_maps,
    ligand_descriptor,
    product_feature,
    protein_descriptor,
)


def test_protein_descriptor_is_normalized_and_length_aware():
    descriptor = protein_descriptor("ACDEFGHIKLMNPQRSTVWY")
    assert descriptor.shape == (27,)
    assert np.isclose(descriptor[:20].sum(), 1.0)
    assert descriptor[-1] > 0


def test_ligand_descriptor_is_finite_and_contains_fingerprint_bits():
    descriptor = ligand_descriptor("CCO", fp_bits=32)
    assert descriptor.shape == (40,)
    assert np.isfinite(descriptor).all()
    assert descriptor[:32].sum() > 0


def test_product_feature_is_outer_product_contract():
    protein = np.asarray([1.0, 2.0])
    ligand = np.asarray([3.0, 5.0, 7.0])
    assert np.array_equal(
        product_feature(protein, ligand),
        np.asarray([3.0, 5.0, 7.0, 6.0, 10.0, 14.0]))


def test_donor_maps_choose_different_protein_group_and_scaffold():
    cells = [
        {"target_id": "t1", "protein_group_40": "g1", "ligand_id": "l1", "scaffold": "s1"},
        {"target_id": "t2", "protein_group_40": "g2", "ligand_id": "l2", "scaffold": "s2"},
        {"target_id": "t3", "protein_group_40": "g3", "ligand_id": "l3", "scaffold": "s3"},
    ]
    protein_donor, ligand_donor = donor_maps(cells)
    target_group = {cell["target_id"]: cell["protein_group_40"] for cell in cells}
    ligand_scaffold = {cell["ligand_id"]: cell["scaffold"] for cell in cells}
    assert all(target_group[target] != target_group[donor]
               for target, donor in protein_donor.items())
    assert all(ligand_scaffold[ligand] != ligand_scaffold[donor]
               for ligand, donor in ligand_donor.items())
