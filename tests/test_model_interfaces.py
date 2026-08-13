import pytest
import torch

from model.encoders import LigandEncoder, ProteinEncoder
from model.bpsf import BipartitePairSectionFormer
from scripts.build_ligand_bank import ATOM_FEAT_DIM, BOND_FEAT_DIM, featurize_smiles


DTYPE = torch.float64


def test_preprocessing_graph_feeds_ligand_encoder_directly():
    graph = featurize_smiles("c1ccccc1O", max_atoms=16)
    encoder = LigandEncoder(d_h=16, n_layers=2, dtype=DTYPE)
    X = torch.as_tensor(graph["X"], dtype=DTYPE).unsqueeze(0)
    A = torch.as_tensor(graph["A"], dtype=DTYPE).unsqueeze(0)
    mask = torch.as_tensor(graph["mask"], dtype=DTYPE).unsqueeze(0)
    pooled, atoms = encoder(X, A, mask)

    assert tuple(pooled.shape) == (1, 16)
    assert tuple(atoms.shape) == (1, 16, 16)
    assert X.shape[-1] == ATOM_FEAT_DIM
    assert A.shape[-1] == BOND_FEAT_DIM


def test_ligand_encoder_rejects_a_zero_atom_sample():
    encoder = LigandEncoder(d_h=8, n_layers=1, dtype=DTYPE)
    atoms = torch.zeros(1, 4, ATOM_FEAT_DIM, dtype=DTYPE)
    bonds = torch.zeros(1, 4, 4, BOND_FEAT_DIM, dtype=DTYPE)
    mask = torch.zeros(1, 4, dtype=DTYPE)
    with pytest.raises(ValueError, match="zero-atom"):
        encoder(atoms, bonds, mask)


def test_protein_encoder_preserves_residue_axis():
    encoder = ProteinEncoder(protein_dim=6, d_h=8, dtype=DTYPE)
    pooled, residues = encoder(
        torch.randn(2, 6, dtype=DTYPE),
        torch.randn(2, 5, 6, dtype=DTYPE),
    )
    assert tuple(pooled.shape) == (2, 8)
    assert tuple(residues.shape) == (2, 5, 8)


def test_pair_trunk_rejects_empty_partners():
    bridge = BipartitePairSectionFormer(4, 2, pair_dim=4, blocks=1,
                                        latent_count=2, heads=1, dtype=DTYPE)
    atoms = torch.zeros(1, 3, 4, dtype=DTYPE)
    residues = torch.zeros(1, 5, 4, dtype=DTYPE)
    atom_mask = torch.zeros(1, 3, dtype=DTYPE)
    residue_mask = torch.ones(1, 5, dtype=DTYPE)
    with pytest.raises(ValueError, match="zero-atom"):
        bridge(atoms, residues, atom_mask, residue_mask)
    atom_mask.fill_(1)
    residue_mask.zero_()
    with pytest.raises(ValueError, match="zero-residue"):
        bridge(atoms, residues, atom_mask, residue_mask)
