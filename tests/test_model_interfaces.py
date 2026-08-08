import torch

from model.encoders import LigandEncoder, ProteinEncoder
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


def test_protein_encoder_preserves_residue_axis():
    encoder = ProteinEncoder(protein_dim=6, d_h=8, dtype=DTYPE)
    pooled, residues = encoder(
        torch.randn(2, 6, dtype=DTYPE),
        torch.randn(2, 5, 6, dtype=DTYPE),
    )
    assert tuple(pooled.shape) == (2, 8)
    assert tuple(residues.shape) == (2, 5, 8)
