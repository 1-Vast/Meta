import numpy as np

from research.e0_identifiability.audit_evidence import _safe_correlation, _teacher_values


def test_teacher_sufficient_statistic_reconstructs_direct_formula():
    rows = [{"ligand_state_key": "l", "active_protein_key": "p", "example_id": 0}]
    atom = np.zeros((2, 40), dtype=np.float32)
    atom[:, 32] = 1
    residue = np.zeros((3, 6), dtype=np.float32)
    residue[:, 0] = 1
    proteins = {"p": {"chemistry": residue}}
    ligands = {"l": {"chemistry": atom}}
    contact = np.full((2, 3), 0.5, dtype=np.float32)
    distance = np.zeros((2, 3, 5), dtype=np.float32)
    distance[..., 1] = 1
    direct, sufficient, _ = _teacher_values(
        rows, proteins, ligands, {0: {"contact": contact, "distance": distance}},
        np.ones((8, 6), dtype=np.float32))
    assert np.allclose(direct, sufficient, atol=1e-7)


def test_constant_correlation_is_reported_as_unavailable():
    result = _safe_correlation(np.ones(4), np.arange(4))
    assert result == {"pearson": None, "spearman": None}
