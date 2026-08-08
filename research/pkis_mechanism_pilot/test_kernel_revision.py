import numpy as np

from .kernel_revision import (
    bounded_z,
    fit_bounds,
    fit_convex_channel_weights,
    fit_separable_krr,
    kernel_ligand_from_smiles,
    predict_separable,
    protein_channel_features,
    source_rbf,
    tanimoto_kernel,
)


def test_channel_ligand_is_five_bounded_bit_vectors():
    record = kernel_ligand_from_smiles("x", "CC(=O)Nc1cc[nH]c1")
    assert record.channel_fingerprints.shape == (5, 1024)
    assert set(np.unique(record.channel_fingerprints)).issubset({0, 1})
    assert record.channel_fingerprints[-1].sum() > 0


def test_tanimoto_has_exact_self_similarity_for_nonempty_rows():
    bits = np.asarray([[1, 0, 1], [0, 1, 0], [0, 0, 0]], dtype=np.uint8)
    kernel = tanimoto_kernel(bits, bits)
    assert np.allclose(kernel, kernel.T)
    assert np.allclose(np.diag(kernel)[:2], 1.0)
    assert kernel[2, 2] == 0.0


def test_protein_channels_have_declared_aligned_dimensions():
    pocket = "A" * 85
    channels = protein_channel_features(pocket)
    assert len(channels) == 5
    assert [len(channel) for channel in channels] == [170, 170, 85, 85, 170]
    kernel, tau = source_rbf(np.stack([channels[-1], channels[-1]]))
    assert tau == 1.0
    assert np.allclose(kernel, 1.0)


def test_separable_solver_reconstructs_low_regularization_training_values():
    kl = np.asarray([[1.0, 0.2], [0.2, 1.0]])
    kp = np.asarray([[1.0, 0.1], [0.1, 1.0]])
    residual = np.asarray([[0.5, -0.2], [-0.4, 0.1]])
    coefficient = fit_separable_krr(kl, kp, residual, 1e-9)
    prediction = predict_separable(kl, coefficient, kp)
    assert np.allclose(prediction, residual, atol=1e-7)


def test_convex_channel_weights_and_bounded_z_contract():
    prediction = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    outcome = np.asarray([0.2, 0.3, 0.5])
    result = fit_convex_channel_weights(prediction, outcome)
    assert np.all(result["weights"] >= 0.0)
    assert result["sum"] <= 1.0 + 1e-10
    contributions = np.zeros((3, 2, 5), dtype=np.float64)
    contributions[..., :2] = prediction[:, None, :]
    z = bounded_z(contributions, fit_bounds(contributions))
    assert z.shape == contributions.shape
    assert np.all((0.0 <= z) & (z <= 1.0))
