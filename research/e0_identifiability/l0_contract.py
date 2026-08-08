"""Frozen contracts for the E-AFF-L0 protein affinity-location Gate.

Pure functions only: no data, no labels, no model. Frozen before execution.
"""
from __future__ import annotations

import numpy as np

from model import bands

# --- z_bio: bounded, permutation-invariant, pair-local coordinates ---------
# phi has axes (atom channel, residue class, radial bin):
#   atom channels : hydrophobe aromatic donor acceptor positive negative halogen other
#   residue class : 0 aliphatic 1 aromatic 2 polar 3 positive 4 negative 5 other
#   radial bins   : centres 2.0 3.5 5.0 6.5 8.0 9.5 Angstrom
APOLAR_ATOMS = (0, 1)
POLAR_ATOMS = (2, 3)
CHARGED_ATOMS = (4, 5)
APOLAR_RESIDUES = (0, 1)
POLAR_RESIDUES = (2,)
CHARGED_RESIDUES = (3, 4)
CLOSE_BINS = (0, 1)
CONTACT_BIN = (0,)
Z_BIO_NAMES = (
    "interface_extent", "short_range_fraction", "apolar_complementarity",
    "polar_complementarity", "charge_complementarity", "steric_proxy",
    "ligand_size_gauge",
)
Z_BIO_GAUGE_COORDINATES = ("ligand_size_gauge",)
MAX_LIGAND_ATOMS = 128.0


def z_bio(phi: np.ndarray, heavy_atoms: np.ndarray) -> np.ndarray:
    """Map the frozen 8x6x6 radial tensor to bounded pair-local coordinates.

    Every coordinate is a ratio of nonnegative masses or a bounded saturation of
    one, so all lie in [0,1]. Every coordinate sums over atoms and residues, so
    all are permutation invariant. No target, assay, task or document identifier
    is used.
    """
    value = np.abs(np.asarray(phi, dtype=np.float64))
    total = value.sum(axis=(1, 2, 3))
    safe = np.where(total > 1e-12, total, 1.0)

    def share(atoms, residues, bins) -> np.ndarray:
        block = value[:, atoms, :, :][:, :, residues, :][:, :, :, bins]
        return block.sum(axis=(1, 2, 3)) / safe

    columns = [
        total / (1.0 + total),
        value[:, :, :, CLOSE_BINS].sum(axis=(1, 2, 3)) / safe,
        share(APOLAR_ATOMS, APOLAR_RESIDUES, slice(None)),
        share(POLAR_ATOMS, POLAR_RESIDUES, CLOSE_BINS),
        share(CHARGED_ATOMS, CHARGED_RESIDUES, slice(None)),
        value[:, :, :, CONTACT_BIN].sum(axis=(1, 2, 3)) / safe,
        np.clip(np.asarray(heavy_atoms, dtype=np.float64) / MAX_LIGAND_ATOMS, 0.0, 1.0),
    ]
    out = np.stack(columns, axis=1)
    if out.shape[1] != len(Z_BIO_NAMES):
        raise RuntimeError("z_bio coordinate count disagrees with the frozen contract")
    return out


def assert_bounded_observable(features: np.ndarray, name: str) -> None:
    if not np.all(np.isfinite(features)):
        raise RuntimeError(f"{name} is not finite")
    if features.min() < -1e-9 or features.max() > 1.0 + 1e-9:
        raise RuntimeError(f"{name} left the declared [0,1] bounds")


def bounded_projection(train: np.ndarray, values: np.ndarray, dim: int) -> np.ndarray:
    """Unsupervised PCA to `dim`, then min-max bounded by the training range.

    Fitted on the training pool only, so every arm receives the same treatment
    and no evaluation label or fold influences the projection.
    """
    centre = train.mean(axis=0)
    _, _, right = np.linalg.svd(train - centre, full_matrices=False)
    basis = right[:dim].T
    projected_train = (train - centre) @ basis
    low = projected_train.min(axis=0)
    high = projected_train.max(axis=0)
    span = np.where(high - low > 1e-12, high - low, 1.0)
    return np.clip(((values - centre) @ basis - low) / span, 0.0, 1.0)


def knn_conditional_cdf(train_features: np.ndarray, train_y: np.ndarray,
                        query_features: np.ndarray, neighbours: int,
                        grid: np.ndarray) -> np.ndarray:
    """k-nearest-neighbour estimate of P(Y <= v | features) on the fixed grid.

    Dimension-robust replacement for the fixed-bandwidth kernel estimate, which
    returns nearly the marginal in seven dimensions. `neighbours` is fixed by the
    registered rule and is identical for every arm.
    """
    train_features = np.asarray(train_features, dtype=np.float64)
    query_features = np.asarray(query_features, dtype=np.float64)
    indicator = (np.asarray(train_y, dtype=np.float64)[:, None]
                 <= np.asarray(grid, dtype=np.float64)[None, :]).astype(np.float64)
    k = int(max(1, min(neighbours, len(train_features))))
    out = np.empty((len(query_features), len(grid)), dtype=np.float64)
    for index, query in enumerate(query_features):
        distance = np.sum(np.square(train_features - query), axis=1)
        chosen = np.argpartition(distance, k - 1)[:k]
        out[index] = np.maximum.accumulate(
            np.clip(indicator[chosen].mean(axis=0), 0.0, 1.0))
    return out


def band_mean_interval(beta: np.ndarray, grid: np.ndarray, a_max: float) -> np.ndarray:
    """Band-induced bounds on E[X]: E = a_max - integral F."""
    lower, upper = bands.split(np.atleast_2d(np.asarray(beta, dtype=np.float64)))
    return np.stack([a_max - np.trapz(upper, grid, axis=-1),
                     a_max - np.trapz(lower, grid, axis=-1)], axis=1)


def pooled_replicate_sigma(groups: list[np.ndarray]) -> tuple[float, int, int]:
    """Frozen sigma_assay estimator: pooled within-cell SD over replicate cells.

    A cell is one (task, ligand) pair, and a task key already fixes one exact
    assay, so within-cell replicates are within-assay replicates. Cells with
    fewer than two measurements are excluded and never treated as zero noise.
    """
    numerator = denominator = 0.0
    cells = 0
    for values in groups:
        if len(values) < 2:
            continue
        cells += 1
        numerator += float(np.sum(np.square(values - values.mean())))
        denominator += len(values) - 1
    if denominator <= 0:
        return float("nan"), 0, 0
    return float(np.sqrt(numerator / denominator)), cells, int(denominator)


def sigma_confidence(sigma: float, degrees: int, alpha: float = 0.05) -> tuple[float, float]:
    from scipy.stats import chi2
    if degrees <= 0 or not np.isfinite(sigma):
        return float("nan"), float("nan")
    lower = sigma * np.sqrt(degrees / chi2.ppf(1.0 - alpha / 2.0, degrees))
    upper = sigma * np.sqrt(degrees / chi2.ppf(alpha / 2.0, degrees))
    return float(lower), float(upper)


def component_macro(values: np.ndarray, components: np.ndarray) -> float:
    unique = np.unique(components)
    return float(np.mean([values[components == name].mean() for name in unique]))


def component_bootstrap(values: np.ndarray, components: np.ndarray,
                        draws: int = 2000, seed: int = 20260808) -> tuple[float, float]:
    unique = np.unique(components)
    cache = {name: values[components == name] for name in unique}
    rng = np.random.default_rng(seed)
    samples = np.empty(draws, dtype=np.float64)
    for index in range(draws):
        chosen = rng.choice(unique, size=len(unique), replace=True)
        samples[index] = np.mean([cache[name].mean() for name in chosen])
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))
