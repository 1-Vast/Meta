"""Five-channel product-kernel statistic frozen in the v2 preregistration."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from .mechanism import CHANNEL_NAMES, SITEALIGN


LAMBDA_GRID = (1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3, 1e4)


@dataclass(frozen=True)
class KernelLigand:
    molid: str
    smiles: str
    canonical_smiles: str
    generic_scaffold: str
    nuisance_features: np.ndarray
    channel_fingerprints: np.ndarray


def _feature_factory():
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def ligand_channel_fingerprints(mol, factory=None, n_bits: int = 1024) -> np.ndarray:
    """Return center-restricted Morgan bits in the five biological channels."""
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit import DataStructs

    if factory is None:
        factory = _feature_factory()
    atom_sets = {name: set() for name in (
        "Donor", "Acceptor", "PosIonizable", "NegIonizable",
        "Aromatic", "Hydrophobe", "LumpedHydrophobe",
    )}
    for feature in factory.GetFeaturesForMol(mol):
        family = feature.GetFamily()
        if family in atom_sets:
            atom_sets[family].update(feature.GetAtomIds())
    atom_sets["Hydrophobe"].update(atom_sets["LumpedHydrophobe"])
    centres = (
        atom_sets["Donor"] | atom_sets["Acceptor"],
        atom_sets["PosIonizable"] | atom_sets["NegIonizable"],
        atom_sets["Aromatic"],
        atom_sets["Hydrophobe"],
        {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1},
    )
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
    out = np.zeros((len(CHANNEL_NAMES), n_bits), dtype=np.uint8)
    for channel, atom_ids in enumerate(centres):
        if not atom_ids:
            continue
        bit_vector = generator.GetFingerprint(mol, fromAtoms=sorted(atom_ids))
        DataStructs.ConvertToNumpyArray(bit_vector, out[channel])
    return out


def kernel_ligand_from_smiles(molid: str, smiles: str, factory=None) -> KernelLigand:
    from rdkit import Chem
    from .mechanism import ligand_from_smiles

    base = ligand_from_smiles(molid, smiles, factory=factory)
    mol = Chem.MolFromSmiles(base.canonical_smiles)
    if mol is None:
        raise ValueError("RDKit could not reconstruct canonical SMILES")
    return KernelLigand(
        molid=base.molid,
        smiles=base.smiles,
        canonical_smiles=base.canonical_smiles,
        generic_scaffold=base.generic_scaffold,
        nuisance_features=base.model_features,
        channel_fingerprints=ligand_channel_fingerprints(mol, factory=factory),
    )


def protein_channel_features(pocket: str) -> tuple[np.ndarray, ...]:
    """Aligned SiteAlign rows for the five preregistered protein channels."""
    if len(pocket) != 85:
        raise ValueError("KLIFS pocket must have exactly 85 positions")
    values = np.asarray([SITEALIGN.get(aa, (0.0,) * 6) for aa in pocket], dtype=np.float64)
    size, donor, acceptor, charge, aromatic, aliphatic = values.T
    positive = np.clip(charge, 0.0, None)
    negative = np.clip(-charge, 0.0, None)
    free_space = np.clip((3.0 - size) / 2.0, 0.0, 1.0)
    return (
        np.stack([donor / 3.0, acceptor / 2.0]).reshape(-1),
        np.stack([positive, negative]).reshape(-1),
        aromatic.reshape(-1),
        aliphatic.reshape(-1),
        np.stack([size / 3.0, free_space]).reshape(-1),
    )


def tanimoto_kernel(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    if left.ndim != 2 or right.ndim != 2 or left.shape[1] != right.shape[1]:
        raise ValueError("Tanimoto inputs must be two matrices with equal width")
    intersection = left @ right.T
    union = left.sum(axis=1, keepdims=True) + right.sum(axis=1)[None, :] - intersection
    return np.divide(intersection, union, out=np.zeros_like(intersection), where=union > 0.0)


def squared_distance(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    value = (
        np.square(left).mean(axis=1, keepdims=True)
        + np.square(right).mean(axis=1)[None, :]
        - 2.0 * (left @ right.T) / left.shape[1]
    )
    return np.clip(value, 0.0, None)


def source_rbf(features: np.ndarray) -> tuple[np.ndarray, float]:
    distance = squared_distance(features, features)
    upper = distance[np.triu_indices(len(distance), 1)]
    nonzero = upper[upper > 1e-12]
    tau = float(np.median(nonzero)) if len(nonzero) else 1.0
    return np.exp(-distance / tau), tau


def cross_rbf(left: np.ndarray, right: np.ndarray, tau: float) -> np.ndarray:
    if tau <= 0.0:
        raise ValueError("RBF scale must be positive")
    return np.exp(-squared_distance(left, right) / float(tau))


def fit_separable_krr(
    ligand_kernel: np.ndarray,
    protein_kernel: np.ndarray,
    residual: np.ndarray,
    regularization: float,
) -> np.ndarray:
    """Solve K_L C K_P + lambda C = R by symmetric eigendecomposition."""
    ligand_kernel = np.asarray(ligand_kernel, dtype=np.float64)
    protein_kernel = np.asarray(protein_kernel, dtype=np.float64)
    residual = np.asarray(residual, dtype=np.float64)
    if residual.shape != (len(ligand_kernel), len(protein_kernel)):
        raise ValueError("residual rectangle and kernel dimensions disagree")
    sl, ul = np.linalg.eigh((ligand_kernel + ligand_kernel.T) / 2.0)
    sp, up = np.linalg.eigh((protein_kernel + protein_kernel.T) / 2.0)
    sl = np.clip(sl, 0.0, None)
    sp = np.clip(sp, 0.0, None)
    rotated = ul.T @ residual @ up
    denominator = sl[:, None] * sp[None, :] + float(regularization)
    return ul @ (rotated / denominator) @ up.T


def predict_separable(
    ligand_cross_kernel: np.ndarray,
    coefficients: np.ndarray,
    protein_cross_kernel: np.ndarray,
) -> np.ndarray:
    return np.asarray(ligand_cross_kernel) @ coefficients @ np.asarray(protein_cross_kernel)


def fit_convex_channel_weights(prediction: np.ndarray, outcome: np.ndarray) -> dict:
    """Least-squares weights with w >= 0 and sum(w) <= 1; zero is feasible."""
    prediction = np.asarray(prediction, dtype=np.float64)
    outcome = np.asarray(outcome, dtype=np.float64).reshape(-1)
    if prediction.ndim != 2 or prediction.shape[0] != len(outcome):
        raise ValueError("prediction/outcome dimensions disagree")
    n_channels = prediction.shape[1]

    def objective(weight):
        error = prediction @ weight - outcome
        return float(error @ error)

    constraints = ({"type": "ineq", "fun": lambda weight: 1.0 - float(np.sum(weight))},)
    fitted = minimize(
        objective,
        np.full(n_channels, 1.0 / (2.0 * n_channels)),
        method="SLSQP",
        bounds=[(0.0, 1.0)] * n_channels,
        constraints=constraints,
        options={"ftol": 1e-12, "maxiter": 2000},
    )
    if not fitted.success:
        raise RuntimeError(f"convex channel fit failed: {fitted.message}")
    weight = np.clip(fitted.x, 0.0, 1.0)
    if weight.sum() > 1.0:
        weight /= weight.sum()
    zero_mse = float(np.mean(np.square(outcome)))
    fitted_mse = float(np.mean(np.square(prediction @ weight - outcome)))
    return {
        "weights": weight,
        "sum": float(weight.sum()),
        "zero_mse": zero_mse,
        "fitted_mse": fitted_mse,
        "mse_reduction": zero_mse - fitted_mse,
        "optimizer_message": str(fitted.message),
    }


def fit_bounds(contributions: np.ndarray) -> dict:
    values = np.asarray(contributions, dtype=np.float64).reshape(-1, len(CHANNEL_NAMES))
    low = np.quantile(values, 0.01, axis=0)
    high = np.quantile(values, 0.99, axis=0)
    span = np.where(high - low > 1e-8, high - low, 1.0)
    return {"low": low.tolist(), "high": high.tolist(), "span": span.tolist()}


def bounded_z(contributions: np.ndarray, bounds: dict) -> np.ndarray:
    values = np.asarray(contributions, dtype=np.float64)
    low = np.asarray(bounds["low"], dtype=np.float64)
    span = np.asarray(bounds["span"], dtype=np.float64)
    transformed = 2.0 * (values - low) / span - 1.0
    return 0.5 * (np.tanh(transformed) + 1.0)
