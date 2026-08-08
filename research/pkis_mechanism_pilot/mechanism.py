"""Pure chemistry/protein feature construction for the PKIS external pilot.

The module intentionally contains no affinity labels, target identifiers in a
feature vector, learned embeddings, or frozen-operator modifications.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable

import numpy as np


CHANNEL_NAMES = (
    "hbond_complementarity",
    "ionic_complementarity",
    "aromatic_packing",
    "hydrophobic_packing",
    "steric_accommodation",
)
SHELL_NAMES = ("core", "one_bond", "two_bonds", "three_or_more_bonds")
LIGAND_FEATURE_NAMES = (
    "donor", "acceptor", "positive", "negative", "aromatic", "hydrophobe", "occupancy"
)

# SiteAlign residue definitions, represented in one-letter form. Values are
# size, HBD, HBA, charge, aromatic, aliphatic. The table is the public KiSSim /
# SiteAlign encoding, copied here so this pilot has no runtime dependency on
# structure-download machinery.
SITEALIGN = {
    "A": (1.0, 0.0, 0.0, 0.0, 0.0, 1.0),
    "R": (3.0, 3.0, 0.0, 1.0, 0.0, 0.0),
    "N": (2.0, 1.0, 1.0, 0.0, 0.0, 0.0),
    "D": (2.0, 0.0, 2.0, -1.0, 0.0, 0.0),
    "C": (1.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    "Q": (2.0, 1.0, 1.0, 0.0, 0.0, 0.0),
    "E": (2.0, 0.0, 2.0, -1.0, 0.0, 0.0),
    "G": (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    "H": (2.0, 1.0, 1.0, 0.0, 1.0, 0.0),
    "I": (2.0, 0.0, 0.0, 0.0, 0.0, 1.0),
    "L": (2.0, 0.0, 0.0, 0.0, 0.0, 1.0),
    "K": (2.0, 1.0, 0.0, 1.0, 0.0, 0.0),
    "M": (2.0, 0.0, 0.0, 0.0, 0.0, 1.0),
    "F": (3.0, 0.0, 0.0, 0.0, 1.0, 0.0),
    "P": (1.0, 0.0, 0.0, 0.0, 0.0, 1.0),
    "S": (1.0, 1.0, 1.0, 0.0, 0.0, 0.0),
    "T": (1.0, 1.0, 1.0, 0.0, 0.0, 1.0),
    "W": (3.0, 1.0, 0.0, 0.0, 1.0, 0.0),
    "Y": (3.0, 1.0, 1.0, 0.0, 1.0, 0.0),
    "V": (1.0, 0.0, 0.0, 0.0, 0.0, 1.0),
}


# Historical assay names that are unambiguous gene synonyms. Assay constructs
# and mutants are filtered before this table is consulted.
TARGET_ALIASES = {
    "ARG": "ABL2",
    "ARK5": "NUAK1",
    "AURORAA": "AURKA",
    "AURORAB": "AURKB",
    "AURORAC": "AURKC",
    "CK1A": "CSNK1A1",
    "CK1G1": "CSNK1G1",
    "CK1G2": "CSNK1G2",
    "CK1G3": "CSNK1G3",
    "CRAF": "RAF1",
    "DCAMKL1": "DCLK1",
    "DCAMKL2": "DCLK2",
    "DCAMKL3": "DCLK3",
    "ERK1": "MAPK3",
    "ERK2": "MAPK1",
    "ERK3": "MAPK6",
    "ERK4": "MAPK4",
    "ERK5": "MAPK7",
    "ERK8": "MAPK15",
    "FMS": "CSF1R",
    "IKKALPHA": "CHUK",
    "IKKBETA": "IKBKB",
    "IKKEPSILON": "IKBKE",
    "IRR": "INSRR",
    "KDR": "KDR",
    "MEK1": "MAP2K1",
    "MEK2": "MAP2K2",
    "MEK3": "MAP2K3",
    "MEK4": "MAP2K4",
    "MEK5": "MAP2K5",
    "MEK6": "MAP2K6",
    "MKK7": "MAP2K7",
    "MLCK": "MYLK",
    "MNK2": "MKNK2",
    "P38ALPHA": "MAPK14",
    "P38BETA": "MAPK11",
    "P38DELTA": "MAPK13",
    "P38GAMMA": "MAPK12",
    "P70S6K1": "RPS6KB1",
    "PAK7": "PAK5",
    "PCTK1": "CDK16",
    "PCTK2": "CDK17",
    "PCTK3": "CDK18",
    "PDPK1": "PDPK1",
    "PFTK1": "CDK14",
    "PKA": "PRKACA",
    "PKACALPHA": "PRKACA",
    "PKACBETA": "PRKACB",
    "PRAK": "MAPKAPK5",
    "PRKR": "EIF2AK2",
    "QSK": "KIAA0999",
    "S6K1": "RPS6KB1",
    "SGK": "SGK1",
    "SNARK": "NUAK2",
    "TAK1": "MAP3K7",
    "TRKA": "NTRK1",
    "TRKB": "NTRK2",
    "TRKC": "NTRK3",
    "VEGFR2": "KDR",
    "YSK1": "STK25",
    "YSK4": "MAP3K19",
}

_CONSTRUCT_PATTERNS = (
    r"(?:^|[_.-])[A-Z]\d+[A-Z](?:$|[_.-])",  # point mutant
    r"CYCLIN", r"PHOSPHORYLATED", r"AUTOINHIBITED", r"PSEUDOKINASE",
    r"KINDOM", r"DOMAIN", r"_P35", r"A1B1G1", r"A2B1G1",
    r"LYNA$", r"LYNB$",
)


def normalize_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def is_construct_or_variant(value: str) -> bool:
    upper = str(value).upper()
    return any(re.search(pattern, upper) is not None for pattern in _CONSTRUCT_PATTERNS)


@dataclass(frozen=True)
class TargetRecord:
    kinase_id: int
    name: str
    hgnc: str
    family: str
    group: str
    uniprot: str
    pocket: str


def load_klifs(path: str | Path) -> tuple[dict[str, TargetRecord], dict[str, list[TargetRecord]]]:
    rows = json.loads(Path(path).read_text(encoding="utf-8"))
    valid = []
    for row in rows:
        pocket = row.get("pocket")
        if row.get("species") != "Human" or not isinstance(pocket, str) or len(pocket) != 85:
            continue
        valid.append(TargetRecord(
            kinase_id=int(row["kinase_ID"]), name=str(row.get("name", "")),
            hgnc=str(row.get("HGNC", "")), family=str(row.get("family", "")),
            group=str(row.get("group", "")), uniprot=str(row.get("uniprot", "")),
            pocket=pocket,
        ))
    by_key: dict[str, list[TargetRecord]] = {}
    for record in valid:
        for key in {normalize_name(record.name), normalize_name(record.hgnc)}:
            if key:
                by_key.setdefault(key, []).append(record)
    unique = {}
    for key, records in by_key.items():
        ids = {record.kinase_id for record in records}
        if len(ids) == 1:
            unique[key] = records[0]
    return unique, by_key


def map_target(value: str, unique_index: dict[str, TargetRecord]) -> tuple[TargetRecord | None, str]:
    if is_construct_or_variant(value):
        return None, "construct_or_variant"
    key = normalize_name(value)
    query = TARGET_ALIASES.get(key, key)
    record = unique_index.get(normalize_name(query))
    if record is None:
        return None, "no_unambiguous_klifs_pocket"
    return record, "included"


def pocket_sitealign(pocket: str) -> np.ndarray:
    if len(pocket) != 85:
        raise ValueError("KLIFS pocket must have exactly 85 aligned positions")
    return np.asarray([SITEALIGN.get(aa, (0.0,) * 6) for aa in pocket], dtype=np.float32)


def protein_pair_properties(pocket: str) -> np.ndarray:
    """Return seven analytic property rows used by the five pair channels."""
    values = pocket_sitealign(pocket)
    size, hbd, hba, charge, aromatic, aliphatic = values.T
    positive = np.clip(charge, 0.0, None)
    negative = np.clip(-charge, 0.0, None)
    inverse_size = np.clip((3.0 - size) / 2.0, 0.0, 1.0)
    return np.stack([hbd, hba, positive, negative, aromatic, aliphatic, inverse_size])


def protein_model_features(pocket: str) -> np.ndarray:
    values = pocket_sitealign(pocket).copy()
    values[:, 0] /= 3.0
    values[:, 1] /= 3.0
    values[:, 2] /= 2.0
    return values.reshape(-1)


@dataclass
class LigandRecord:
    molid: str
    smiles: str
    canonical_smiles: str
    generic_scaffold: str
    model_features: np.ndarray
    pharmacophore_shells: np.ndarray


def _feature_factory():
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def _murcko(mol):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        return scaffold, "ACYCLIC"
    generic = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return scaffold, Chem.MolToSmiles(generic, canonical=True)


def _shells(mol, scaffold) -> np.ndarray:
    from rdkit import Chem
    n = mol.GetNumAtoms()
    if scaffold.GetNumAtoms() == 0:
        return np.zeros(n, dtype=np.int64)
    core = mol.GetSubstructMatch(scaffold)
    if not core:
        raise ValueError("Murcko scaffold did not map back to its parent molecule")
    distance = np.asarray(Chem.GetDistanceMatrix(mol), dtype=np.float64)
    return np.minimum(distance[:, list(core)].min(axis=1).astype(np.int64), 3)


def ligand_from_smiles(molid: str, smiles: str, factory=None) -> LigandRecord:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    scaffold, generic_scaffold = _murcko(mol)
    shell = _shells(mol, scaffold)
    heavy = max(1, mol.GetNumHeavyAtoms())

    if factory is None:
        factory = _feature_factory()
    atom_sets = {name: set() for name in ("Donor", "Acceptor", "PosIonizable", "NegIonizable",
                                           "Aromatic", "Hydrophobe", "LumpedHydrophobe")}
    for feature in factory.GetFeaturesForMol(mol):
        if feature.GetFamily() in atom_sets:
            atom_sets[feature.GetFamily()].update(feature.GetAtomIds())
    atom_sets["Hydrophobe"].update(atom_sets["LumpedHydrophobe"])
    families = ("Donor", "Acceptor", "PosIonizable", "NegIonizable", "Aromatic", "Hydrophobe")
    pharm = np.zeros((7, 4), dtype=np.float32)
    for index, family in enumerate(families):
        for atom_id in atom_sets[family]:
            pharm[index, shell[atom_id]] += 1.0 / heavy
    for atom_id in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(atom_id).GetAtomicNum() > 1:
            pharm[6, shell[atom_id]] += 1.0 / heavy

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    bitvect = generator.GetFingerprint(mol)
    fingerprint = np.zeros(1024, dtype=np.float32)
    DataStructs.ConvertToNumpyArray(bitvect, fingerprint)
    aromatic_atoms = sum(int(atom.GetIsAromatic()) for atom in mol.GetAtoms())
    carbon_atoms = sum(int(atom.GetAtomicNum() == 6) for atom in mol.GetAtoms())
    formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    scalars = np.asarray([
        min(Descriptors.MolWt(mol) / 800.0, 1.5),
        np.clip((Crippen.MolLogP(mol) + 5.0) / 12.0, 0.0, 1.5),
        min(Descriptors.TPSA(mol) / 250.0, 1.5),
        min(Lipinski.NumHDonors(mol) / 10.0, 1.5),
        min(Lipinski.NumHAcceptors(mol) / 20.0, 1.5),
        min(Lipinski.NumRotatableBonds(mol) / 20.0, 1.5),
        min(Lipinski.RingCount(mol) / 10.0, 1.5),
        min(heavy / 100.0, 1.5),
        carbon_atoms / heavy,
        aromatic_atoms / heavy,
        Descriptors.FractionCSP3(mol),
        np.clip(formal_charge / 5.0, -1.0, 1.0),
    ], dtype=np.float32)
    return LigandRecord(
        molid=str(molid), smiles=str(smiles), canonical_smiles=canonical,
        generic_scaffold=generic_scaffold,
        model_features=np.concatenate([fingerprint, scalars]),
        pharmacophore_shells=pharm,
    )


def load_smiles(path: str | Path) -> dict[str, str]:
    out = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            smiles, molid = line.rsplit(maxsplit=1)
            out[str(molid)] = smiles
    return out


def pair_feature_matrix(ligand_shells: np.ndarray, protein_properties: np.ndarray) -> np.ndarray:
    """Return `(ligand, target, 5*4*85)` analytic pair features."""
    ligand_shells = np.asarray(ligand_shells, dtype=np.float32)
    protein_properties = np.asarray(protein_properties, dtype=np.float32)
    if ligand_shells.ndim != 3 or ligand_shells.shape[1:] != (7, 4):
        raise ValueError("ligand shells must have shape (n_ligand, 7, 4)")
    if protein_properties.ndim != 3 or protein_properties.shape[1:] != (7, 85):
        raise ValueError("protein properties must have shape (n_target, 7, 85)")
    donor, acceptor, pos, neg, aromatic, hydrophobe, occupancy = (
        ligand_shells[:, index] for index in range(7)
    )
    p_hbd, p_hba, p_pos, p_neg, p_arom, p_aliph, p_space = (
        protein_properties[:, index] for index in range(7)
    )

    def outer(ligand_value, protein_value):
        return np.einsum("lb,tj->ltbj", ligand_value, protein_value, optimize=True)

    channels = [
        outer(acceptor, p_hbd) + outer(donor, p_hba),
        outer(neg, p_pos) + outer(pos, p_neg),
        outer(aromatic, p_arom),
        outer(hydrophobe, p_aliph),
        outer(occupancy, p_space),
    ]
    tensor = np.stack(channels, axis=2)
    return np.ascontiguousarray(tensor.reshape(len(ligand_shells), len(protein_properties), -1))


def channel_contributions(features: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    features = np.asarray(features, dtype=np.float64)
    coefficients = np.asarray(coefficients, dtype=np.float64)
    expected = len(CHANNEL_NAMES) * len(SHELL_NAMES) * 85
    if features.shape[-1] != expected or coefficients.shape != (expected,):
        raise ValueError("pair feature/coefficient dimensions disagree with the frozen pilot tensor")
    shape = features.shape[:-1]
    block = len(SHELL_NAMES) * 85
    values = features.reshape(*shape, len(CHANNEL_NAMES), block)
    weights = coefficients.reshape(len(CHANNEL_NAMES), block)
    return np.einsum("...cb,cb->...c", values, weights, optimize=True)


def fit_channel_bounds(source_contributions: np.ndarray) -> dict:
    values = np.asarray(source_contributions, dtype=np.float64).reshape(-1, len(CHANNEL_NAMES))
    low = np.quantile(values, 0.01, axis=0)
    high = np.quantile(values, 0.99, axis=0)
    span = np.where(high - low > 1e-8, high - low, 1.0)
    return {"low": low.tolist(), "high": high.tolist(), "span": span.tolist()}


def bounded_biological_z(contributions: np.ndarray, bounds: dict) -> np.ndarray:
    values = np.asarray(contributions, dtype=np.float64)
    low = np.asarray(bounds["low"], dtype=np.float64)
    span = np.asarray(bounds["span"], dtype=np.float64)
    # Smoothly bounded; 1st and 99th source quantiles map near 0.12 and 0.88.
    centre_scaled = 2.0 * ((values - low) / span) - 1.0
    return 0.5 * (np.tanh(centre_scaled) + 1.0)


def ordered_anchor_simplex(location: np.ndarray, uncertainty: np.ndarray | float = 0.0) -> np.ndarray:
    """Map a bounded location to the existing six ordered anchors plus width.

    Column zero is reserved for the population band by the frozen operator.
    Columns one through six are the ordered logistic anchors; column seven is
    the broad uniform/abstention anchor.
    """
    centres = np.asarray([0.15, 0.30, 0.45, 0.60, 0.75, 0.90], dtype=np.float64)
    value = np.clip(np.asarray(location, dtype=np.float64), 0.0, 1.0).reshape(-1)
    width = np.broadcast_to(np.clip(np.asarray(uncertainty, dtype=np.float64), 0.0, 1.0), value.shape)
    ordered = np.zeros((len(value), len(centres)), dtype=np.float64)
    for row, item in enumerate(value):
        if item <= centres[0]:
            ordered[row, 0] = 1.0
        elif item >= centres[-1]:
            ordered[row, -1] = 1.0
        else:
            upper = int(np.searchsorted(centres, item, side="right"))
            lower = upper - 1
            fraction = (item - centres[lower]) / (centres[upper] - centres[lower])
            ordered[row, lower] = 1.0 - fraction
            ordered[row, upper] = fraction
    out = np.zeros((len(value), 8), dtype=np.float64)
    out[:, 1:7] = ordered * (1.0 - width[:, None])
    out[:, 7] = width
    return out


def stable_fold(value: str, n_folds: int, seed: int = 20260808) -> int:
    digest = hashlib.sha256(f"{seed}:{value}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % int(n_folds)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def double_center(matrix: np.ndarray) -> np.ndarray:
    value = np.asarray(matrix, dtype=np.float64)
    return value - value.mean(axis=1, keepdims=True) - value.mean(axis=0, keepdims=True) + value.mean()


def deterministic_derangement(records: Iterable[TargetRecord]) -> np.ndarray:
    records = list(records)
    if len(records) < 2:
        raise ValueError("derangement needs at least two targets")
    by_group: dict[str, list[int]] = {}
    for index, record in enumerate(records):
        by_group.setdefault(record.group, []).append(index)
    mapping = np.full(len(records), -1, dtype=np.int64)
    singleton = []
    for indices in by_group.values():
        if len(indices) == 1:
            singleton.extend(indices)
        else:
            ordered = sorted(indices, key=lambda index: (records[index].hgnc, records[index].kinase_id))
            for left, right in zip(ordered, ordered[1:] + ordered[:1]):
                mapping[left] = right
    if singleton:
        if len(singleton) == 1:
            source = singleton[0]
            mapping[source] = (source + 1) % len(records)
        else:
            ordered = sorted(singleton, key=lambda index: (records[index].hgnc, records[index].kinase_id))
            for left, right in zip(ordered, ordered[1:] + ordered[:1]):
                mapping[left] = right
    if np.any(mapping < 0) or np.any(mapping == np.arange(len(records))):
        raise RuntimeError("deterministic mapping is not a complete derangement")
    return mapping
