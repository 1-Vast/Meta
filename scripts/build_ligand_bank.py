"""Build fixed-schema RDKit ligand graph shards from canonical DTA rows."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM, GRAPH_SCHEMA, MAX_ATOMS
from scripts.data_contract import read_jsonl


ELEMENTS = ("B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "Si")


def _one_hot(index: int, size: int) -> np.ndarray:
    value = np.zeros(size, dtype=np.float32)
    value[min(max(index, 0), size - 1)] = 1.0
    return value


def featurize_molecule(molecule, max_atoms: int = MAX_ATOMS,
                       atom_indices: list[int] | None = None) -> dict[str, np.ndarray]:
    from rdkit import Chem
    indices = list(range(molecule.GetNumAtoms())) if atom_indices is None else list(atom_indices)
    if len(indices) != len(set(indices)) or any(index < 0 or index >= molecule.GetNumAtoms()
                                                for index in indices):
        raise ValueError("atom_indices must be a valid one-to-one molecule ordering")
    count = len(indices)
    if count > max_atoms:
        raise ValueError(f"ligand has {count} atoms; max_atoms={max_atoms}")
    position = {atom_index: output_index for output_index, atom_index in enumerate(indices)}
    X = np.zeros((max_atoms, ATOM_FEAT_DIM), dtype=np.float32)
    A = np.zeros((max_atoms, max_atoms, BOND_FEAT_DIM), dtype=np.float32)
    mask = np.zeros(max_atoms, dtype=np.float32)
    hybrid = {Chem.rdchem.HybridizationType.SP: 0, Chem.rdchem.HybridizationType.SP2: 1,
              Chem.rdchem.HybridizationType.SP3: 2}
    chiral = {Chem.rdchem.ChiralType.CHI_UNSPECIFIED: 0,
              Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW: 1,
              Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW: 2}
    for atom_index in indices:
        atom = molecule.GetAtomWithIdx(atom_index)
        index = position[atom_index]
        symbol = atom.GetSymbol()
        element = ELEMENTS.index(symbol) if symbol in ELEMENTS else len(ELEMENTS)
        charge = {-2: 0, -1: 1, 0: 2, 1: 3}.get(atom.GetFormalCharge(), 4)
        graph_degree = sum(neighbor.GetIdx() in position for neighbor in atom.GetNeighbors())
        features = np.concatenate((
            _one_hot(element, 12), _one_hot(graph_degree, 6), _one_hot(charge, 5),
            _one_hot(hybrid.get(atom.GetHybridization(), 3), 4),
            np.asarray([atom.GetIsAromatic(), atom.IsInRing()], dtype=np.float32),
            _one_hot(chiral.get(atom.GetChiralTag(), 0), 3),
        ))
        X[index] = features
        mask[index] = 1.0
    bond_type = {Chem.rdchem.BondType.SINGLE: 0, Chem.rdchem.BondType.DOUBLE: 1,
                 Chem.rdchem.BondType.TRIPLE: 2, Chem.rdchem.BondType.AROMATIC: 3}
    stereo = {Chem.rdchem.BondStereo.STEREONONE: 0, Chem.rdchem.BondStereo.STEREOZ: 1,
              Chem.rdchem.BondStereo.STEREOE: 2}
    direction = {Chem.rdchem.BondDir.NONE: 0, Chem.rdchem.BondDir.ENDUPRIGHT: 1}
    for bond in molecule.GetBonds():
        raw_left, raw_right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if raw_left not in position or raw_right not in position:
            continue
        features = np.concatenate((
            _one_hot(bond_type.get(bond.GetBondType(), 0), 4),
            _one_hot(stereo.get(bond.GetStereo(), 3), 4),
            np.asarray([bond.GetIsConjugated(), bond.IsInRing()], dtype=np.float32),
            _one_hot(direction.get(bond.GetBondDir(), 1), 2),
        ))
        left, right = position[raw_left], position[raw_right]
        A[left, right] = features
        A[right, left] = features
    return {"X": X, "A": A, "mask": mask}


def featurize_smiles(smiles: str, max_atoms: int = MAX_ATOMS) -> dict[str, np.ndarray]:
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError("RDKit is required to build the ligand bank") from error
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("invalid canonical SMILES")
    return featurize_molecule(molecule, max_atoms=max_atoms)


def build_ligand_bank(rows_path: str | Path, output_dir: str | Path, *, shard_size: int = 2048,
                      max_atoms: int = 128) -> dict:
    rows = read_jsonl(rows_path)
    molecules = {row["drug_key"]: row["smiles"] for row in rows}
    items = sorted(molecules.items())
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    shards = []
    for start in range(0, len(items), shard_size):
        batch = items[start:start + shard_size]
        graphs = [featurize_smiles(smiles, max_atoms=max_atoms) for _, smiles in batch]
        filename = f"shard_{start:08d}_{start + len(batch):08d}.npz"
        np.savez_compressed(root / filename, keys=np.asarray([key for key, _ in batch]),
                            X=np.stack([graph["X"] for graph in graphs]),
                            A=np.stack([graph["A"] for graph in graphs]),
                            mask=np.stack([graph["mask"] for graph in graphs]))
        shards.append(filename)
    manifest = {"schema": GRAPH_SCHEMA, "atom_feature_dim": ATOM_FEAT_DIM,
                "bond_feature_dim": BOND_FEAT_DIM,
                "max_atoms": max_atoms, "molecules": len(items), "shards": shards}
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def load_ligand_bank(directory: str | Path) -> dict[str, dict[str, np.ndarray]]:
    """Load and validate a prebuilt graph bank without importing RDKit."""
    root = Path(directory)
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {"schema": GRAPH_SCHEMA, "atom_feature_dim": ATOM_FEAT_DIM,
                "bond_feature_dim": BOND_FEAT_DIM}
    if any(manifest.get(key) != value for key, value in expected.items()):
        raise ValueError("ligand bank graph schema does not match the runtime contract")
    if int(manifest.get("max_atoms", -1)) != MAX_ATOMS:
        raise ValueError("ligand bank max_atoms does not match the runtime contract")
    graphs = {}
    for filename in manifest.get("shards", []):
        with np.load(root / filename, allow_pickle=False) as shard:
            keys, X, A, mask = shard["keys"], shard["X"], shard["A"], shard["mask"]
            if X.shape[1:] != (MAX_ATOMS, ATOM_FEAT_DIM) or \
                    A.shape[1:] != (MAX_ATOMS, MAX_ATOMS, BOND_FEAT_DIM) or \
                    mask.shape[1:] != (MAX_ATOMS,):
                raise ValueError("ligand bank shard tensor shape does not match its schema")
            for index, key in enumerate(keys):
                key = str(key)
                if key in graphs:
                    raise ValueError(f"duplicate ligand graph key {key}")
                graphs[key] = {"X": X[index], "A": A[index], "mask": mask[index]}
    if len(graphs) != int(manifest.get("molecules", -1)):
        raise ValueError("ligand bank manifest molecule count does not match its shards")
    return graphs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("rows")
    parser.add_argument("output")
    parser.add_argument("--shard-size", type=int, default=2048)
    parser.add_argument("--max-atoms", type=int, default=128)
    args = parser.parse_args()
    print(json.dumps(build_ligand_bank(args.rows, args.output, shard_size=args.shard_size,
                                       max_atoms=args.max_atoms), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
