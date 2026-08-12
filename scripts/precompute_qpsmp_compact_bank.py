"""Pack padded ligand graphs once so QPSMP training is not NPZ-I/O bound."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def pack_shard(source: Path, target: Path) -> dict:
    with np.load(source, allow_pickle=False) as stored:
        keys = stored["keys"].astype(str)
        mask = stored["mask"]
        sizes = mask.sum(axis=1).astype(np.int32)
        atom_offsets = np.concatenate(([0], np.cumsum(sizes, dtype=np.int64)))
        bond_counts = sizes.astype(np.int64) ** 2
        bond_offsets = np.concatenate(([0], np.cumsum(bond_counts, dtype=np.int64)))
        atoms = stored["X"]
        bonds = stored["A"]
        packed_atoms = np.concatenate(
            [atoms[index, :size] for index, size in enumerate(sizes)], axis=0)
        packed_bonds = np.concatenate(
            [bonds[index, :size, :size].reshape(size * size, bonds.shape[-1])
             for index, size in enumerate(sizes)], axis=0)
    np.savez_compressed(
        target, keys=keys, sizes=sizes, atom_offsets=atom_offsets,
        bond_offsets=bond_offsets, X=packed_atoms, A=packed_bonds)
    return {"path": target.name, "records": int(len(keys))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    source_manifest = json.loads(
        (args.source / "manifest.json").read_text(encoding="utf-8"))
    args.output.mkdir(parents=True)
    shards = []
    for item in source_manifest["shards"]:
        name = item["path"] if isinstance(item, dict) else item
        shards.append(pack_shard(args.source / name, args.output / name))
        print(json.dumps(shards[-1]), flush=True)
    manifest = {
        "schema": "MetaSieve.QPSMPCompactLigandBank.v1",
        "source_manifest": str((args.source / "manifest.json").resolve()),
        "source_schema": source_manifest["schema"],
        "records": sum(item["records"] for item in shards),
        "atom_feature_dim": source_manifest["atom_feature_dim"],
        "bond_feature_dim": source_manifest["bond_feature_dim"],
        "shards": shards,
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
