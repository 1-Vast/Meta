"""Build frozen-contract ligand graphs for the chemistry-novel R0-C panel."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import torch

from scripts.build_holo_complex_index import _ccd_molecule
from scripts.build_ligand_bank import featurize_molecule
from scripts.build_structure_supervision import _sha256_json
from scripts.data_contract import read_jsonl
from scripts.structure_sources.rcsb import sha256_file


FROZEN_PANEL_SHA256 = "fb7232ae70c974fc5c4e6ff8ef5517cd5a982cebf10cbcb80e651e8e043b48b9"


def _graph(record: dict) -> dict:
    ccd_path = Path(record["ccd_path"])
    if sha256_file(ccd_path) != record["ccd_sha256"]:
        raise ValueError("R0-C CCD hash differs from its governed record")
    chemistry = _ccd_molecule(ccd_path)
    molecule = chemistry["molecule"]
    index_by_name = {
        atom.GetProp("_CCDAtomName"): atom.GetIdx() for atom in molecule.GetAtoms()
    }
    names = list(chemistry["heavy_atom_names"])
    indices = [index_by_name[name] for name in names]
    features = featurize_molecule(molecule, atom_indices=indices)
    rows, columns = np.nonzero(features["A"].sum(axis=-1))
    return {
        "X": torch.from_numpy(features["X"]),
        "mask": torch.from_numpy(features["mask"]),
        "edge_index": torch.from_numpy(np.stack((rows, columns))).long(),
        "edge_attr": torch.from_numpy(features["A"][rows, columns]),
        "atom_mapping_hash": _sha256_json(names),
        "ligand_comp_id": record["ligand_comp_id"],
    }


def build(panel_path: str | Path, output_path: str | Path) -> dict:
    panel_file, output = Path(panel_path), Path(output_path)
    manifest_path = output.with_suffix(output.suffix + ".manifest.json")
    if output.exists() or manifest_path.exists():
        raise FileExistsError(f"R0-C ligand bank output already exists: {output}")
    if sha256_file(panel_file) != FROZEN_PANEL_SHA256:
        raise ValueError("R0-C ligand panel does not match the frozen SHA256")
    records = read_jsonl(panel_file)
    if len(records) < 120:
        raise ValueError("R0-C ligand panel is below the frozen component Gate")
    bank = {}
    for record in sorted(records, key=lambda row: row["ccd_sha256"]):
        key = str(record["ccd_sha256"])
        if key in bank:
            raise ValueError("R0-C final panel contains a duplicate CCD hash")
        graph = _graph(record)
        if int(graph["mask"].sum()) != int(record["ligand_heavy_atoms"]):
            raise ValueError("R0-C ligand graph atom count differs from governed record")
        bank[key] = graph
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    torch.save(bank, temporary)
    temporary.replace(output)
    manifest = {
        "schema": "MetaSieve.R0C.LigandGraphBank.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records": len(records),
        "graphs": len(bank),
        "panel_path": str(panel_file.resolve()),
        "panel_sha256": sha256_file(panel_file),
        "bank_path": str(output.resolve()),
        "bank_sha256": sha256_file(output),
        "graph_contract": "frozen_P1B_GINE_X32_edge12_max_atoms128_CCD_canonical_heavy_order",
        "affinity_labels_used": False,
        "distance_values_used": False,
        "trainable_parameters": 0,
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("panel")
    parser.add_argument("output")
    args = parser.parse_args()
    print(json.dumps(build(args.panel, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
