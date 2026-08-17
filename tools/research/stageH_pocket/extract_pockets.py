"""Stage H pocket-prior extractor (external data lane).

For every governed meta_train/meta_val target, find homologous holo complexes
in the local pilot20k corpus (MMseqs2, >=30% identity, >=50% query coverage),
and extract pocket descriptors from the holo structure: protein heavy atoms
within 6.0 A of the holo ligand. Descriptors are structure-derived protein
features reported as external data; no BindingDB labels are touched, and
meta_test targets are never processed in this stage.

Output: tools/research/stageH_pocket/pocket_descriptors.npz + manifest.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import gemmi  # noqa: E402

from scripts.qpsmp_data import QPSMPData  # noqa: E402
from scripts.train_qpsmp import (  # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
BASE = Path(__file__).resolve().parent
HOLO = ROOT / "dataset/processed/open_structures/pilot20k_holo_governed_v2"
RADIUS = 6.0
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(AA)}
ELEMENTS = ("C", "N", "O", "S")


def parse_row(structure_path, protein_asym, ligand_comp, ligand_auth_seq):
    path = Path(structure_path)
    if not path.exists():
        return None
    try:
        structure = gemmi.read_structure(str(path))
    except Exception:  # noqa: BLE001
        return None
    model = structure[0]
    protein = model.find_chain(protein_asym)
    if protein is None:
        protein = max(model, key=len, default=None)
    if protein is None:
        return None
    # The holo ligand may sit inside the protein chain or in its own chain;
    # locate it by compound id (and auth seqid when parseable) across chains.
    ligand_atoms = []
    try:
        want_seq = int(str(ligand_auth_seq).split(".")[0])
    except ValueError:
        want_seq = None
    best = None
    for chain in model:
        for residue in chain:
            if residue.het_flag != "H" or residue.name != ligand_comp:
                continue
            score = 0
            if want_seq is not None and residue.seqid.num == want_seq:
                score += 1
            if best is None or score > best[0]:
                best = (score, residue)
    if best is None:
        return None
    for atom in best[1]:
        if atom.element.name != "H" and atom.pos.x == atom.pos.x:
            ligand_atoms.append(atom.pos)
    if len(ligand_atoms) < 4:
        return None
    protein_atoms = []
    for residue in protein:
        if residue.het_flag == "A":
            continue
        aa = AA_INDEX.get(residue.name, -1)
        for atom in residue:
            if atom.element.name != "H" and atom.pos.x == atom.pos.x:
                protein_atoms.append((atom, aa))
    if not protein_atoms:
        return None
    positions = np.asarray(
        [[atom.pos.x, atom.pos.y, atom.pos.z] for atom, _ in protein_atoms],
        dtype=np.float32)
    ligand_positions = np.asarray(
        [[pos.x, pos.y, pos.z] for pos in ligand_atoms], dtype=np.float32)
    delta = positions[:, None, :] - ligand_positions[None, :, :]
    distance = np.sqrt((delta ** 2).sum(-1))      # [N_protein, N_ligand]
    mask = (distance <= RADIUS).any(-1)
    indices = np.flatnonzero(mask)
    if not len(indices):
        return None
    aa_counts = np.zeros(len(AA), dtype=np.float32)
    elem_counts = {e: 0.0 for e in ELEMENTS}
    coords = []
    for i in indices:
        atom, aa = protein_atoms[int(i)]
        coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
        if aa >= 0:
            aa_counts[aa] += 1.0
        for e in ELEMENTS:
            if atom.element.name == e:
                elem_counts[e] += 1.0
                break
    coords = np.asarray(coords, dtype=np.float32)
    span = coords.max(0) - coords.min(0)
    volume = float(np.prod(span + 1e-3))
    n = len(coords)
    return {
        "pocket_atoms": n,
        "pocket_volume_a3": volume,
        "aa_fraction": (aa_counts / n).tolist(),
        "elem_fraction": [elem_counts[e] / n for e in ELEMENTS],
        "holo_ligand_heavy_atoms": len(ligand_atoms),
    }


def main():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    targets = {}
    for split in ("meta_train", "meta_val"):
        for target in sorted(data.tasks[split]):
            targets[target] = split
    hits = {}
    with (BASE / "target_vs_holo.tsv").open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip().split()
            if len(parts) < 6:
                continue
            query, holo_key, pident, _, qcov, _ = parts[:6]
            pident, qcov = float(pident), float(qcov)
            if query not in targets:
                continue
            if pident < 30.0 or qcov < 0.5:
                continue
            hits.setdefault(query, []).append((pident, qcov, holo_key))
    for query in hits:
        hits[query].sort(key=lambda item: (-item[0], -item[1]))
        hits[query] = hits[query][:5]

    holo_rows = {}
    with (HOLO / "complexes.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if not row.get("biologically_relevant") or not row.get("noncovalent"):
                continue
            if int(row.get("ligand_heavy_atoms", 0)) < 4:
                continue
            holo_rows.setdefault(row["sequence_sha256"], []).append(row)

    keys, vectors, meta = [], [], []
    covered = 0
    for target, split in sorted(targets.items()):
        candidates = hits.get(target, [])
        descriptors = []
        used = 0
        for pident, qcov, holo_key in candidates:
            for row in holo_rows.get(holo_key, [])[:5]:
                desc = parse_row(row["structure_path"], row["protein_asym_id"],
                                row["ligand_comp_id"], row["ligand_auth_seq_id"])
                if desc is None:
                    continue
                descriptors.append(desc)
                used += 1
                if used >= 5:
                    break
            if used >= 5:
                break
        if not descriptors:
            continue
        covered += 1
        stack = np.asarray([
            np.concatenate([
                np.asarray([d["pocket_atoms"], d["pocket_volume_a3"],
                            d["holo_ligand_heavy_atoms"]], dtype=np.float32),
                np.asarray(d["aa_fraction"], dtype=np.float32),
                np.asarray(d["elem_fraction"], dtype=np.float32),
            ]) for d in descriptors])
        vector = np.concatenate([
            stack.mean(0),
            np.asarray([candidates[0][0], candidates[0][1], float(used)],
                       dtype=np.float32),
        ])
        keys.append(target)
        vectors.append(vector)
        meta.append({"target": target, "split": split,
                     "best_pident": candidates[0][0],
                     "best_qcov": candidates[0][1],
                     "n_complexes": used})
    vectors = np.stack(vectors)
    np.savez_compressed(BASE / "pocket_descriptors.npz",
                        keys=np.asarray(keys), vectors=vectors)
    order = ("pocket_atoms, pocket_volume_a3, holo_ligand_heavy_atoms, "
             "20x aa_fraction, 4x elem_fraction(C,N,O,S), best_pident, "
             "best_qcov, n_complexes")
    (BASE / "pocket_descriptors.manifest.json").write_text(json.dumps({
        "schema": "MetaSieve.StageH.PocketDescriptors.v1",
        "source": "pilot20k_holo_governed_v2 (local BioLiP2-derived corpus)",
        "radius_angstrom": RADIUS,
        "min_pident": 30.0, "min_qcov": 0.5,
        "targets_total": len(targets), "targets_covered": covered,
        "vector_dim": int(vectors.shape[1]),
        "descriptor_order": order,
        "meta_test": data.seal_record(),
    }, indent=1), encoding="utf-8")
    print("covered", covered, "of", len(targets), "targets; vector dim",
          vectors.shape[1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
