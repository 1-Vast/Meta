"""Stage R step 1: frozen A0 features and zero-shot endpoints, per (target, ligand).

Features depend only on `(protein, ligand)`, never on the episode, so they are
extracted once per pair and episodes are formed afterwards as index subsets.
That makes the number of training episodes free rather than a compute budget,
which matters because the exact A2 operator needs *episodes* — a support set
and a query panel — not the ligand pairs the v2 probe used.

Stored per (target, ligand), for each protein condition:

* `embed` (96) — the A2 plan's `e0`;
* `max_state` (192) — the representation Phase 3 measured as retaining the most
  protein-differential of any stage;
* `ligand` (192) — the protein-blind encoder output, the negative reference;
* `f0` — the frozen zero-shot endpoint in pK, which defines the residual
  `r_i = y_i − f0(P, L_i)` the operator consumes;
* the Morgan fingerprint, for the Tanimoto comparator.

Protein conditions: `correct`, `wrong` (nearest legal cross-component donor
from the same split, whitened on meta_train), `reference` (one fixed meta_train
protein for every target).

No training. `meta_test` unreachable. Query labels are stored as loss targets
and enter no encoder, selector, donor rule or normalisation statistic.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    training_label_scale,
)
from tools.research.a2_readiness_v2 import _frozen                    # noqa: E402
from tools.research.a2_readiness_v2._arms import trained_arm          # noqa: E402
from tools.research.a2_readiness_v2._donors import stratified_donors  # noqa: E402
from tools.research.a2_readiness_v2._features import Capture          # noqa: E402

CONDITIONS = ("correct", "wrong", "reference")
KEEP = ("embed", "max_state", "ligand")
MAX_LIGANDS_PER_TARGET = 24          # k=5 support + 16 query + headroom


def protein_inputs(data, target: str, device: str, dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return [pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0)]


def ligand_tensors(data, indices, device, dtype):
    """Pad one target's ligand graphs to a common width and batch them."""
    graphs = []
    for index in indices:
        ligand = data.cells[int(index)]["ligand_id"]
        atoms, bonds, mask = data.ligand_bank.get(ligand)
        graphs.append((torch.from_numpy(atoms.copy()),
                       torch.from_numpy(bonds.copy()),
                       torch.from_numpy(mask.copy())))
    width = max(item[0].shape[0] for item in graphs)
    atoms, bonds, masks = [], [], []
    for atom, bond, mask in graphs:
        missing = width - atom.shape[0]
        atoms.append(torch.nn.functional.pad(atom, (0, 0, 0, missing)))
        bonds.append(torch.nn.functional.pad(bond, (0, 0, 0, missing, 0, missing)))
        masks.append(torch.nn.functional.pad(mask, (0, missing)))
    return (torch.stack(atoms).to(device, dtype).unsqueeze(0),
            torch.stack(bonds).to(device, dtype).unsqueeze(0),
            torch.stack(masks).to(device, dtype).unsqueeze(0))


def encode(model, protein_parts, atoms, bonds, mask, fingerprint, scale):
    """Frozen representations plus the zero-shot endpoint, one row per ligand.

    `f0` is returned **in pK**. The trunk predicts in normalised label units
    (`LabelScale`, fitted on `meta_train`), and the residual the operator
    consumes is `y − f0` with `y` in pK, so the conversion has to happen here.
    Getting it wrong is silent: the operator still trains, it just regresses
    onto a residual inflated by the label scale.
    """
    pooled, tokens, protein_mask, chemistry = protein_parts
    capture = Capture(model)
    try:
        with torch.no_grad():
            output = model(
                pooled, tokens, protein_mask,
                atoms[:, :0], bonds[:, :0], mask[:, :0],
                torch.zeros(1, 0, device=pooled.device, dtype=pooled.dtype),
                atoms, bonds, mask, adapt=False,
                protein_chemistry=chemistry,
                support_fingerprint=fingerprint[:, :0],
                query_fingerprint=fingerprint)
    finally:
        values = dict(capture.values)
        capture.close()
    features = {name: values[name].float().cpu().numpy() for name in KEEP}
    raw = output.zero_shot.squeeze(0).detach().float().cpu().numpy()
    features["f0"] = raw * scale.scale + scale.mean
    return features


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", required=True,
                        choices=("meta_train", "meta_val"))
    parser.add_argument("--checkpoint", type=Path,
                        default=_frozen.A0_CHECKPOINTS[0])
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=_frozen.SPLIT_DIRECTORY)
    scale = training_label_scale(data)
    donors = stratified_donors(data, arguments.split, arguments.split,
                               _frozen.WHITENING_POOL)
    reference_target = sorted(data.tasks["meta_train"])[0]

    model, _, seed = trained_arm(arguments.checkpoint, data, arguments.device)
    model.eval()
    dtype = next(model.parameters()).dtype

    # One deterministic ligand order per target, label-blind: dedupe by ligand
    # identity in the order the corpus lists them, then cap. Episodes are drawn
    # from this order later, so the cap never depends on any label.
    targets, components, rows = [], [], []
    store = {f"{condition}__{name}": []
             for condition in CONDITIONS for name in (*KEEP, "f0")}
    labels, fingerprints, owner, ligand_ids = [], [], [], []

    for scanned, target in enumerate(sorted(data.tasks[arguments.split])):
        # `order` counts *kept* targets, not scanned ones: targets with too few
        # unique ligands are skipped, and using the scan index would leave
        # `target_index` pointing past the end of `component_of_target`.
        order = len(targets)
        indices, seen = [], set()
        for index in data.tasks[arguments.split][target]:
            ligand = data.cells[int(index)]["ligand_id"]
            if ligand in seen:
                continue
            seen.add(ligand)
            indices.append(int(index))
            if len(indices) >= MAX_LIGANDS_PER_TARGET:
                break
        if len(indices) < 6:                 # k=5 support plus one query
            continue
        atoms, bonds, mask = ligand_tensors(data, indices, arguments.device, dtype)
        fingerprint = data.fingerprint_rows(tuple(indices))
        sources = {"correct": target,
                   "wrong": donors[target]["nearest"][0],
                   "reference": reference_target}
        for condition, source in sources.items():
            features = encode(
                model, protein_inputs(data, source, arguments.device, dtype),
                atoms, bonds, mask,
                fingerprint.to(arguments.device, dtype).unsqueeze(0), scale)
            for name, value in features.items():
                store[f"{condition}__{name}"].append(value.astype(np.float32))
        labels.append(np.asarray([data.cells[i]["pK"] for i in indices],
                                 dtype=np.float32))
        fingerprints.append(fingerprint.numpy().astype(np.float32))
        owner.append(np.full(len(indices), order, dtype=np.int64))
        # Stable biological identities. Episode seeding keys on the *target
        # name*, not on the positional index, so a change to the ligand cap or
        # the eligibility filter cannot silently renumber the banks; and the
        # reproducibility test compares ligand identities rather than row
        # offsets.
        ligand_ids.extend(data.cells[i]["ligand_id"] for i in indices)
        targets.append(target)
        components.append(data.cells[indices[0]]["protein_group_40"])
        rows.append(len(indices))
        if (scanned + 1) % 50 == 0:
            print(f"  scanned {scanned + 1}, kept {len(targets)}")

    unique_components = sorted(set(components))
    payload = {key: np.concatenate(value, 0) for key, value in store.items()}
    payload["y"] = np.concatenate(labels, 0)
    payload["fingerprint"] = np.concatenate(fingerprints, 0)
    payload["target_index"] = np.concatenate(owner, 0)
    payload["component_of_target"] = np.asarray(
        [unique_components.index(c) for c in components], dtype=np.int64)
    payload["ligands_per_target"] = np.asarray(rows, dtype=np.int64)
    # Fixed-width unicode, not object dtype: these must load under
    # `allow_pickle=False` like every other governed bank in the project.
    payload["target_names"] = np.asarray(targets)
    payload["ligand_ids"] = np.asarray(ligand_ids)
    payload["protein_pooled"] = np.stack([
        np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
        for t in targets])

    arguments.output.mkdir(parents=True, exist_ok=True)
    path = arguments.output / f"{arguments.split}.npz"
    np.savez_compressed(path, **payload)
    (arguments.output / f"{arguments.split}.meta.json").write_text(json.dumps({
        "schema": "MetaSieve.A2ExactProbe.LigandFeatures.v1",
        "split": arguments.split,
        "checkpoint": str(Path(arguments.checkpoint).relative_to(ROOT)),
        "checkpoint_sha256": _frozen.sha256(Path(arguments.checkpoint)),
        "checkpoint_seed": seed,
        "targets": len(targets), "components": len(unique_components),
        "ligand_rows": int(payload["y"].shape[0]),
        "max_ligands_per_target": MAX_LIGANDS_PER_TARGET,
        "min_ligands_required": 6,
        "conditions": list(CONDITIONS),
        "representations": list(KEEP),
        "reference_protein_target": reference_target,
        "label_scale": {"mean": float(scale.mean), "scale": float(scale.scale),
                        "fitted_on": "meta_train"},
        "labels_are_pk": True,
        "meta_test": data.seal_record(),
        "frozen_design": _frozen.frozen_manifest(),
    }, indent=1), encoding="utf-8")
    print(f"wrote {path}  ({len(targets)} targets, "
          f"{payload['y'].shape[0]} ligand rows, "
          f"{len(unique_components)} components)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
