"""Extract frozen A0 representations for every panel, under both protein conditions.

This is the expensive half of the shared discriminator and it is pure inference:
no gradient, no training, no label read. It writes one cache file per
(checkpoint seed, protein condition) so the analysis can stream them one at a
time instead of holding every representation in memory at once.

Two conditions are extracted:

* `correct` — the panel's own protein;
* `wrong` — the nearest cross-component protein (`_donors.stratified_donors`),
  which is the hardest legal substitution and therefore the one a specificity
  claim has to beat.

The protein-permutation control is deliberately *not* extracted here. It is a
confirmatory control, and spending three more checkpoint sweeps confirming the
specificity of an effect that has not yet cleared its primary contrast would be
compute spent on a conclusion that is not in evidence. It runs only if the gate
passes.

Duplicate representations are not stored: `ContactGrammar` applies a 0/1 atom
mask to `state` and then pools it, so the captured `state_mean` and `state_max`
are bitwise the model's own `mean_state` and `max_state`. Only `state_rms` is
new post-fusion information, and a test pins that identity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)
from tools.research.a2_readiness_v2._donors import stratified_donors   # noqa: E402
from tools.research.method_ladder._shared import _frozen               # noqa: E402
from tools.research.method_ladder._shared.capture import extract       # noqa: E402
from tools.research.method_ladder._shared.panels import (              # noqa: E402
    build_panels, panel_inputs, protein_parts,
)

# Stored representations. `state_mean`/`state_max` are omitted as exact
# duplicates of `mean_state`/`max_state`; see the module docstring.
STORED = ("ligand", "context_mean", "context_max", "context_rms",
          "mean_state", "max_state", "state_rms",
          "occupancy", "embed", "section", "interaction")

CACHE = Path(__file__).resolve().parent / "cache"


def open_data() -> QPSMPData:
    """Fail-closed: `meta_test` is excluded by the default and stays that way."""
    return QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=_frozen.SPLIT_DIRECTORY)


def sweep(model, data, panels, condition: str, device: str) -> dict:
    """One forward pass per panel; returns stacked features and panel metadata."""
    blocks: dict[str, list[np.ndarray]] = {name: [] for name in STORED}
    labels, targets, components, counts = [], [], [], []
    for panel in panels:
        atoms, bonds, mask, fingerprint = panel_inputs(data, panel, device)
        source = panel.target if condition == "correct" else panel.donor
        features = extract(model, protein_parts(data, source, device),
                           atoms, bonds, mask, fingerprint)
        for name in STORED:
            blocks[name].append(features[name].astype(np.float32))
        labels.append(panel.labels.astype(np.float32))
        targets.append(panel.target)
        components.append(panel.component)
        counts.append(len(panel.cells))
    payload = {name: np.concatenate(values, axis=0)
               for name, values in blocks.items()}
    payload["_labels"] = np.concatenate(labels, axis=0)
    payload["_counts"] = np.asarray(counts, dtype=np.int64)
    payload["_targets"] = np.asarray(targets, dtype="U64")
    payload["_components"] = np.asarray(components, dtype="U64")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--splits", nargs="+", default=["meta_train", "meta_val"])
    parser.add_argument("--conditions", nargs="+", default=["correct", "wrong"])
    arguments = parser.parse_args()

    CACHE.mkdir(parents=True, exist_ok=True)
    data = open_data()
    manifest = {"schema": "MethodLadder.FeatureCache.v1",
                "frozen": _frozen.frozen_manifest(),
                "stored": list(STORED), "panels": {}, "files": {},
                "meta_test": data.seal_record()}

    donors = {split: stratified_donors(data, split, "meta_val", "meta_train")
              for split in arguments.splits}
    panels = {split: build_panels(data, split, donors[split])
              for split in arguments.splits}
    for split, items in panels.items():
        manifest["panels"][split] = {
            "targets": len(items),
            "components": len({p.component for p in items}),
            "ligands": int(sum(len(p.cells) for p in items))}
        print(f"{split}: {len(items)} panels, "
              f"{len({p.component for p in items})} components, "
              f"{sum(len(p.cells) for p in items)} ligands")

    for checkpoint in _frozen.A0_CHECKPOINTS:
        if not checkpoint.exists():
            raise FileNotFoundError(f"frozen A0 checkpoint missing: {checkpoint}")
        seed = checkpoint.parent.name.replace("A0_incumbent_seed", "")
        model, _kind, _seed = load_arm(checkpoint, data, arguments.device)
        model.eval()
        for split in arguments.splits:
            for condition in arguments.conditions:
                name = f"{split}.{condition}.seed{seed}.npz"
                path = CACHE / name
                if path.exists():
                    print(f"  {name}: cached")
                    continue
                payload = sweep(model, data, panels[split], condition,
                                arguments.device)
                np.savez_compressed(path, **payload)
                manifest["files"][name] = _frozen.sha256(path)
                print(f"  {name}: {payload['_labels'].shape[0]} rows")
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    (CACHE / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=1, sort_keys=True), encoding="utf-8")
    print(f"\nwrote {CACHE / 'MANIFEST.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
