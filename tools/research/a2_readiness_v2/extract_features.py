"""Phase 2, step 1: cache frozen A0 internal representations for both splits.

Writes one `.npz` per (split, protein condition). No training, no gradient,
`meta_test` unreachable. Query labels are stored as *targets only* — nothing in
this file lets a label reach an encoder, a donor rule or a normalisation
statistic.

Protein conditions:

``correct``     the recipient's own protein — the deployed condition;
``wrong``       the nearest legal cross-component donor, whitened on meta_train;
``reference``   one fixed meta_train protein for *every* target. This is the
                ligand-only control done properly: the feature space and the
                trunk are identical, only the recipient-specific protein
                information is removed. Comparing `correct` against `reference`
                isolates exactly what knowing the right protein buys.

Run:
```
conda run -n drug python -m tools.research.a2_readiness_v2.extract_features \
  --split meta_train --output tools/research/a2_readiness_v2/features
```
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
from scripts.train_level_shape import normalized                      # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)
from tools.research.a2_readiness_v2 import _frozen                    # noqa: E402
from tools.research.a2_readiness_v2._arms import trained_arm          # noqa: E402
from tools.research.a2_readiness_v2._donors import stratified_donors  # noqa: E402
from tools.research.a2_readiness_v2._features import extract          # noqa: E402

CONDITIONS = ("correct", "wrong", "reference")


def protein_inputs(data, target: str, device: str, dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return [pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0)]


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
    specs = data.fixed_nested_episode_banks(
        arguments.split, (0,), _frozen.QUERY_SIZE, 1,
        _frozen.EVALUATION_SEED, None)[0]

    # Donors for `meta_train` come from `meta_train`, for `meta_val` from
    # `meta_val`; whitening statistics are always `meta_train`. A donor drawn
    # from the other split would confound wrong-identity with seen-versus-unseen.
    donors = stratified_donors(data, arguments.split, arguments.split,
                               _frozen.WHITENING_POOL)
    # One fixed reference protein, chosen deterministically as the first
    # meta_train target in sorted order. It is a *choice made before any result
    # is observed*, not a search.
    reference_target = sorted(data.tasks["meta_train"])[0]

    model, _, seed = trained_arm(arguments.checkpoint, data, arguments.device)
    model.eval()
    dtype = next(model.parameters()).dtype

    store: dict[str, list] = {f"{condition}__{name}": []
                              for condition in CONDITIONS
                              for name in ("occupancy", "mean_state",
                                           "max_state", "embed", "section",
                                           "ligand", "interaction",
                                           "_endpoint", "_protein_value")}
    labels, target_index = [], []
    targets, components = [], []
    protein_pooled, wrong_protein_pooled, reference_protein_pooled = [], [], []

    for order, spec in enumerate(specs):
        episode = compact_episode(normalized(data.materialize(spec), scale))
        query = (episode.query_atoms.to(arguments.device, dtype).unsqueeze(0),
                 episode.query_bonds.to(arguments.device, dtype).unsqueeze(0),
                 episode.query_mask.to(arguments.device, dtype).unsqueeze(0),
                 episode.query_fingerprint.to(arguments.device, dtype).unsqueeze(0))
        sources = {
            "correct": spec.target,
            "wrong": donors[spec.target]["nearest"][0],
            "reference": reference_target,
        }
        for condition, source in sources.items():
            features = extract(
                model, protein_inputs(data, source, arguments.device, dtype),
                *query)
            for name, value in features.items():
                key = f"{condition}__{name}"
                if key in store:
                    store[key].append(value.astype(np.float32))
        count = int(episode.query_y.numel())
        labels.append(episode.query_y.numpy() * scale.scale + scale.mean)
        target_index.append(np.full(count, order, dtype=np.int64))
        targets.append(spec.target)
        components.append(spec.component)
        # The pooled protein vector for each condition, so a wrong-protein arm
        # can replace the protein *everywhere* — both the features and the
        # protein-side conditioning — rather than only half of it.
        for store_list, source in ((protein_pooled, sources["correct"]),
                                   (wrong_protein_pooled, sources["wrong"]),
                                   (reference_protein_pooled, sources["reference"])):
            store_list.append(
                np.asarray(data.protein_for_target(source)[0], dtype=np.float32))
        if (order + 1) % 50 == 0:
            print(f"  {order + 1}/{len(specs)} targets")

    unique_components = sorted(set(components))
    component_of_target = np.asarray(
        [unique_components.index(c) for c in components], dtype=np.int64)

    payload = {key: np.concatenate(value, 0) for key, value in store.items()
               if value}
    payload["y"] = np.concatenate(labels, 0).astype(np.float32)
    payload["target_index"] = np.concatenate(target_index, 0)
    payload["component_of_target"] = component_of_target
    payload["protein_pooled"] = np.stack(protein_pooled)
    payload["wrong_protein_pooled"] = np.stack(wrong_protein_pooled)
    payload["reference_protein_pooled"] = np.stack(reference_protein_pooled)

    arguments.output.mkdir(parents=True, exist_ok=True)
    path = arguments.output / f"{arguments.split}.npz"
    np.savez_compressed(path, **payload)
    meta = {
        "schema": "MetaSieve.A2ReadinessV2.Features.v1",
        "split": arguments.split,
        "checkpoint": str(Path(arguments.checkpoint).relative_to(ROOT)),
        "checkpoint_sha256": _frozen.sha256(Path(arguments.checkpoint)),
        "checkpoint_seed": seed,
        "targets": len(specs), "components": len(unique_components),
        "rows": int(payload["y"].shape[0]),
        "conditions": list(CONDITIONS),
        "reference_protein_target": reference_target,
        "reference_protein_rule": "first meta_train target in sorted order",
        "widths": {name: int(payload[f"correct__{name}"].shape[1])
                   for name in ("occupancy", "mean_state", "max_state",
                                "embed", "section", "ligand", "interaction")},
        "label_scale": {"mean": float(scale.mean), "scale": float(scale.scale),
                        "fitted_on": "meta_train"},
        "meta_test": data.seal_record(),
        "frozen_design": _frozen.frozen_manifest(),
    }
    (arguments.output / f"{arguments.split}.meta.json").write_text(
        json.dumps(meta, indent=1), encoding="utf-8")
    print(f"wrote {path}  ({payload['y'].shape[0]} rows, "
          f"{len(specs)} targets, {len(unique_components)} components)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
