"""Stage E level-head audit: panel-covariate shuffle.

For the frozen LSP checkpoint, per meta_val episode we substitute the query
ligand set with the query ligands of a DIFFERENT episode (same split, cyclic
offset), keep the protein and support unchanged, and measure the resulting
level-head output and zero-shot MSE. If the level head uses real panel
information, the shuffled panel must hurt the level alignment. This is an
audit of the trained mechanism, not a new training run.
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

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageD_level_panel.evaluate_staged import build_model

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = build_model(payload, data)
    model.load_state_dict(payload["model_state"])
    model.to(args.device).eval()

    specs = data.fixed_nested_episode_banks(
        "meta_val", (0,), 16, 2, 73101, None)[0]
    episodes = [compact_episode(normalized_episode(data.materialize(s), label_scale))
                for s in specs]
    rows = []
    with torch.no_grad():
        for index, episode in enumerate(episodes):
            from dataclasses import replace
            donor = episodes[(index + 1) % len(episodes)]
            shuffled = replace(
                episode,
                query_atoms=donor.query_atoms, query_bonds=donor.query_bonds,
                query_mask=donor.query_mask, query_y=donor.query_y,
                query_fingerprint=donor.query_fingerprint)
            from scripts.train_qpsmp import forward
            truth = episode.query_y.numpy() * label_scale.scale + label_scale.mean
            own = forward(model, episode)
            shuf = forward(model, shuffled)
            own_level = (own.zero_shot.detach().cpu().numpy() * label_scale.scale
                         + label_scale.mean).mean()
            shuf_level = (shuf.zero_shot.detach().cpu().numpy() * label_scale.scale
                          + label_scale.mean).mean()
            own_mse = float(((own.prediction.detach().cpu().numpy()
                              * label_scale.scale + label_scale.mean
                              - truth) ** 2).mean())
            rows.append({"component": episode.spec.component,
                         "target": episode.spec.target,
                         "truth_mean": float(truth.mean()),
                         "own_level": float(own_level),
                         "shuffled_level": float(shuf_level),
                         "own_k0_mse": own_mse,
                         "level_move": float(abs(own_level - shuf_level))})
    payload_out = {
        "schema": "MetaSieve.StageE.PanelShuffleAudit.v1",
        "checkpoint": str(args.checkpoint.resolve()),
        "n_episodes": len(rows),
        "mean_level_move_pk": float(np.mean([r["level_move"] for r in rows])),
        "own_level_error_pk2": float(np.mean([
            (r["own_level"] - r["truth_mean"]) ** 2 for r in rows])),
        "shuffled_level_error_pk2": float(np.mean([
            (r["shuffled_level"] - r["truth_mean"]) ** 2 for r in rows])),
        "meta_test": data.seal_record(),
    }
    args.output.write_text(json.dumps(payload_out, indent=1), encoding="utf-8")
    print(json.dumps(payload_out, indent=1))
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
