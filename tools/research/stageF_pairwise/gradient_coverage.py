"""Post-hoc gradient-coverage audit for a Stage F checkpoint (one episode)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, forward, normalized_episode, training_label_scale,
)
from tools.research.stageB_complementary.train_stageb import (
    draw_fit_episode, partition_components,
)
from tools.research.stageF_pairwise.model import PairwiseTransportModel

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
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    config = TrainConfig(**{k: v for k, v in payload["config"].items()
                            if k in {f.name for f in TrainConfig.__dataclass_fields__.values()}})
    model = PairwiseTransportModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks, adapter_scale=config.adapter_scale,
        use_learned_key=False, dtype=torch.float32)
    model.load_state_dict(payload["model_state"])
    model.to(args.device).train()
    torch.manual_seed(0)
    rng = np.random.default_rng(0)
    fit, _ = partition_components(data)
    label_scale = training_label_scale(data)
    spec = draw_fit_episode(data, fit, 3, 12, rng)
    episode = compact_episode(normalized_episode(data.materialize(spec), label_scale))
    output = forward(model, episode)
    query_y = episode.query_y.to(device=output.prediction.device,
                                 dtype=output.prediction.dtype)
    loss = F.smooth_l1_loss(output.prediction, query_y)
    model.zero_grad(set_to_none=True)
    loss.backward()
    groups = {}
    for name, param in model.named_parameters():
        if param.grad is None:
            continue
        prefix = name.split(".")[0]
        groups.setdefault(prefix, []).append(float(param.grad.norm()))
    coverage = {
        "schema": "MetaSieve.StageF.GradientCoverage.v1",
        "checkpoint": str(args.checkpoint.resolve()),
        "per_module_mean_grad_norm": {
            key: float(np.mean(values)) for key, values in groups.items()},
        "modules_with_gradient": sorted(groups),
    }
    args.output.write_text(json.dumps(coverage, indent=1), encoding="utf-8")
    print(json.dumps(coverage, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
