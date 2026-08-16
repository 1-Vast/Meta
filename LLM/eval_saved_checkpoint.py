from pathlib import Path
import json
import sys
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from model.qpsmp_meta import QPSMPBioModel
from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    TrainConfig, evaluate, training_label_scale,
)

run = Path("report/meta_fewshot/lcipf_elmt_k3_performance_80step_rerun_20260814")
payload = torch.load(run / "checkpoint.pt", map_location="cpu", weights_only=False)
config = TrainConfig(**payload["config"])
data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
model = QPSMPBioModel(
    protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
    hidden_dim=config.hidden_dim, task_dim=config.task_dim,
    ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
    pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
    pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
    support_hidden_dim=config.support_hidden_dim,
    support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
    adaptive_blocks=config.adaptive_blocks, adapter_scale=config.adapter_scale,
    use_cartesian=config.use_cartesian)
model.load_state_dict(payload["model_state"])
model.to(config.device)
scale = training_label_scale(data)
bank = data.fixed_episode_bank(
    "meta_test", config.support_size, config.query_size,
    config.test_draws_per_target, config.seed,
    config.eval_targets_per_component)
metrics = evaluate(model, data, bank, controls=True, label_scale=scale)
(run / "EVALUATION.json").write_text(
    json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(metrics, indent=2, sort_keys=True))
