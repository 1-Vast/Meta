"""Stage L frozen evaluation on meta_val (gate: level head active at k=0 only)."""
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
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.evaluate import (
    concordance, pearson, r_squared, spearman,
)
from tools.research.stageJ_assay.model import AssayLevelModel
from tools.research.stageJ_assay.train_stagej import episode_journal_ids
from tools.research.stageL_gated.train_stagel import forward_gated

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)


def tanimoto_matrix(fingerprint):
    intersection = fingerprint @ fingerprint.T
    total = fingerprint.sum(1)
    union = total[:, None] + total[None, :] - intersection
    return intersection / np.maximum(union, 1e-9)


def cliff_sign(prediction, truth, similarity):
    rows, cols = np.triu_indices(len(truth), 1)
    keep = (similarity[rows, cols] >= 0.6) & (np.abs(truth[rows] - truth[cols]) >= 1.0)
    if not keep.any():
        return float("nan")
    delta = truth[rows][keep] - truth[cols][keep]
    gap = prediction[rows][keep] - prediction[cols][keep]
    return float((np.sign(gap) == np.sign(delta)).mean())


def episode_metrics(values, truth, similarity):
    return {
        "mse_pk": float(((values - truth) ** 2).mean()),
        "rmse_pk": float(np.sqrt(((values - truth) ** 2).mean())),
        "level_squared": float((values.mean() - truth.mean()) ** 2),
        "centered_mse_pk": float((((values - values.mean())
                                  - (truth - truth.mean())) ** 2).mean()),
        "spearman": float(spearman(values, truth)),
        "pearson": float(pearson(values, truth)),
        "r_squared": float(r_squared(values, truth)),
        "ci": float(concordance(values, truth)),
        "cliff_sign": cliff_sign(values, truth, similarity),
    }


def main():
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
    from dataclasses import fields
    config = TrainConfig(**{k: v for k, v in payload["config"].items()
                            if k in {f.name for f in fields(TrainConfig)}})
    vocab = payload["journal_vocab"]
    model = AssayLevelModel(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks,
        adapter_scale=config.adapter_scale,
        journal_vocab=len(vocab), use_learned_key=False, dtype=torch.float32)
    model.load_state_dict(payload["model_state"])
    model.to(args.device)
    model.eval()
    index = {code: i for i, code in enumerate(vocab)}

    train_ligands = [data.cells[i]["ligand_id"] for i in range(len(data.cells))
                     if data.cells[i]["split"] == "meta_train"]
    fp_table = data.fingerprints
    train_fp_matrix = np.stack([fp_table[ligand].numpy()
                                for ligand in train_ligands]).astype(np.float32)
    train_fp_sum = train_fp_matrix.sum(1)

    specs = data.fixed_nested_episode_banks("meta_val", SUPPORT_SIZES, 16, 2,
                                            73101, None)
    scale, mean = label_scale.scale, label_scale.mean
    rows = []
    draw_counts = {}
    for k, bank in specs.items():
        gate = 1.0 if k == 0 else 0.0
        for spec in bank:
            draw = draw_counts.get((k, spec.target), 0)
            draw_counts[(k, spec.target)] = draw + 1
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            ids = episode_journal_ids(data, spec, index)
            truth = episode.query_y.numpy() * scale + mean
            sim = tanimoto_matrix(episode.query_fingerprint.numpy())
            nov = []
            for i in spec.query:
                fp = fp_table[data.cells[i]["ligand_id"]].numpy()
                inter = train_fp_matrix @ fp
                union = train_fp_sum + float(fp.sum()) - inter
                nov.append(float((inter / np.maximum(union, 1e-9)).max()))

            def record(condition, prediction, **extra):
                values = prediction.detach().cpu().numpy().ravel() * scale + mean
                rows.append({"k": k, "component": spec.component,
                             "target": spec.target, "draw": draw,
                             "condition": condition,
                             "max_train_tanimoto": float(max(nov)),
                             "mean_train_tanimoto": float(np.mean(nov)),
                             **episode_metrics(values, truth, sim), **extra})

            with torch.no_grad():
                output = forward_gated(model, episode, ids, gate)
                record("correct", output.prediction)
                record("zero_shot", output.zero_shot)
                if k == 0:
                    continue
                raw_residual = (output.support_residual_quotient
                                + output.level_adjustment[..., :1])
                wrong_labels = episode.support_y.numpy() - 2.0 * raw_residual.detach().cpu().numpy().ravel()
                from dataclasses import replace
                wrong_ep = replace(episode, support_y=torch.as_tensor(
                    wrong_labels, dtype=episode.support_y.dtype))
                with torch.no_grad():
                    wrong_out = forward_gated(model, wrong_ep, ids, gate)
                record("matched_wrong_support", wrong_out.prediction)
                if k > 1:
                    rolled = replace(episode, support_y=torch.roll(
                        episode.support_y, 1, dims=-1))
                    with torch.no_grad():
                        rolled_out = forward_gated(model, rolled, ids, gate)
                    record("permuted_support", rolled_out.prediction)
                pooled, tokens, mask = data.protein_for_target(spec.donor_target)
                chemistry = data.protein_chemistry_for_target(spec.donor_target)
                wrong_ep = replace(episode, protein_pooled=pooled,
                                  protein_tokens=tokens, protein_mask=mask,
                                  protein_chemistry=chemistry)
                with torch.no_grad():
                    wrong_protein = forward_gated(model, wrong_ep, ids, gate)
                record("wrong_protein", wrong_protein.prediction)

    rows.sort(key=lambda r: (r["k"], r["condition"], r["target"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    from scripts.stageR0_retrieval_falsification import component_target_mean
    aggregates = {}
    for k in SUPPORT_SIZES:
        aggregates[str(k)] = {}
        for condition in sorted({r["condition"] for r in rows if r["k"] == k}):
            metrics = {}
            for field in ("mse_pk", "rmse_pk", "level_squared", "centered_mse_pk",
                          "spearman", "pearson", "r_squared", "ci", "cliff_sign"):
                metrics[field] = component_target_mean(
                    [(r["component"], r["target"], r[field]) for r in rows
                     if r["k"] == k and r["condition"] == condition])
            aggregates[str(k)][condition] = metrics
    result = {"schema": "MetaSieve.StageL.Evaluation.v1",
              "checkpoint": str(args.checkpoint.resolve()),
              "arm": payload["arm"], "split": "meta_val",
              "aggregates": aggregates, "meta_test": data.seal_record()}
    result_path = args.output.with_suffix(".summary.json")
    result_path.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(json.dumps({k: {c: v["mse_pk"] for c, v in aggregates[str(k)].items()}
                      for k in SUPPORT_SIZES}, indent=1))
    print(f"wrote {args.output} and {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
