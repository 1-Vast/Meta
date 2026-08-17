"""Stage E frozen evaluation on meta_val (read once, after freezing).

Metrics per episode per condition: MSE (restored pK^2), RMSE, level^2,
centered MSE, Spearman, Pearson, R^2, concordance index, activity-cliff sign,
and ligand-novelty strata (max/mean Tanimoto to meta_train ligands).

Conditions: correct (full), zero_shot (adapt=False), permuted_support
(k>=2), matched_wrong_support (k=1), wrong_protein (full-system perturbation).

No training; meta_test never constructed.
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
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, forward, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.evaluate import (                # noqa: E402
    concordance, pearson, r_squared, spearman,
)
from tools.research.stageD_level_panel.model import PanelLevelShapeModel  # noqa: E402

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
EVALUATION_SEED = 73101
QUERY_SIZE = 16
DRAWS = 2
CLIFF_SIMILARITY = 0.6
CLIFF_GAP = 1.0


def tanimoto_matrix(fingerprint: np.ndarray) -> np.ndarray:
    intersection = fingerprint @ fingerprint.T
    total = fingerprint.sum(1)
    union = total[:, None] + total[None, :] - intersection
    return intersection / np.maximum(union, 1e-9)


def cliff_sign(prediction: np.ndarray, truth: np.ndarray,
               similarity: np.ndarray) -> float:
    rows, cols = np.triu_indices(len(truth), 1)
    keep = (similarity[rows, cols] >= CLIFF_SIMILARITY) &            (np.abs(truth[rows] - truth[cols]) >= CLIFF_GAP)
    if not keep.any():
        return float("nan")
    delta = truth[rows][keep] - truth[cols][keep]
    gap = prediction[rows][keep] - prediction[cols][keep]
    return float((np.sign(gap) == np.sign(delta)).mean())


def episode_metrics(values: np.ndarray, truth: np.ndarray,
                    similarity: np.ndarray) -> dict:
    out = {
        "mse_pk": float(((values - truth) ** 2).mean()),
        "rmse_pk": float(np.sqrt(((values - truth) ** 2).mean())),
        "level_squared": float((values.mean() - truth.mean()) ** 2),
        "centered_mse_pk": float(
            (((values - values.mean()) - (truth - truth.mean())) ** 2).mean()),
        "spearman": float(spearman(values, truth)),
        "pearson": float(pearson(values, truth)),
        "r_squared": float(r_squared(values, truth)),
        "ci": float(concordance(values, truth)),
        "cliff_sign": cliff_sign(values, truth, similarity),
    }
    return out


def build_model(payload: dict, data: QPSMPData):
    from dataclasses import fields
    config = TrainConfig(**{k: v for k, v in payload["config"].items()
                            if k in {f.name for f in fields(TrainConfig)}})
    if payload["arm"] in ("K", "K-REG"):
        from tools.research.stageK_contrastive.model import CoembedModel
        return CoembedModel(
            protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
            hidden_dim=config.hidden_dim, task_dim=config.task_dim,
            ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
            pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
            pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
            support_hidden_dim=config.support_hidden_dim,
            support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
            adaptive_blocks=config.adaptive_blocks,
            adapter_scale=config.adapter_scale,
            coembed_dim=128, dtype=torch.float32)
    if payload["arm"] in ("F", "F-ABS"):
        from tools.research.stageF_pairwise.model import PairwiseTransportModel
        return PairwiseTransportModel(
            protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
            hidden_dim=config.hidden_dim, task_dim=config.task_dim,
            ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
            pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
            pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
            support_hidden_dim=config.support_hidden_dim,
            support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
            adaptive_blocks=config.adaptive_blocks,
            adapter_scale=config.adapter_scale,
            use_learned_key=False, dtype=torch.float32)
    if payload["arm"] in ("LSP", "LSP-NOROUTE"):
        return PanelLevelShapeModel(
            protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
            hidden_dim=config.hidden_dim, task_dim=config.task_dim,
            ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
            pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
            pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
            support_hidden_dim=config.support_hidden_dim,
            support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
            adaptive_blocks=config.adaptive_blocks,
            adapter_scale=config.adapter_scale,
            use_learned_key=False, dtype=torch.float32)
    from scripts.train_qpsmp import resolve_architecture
    return resolve_architecture(config.arch)(
        protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
        hidden_dim=config.hidden_dim, task_dim=config.task_dim,
        ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
        pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
        pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
        support_hidden_dim=config.support_hidden_dim,
        support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
        adaptive_blocks=config.adaptive_blocks, adapter_scale=config.adapter_scale)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protein-bank", type=Path, default=None,
                        help="override the protein bank (external-representation lanes)")
    parser.add_argument("--split", default="meta_val",
                        choices=("meta_train", "meta_val", "meta_test"))
    parser.add_argument("--include-meta-test", action="store_true",
                        help="physically unseal meta_test cells (fail-closed)")
    parser.add_argument("--meta-test-authorization", default=None,
                        help="written reason, mandatory with --include-meta-test")
    parser.add_argument("--live-esm", type=Path, default=None,
                        help="evaluate a Stage I checkpoint with the live LoRA encoder")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    if args.split == "meta_test" and not args.include_meta_test:
        parser.error("evaluating meta_test requires --include-meta-test")
    if args.include_meta_test and not (args.meta_test_authorization or "").strip():
        parser.error("--include-meta-test requires a written --meta-test-authorization")
    protein_bank = PROTEIN_BANK if args.protein_bank is None else args.protein_bank
    data = QPSMPData(CORPUS, protein_bank, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT,
                     include_meta_test=args.include_meta_test,
                     meta_test_authorization=args.meta_test_authorization)
    label_scale = training_label_scale(data)
    payload = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    model = build_model(payload, data)
    model.load_state_dict(payload["model_state"])
    model.to(args.device)
    model.eval()

    live_encoder = None
    if args.live_esm is not None:
        from tools.research.stageI_lm.lora_esm import LiveESMProteinEncoder
        live_encoder = LiveESMProteinEncoder(str(args.live_esm), args.device)
        if payload.get("lora_state"):
            live_encoder.load_lora_state(payload["lora_state"])

    train_ligand_fps = data.fingerprints  # corpus table, used for novelty strata
    train_ligands = [data.cells[i]["ligand_id"]
                     for i in range(len(data.cells))
                     if data.cells[i]["split"] == "meta_train"]
    train_fp_matrix = np.stack(
        [train_ligand_fps[ligand].numpy() for ligand in train_ligands]
    ).astype(np.float32)
    train_fp_sum = train_fp_matrix.sum(1)

    specs = data.fixed_nested_episode_banks(
        args.split, SUPPORT_SIZES, QUERY_SIZE, DRAWS, EVALUATION_SEED, None)
    scale, mean = label_scale.scale, label_scale.mean
    rows = []
    draw_counts = {}
    for k, bank in specs.items():
        for spec in bank:
            draw = draw_counts.get((k, spec.target), 0)
            draw_counts[(k, spec.target)] = draw + 1
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            if live_encoder is not None:
                from dataclasses import replace
                sequence = data._protein_sequences[spec.target]
                with torch.no_grad():
                    pooled, slots, mask = live_encoder.encode(sequence)
                episode = replace(episode, protein_pooled=pooled,
                                  protein_tokens=slots, protein_mask=mask)
            truth = episode.query_y.numpy() * scale + mean
            sim = tanimoto_matrix(episode.query_fingerprint.numpy())
            # ligand novelty: max/mean Tanimoto to any meta_train ligand
            nov = []
            for i in spec.query:
                ligand = data.cells[i]["ligand_id"]
                fp = train_ligand_fps[ligand].numpy()
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
                output = forward(model, episode)
                record("correct", output.prediction,
                       level_abs=float(output.level_adjustment.abs().mean()),
                       sar_abs=float(output.sar_adaptation.abs().mean()))
                record("zero_shot", output.zero_shot)
                if k == 0:
                    continue
                frozen = forward(model, episode, adapt=False)
                record("no_adaptation", frozen.prediction)
                # matched-wrong support (k=1): flip residual magnitude
                raw_residual = (output.support_residual_quotient
                                + output.level_adjustment[..., :1])
                wrong_labels = episode.support_y.numpy() - 2.0 *                     raw_residual.detach().cpu().numpy().ravel()
                from dataclasses import replace
                wrong_ep = replace(episode, support_y=torch.as_tensor(
                    wrong_labels, dtype=episode.support_y.dtype))
                with torch.no_grad():
                    wrong_out = forward(model, wrong_ep)
                record("matched_wrong_support", wrong_out.prediction)
                if k > 1:
                    rolled = replace(episode, support_y=torch.roll(
                        episode.support_y, 1, dims=-1))
                    with torch.no_grad():
                        rolled_out = forward(model, rolled)
                    record("permuted_support", rolled_out.prediction)
                # wrong protein: full-system perturbation
                pooled, tokens, mask = data.protein_for_target(spec.donor_target)
                chemistry = data.protein_chemistry_for_target(spec.donor_target)
                wrong_ep = replace(
                    episode, protein_pooled=pooled, protein_tokens=tokens,
                    protein_mask=mask, protein_chemistry=chemistry)
                with torch.no_grad():
                    wrong_protein = forward(model, wrong_ep)
                record("wrong_protein", wrong_protein.prediction)

    rows.sort(key=lambda r: (r["k"], r["condition"], r["target"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    # aggregate per k and condition with component_target_mean weighting
    from scripts.stageR0_retrieval_falsification import component_target_mean
    aggregates = {}
    for k in SUPPORT_SIZES:
        aggregates[str(k)] = {}
        for condition in sorted({r["condition"] for r in rows if r["k"] == k}):
            metrics = {}
            for field in ("mse_pk", "rmse_pk", "level_squared", "centered_mse_pk",
                          "spearman", "pearson", "r_squared", "ci", "cliff_sign"):
                metrics[field] = component_target_mean(
                    [(r["component"], r["target"], r[field])
                     for r in rows if r["k"] == k and r["condition"] == condition])
            aggregates[str(k)][condition] = metrics
    result = {"schema": "MetaSieve.StageE.Evaluation.v1",
              "checkpoint": str(args.checkpoint.resolve()),
              "arm": payload["arm"],
              "split": args.split,
              "episodes_per_k": {str(k): len([r for r in rows if r["k"] == k])
                                 for k in SUPPORT_SIZES},
              "aggregates": aggregates,
              "meta_test": data.seal_record()}
    result_path = args.output.with_suffix(".summary.json")
    result_path.write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(json.dumps({k: {c: v["mse_pk"] for c, v in aggregates[str(k)].items()}
                      for k in SUPPORT_SIZES}, indent=1))
    print(f"wrote {args.output} and {result_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
