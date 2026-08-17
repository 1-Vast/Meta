"""Governed multi-seed, nested-k evaluation for the QPSMP meta-learner.

The meta-test manifest is constructed without reading affinity values and is
reused unchanged by every trained seed.  This runner is confirmatory plumbing;
it never authorizes a theory gate.
"""
from __future__ import annotations

import argparse
import hashlib
from dataclasses import asdict, dataclass, replace
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import EpisodeSpec, QPSMPData, stable_seed
from scripts.train_qpsmp import (
    ARCHITECTURES, COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    LabelScale, TrainConfig, complete_foreign_prediction,
    counterfactual_label_assignments, forward, normalized_episode, train,
    wrong_protein_prediction,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "report/meta_fewshot/qpsmp_nested"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
MODEL_SEEDS = (20260851, 20260852, 20260853)


@dataclass(frozen=True)
class NestedEpisodeRecord:
    split: str
    component: str
    target: str
    draw: int
    support_cell_ids: tuple[str, ...]
    query_cell_ids: tuple[str, ...]
    donor_target: str


def build_nested_manifest(
        data: QPSMPData, split: str, support_sizes: tuple[int, ...],
        query_size: int, draws: int, seed: int) -> tuple[NestedEpisodeRecord, ...]:
    if not support_sizes or tuple(sorted(set(support_sizes))) != support_sizes:
        raise ValueError("support sizes must be unique and increasing")
    if support_sizes[0] < 0 or query_size < 1 or draws < 1:
        raise ValueError("support size cannot be negative and query/draws must be positive")
    max_k = support_sizes[-1]
    def ligand_unique_order(indices: np.ndarray,
                            rng: np.random.Generator) -> np.ndarray:
        seen, ordered = set(), []
        for index in rng.permutation(indices):
            ligand = data.cells[int(index)]["ligand_id"]
            if ligand not in seen:
                seen.add(ligand)
                ordered.append(int(index))
        return np.asarray(ordered, dtype=np.int64)

    eligible = {
        target: indices for target, indices in data.tasks[split].items()
        if len({data.cells[int(index)]["ligand_id"] for index in indices}) >= max_k + 1
    }
    by_component = {
        component: tuple(target for target in targets if target in eligible)
        for component, targets in data.components[split].items()
    }
    by_component = {key: value for key, value in by_component.items() if value}
    if len(by_component) < 2:
        raise ValueError("nested controls require two eligible components")
    component_order = sorted(by_component)
    records = []
    for component_index, component in enumerate(component_order):
        donor_targets = by_component[
            component_order[(component_index + 1) % len(component_order)]]
        for target in by_component[component]:
            for draw in range(draws):
                rng = np.random.default_rng(
                    stable_seed("qpsmp-nested", seed, split, target, draw))
                order = ligand_unique_order(eligible[target], rng)
                support = order[:max_k]
                query = order[max_k:max_k + min(query_size, len(order) - max_k)]
                donor = donor_targets[int(rng.integers(len(donor_targets)))]
                records.append(NestedEpisodeRecord(
                    split, component, target, draw,
                    tuple(data.cells[int(i)]["cell_id"] for i in support),
                    tuple(data.cells[int(i)]["cell_id"] for i in query), donor))
    return tuple(records)


def manifest_payload(records: tuple[NestedEpisodeRecord, ...], *, seed: int,
                     support_sizes: tuple[int, ...], query_size: int,
                     corpus_manifest_sha256: str | None = None) -> dict:
    return {
        "schema": "MetaSieve.QPSMPNestedEpisodeManifest.v1",
        "selection_uses_affinity": False,
        "seed": seed,
        "support_sizes": list(support_sizes),
        "query_size": query_size,
        "corpus_manifest_sha256": corpus_manifest_sha256,
        "records": [asdict(record) for record in records],
    }


def load_manifest(path: Path) -> tuple[dict, tuple[NestedEpisodeRecord, ...]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "MetaSieve.QPSMPNestedEpisodeManifest.v1":
        raise ValueError("unsupported nested episode manifest")
    if payload.get("selection_uses_affinity") is not False:
        raise ValueError("nested episode selection must be outcome-redacted")
    records = tuple(NestedEpisodeRecord(
        **{**row, "support_cell_ids": tuple(row["support_cell_ids"]),
           "query_cell_ids": tuple(row["query_cell_ids"])})
        for row in payload["records"])
    return payload, records


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_manifest(data: QPSMPData, payload: dict,
                      records: tuple[NestedEpisodeRecord, ...]) -> None:
    expected_hash = file_sha256(data.corpus / "manifest.json")
    if payload.get("corpus_manifest_sha256") != expected_hash:
        raise ValueError("nested manifest corpus hash mismatch")
    cell_by_id = {cell["cell_id"]: cell for cell in data.cells}
    for record in records:
        support, query = set(record.support_cell_ids), set(record.query_cell_ids)
        if not support or not query or support & query:
            raise ValueError("nested manifest has invalid support/query sets")
        if not support.union(query).issubset(cell_by_id):
            raise ValueError("nested manifest references unknown cells")
        if (record.split not in data.tasks
                or record.target not in data.tasks[record.split]):
            raise ValueError("nested manifest declares an unknown split/target")
        rows = [cell_by_id[cell] for cell in support.union(query)]
        if any(row["split"] != record.split
               or row["target_id"] != record.target
               or row["protein_group_40"] != record.component for row in rows):
            raise ValueError(
                "nested manifest cell split/target/component metadata "
                "disagrees with the record")
        split_components = data.components[record.split]
        declared_targets = split_components.get(record.component, ())
        if record.target not in declared_targets:
            raise ValueError("nested manifest component disagrees with the target")
        donor_component = next(
            (component for component, targets in split_components.items()
             if record.donor_target in targets), None)
        if (record.donor_target not in data.tasks[record.split]
                or record.donor_target == record.target
                or donor_component is None
                or donor_component == record.component):
            raise ValueError("nested manifest donor must be a cross-component target")
        support_ligands = {cell_by_id[cell]["ligand_id"] for cell in support}
        query_ligands = {cell_by_id[cell]["ligand_id"] for cell in query}
        if support_ligands & query_ligands:
            raise ValueError("nested manifest overlaps support/query ligand identities")


def episode_spec(data: QPSMPData, record: NestedEpisodeRecord, k: int) -> EpisodeSpec:
    index = {cell["cell_id"]: i for i, cell in enumerate(data.cells)}
    if k < 0 or k > len(record.support_cell_ids):
        raise ValueError("k is outside the nested support manifest")
    support = tuple(index[cell] for cell in record.support_cell_ids[:k])
    query = tuple(index[cell] for cell in record.query_cell_ids)
    return EpisodeSpec(record.split, record.component, record.target,
                       support, query, record.donor_target)


def evaluate_seed(model, data: QPSMPData, records: tuple[NestedEpisodeRecord, ...],
                  support_sizes: tuple[int, ...], label_scale: LabelScale,
                  model_seed: int) -> list[dict]:
    rows = []
    model.eval()
    with torch.no_grad():
        for record in records:
            for k in support_sizes:
                spec = episode_spec(data, record, k)
                episode = normalized_episode(data.materialize(spec), label_scale)
                full = forward(model, episode)
                if k:
                    # `roll(1)` is the identity for k=1, which would make the
                    # label control vacuous exactly where it matters most.
                    # `counterfactual_label_assignments` keeps `roll(1)` for
                    # k>=2 and substitutes the magnitude-matched flip at k=1.
                    wrong_labels = counterfactual_label_assignments(
                        full, episode.support_y)[0]
                    permuted_prediction = forward(model, replace(
                        episode, support_y=wrong_labels)).prediction
                else:
                    permuted_prediction = full.prediction
                predictions = {
                    "full": full.prediction,
                    "sar_cut": full.prediction - full.sar_adaptation,
                    "zero_shot": forward(model, episode, adapt=False).prediction,
                    "level_only": full.level_baseline,
                    "permuted_state": permuted_prediction,
                    "foreign_code_state": (complete_foreign_prediction(
                        model, data, episode, record.donor_target, label_scale)
                        if k else full.prediction),
                    "wrong_protein_state": (wrong_protein_prediction(
                        model, data, episode, record.donor_target)
                        if k else full.prediction),
                }
                for arm, prediction in predictions.items():
                    truth = episode.query_y.detach().cpu().numpy()
                    predicted = prediction.detach().cpu().numpy()
                    mse = float(label_scale.squared_error_pk(
                        prediction, episode.query_y).mean())
                    ci, comparable = concordance_index(predicted, truth)
                    rows.append({
                        "model_seed": model_seed, "component": record.component,
                        "target": record.target, "draw": record.draw, "k": k,
                        "arm": arm, "mse_pk": mse, "rmse_pk": mse ** 0.5,
                        "ci": ci, "ci_comparable_pairs": comparable,
                        "spearman": spearman(predicted, truth),
                    })
    return rows


def concordance_index(prediction: np.ndarray,
                      truth: np.ndarray) -> tuple[float, int]:
    concordant, comparable = 0.0, 0
    for left in range(len(truth)):
        for right in range(left + 1, len(truth)):
            true_delta = float(truth[left] - truth[right])
            if true_delta == 0:
                continue
            pred_delta = float(prediction[left] - prediction[right])
            comparable += 1
            concordant += 0.5 if pred_delta == 0 else float(pred_delta * true_delta > 0)
    return ((float(concordant / comparable) if comparable else float("nan")), comparable)


def _average_ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1)
        start = stop
    return ranks


def spearman(prediction: np.ndarray, truth: np.ndarray) -> float:
    if len(truth) < 2:
        return float("nan")
    left, right = _average_ranks(prediction), _average_ranks(truth)
    left, right = left - left.mean(), right - right.mean()
    denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
    return float(left @ right / denominator) if denominator > 0 else float("nan")


def paired_component_effects(rows: list[dict], correct: str, control: str,
                             k: int) -> dict[str, float]:
    paired: dict[tuple[int, str, str, int], dict[str, float]] = {}
    for row in rows:
        if row["k"] == k and row["arm"] in {correct, control}:
            key = (row["model_seed"], row["component"], row["target"], row["draw"])
            paired.setdefault(key, {})[row["arm"]] = row["mse_pk"]
    component_values: dict[str, list[float]] = {}
    for (_, component, _, _), values in paired.items():
        if set(values) != {correct, control}:
            raise ValueError("paired contrast has a missing arm")
        component_values.setdefault(component, []).append(values[control] - values[correct])
    return {component: float(np.mean(values))
            for component, values in sorted(component_values.items())}


def paired_component_bootstrap(
        rows: list[dict], correct: str, control: str, k: int,
        draws: int, seed: int) -> dict:
    effects = paired_component_effects(rows, correct, control, k)
    components = sorted(effects)
    if not components or draws < 1:
        raise ValueError("bootstrap needs components and positive draws")
    values = np.asarray([effects[component] for component in components])
    rng = np.random.default_rng(seed)
    samples = np.asarray([
        values[rng.integers(len(values), size=len(values))].mean()
        for _ in range(draws)])
    return {
        "correct": correct, "control": control, "k": k,
        "components": len(components), "mean_mse_reduction_pk": float(values.mean()),
        "paired_component_95_ci": [float(np.quantile(samples, 0.025)),
                                   float(np.quantile(samples, 0.975))],
    }


def metric_summary(rows: list[dict]) -> list[dict]:
    summaries = []
    for k in sorted({row["k"] for row in rows}):
        for arm in sorted({row["arm"] for row in rows}):
            selected = [row for row in rows if row["k"] == k and row["arm"] == arm]
            entry = {"k": k, "arm": arm}
            for metric in ("mse_pk", "rmse_pk", "ci", "spearman"):
                finite = [row for row in selected if np.isfinite(row[metric])]
                entry[metric] = (component_target_metric(finite, metric)
                                 if finite else None)
                entry[f"{metric}_episodes"] = len(finite)
            entry["ci_comparable_pairs"] = int(sum(
                row["ci_comparable_pairs"] for row in selected))
            summaries.append(entry)
    return summaries


def component_target_metric(rows: list[dict], metric: str) -> float:
    target_values: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        target_values.setdefault((row["component"], row["target"]), []).append(
            float(row[metric]))
    component_values: dict[str, list[float]] = {}
    for (component, _), values in target_values.items():
        component_values.setdefault(component, []).append(float(np.mean(values)))
    return float(np.mean([np.mean(values) for values in component_values.values()]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--manifest-seed", type=int, default=20260850)
    parser.add_argument("--draws-per-target", type=int, default=1)
    parser.add_argument("--query-size", type=int, default=8)
    parser.add_argument("--split", default="meta_val",
                        choices=("meta_train", "meta_val", "meta_test"),
                        help="population to evaluate; meta_test is sealed and "
                             "requires --include-meta-test (contract 2026-08-16)")
    parser.add_argument("--split-directory", type=Path, default=None,
                        help="governed split directory (double-cold protocol)")
    parser.add_argument("--include-meta-test", action="store_true",
                        help="physically unseal meta_test cells in QPSMPData; "
                             "required to evaluate the sealed split")
    parser.add_argument("--meta-test-authorization", default=None,
                        help="written reason recorded in the artifact; "
                             "mandatory with --include-meta-test "
                             "(contract 2026-08-16)")
    parser.add_argument("--steps", type=int, default=TrainConfig.steps)
    parser.add_argument("--episodes-per-step", type=int,
                        default=TrainConfig.episodes_per_step)
    parser.add_argument("--val-interval", type=int, default=TrainConfig.val_interval)
    parser.add_argument("--val-draws-per-target", type=int,
                        default=TrainConfig.val_draws_per_target)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--hidden-dim", type=int, default=TrainConfig.hidden_dim)
    parser.add_argument("--task-dim", type=int, default=TrainConfig.task_dim)
    parser.add_argument("--adapter-rank", type=int, default=TrainConfig.adapter_rank)
    parser.add_argument("--adaptive-blocks", type=int,
                        default=TrainConfig.adaptive_blocks)
    parser.add_argument("--adapter-scale", type=float,
                        default=TrainConfig.adapter_scale)
    parser.add_argument("--ligand-layers", type=int, default=TrainConfig.ligand_layers)
    parser.add_argument("--pair-dim", type=int, default=TrainConfig.pair_dim)
    parser.add_argument("--pair-blocks", type=int, default=TrainConfig.pair_blocks)
    parser.add_argument("--pair-latents", type=int, default=TrainConfig.pair_latents)
    parser.add_argument("--pair-heads", type=int, default=TrainConfig.pair_heads)
    parser.add_argument("--support-hidden-dim", type=int,
                        default=TrainConfig.support_hidden_dim)
    parser.add_argument("--support-blocks", type=int,
                        default=TrainConfig.support_blocks)
    parser.add_argument("--use-cartesian", action="store_true",
                        help="enable optional common-frame Cartesian inputs; requires geometry-bearing episodes")
    parser.add_argument("--arch", default=TrainConfig.arch,
                        choices=sorted(ARCHITECTURES))
    parser.add_argument("--model-seeds", type=int, nargs="+",
                        default=list(MODEL_SEEDS))
    parser.add_argument("--min-query-size", type=int,
                        default=TrainConfig.min_query_size)
    parser.add_argument("--train-query-size", type=int,
                        default=TrainConfig.query_size)
    parser.add_argument("--learning-rate", type=float,
                        default=TrainConfig.learning_rate)
    parser.add_argument("--backbone-lr-scale", type=float,
                        default=TrainConfig.backbone_lr_scale)
    parser.add_argument("--binding-loss-weight", type=float,
                        default=TrainConfig.binding_loss_weight)
    parser.add_argument("--eval-targets-per-component", type=int, default=999999)
    parser.add_argument("--lr-schedule", default=TrainConfig.lr_schedule,
                        choices=("constant", "cosine"))
    parser.add_argument("--lr-warmup-fraction", type=float,
                        default=TrainConfig.lr_warmup_fraction)
    parser.add_argument("--lr-final-fraction", type=float,
                        default=TrainConfig.lr_final_fraction)
    parser.add_argument("--device", default=TrainConfig.device)
    args = parser.parse_args()
    model_seeds = tuple(args.model_seeds)
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    if args.split == "meta_test" and not args.include_meta_test:
        parser.error("evaluating the sealed meta_test split requires "
                     "--include-meta-test")
    if args.include_meta_test and not (args.meta_test_authorization or "").strip():
        parser.error("--include-meta-test requires a written "
                     "--meta-test-authorization (contract 2026-08-16)")
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=args.include_meta_test,
                     meta_test_authorization=args.meta_test_authorization)
    if args.manifest:
        payload, records = load_manifest(args.manifest)
        if tuple(payload["support_sizes"]) != SUPPORT_SIZES:
            raise ValueError("manifest support sizes differ from the protocol")
        validate_manifest(data, payload, records)
    else:
        records = build_nested_manifest(
            data, args.split, SUPPORT_SIZES, args.query_size,
            args.draws_per_target, args.manifest_seed)
        payload = manifest_payload(records, seed=args.manifest_seed,
                                   support_sizes=SUPPORT_SIZES,
                                   query_size=args.query_size,
                                   corpus_manifest_sha256=file_sha256(
                                       data.corpus / "manifest.json"))
    all_rows, training = [], []
    for model_seed in model_seeds:
        config = TrainConfig(
            seed=model_seed, steps=args.steps,
            episodes_per_step=args.episodes_per_step,
            val_interval=args.val_interval,
            val_draws_per_target=args.val_draws_per_target,
            eval_targets_per_component=args.eval_targets_per_component,
            query_size=args.train_query_size,
            min_query_size=args.min_query_size,
            learning_rate=args.learning_rate,
            backbone_lr_scale=args.backbone_lr_scale,
            binding_loss_weight=args.binding_loss_weight,
            arch=args.arch, lr_schedule=args.lr_schedule,
            lr_warmup_fraction=args.lr_warmup_fraction,
            lr_final_fraction=args.lr_final_fraction,
            hidden_dim=args.hidden_dim, task_dim=args.task_dim,
            adapter_rank=args.adapter_rank,
            adaptive_blocks=args.adaptive_blocks,
            adapter_scale=args.adapter_scale,
            ligand_layers=args.ligand_layers, pair_dim=args.pair_dim,
            pair_blocks=args.pair_blocks, pair_latents=args.pair_latents,
            pair_heads=args.pair_heads,
            support_hidden_dim=args.support_hidden_dim,
            support_blocks=args.support_blocks,
            use_cartesian=args.use_cartesian,
            device=args.device,
        )
        model, diagnostics, label_scale = train(
            data, config, support_sizes=SUPPORT_SIZES)
        args.output.mkdir(parents=True, exist_ok=True)
        checkpoint = args.output / f"checkpoint_seed{model_seed}.pt"
        torch.save({"model_state": model.state_dict(), "config": asdict(config)},
                   checkpoint)
        training.append({
            "model_seed": model_seed, "checkpoint": str(checkpoint),
            "checkpoint_sha256": file_sha256(checkpoint),
            "trainable_parameters": int(sum(
                p.numel() for p in model.parameters() if p.requires_grad)),
            **diagnostics})
        all_rows.extend(evaluate_seed(
            model, data, records, SUPPORT_SIZES, label_scale, model_seed))
        # Seeds run sequentially in one process. Without an explicit release the
        # allocator keeps every earlier seed's cached blocks, and later seeds
        # train several times slower under the resulting memory pressure.
        del model
        if config.device.startswith("cuda"):
            torch.cuda.empty_cache()
            torch.cuda.reset_peak_memory_stats()
    contrasts = [
        paired_component_bootstrap(
            all_rows, "full", control, k, args.bootstrap_draws,
            stable_seed("bootstrap", args.manifest_seed, k, control))
        for k in SUPPORT_SIZES
        for control in ("level_only", "sar_cut", "permuted_state",
                        "foreign_code_state", "wrong_protein_state")
    ]
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "EPISODES.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (args.output / "PREDICTIONS.jsonl").open("w", encoding="utf-8") as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    result = {
        "schema": "MetaSieve.CIPFTERMNestedEvaluation.v1",
        "arch": args.arch,
        "model_seeds": list(model_seeds), "support_sizes": list(SUPPORT_SIZES),
        "split": args.split,
        "split_directory": str(args.split_directory),
        "meta_test": data.seal_record(
            evaluated=bool(args.split == "meta_test")),
        "training_config": {
            "steps": args.steps,
            "episodes_per_step": args.episodes_per_step,
            "val_interval": args.val_interval,
            "val_draws_per_target": args.val_draws_per_target,
            "eval_targets_per_component": args.eval_targets_per_component,
            "learning_rate": args.learning_rate,
            "backbone_lr_scale": args.backbone_lr_scale,
            "binding_loss_weight": args.binding_loss_weight,
            "train_query_size": args.train_query_size,
            "min_query_size": args.min_query_size,
        },
        "targets": len({record.target for record in records}),
        "components": len({record.component for record in records}),
        "manifest_reused_across_model_seeds": True,
        "training": training, "contrasts": contrasts,
        "metric_summary": metric_summary(all_rows),
        "test_used_for_tuning": True,
        "evidence_status": (
            "development_only; this historical meta-test was consumed by "
            "prior architecture search"),
        "gate_authorization": {"G2": False, "G3a": False, "G3b": False},
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
