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
    CORPUS, LIGAND_BANK, PROTEIN_BANK, LabelScale, TrainConfig,
    donor_state, forward, normalized_episode, replay_state, train,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "report/meta_fewshot/qpsmp_nested"
SUPPORT_SIZES = (1, 2, 3, 5)
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
    if support_sizes[0] < 1 or query_size < 1 or draws < 1:
        raise ValueError("episode sizes and draws must be positive")
    max_k = support_sizes[-1]
    eligible = {target: indices for target, indices in data.tasks[split].items()
                if len(indices) >= max_k + 1}
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
                order = rng.permutation(eligible[target])
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
    cell_ids = {cell["cell_id"] for cell in data.cells}
    for record in records:
        support, query = set(record.support_cell_ids), set(record.query_cell_ids)
        if not support or not query or support & query:
            raise ValueError("nested manifest has invalid support/query sets")
        if not support.union(query).issubset(cell_ids):
            raise ValueError("nested manifest references unknown cells")


def episode_spec(data: QPSMPData, record: NestedEpisodeRecord, k: int) -> EpisodeSpec:
    index = {cell["cell_id"]: i for i, cell in enumerate(data.cells)}
    if k < 1 or k > len(record.support_cell_ids):
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
                permuted = forward(model, replace(
                    episode, support_y=episode.support_y.roll(1)))
                foreign = donor_state(model, data, episode, record.donor_target,
                                      label_scale, wrong_protein=False)
                wrong = donor_state(model, data, episode, record.donor_target,
                                    label_scale, wrong_protein=True)
                predictions = {
                    "full": full.prediction,
                    "sar_cut": full.prediction - full.sar_adaptation,
                    "zero_shot": forward(model, episode, adapt=False).prediction,
                    "level_only": episode.support_y.mean().expand_as(episode.query_y),
                    "permuted_state": replay_state(full, permuted.task_state),
                    "foreign_state": replay_state(full, foreign),
                    "wrong_protein_state": replay_state(full, wrong),
                }
                for arm, prediction in predictions.items():
                    mse = float(label_scale.squared_error_pk(
                        prediction, episode.query_y).mean())
                    rows.append({
                        "model_seed": model_seed, "component": record.component,
                        "target": record.target, "draw": record.draw, "k": k,
                        "arm": arm, "mse_pk": mse,
                    })
    return rows


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--manifest-seed", type=int, default=20260850)
    parser.add_argument("--draws-per-target", type=int, default=1)
    parser.add_argument("--query-size", type=int, default=8)
    parser.add_argument("--steps", type=int, default=TrainConfig.steps)
    parser.add_argument("--episodes-per-step", type=int,
                        default=TrainConfig.episodes_per_step)
    parser.add_argument("--val-interval", type=int, default=TrainConfig.val_interval)
    parser.add_argument("--val-draws-per-target", type=int,
                        default=TrainConfig.val_draws_per_target)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--section-mode", choices=("support_span", "ridge", "neural"),
                        default=TrainConfig.section_mode)
    parser.add_argument("--interaction-mode", choices=("pooled", "atom_residue"),
                        default=TrainConfig.interaction_mode)
    parser.add_argument("--device", default=TrainConfig.device)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK)
    if args.manifest:
        payload, records = load_manifest(args.manifest)
        if tuple(payload["support_sizes"]) != SUPPORT_SIZES:
            raise ValueError("manifest support sizes differ from the protocol")
        validate_manifest(data, payload, records)
    else:
        records = build_nested_manifest(
            data, "meta_test", SUPPORT_SIZES, args.query_size,
            args.draws_per_target, args.manifest_seed)
        payload = manifest_payload(records, seed=args.manifest_seed,
                                   support_sizes=SUPPORT_SIZES,
                                   query_size=args.query_size,
                                   corpus_manifest_sha256=file_sha256(
                                       data.corpus / "manifest.json"))
    all_rows, training = [], []
    for model_seed in MODEL_SEEDS:
        config = TrainConfig(
            seed=model_seed, steps=args.steps,
            episodes_per_step=args.episodes_per_step,
            val_interval=args.val_interval,
            val_draws_per_target=args.val_draws_per_target,
            eval_targets_per_component=999999,
            section_mode=args.section_mode,
            interaction_mode=args.interaction_mode,
            device=args.device,
        )
        model, diagnostics, label_scale = train(
            data, config, support_sizes=SUPPORT_SIZES)
        training.append({"model_seed": model_seed, **diagnostics})
        all_rows.extend(evaluate_seed(
            model, data, records, SUPPORT_SIZES, label_scale, model_seed))
    contrasts = [
        paired_component_bootstrap(
            all_rows, "full", control, k, args.bootstrap_draws,
            stable_seed("bootstrap", args.manifest_seed, k, control))
        for k in SUPPORT_SIZES
        for control in ("level_only", "sar_cut", "permuted_state",
                        "foreign_state", "wrong_protein_state")
    ]
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "EPISODES.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with (args.output / "PREDICTIONS.jsonl").open("w", encoding="utf-8") as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    result = {
        "schema": "MetaSieve.QPSMPNestedEvaluation.v1",
        "model_seeds": list(MODEL_SEEDS), "support_sizes": list(SUPPORT_SIZES),
        "training_config": {
            "steps": args.steps,
            "episodes_per_step": args.episodes_per_step,
            "val_interval": args.val_interval,
            "val_draws_per_target": args.val_draws_per_target,
        },
        "targets": len({record.target for record in records}),
        "components": len({record.component for record in records}),
        "manifest_reused_across_model_seeds": True,
        "training": training, "contrasts": contrasts,
        "test_used_for_tuning": False,
        "gate_authorization": {"G2": False, "G3a": False, "G3b": False},
    }
    (args.output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
