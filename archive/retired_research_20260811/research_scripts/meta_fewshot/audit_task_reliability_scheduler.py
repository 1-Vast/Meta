"""Development-only AdaMBind-inspired task scheduler falsification audit."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import torch

from model.metasieve_v1 import MetaSieveV1
from research.meta_fewshot.task_reliability_scheduler import (
    apply_cross_fitted,
    component_bootstrap,
    component_macro,
    cross_fitted_ridge,
    permute_informative_rows,
    residualize,
    spearman,
)
from research.meta_fewshot.train_main_v0 import CORPUS, FEATURES, load_data, sha256
from research.meta_fewshot.train_main_v1 import (
    _flat_batched_gradients,
    build_tasks,
    cluster_tasks,
    draw_episode,
    gradient_cosine_rows,
    pack_episode_indices,
)


ROOT = Path(__file__).resolve().parents[2]
CHECKPOINT_DIR = (
    ROOT / "report/meta_fewshot/main_v1_support_only_mlp64_d4_final/checkpoints")
OUT = ROOT / "report/meta_fewshot/task_reliability_scheduler_v1"
SEEDS = (20260831, 20260832, 20260833)
SUPPORT_SIZES = (1, 2, 3, 5)
CONTROLS = (
    "protein_shuffle", "wrong_support", "label_permutation",
    "ligand_only", "intercept_only")


def stable_seed(seed: int, *parts: object) -> int:
    text = "|".join((str(seed), *(str(part) for part in parts)))
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def load_model(checkpoint_path: Path, device: str) -> tuple[MetaSieveV1, dict]:
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    config = checkpoint["config"]
    model = MetaSieveV1(
        input_dim=288, section_dim=config["section_dim"], ridge=config["ridge"],
        support_only_section=config["support_only_section"],
        population_hidden_dim=config["population_hidden_dim"],
        pair_hidden_dim=config["pair_hidden_dim"],
    ).to(device)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    return model, checkpoint


def task_covariates(cells: list[dict], tasks: dict[str, np.ndarray],
                    validation_scaffolds: set[str],
                    scaffold_by_ligand: dict[str, str]
                    ) -> tuple[np.ndarray, list[str], dict]:
    names = [
        "log_task_size", "label_spread",
        "scaffold_overlap_meta_val_nonempty", "scaffold_missing_fraction",
        "mean_log1p_replicate_count", "multi_replicate_fraction",
        "mean_log1p_panel_count", "multi_panel_fraction",
    ]
    rows = []
    for target in sorted(tasks):
        values = [cells[int(index)] for index in tasks[target]]
        replicates = np.asarray([row["replicate_count"] for row in values], dtype=float)
        panels = np.asarray([row["panel_count"] for row in values], dtype=float)
        scaffolds = [scaffold_by_ligand[row["ligand_id"]] for row in values]
        known_scaffolds = [scaffold for scaffold in scaffolds if scaffold]
        rows.append([
            np.log1p(len(values)), np.std([row["pK"] for row in values]),
            (np.mean([scaffold in validation_scaffolds for scaffold in known_scaffolds])
             if known_scaffolds else np.nan),
            1.0 - len(known_scaffolds) / len(scaffolds),
            np.mean(np.log1p(replicates)), np.mean(replicates > 1),
            np.mean(np.log1p(panels)), np.mean(panels > 1),
        ])
    result = np.asarray(rows, dtype=np.float64)
    missing = ~np.isfinite(result[:, 2])
    if missing.all():
        raise ValueError("scaffold familiarity is unavailable for every task")
    result[missing, 2] = np.mean(result[~missing, 2])
    metadata = {
        "empty_scaffolds_excluded_from_overlap": True,
        "all_scaffolds_missing_tasks": int(missing.sum()),
        "all_missing_task_overlap_handling": (
            "mean_imputation_for_nuisance_residualization_plus_missing_fraction"),
        "enters_scheduler_scorer": False,
    }
    return result, names, metadata


def episodes_for_tasks(tasks: dict[str, np.ndarray], targets: list[str], *,
                       seed: int, k: int, max_query: int = 32
                       ) -> list[tuple[np.ndarray, np.ndarray]]:
    return [draw_episode(
        tasks[target], np.random.default_rng(stable_seed(seed, k, target)),
        k, max_query) for target in targets]


def different_group_donors(cells: list[dict], targets: list[str]) -> dict[str, str]:
    group = {row["target_id"]: row["protein_group_40"] for row in cells}
    result = {}
    for index, target in enumerate(targets):
        for offset in range(1, len(targets)):
            donor = targets[(index + offset) % len(targets)]
            if group[donor] != group[target]:
                result[target] = donor
                break
        if target not in result:
            raise ValueError("no different-component support donor")
    return result


def _final_population_bias(model: MetaSieveV1) -> torch.Tensor:
    if isinstance(model.population, torch.nn.Sequential):
        return model.population[-1].bias.squeeze()
    return model.population.bias.squeeze()


def episode_view(
        model: MetaSieveV1, tensors: dict,
        episodes: list[tuple[np.ndarray, np.ndarray]], *, family: str,
        mode: str = "full", support_override: torch.Tensor | None = None,
        query_override: torch.Tensor | None = None,
        chunk_size: int = 32) -> tuple[np.ndarray, np.ndarray, torch.Tensor]:
    """Return difficulty, support/query alignment and per-task query gradients."""
    parameters = tuple(model.parameters())
    difficulty, alignment, query_gradients = [], [], []
    for start in range(0, len(episodes), chunk_size):
        stop = min(start + chunk_size, len(episodes))
        support_index, query_index, query_mask = pack_episode_indices(
            episodes[start:stop], tensors["y"].device)
        support_y = (tensors["y"][support_index] if support_override is None
                     else support_override[start:stop])
        query_y = (tensors["y"][query_index] if query_override is None
                   else query_override[start:stop, :query_index.shape[1]])
        if mode == "full":
            prediction, support_prediction = model.batched_episode(
                tensors["ligand"][support_index], tensors[family][support_index],
                support_y, tensors["ligand"][query_index],
                tensors[family][query_index])
        elif mode == "ligand_only":
            support_prediction = model.population(
                tensors["ligand"][support_index]).squeeze(-1)
            prediction = model.population(
                tensors["ligand"][query_index]).squeeze(-1)
        elif mode == "intercept_only":
            bias = _final_population_bias(model)
            support_prediction = bias.expand_as(support_y)
            prediction = bias.expand_as(query_y)
        else:
            raise ValueError(f"unknown audit view: {mode}")
        squared = (prediction - query_y).square()
        query_loss = ((squared * query_mask).sum(dim=1)
                      / query_mask.sum(dim=1).to(squared.dtype))
        support_loss = (support_prediction - support_y).square().mean(dim=1)
        support_gradient = _flat_batched_gradients(
            support_loss, parameters, retain_graph=True).detach()
        query_gradient = _flat_batched_gradients(
            query_loss, parameters, retain_graph=False).detach()
        difficulty.append(torch.log1p(query_loss.detach()).cpu().numpy())
        alignment.append(gradient_cosine_rows(
            support_gradient, query_gradient).cpu().numpy())
        query_gradients.append(query_gradient)
    return (
        np.concatenate(difficulty), np.concatenate(alignment),
        torch.cat(query_gradients),
    )


def validation_reference(gradients: torch.Tensor, groups: np.ndarray) -> torch.Tensor:
    references = []
    for group in np.unique(groups):
        rows = torch.as_tensor(groups == group, device=gradients.device)
        reference = gradients[rows].mean(dim=0)
        references.append(reference / reference.norm().clamp_min(1e-12))
    return torch.stack(references)


def transfer_utility(gradients: torch.Tensor, references: torch.Tensor) -> np.ndarray:
    normalized = gradients / gradients.norm(dim=1, keepdim=True).clamp_min(1e-12)
    return (normalized @ references.T).mean(dim=1).detach().cpu().numpy()


def permuted_labels(episodes: list[tuple[np.ndarray, np.ndarray]], tensors: dict
                    ) -> tuple[torch.Tensor, torch.Tensor]:
    support_index, query_index, query_mask = pack_episode_indices(
        episodes, tensors["y"].device)
    support = tensors["y"][support_index]
    support = torch.roll(support.flatten(), 1).reshape_as(support)
    query = tensors["y"][query_index].clone()
    valid = query[query_mask]
    query[query_mask] = torch.roll(valid, 1)
    return support, query


def simple_metrics(score: np.ndarray, utility: np.ndarray, null_score: np.ndarray,
                   transformed_scores: dict[str, np.ndarray],
                   transformed_utilities: dict[str, np.ndarray],
                   nuisance: np.ndarray, groups: np.ndarray) -> dict:
    _, joined = component_macro(np.column_stack((
        score, utility, null_score,
        *[transformed_scores[name] for name in CONTROLS],
        *[transformed_utilities[name] for name in CONTROLS], nuisance)), groups)
    clean, actual, null = joined[:, 0], joined[:, 1], joined[:, 2]
    score_controls = {name: joined[:, 3 + index]
                      for index, name in enumerate(CONTROLS)}
    offset = 3 + len(CONTROLS)
    utility_controls = {name: joined[:, offset + index]
                        for index, name in enumerate(CONTROLS)}
    nuisance_m = joined[:, offset + len(CONTROLS):]
    correlation = spearman(clean, actual)
    return {
        "clean_correlation": correlation,
        "matched_null_correlation": spearman(null, actual),
        "matched_null_advantage": correlation - spearman(null, actual),
        "score_correlation_losses": {
            name: correlation - spearman(values, actual)
            for name, values in score_controls.items()},
        "utility_deltas": {
            name: float(np.mean(actual - values))
            for name, values in utility_controls.items()},
        "nuisance_residual_advantage": (
            spearman(residualize(clean, nuisance_m), residualize(actual, nuisance_m))
            - spearman(residualize(null, nuisance_m), residualize(actual, nuisance_m))),
    }


def audit_seed_k(model: MetaSieveV1, cells: list[dict], tensors: dict, *,
                 seed: int, k: int, covariates: np.ndarray) -> dict:
    source_tasks = build_tasks(cells, "meta_train", k, 3)
    validation_tasks = build_tasks(cells, "meta_val", k, 3)
    source_targets = sorted(source_tasks)
    validation_targets = sorted(validation_tasks)
    source_episodes = episodes_for_tasks(
        source_tasks, source_targets, seed=seed, k=k)
    validation_episodes = episodes_for_tasks(
        validation_tasks, validation_targets, seed=seed + 700_001, k=k)
    _, _, validation_gradients = episode_view(
        model, tensors, validation_episodes, family="correct")
    group_by_target = {row["target_id"]: row["protein_group_40"] for row in cells}
    validation_groups = np.asarray([
        group_by_target[target] for target in validation_targets])
    references = validation_reference(validation_gradients, validation_groups)

    clean_difficulty, clean_alignment, clean_gradient = episode_view(
        model, tensors, source_episodes, family="correct")
    clean_utility = transfer_utility(clean_gradient, references)
    clean_features = np.column_stack((clean_difficulty, clean_alignment))
    source_groups = np.asarray([group_by_target[target] for target in source_targets])
    score, fitted = cross_fitted_ridge(
        clean_features, clean_utility, source_groups)
    null_features = permute_informative_rows(
        clean_features, covariates[:, 0], seed=stable_seed(seed, k, "null"))
    null_score, _ = cross_fitted_ridge(
        null_features, clean_utility, source_groups)

    transformed_stats: dict[str, tuple[np.ndarray, np.ndarray, torch.Tensor]] = {}
    transformed_stats["protein_shuffle"] = episode_view(
        model, tensors, source_episodes, family="wrong")
    donors = different_group_donors(cells, source_targets)
    donor_episodes = episodes_for_tasks(
        source_tasks, [donors[target] for target in source_targets],
        seed=seed + 900_001, k=k)
    wrong_support_episodes = [
        (donor[0], recipient[1])
        for donor, recipient in zip(donor_episodes, source_episodes)]
    transformed_stats["wrong_support"] = episode_view(
        model, tensors, wrong_support_episodes, family="correct")
    support_y, query_y = permuted_labels(source_episodes, tensors)
    transformed_stats["label_permutation"] = episode_view(
        model, tensors, source_episodes, family="correct",
        support_override=support_y, query_override=query_y)
    transformed_stats["ligand_only"] = episode_view(
        model, tensors, source_episodes, family="correct", mode="ligand_only")
    transformed_stats["intercept_only"] = episode_view(
        model, tensors, source_episodes, family="correct", mode="intercept_only")

    transformed_scores, transformed_utilities = {}, {}
    for name, (difficulty, alignment, gradient) in transformed_stats.items():
        features = np.column_stack((difficulty, alignment))
        transformed_scores[name] = apply_cross_fitted(fitted, features)
        transformed_utilities[name] = transfer_utility(gradient, references)
    return {
        "targets": source_targets,
        "groups": source_groups,
        "score": score,
        "utility": clean_utility,
        "null_score": null_score,
        "transformed_scores": transformed_scores,
        "transformed_utilities": transformed_utilities,
        "nuisance": covariates,
        "metrics": simple_metrics(
            score, clean_utility, null_score, transformed_scores,
            transformed_utilities, covariates, source_groups),
        "selected_alphas": [fold.alpha for fold in fitted],
    }


def gate_decision(bootstrap: dict, seed_metrics: list[dict]) -> dict:
    observed, lower = bootstrap["observed"], bootstrap["lcb95"]
    score_loss = observed["score_correlation_losses"]
    utility_delta = observed["utility_deltas"]
    criteria = {
        "matched_null_superiority": (
            observed["matched_null_advantage"] >= 0.05
            and lower["matched_null_advantage"] > 0),
        "absolute_utility_tracking": (
            observed["clean_correlation"] >= 0.20
            and lower["clean_correlation"] > 0),
        "partner_support_necessity": all(
            utility_delta[name] > 0
            and lower[f"utility_delta:{name}"] > 0
            and score_loss[name] >= 0.05
            and lower[f"score_loss:{name}"] > 0
            for name in ("protein_shuffle", "wrong_support")),
        "label_necessity": (
            score_loss["label_permutation"] >= 0.05
            and lower["score_loss:label_permutation"] > 0),
        "shortcut_controls": all(
            score_loss[name] >= 0.05 and lower[f"score_loss:{name}"] > 0
            for name in ("ligand_only", "intercept_only")),
        "nuisance_survival": (
            observed["nuisance_residual_advantage"] > 0
            and lower["nuisance_residual_advantage"] > 0),
    }
    seed_direction = all(
        metric["matched_null_advantage"] > 0
        and metric["clean_correlation"] > 0
        and metric["nuisance_residual_advantage"] > 0
        and all(metric["score_correlation_losses"][name] > 0 for name in CONTROLS)
        and all(metric["utility_deltas"][name] > 0
                for name in ("protein_shuffle", "wrong_support"))
        for metric in seed_metrics)
    criteria["all_seed_directions"] = seed_direction
    return {"criteria": criteria, "pass": all(criteria.values())}


def availability(cells: list[dict]) -> dict:
    source = build_tasks(cells, "meta_train", 5, 3)
    validation = build_tasks(cells, "meta_val", 5, 3)
    components = cluster_tasks(cells, validation)
    largest = max(map(len, components.values())) / len(validation)
    return {
        "source_targets_k5": len(source),
        "source_components_k5": len(cluster_tasks(cells, source)),
        "meta_val_targets_k5": len(validation),
        "meta_val_components_k5": len(components),
        "largest_meta_val_component_share": largest,
        "empirical_replicate_disagreement": "NA",
        "continuous_protein_familiarity": "NA",
        "scaffold_familiarity": "available_from_governed_ligand_records",
        "pass": len(source) >= 100 and len(components) >= 8 and largest <= 0.35,
    }


def main() -> None:
    started = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    parser.add_argument("--support-sizes", type=int, nargs="+",
                        default=list(SUPPORT_SIZES))
    args = parser.parse_args()
    if args.device != "cuda" or not torch.cuda.is_available():
        raise RuntimeError("formal scheduler audit requires CUDA")
    cells, tensors, _, normalization = load_data(args.device)
    gate0 = availability(cells)
    args.out.mkdir(parents=True, exist_ok=True)
    if not gate0["pass"]:
        (args.out / "RESULT.json").write_text(json.dumps({
            "schema": "MetaSieve.TaskReliabilitySchedulerAudit.v1",
            "development_only": True, "gate0": gate0,
            "terminal_verdict": "GATE0_FAIL_CLOSED",
        }, indent=2), encoding="utf-8")
        return

    ligands = {
        row["drug_key"]: row["scaffold"]
        for row in map(json.loads, (CORPUS / "ligands.jsonl").read_text().splitlines())}
    validation_scaffolds = {
        ligands[row["ligand_id"]] for row in cells
        if row["split"] == "meta_val" and ligands[row["ligand_id"]]}
    raw: dict[int, list[dict]] = defaultdict(list)
    checkpoint_records = []
    requested_seeds = tuple(args.seeds)
    requested_support_sizes = tuple(args.support_sizes)
    if not set(requested_support_sizes).issubset(SUPPORT_SIZES):
        raise ValueError("unsupported few-shot size")
    formal = (requested_seeds == SEEDS
              and requested_support_sizes == SUPPORT_SIZES
              and args.bootstrap_draws == 9999)
    for seed in requested_seeds:
        path = CHECKPOINT_DIR / f"uniform_clean_seed{seed}.pt"
        model, checkpoint = load_model(path, args.device)
        checkpoint_records.append({"seed": seed, "path": str(path.relative_to(ROOT)),
                                   "sha256": sha256(path)})
        for k in requested_support_sizes:
            tasks = build_tasks(cells, "meta_train", k, 3)
            covariates, covariate_names, scaffold_handling = task_covariates(
                cells, tasks, validation_scaffolds, ligands)
            raw[k].append(audit_seed_k(
                model, cells, tensors, seed=seed, k=k, covariates=covariates))

    by_k = {}
    all_pass = True
    for k in requested_support_sizes:
        first = raw[k][0]
        aggregate = {
            "score": np.mean([item["score"] for item in raw[k]], axis=0),
            "utility": np.mean([item["utility"] for item in raw[k]], axis=0),
            "null_score": np.mean([item["null_score"] for item in raw[k]], axis=0),
            "transformed_scores": {
                name: np.mean([item["transformed_scores"][name]
                               for item in raw[k]], axis=0) for name in CONTROLS},
            "transformed_utilities": {
                name: np.mean([item["transformed_utilities"][name]
                               for item in raw[k]], axis=0) for name in CONTROLS},
            "nuisance": first["nuisance"],
        }
        bootstrap = component_bootstrap(
            aggregate["score"], aggregate["utility"], aggregate["null_score"],
            aggregate["transformed_scores"], aggregate["transformed_utilities"],
            aggregate["nuisance"], first["groups"],
            seed=stable_seed(20260811, k, "bootstrap"), draws=args.bootstrap_draws)
        decision = gate_decision(bootstrap, [item["metrics"] for item in raw[k]])
        all_pass = all_pass and decision["pass"]
        by_k[str(k)] = {
            "targets": len(first["targets"]), "components": bootstrap["components"],
            "seed_metrics": [item["metrics"] for item in raw[k]],
            "selected_alphas": [item["selected_alphas"] for item in raw[k]],
            "bootstrap": bootstrap, "gate": decision,
        }
    result = {
        "schema": "MetaSieve.TaskReliabilitySchedulerAudit.v1",
        "development_only": True,
        "scientific_axis": "task_reliability_transferability_scheduler",
        "device": args.device,
        "preregistration": str((ROOT / "research/meta_fewshot/PREREG_TASK_RELIABILITY_SCHEDULER_V1.md").relative_to(ROOT)),
        "data_firewall": {"used_for_fit_or_score": ["meta_train", "meta_val"],
                          "not_indexed_scored_or_selected": [
                              "meta_test", "fresh_confirmation"],
                          "meta_test_labels_used": False},
        "gate0": gate0,
        "covariates": covariate_names,
        "scheduler_input_columns": [
            "log1p_clean_query_mse", "support_query_gradient_cosine"],
        "covariates_enter_scheduler_scorer": False,
        "scaffold_missingness": scaffold_handling,
        "covariate_semantics": {
            "difficulty": "log1p_clean_query_mse",
            "transferability": "support_query_gradient_cosine_and_meta_val_first_order_utility",
            "familiarity": "exact_Bemis_Murcko_overlap_only; continuous_protein_similarity_NA",
            "reliability": "NA; provenance density is not replicate disagreement",
            "biological_controls": "offline_only_never_scheduler_supervision",
        },
        "corpus_sha256": sha256(CORPUS / "cells.jsonl.gz"),
        "features_sha256": sha256(FEATURES),
        "checkpoints": checkpoint_records,
        "provenance": {
            "preregistration_sha256": sha256(
                ROOT / "research/meta_fewshot/PREREG_TASK_RELIABILITY_SCHEDULER_V1.md"),
            "runner_sha256": sha256(Path(__file__)),
            "scorer_sha256": sha256(
                ROOT / "research/meta_fewshot/task_reliability_scheduler.py"),
            "model_sha256": sha256(ROOT / "model/metasieve_v1.py"),
            "git_head": subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
                capture_output=True, text=True).stdout.strip(),
            "torch": torch.__version__,
            "cuda_device": torch.cuda.get_device_name(),
            "elapsed_seconds": time.perf_counter() - started,
        },
        "by_k": by_k,
        "run_scope": "formal" if formal else "smoke_only",
        "gate1_all_k_pass": all_pass if formal else False,
        "short_cuda_training_authorized": formal and all_pass,
        "training_code_integration_authorized": formal and all_pass,
        "terminal_verdict": (
            ("ADOPT_WITH_MODIFICATION" if all_pass
             else "REJECT_TASK_SCHEDULER_GATE1_FAIL_CLOSED")
            if formal else "SMOKE_ONLY_NO_SCIENTIFIC_DECISION"),
        "biology_admission_authorized": False,
        "production_migration_authorized": False,
    }
    (args.out / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "scope": result["run_scope"], "gate0": gate0["pass"],
        "gate1_all_k": result["gate1_all_k_pass"],
        "verdict": result["terminal_verdict"],
        "out": str(args.out / "RESULT.json")}, indent=2))


if __name__ == "__main__":
    main()
