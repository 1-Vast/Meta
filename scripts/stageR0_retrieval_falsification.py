"""Stage R0: falsify the train-only retrieval prior at k=0.

Gates and arms are fixed in
`report/meta_fewshot/stageR0_retrieval_falsification_20260815/PREREGISTRATION.md`.

No training. At k=0 the support transport is inactive, so `w = 0` is exactly the
checkpoint's zero-shot endpoint and the Stage 10 `transport-beta` defect cannot
apply here. Selection is separated from inference by leave-one-protein-component
-out folds: hyperparameters are chosen on nine components and applied once to the
tenth, and only outer-fold predictions are analysed.

One row is emitted per (checkpoint, query cell). Aggregation is equal-component
then equal-target, so the three seeds are averaged inside a target before
components are resampled: **every interval is conditional on those three trained
seeds**.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import fields
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData, stable_seed
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK, TrainConfig,
    compact_episode, normalized_episode, resolve_architecture, training_label_scale,
)

BETAS = (8.0, 16.0, 24.0, 32.0)
WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)
CONFIDENCE = ("fixed", "novelty_gated")
SELECTABLE_SOURCES = ("ligand", "dual_correct", "dual_centered",
                      "blend_correct", "blend_centered")
PROTEIN_ARMS = ("correct", "centered", "shuffled", "random16", "matched")
FALSIFIER_ARMS = ("shuffled", "random16", "matched")
STRATA = ("all", "exact_free", "exact_overlap", "scaffold_disjoint",
          "scaffold_overlap", "tanimoto_lt40", "tanimoto_40_60",
          "tanimoto_60_80", "near_duplicate")
CLIFF_TANIMOTO = 0.6
CLIFF_GAP_PK = 1.0


def tanimoto_rows(query: np.ndarray, bank: np.ndarray) -> np.ndarray:
    inter = query @ bank.T
    union = query.sum(-1)[:, None] + bank.sum(-1)[None, :] - inter
    return inter / np.maximum(union, 1e-9)


def softmax_rows(scores: np.ndarray, beta: float) -> np.ndarray:
    logits = beta * scores
    shifted = np.exp(logits - logits.max(-1, keepdims=True))
    return shifted / shifted.sum(-1, keepdims=True)


def murcko_scaffolds(smiles: dict[str, str | None]) -> dict[str, str]:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    out: dict[str, str] = {}
    for key, value in smiles.items():
        molecule = Chem.MolFromSmiles(value) if value else None
        if molecule is None:
            out[key] = ""
            continue
        try:
            core = MurckoScaffold.GetScaffoldForMol(molecule)
            out[key] = Chem.MolToSmiles(core) if core is not None else ""
        except Exception:                                   # noqa: BLE001
            out[key] = ""
    return out


def component_target_mean(values) -> float:
    """Equal component, then equal target. `values` are (component, target, x)."""
    by_target: dict[tuple[str, str], list[float]] = defaultdict(list)
    for component, target, value in values:
        if value is not None and np.isfinite(value):
            by_target[(component, target)].append(float(value))
    by_component: dict[str, list[float]] = defaultdict(list)
    for (component, _), items in by_target.items():
        by_component[component].append(float(np.mean(items)))
    if not by_component:
        return float("nan")
    return float(np.mean([np.mean(v) for v in by_component.values()]))


def component_bootstrap(values, draws: int, seed: int) -> dict:
    by_target: dict[tuple[str, str], list[float]] = defaultdict(list)
    for component, target, value in values:
        if value is not None and np.isfinite(value):
            by_target[(component, target)].append(float(value))
    by_component: dict[str, list[float]] = defaultdict(list)
    for (component, _), items in by_target.items():
        by_component[component].append(float(np.mean(items)))
    keys = sorted(by_component)
    if not keys:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "components": 0}
    per_component = np.asarray([float(np.mean(by_component[k])) for k in keys])
    rng = np.random.default_rng(seed)
    index = rng.integers(0, len(keys), size=(draws, len(keys)))
    samples = per_component[index].mean(-1)
    return {"mean": float(per_component.mean()),
            "lo": float(np.quantile(samples, 0.025)),
            "hi": float(np.quantile(samples, 0.975)),
            "components": len(keys)}


def zero_shot_endpoints(checkpoints: list[str], data: QPSMPData, split: str,
                        evaluation_seed: int, query_size: int, device: str | None,
                        support_sizes: tuple[int, ...] = (0,)):
    """Zero-shot endpoint in pK for every k=0 episode, per checkpoint."""
    scale = training_label_scale(data)
    specs = list(data.fixed_nested_episode_banks(
        split, support_sizes, query_size, 1, evaluation_seed, None)[0])
    values: dict[str, dict[str, np.ndarray]] = {}
    for path in checkpoints:
        payload = torch.load(Path(path), map_location="cpu", weights_only=False)
        valid = {f.name for f in fields(TrainConfig)}
        config_values = {k: v for k, v in payload["config"].items() if k in valid}
        if device is not None:
            config_values["device"] = device
        config = TrainConfig(**config_values)
        model = resolve_architecture(config.arch)(
            protein_dim=int(data.protein_bank.manifest["hidden_dim"]),
            hidden_dim=config.hidden_dim, task_dim=config.task_dim,
            ligand_layers=config.ligand_layers, pair_dim=config.pair_dim,
            pair_blocks=config.pair_blocks, pair_latents=config.pair_latents,
            pair_heads=config.pair_heads, pair_chunk_size=config.pair_chunk_size,
            support_hidden_dim=config.support_hidden_dim,
            support_blocks=config.support_blocks, adapter_rank=config.adapter_rank,
            adaptive_blocks=config.adaptive_blocks,
            adapter_scale=config.adapter_scale, use_cartesian=config.use_cartesian)
        model.load_state_dict(payload["model_state"])
        model.to(config.device).eval()
        dtype = next(model.parameters()).dtype
        name = Path(path).parent.name
        values[name] = {}
        with torch.no_grad():
            for spec in specs:
                episode = compact_episode(normalized_episode(
                    data.materialize(spec), scale))

                def cast(value):
                    return value.unsqueeze(0).to(config.device, dtype)

                endpoint = model.encode(
                    cast(episode.protein_pooled), cast(episode.protein_tokens),
                    cast(episode.protein_mask), cast(episode.query_atoms),
                    cast(episode.query_bonds), cast(episode.query_mask),
                    cast(episode.protein_chemistry))[0][0].float().cpu().numpy()
                values[name][spec.target] = endpoint * scale.scale + scale.mean
        del model
        if str(config.device).startswith("cuda"):
            torch.cuda.empty_cache()
    return values, specs


def blended(row: dict, source: str, beta: float, w: float,
            confidence: str) -> float:
    weight = w
    if confidence == "novelty_gated":
        weight = w * float(np.clip((row["novelty"] - 0.2) / 0.4, 0.0, 1.0))
    return (1.0 - weight) * row["f0"] + weight * row[f"{source}|{beta:g}"]


def in_stratum(row: dict, name: str) -> bool:
    novelty = row["novelty"]
    return {
        "all": True,
        "exact_overlap": row["exact_overlap"],
        "exact_free": not row["exact_overlap"],
        "scaffold_overlap": row["scaffold_overlap"],
        "scaffold_disjoint": not row["scaffold_overlap"],
        "tanimoto_lt40": novelty < 0.4,
        "tanimoto_40_60": 0.4 <= novelty < 0.6,
        "tanimoto_60_80": 0.6 <= novelty < 0.8,
        "near_duplicate": novelty >= 0.8,
    }[name]


def field_metrics(selected: list[dict], field: str) -> dict:
    squared = [(r["component"], r["target"], (r[field] - r["truth"]) ** 2)
               for r in selected]
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in selected:
        grouped[(row["component"], row["target"], row["checkpoint"])].append(row)
    ci_values, rho_values = [], []
    for (component, target, _), items in grouped.items():
        if len(items) < 2:
            continue
        prediction = np.asarray([r[field] for r in items])
        truth = np.asarray([r["truth"] for r in items])
        value, comparable = concordance_index(prediction, truth)
        if comparable:
            ci_values.append((component, target, value))
        rho = spearman(prediction, truth)
        if rho is not None and np.isfinite(rho):
            rho_values.append((component, target, float(rho)))
    return {"mse_pk": component_target_mean(squared),
            "ci": component_target_mean(ci_values),
            "spearman": component_target_mean(rho_values),
            "rows": len(selected),
            "targets": len({(r["component"], r["target"]) for r in selected})}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", action="append", required=True)
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--support-sizes", default="0",
                        help="nested bank widths; use 0,1,2,3,5 to reproduce the "
                             "exact Stage 10 k=0 population")
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=20)
    parser.add_argument("--protein-beta", type=float, default=16.0)
    parser.add_argument("--neighbors", type=int, default=16)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--device", default=None)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    fingerprints = data.fingerprints
    scaffolds = murcko_scaffolds(data._ligand_smiles)

    train_cells = [c for c in data.cells if c["split"] == "meta_train"]
    ligand_values: dict[str, list[float]] = defaultdict(list)
    by_target: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for cell in train_cells:
        ligand_values[cell["ligand_id"]].append(float(cell["pK"]))
        by_target[cell["target_id"]].append((cell["ligand_id"], float(cell["pK"])))
    train_ligand_ids = set(ligand_values)
    train_scaffolds = {scaffolds[k] for k in train_ligand_ids if scaffolds[k]}
    train_ligands = sorted(ligand_values)
    train_fp = np.stack([fingerprints[k] for k in train_ligands])
    train_mean = np.asarray([float(np.mean(ligand_values[k])) for k in train_ligands])
    train_targets = sorted(by_target)
    target_fp = {t: np.stack([fingerprints[k] for k, _ in by_target[t]])
                 for t in train_targets}
    target_y = {t: np.asarray([v for _, v in by_target[t]]) for t in train_targets}

    pooled = {t: np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
              for t in set(train_targets) | set(data.tasks[args.split])}
    train_pooled = np.stack([pooled[t] for t in train_targets])
    center = train_pooled.mean(0, keepdims=True)      # meta_train-only statistic

    def unit(matrix: np.ndarray) -> np.ndarray:
        return matrix / np.maximum(np.linalg.norm(matrix, axis=-1, keepdims=True), 1e-9)

    train_raw_unit, train_centered_unit = unit(train_pooled), unit(train_pooled - center)

    support_sizes = tuple(int(v) for v in args.support_sizes.split(",") if v != "")
    endpoints, specs = zero_shot_endpoints(
        args.checkpoint, data, args.split, args.evaluation_seed,
        args.query_size, args.device, support_sizes)
    checkpoint_names = sorted(endpoints)

    rows: list[dict] = []
    for spec in specs:
        query_ids = [data.cells[i]["ligand_id"] for i in spec.query]
        truth = np.asarray([data.cells[i]["pK"] for i in spec.query])
        query_fp = np.stack([fingerprints[k] for k in query_ids])
        ligand_sim = tanimoto_rows(query_fp, train_fp)
        novelty = ligand_sim.max(-1)

        vector = pooled[spec.target]
        raw_sim = train_raw_unit @ (vector / max(float(np.linalg.norm(vector)), 1e-9))
        deviation = vector - center[0]
        centered_sim = train_centered_unit @ (
            deviation / max(float(np.linalg.norm(deviation)), 1e-9))

        rng = np.random.default_rng(stable_seed("stageR0", args.split, spec.target))
        order = np.argsort(-raw_sim)
        correct_top = order[:args.neighbors]
        centered_order = np.argsort(-centered_sim)[:args.neighbors]
        shuffled_sim = rng.permutation(raw_sim)
        shuffled_top = np.argsort(-shuffled_sim)[:args.neighbors]
        random_top = rng.choice(len(train_targets), size=args.neighbors, replace=False)
        outside = order[args.neighbors:]
        outside_sim = raw_sim[outside]
        matched, used = [], set()
        for value in raw_sim[correct_top]:
            for candidate in np.argsort(np.abs(outside_sim - value)):
                if int(candidate) not in used:
                    used.add(int(candidate))
                    matched.append(int(outside[candidate]))
                    break
        matched = np.asarray(matched, dtype=np.int64)
        arm_neighbors = {
            "correct": (correct_top, raw_sim[correct_top]),
            "centered": (centered_order, centered_sim[centered_order]),
            "shuffled": (shuffled_top, shuffled_sim[shuffled_top]),
            "random16": (random_top, np.zeros(len(random_top))),
            "matched": (matched, raw_sim[matched]),
        }

        predictions: dict[str, np.ndarray] = {}
        for beta in BETAS:
            ligand_prediction = (softmax_rows(ligand_sim, beta)
                                 * train_mean[None, :]).sum(-1)
            predictions[f"ligand|{beta:g}"] = ligand_prediction
            for arm in PROTEIN_ARMS:
                index, similarity = arm_neighbors[arm]
                share = (np.ones(len(index)) if arm == "random16"
                         else np.exp(args.protein_beta
                                     * (similarity - similarity.max())))
                numerator = np.zeros(len(query_ids))
                denominator = 1e-9
                for local, position in enumerate(index):
                    target = train_targets[int(position)]
                    weight = softmax_rows(
                        tanimoto_rows(query_fp, target_fp[target]), beta)
                    numerator += float(share[local]) * (
                        weight * target_y[target][None, :]).sum(-1)
                    denominator += float(share[local])
                dual = numerator / denominator
                predictions[f"dual_{arm}|{beta:g}"] = dual
                predictions[f"blend_{arm}|{beta:g}"] = 0.5 * (
                    ligand_prediction + dual)

        for name in checkpoint_names:
            f0 = endpoints[name][spec.target]
            for local, ligand in enumerate(query_ids):
                rows.append({
                    "checkpoint": name, "component": spec.component,
                    "target": spec.target, "ligand_id": ligand,
                    "query_index": local, "truth": float(truth[local]),
                    "f0": float(f0[local]), "novelty": float(novelty[local]),
                    "exact_overlap": ligand in train_ligand_ids,
                    "scaffold_overlap": bool(
                        scaffolds[ligand] and scaffolds[ligand] in train_scaffolds),
                    **{key: float(value[local])
                       for key, value in predictions.items()},
                })

    # ---------------- nested leave-one-component-out selection -------------
    components = sorted({row["component"] for row in rows})
    folds: list[dict] = []
    for held_out in components:
        inner = [r for r in rows
                 if r["component"] != held_out and not r["exact_overlap"]]
        best, best_score = None, float("inf")
        for source in SELECTABLE_SOURCES:
            for beta in BETAS:
                for w in WEIGHTS:
                    for confidence in CONFIDENCE:
                        score = component_target_mean(
                            (r["component"], r["target"],
                             (blended(r, source, beta, w, confidence) - r["truth"]) ** 2)
                            for r in inner)
                        if score < best_score:
                            best, best_score = (source, beta, w, confidence), score
        folds.append({"component": held_out, "source": best[0], "beta": best[1],
                      "weight": best[2], "confidence": best[3],
                      "inner_exact_free_mse": best_score})
    selection = {fold["component"]: fold for fold in folds}

    # Selection bias, measured rather than assumed: the single configuration
    # that minimises exact-free MSE on the *whole* population. This is the
    # development-grade number a non-nested analysis would have reported, and
    # the gap against the outer-fold estimate is the cost of tuning on the data
    # you then infer from.
    free_rows = [r for r in rows if not r["exact_overlap"]]
    tuned, tuned_score = None, float("inf")
    for source in SELECTABLE_SOURCES:
        for beta in BETAS:
            for w in WEIGHTS:
                for confidence in CONFIDENCE:
                    score = component_target_mean(
                        (r["component"], r["target"],
                         (blended(r, source, beta, w, confidence) - r["truth"]) ** 2)
                        for r in free_rows)
                    if score < tuned_score:
                        tuned, tuned_score = (source, beta, w, confidence), score
    oracle_tuned = {"source": tuned[0], "beta": tuned[1], "weight": tuned[2],
                    "confidence": tuned[3], "exact_free_mse": tuned_score,
                    "grade": "development: selected on the population it is "
                             "reported on; shown only to quantify selection bias"}

    for row in rows:
        fold = selection[row["component"]]
        row["selected"] = blended(row, fold["source"], fold["beta"],
                                  fold["weight"], fold["confidence"])
        family = fold["source"].split("_")[0]
        for arm in FALSIFIER_ARMS:
            # A protein-blind selection cannot exhibit protein specificity, so
            # its counterfactual is itself and the G4 contrast is exactly zero.
            row[f"cf_{arm}"] = (row["selected"] if family == "ligand" else
                                blended(row, f"{family}_{arm}", fold["beta"],
                                        fold["weight"], fold["confidence"]))

    # ---------------- strata ------------------------------------------------
    summary: dict[str, dict] = {}
    for name in STRATA:
        selected = [r for r in rows if in_stratum(r, name)]
        if not selected:
            continue
        summary[name] = {
            "f0": field_metrics(selected, "f0"),
            "selected": field_metrics(selected, "selected"),
            **{f"cf_{arm}": field_metrics(selected, f"cf_{arm}")
               for arm in FALSIFIER_ARMS},
            "improvement_bootstrap": component_bootstrap(
                [(r["component"], r["target"],
                  (r["f0"] - r["truth"]) ** 2 - (r["selected"] - r["truth"]) ** 2)
                 for r in selected], args.bootstrap_draws, 20260815),
            **{f"vs_{arm}_bootstrap": component_bootstrap(
                [(r["component"], r["target"],
                  (r[f"cf_{arm}"] - r["truth"]) ** 2
                  - (r["selected"] - r["truth"]) ** 2)
                 for r in selected], args.bootstrap_draws, 20260815)
               for arm in FALSIFIER_ARMS},
            "per_component": {
                component: {
                    "rows": len([r for r in selected
                                 if r["component"] == component]),
                    "f0_mse": component_target_mean(
                        (r["component"], r["target"], (r["f0"] - r["truth"]) ** 2)
                        for r in selected if r["component"] == component),
                    "selected_mse": component_target_mean(
                        (r["component"], r["target"],
                         (r["selected"] - r["truth"]) ** 2)
                        for r in selected if r["component"] == component),
                }
                for component in sorted({r["component"] for r in selected})},
        }

    # ---------------- activity-cliff ordering ------------------------------
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["target"], row["checkpoint"])].append(row)
    cliff_scores: dict[str, list] = {f"{field}|{scope}": []
                                     for field in ("f0", "selected")
                                     for scope in ("all", "exact_free")}
    cliff_pairs = {"all": 0, "exact_free": 0}
    for (target, checkpoint), items in grouped.items():
        items = sorted(items, key=lambda r: r["query_index"])
        matrix = tanimoto_rows(
            np.stack([fingerprints[r["ligand_id"]] for r in items]),
            np.stack([fingerprints[r["ligand_id"]] for r in items]))
        for i in range(len(items)):
            for j in range(i + 1, len(items)):
                gap = items[i]["truth"] - items[j]["truth"]
                if matrix[i, j] < CLIFF_TANIMOTO or abs(gap) < CLIFF_GAP_PK:
                    continue
                free = not (items[i]["exact_overlap"] or items[j]["exact_overlap"])
                if checkpoint == checkpoint_names[0]:
                    cliff_pairs["all"] += 1
                    cliff_pairs["exact_free"] += int(free)
                for field in ("f0", "selected"):
                    correct = float(np.sign(items[i][field] - items[j][field])
                                    == np.sign(gap))
                    cliff_scores[f"{field}|all"].append(
                        (items[i]["component"], target, correct))
                    if free:
                        cliff_scores[f"{field}|exact_free"].append(
                            (items[i]["component"], target, correct))
    cliffs = {"pairs": cliff_pairs,
              **{key: component_target_mean(value)
                 for key, value in cliff_scores.items()}}

    # ---------------- gates -------------------------------------------------
    exact_free, all_cells = summary["exact_free"], summary["all"]
    low = summary.get("tanimoto_lt40")
    gates = {
        "G1_exact_free_lower_bound": {
            "bootstrap": exact_free["improvement_bootstrap"],
            "pass": bool(exact_free["improvement_bootstrap"]["lo"] > 0.0)},
        "G2_low_tanimoto_not_worse": {
            "f0_mse": low["f0"]["mse_pk"] if low else None,
            "selected_mse": low["selected"]["mse_pk"] if low else None,
            "pass": bool(low and low["selected"]["mse_pk"] <= low["f0"]["mse_pk"])},
        "G3_ranking_not_regressed": {
            "d_ci": exact_free["selected"]["ci"] - exact_free["f0"]["ci"],
            "d_spearman": (exact_free["selected"]["spearman"]
                           - exact_free["f0"]["spearman"]),
            "pass": bool(
                exact_free["selected"]["ci"] - exact_free["f0"]["ci"] >= -0.01
                and exact_free["selected"]["spearman"]
                - exact_free["f0"]["spearman"] >= -0.01)},
        "G4_protein_specificity": {
            arm: exact_free[f"vs_{arm}_bootstrap"] for arm in FALSIFIER_ARMS} | {
            "pass": bool(
                exact_free["vs_shuffled_bootstrap"]["lo"] > 0.0
                and exact_free["vs_random16_bootstrap"]["mean"] > 0.0
                and exact_free["vs_matched_bootstrap"]["mean"] > 0.0)},
        "G5_not_only_exact_recall": {
            "exact_free_gain": exact_free["improvement_bootstrap"]["mean"],
            "all_gain": all_cells["improvement_bootstrap"]["mean"],
            "ratio": (exact_free["improvement_bootstrap"]["mean"]
                      / all_cells["improvement_bootstrap"]["mean"]
                      if all_cells["improvement_bootstrap"]["mean"] else float("nan")),
            "pass": bool(all_cells["improvement_bootstrap"]["mean"] > 0
                         and exact_free["improvement_bootstrap"]["mean"]
                         >= 0.4 * all_cells["improvement_bootstrap"]["mean"])},
    }

    payload = {
        "schema": "MetaSieve.StageR0RetrievalFalsification.v1",
        "split": args.split, "episodes": len(specs), "rows": len(rows),
        "query_cells": len(rows) // max(len(checkpoint_names), 1),
        "checkpoints": checkpoint_names,
        "interval_semantics": "conditional on the three trained seeds; seeds are "
                              "averaged inside a target before components are "
                              "resampled",
        "support_sizes": list(support_sizes),
        "selection": folds, "globally_tuned_reference": oracle_tuned,
        "summary": summary, "activity_cliffs": cliffs, "gates": gates,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.with_suffix(".rows.jsonl").write_text(
        "\n".join(json.dumps(r, sort_keys=True) for r in rows) + "\n",
        encoding="utf-8")
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")

    print(f"{args.split}: {len(specs)} episodes, "
          f"{len(rows) // len(checkpoint_names)} query cells, "
          f"{len(checkpoint_names)} checkpoints")
    print("\nselected per outer fold")
    for fold in folds:
        print("  %-10s %-16s beta=%2g w=%-5g %s" % (
            fold["component"][:8], fold["source"], fold["beta"],
            fold["weight"], fold["confidence"]))
    print("\n%-18s %6s %8s %8s %8s %8s" % (
        "stratum", "rows", "f0 MSE", "sel MSE", "dCI", "drho"))
    for name in STRATA:
        if name not in summary:
            continue
        entry = summary[name]
        print("%-18s %6d %8.4f %8.4f %+8.4f %+8.4f" % (
            name, entry["selected"]["rows"], entry["f0"]["mse_pk"],
            entry["selected"]["mse_pk"],
            entry["selected"]["ci"] - entry["f0"]["ci"],
            entry["selected"]["spearman"] - entry["f0"]["spearman"]))
    print("\nexact-free, per component (outer fold)")
    for component, entry in summary["exact_free"]["per_component"].items():
        print("  %-10s n=%4d  f0 %7.4f -> selected %7.4f  (%+7.4f)" % (
            component[:8], entry["rows"], entry["f0_mse"], entry["selected_mse"],
            entry["f0_mse"] - entry["selected_mse"]))
    print("\nglobally tuned reference (development grade, selection on the same "
          "data): %s beta=%g w=%g %s -> exact-free MSE %.4f against f0 %.4f" % (
              oracle_tuned["source"], oracle_tuned["beta"], oracle_tuned["weight"],
              oracle_tuned["confidence"], oracle_tuned["exact_free_mse"],
              summary["exact_free"]["f0"]["mse_pk"]))
    print("\nactivity cliffs: %d pairs (%d exact-free); ordering accuracy "
          "f0 %.4f -> selected %.4f (exact-free %.4f -> %.4f)" % (
              cliffs["pairs"]["all"], cliffs["pairs"]["exact_free"],
              cliffs["f0|all"], cliffs["selected|all"],
              cliffs["f0|exact_free"], cliffs["selected|exact_free"]))
    print("\ngates")
    for name, entry in gates.items():
        print(f"  {name:32s} {'PASS' if entry['pass'] else 'FAIL'}")


if __name__ == "__main__":
    main()
