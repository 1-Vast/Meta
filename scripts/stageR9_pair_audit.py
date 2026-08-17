"""Stage R9 pair-level audit: where does B1's global CI loss come from?

No training. Loads the frozen A0 and R8 B1 checkpoints, runs the identical
fixed meta_val bank, and decomposes every comparable within-target pair into
strata. The R8 result posed the question: B1 improves the shape term
(0.913 -> 0.896) and k=5 activity-cliff sign accuracy (0.675 -> 0.768) while
the global CI drops (0.580 -> 0.535 at k=0). This audit answers which pair
stratum carries the loss, with component bootstraps on the B1-A0
differences. `meta_test` is never read.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData, stable_seed
from scripts.stageR6_compare_arms import (
    load_arm, predict_grammar, predict_reltransport, predict_level_shape,
)
from scripts.stageR0_retrieval_falsification import component_bootstrap, tanimoto_rows
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)
from scripts.train_level_shape import matched_donors, normalized
from scripts.stageR6_compare_arms import SUPPORT_SIZES

CLIFF_SIM = 0.6
CLIFF_GAP = 1.0
LOW_SIM = 0.4


def pair_stats(similarity: np.ndarray, prediction: np.ndarray,
               truth: np.ndarray) -> dict[str, dict]:
    """Per-stratum counts, sign accuracy, margins and RankNet loss."""
    count = len(truth)
    rows, cols = np.triu_indices(count, 1)
    comparable = np.abs(truth[rows] - truth[cols]) > 0
    delta_y = truth[rows] - truth[cols]
    delta_p = prediction[rows] - prediction[cols]
    sim = similarity[rows, cols]
    abs_dy = np.abs(delta_y)
    top = np.argsort(-truth)[:5]
    top_pair = np.isin(rows, top) & np.isin(cols, top)
    strata = {
        "all": np.ones(len(rows), dtype=bool),
        "cliff": (sim >= CLIFF_SIM) & (abs_dy >= CLIFF_GAP),
        "similar_small_gap": (sim >= CLIFF_SIM) & (abs_dy < CLIFF_GAP),
        "high_sim": sim >= CLIFF_SIM,
        "mid_sim": (sim >= LOW_SIM) & (sim < CLIFF_SIM),
        "low_sim": sim < LOW_SIM,
        "large_gap": abs_dy >= CLIFF_GAP,
        "small_gap": abs_dy < 0.5,
        "mid_gap": (abs_dy >= 0.5) & (abs_dy < CLIFF_GAP),
        "top5_pairs": top_pair,
        "non_cliff": ~((sim >= CLIFF_SIM) & (abs_dy >= CLIFF_GAP)),
    }
    stats = {}
    for name, mask in strata.items():
        sel = comparable & mask
        if not bool(sel.any()):
            stats[name] = {"pairs": 0}
            continue
        signed = np.sign(delta_y[sel]) * delta_p[sel]
        concordant = (signed > 0).astype(float) + 0.5 * (signed == 0).astype(float)
        loss = np.log1p(np.exp(-signed))          # RankNet softplus
        stats[name] = {
            "pairs": int(sel.sum()),
            "sign_accuracy": float(concordant.mean()),
            "mean_margin_pk": float(delta_p[sel].mean()),
            "mean_abs_margin_pk": float(np.abs(delta_p[sel]).mean()),
            "mean_ranktnet_loss": float(loss.mean()),
            "mean_abs_dy_pk": float(np.abs(delta_y[sel]).mean()),
        }
    return stats


def predictions_for(model, kind, data, spec, episode, donor, device, dtype):
    if kind == "reltransport":
        values = predict_reltransport(
            model, data, spec, episode, donor, device, dtype)
    elif kind == "level_shape":
        values = predict_level_shape(
            model, data, spec, episode, donor, device, dtype)
    else:
        values = predict_grammar(
            model, data, spec, episode, donor, device, dtype)
    return values["full"].squeeze(0).float().cpu().numpy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", action="append", required=True,
                        help="name=path/to/checkpoint.pt (may repeat; same "
                             "name aggregates seeds)")
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=16)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=False)
    scale = training_label_scale(data)
    donors = matched_donors(data, "meta_val", donor_pool="meta_val",
                            whitening_pool="meta_train")
    fingerprints = data.fingerprints
    banks = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, args.query_size, 1, args.evaluation_seed, None)

    models = {}
    for item in args.arm:
        name, _, path = item.partition("=")
        model, kind, seed = load_arm(Path(path), data, args.device)
        models.setdefault(name, []).append((model, kind, seed))

    # records[arm][k] = list of (target, component, similarity, prediction, truth)
    records: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    for name, entries in models.items():
        for model, kind, seed in entries:
            dtype = next(model.parameters()).dtype
            with torch.no_grad():
                for k, specs in banks.items():
                    for spec in specs:
                        episode = compact_episode(
                            normalized(data.materialize(spec), scale))
                        prediction = predictions_for(
                            model, kind, data, spec, episode, donors[spec.target],
                            args.device, dtype)
                        truth = (episode.query_y.numpy()
                                 * scale.scale + scale.mean)
                        prediction = prediction * scale.scale + scale.mean
                        query_fp = np.stack([
                            fingerprints[data.cells[i]["ligand_id"]].numpy()
                            for i in spec.query])
                        similarity = tanimoto_rows(query_fp, query_fp)
                        records[name][k].append(
                            (spec.target, spec.component, similarity,
                             prediction, truth))
            del model
            torch.cuda.empty_cache()

    # Aggregate: per (stratum, k, arm): pooled stats + per-target sign accuracy
    # + component bootstrap on the B1-A0 per-target differences.
    names = sorted(records)
    audit: dict = {"population": {"targets": len(banks[0]), "components": 19},
                   "strata": {}}
    strata_names = ("all", "cliff", "similar_small_gap", "high_sim", "mid_sim",
                    "low_sim", "large_gap", "small_gap", "mid_gap", "top5_pairs",
                    "non_cliff")
    for k in SUPPORT_SIZES:
        for stratum in strata_names:
            cell = {}
            per_target: dict[str, dict] = {name: {} for name in names}
            for name in names:
                pooled = {field: [] for field in (
                    "sign_accuracy", "mean_margin_pk", "mean_abs_margin_pk",
                    "mean_ranktnet_loss")}
                counts = []
                for target, component, sim, pred, truth in records[name][k]:
                    stats = pair_stats(sim, pred, truth).get(stratum, {})
                    if not stats.get("pairs", 0):
                        continue
                    counts.append(stats["pairs"])
                    per_target[name][(component, target)] = stats["sign_accuracy"]
                    for field in pooled:
                        pooled[field].append(stats[field])
                cell[name] = {
                    "pairs": int(sum(counts)),
                    "targets": len(counts),
                }
                for field in pooled:
                    cell[name][field] = (
                        float(np.mean(pooled[field])) if pooled[field] else None)
            pairs_targets = set(per_target[names[0]])
            for name in names[1:]:
                pairs_targets &= set(per_target[name])
            for name in names[1:]:
                values = []
                for key in pairs_targets:
                    component, target = key
                    values.append((component, target,
                                   per_target[names[0]][key]
                                   - per_target[name][key]))
                cell[f"{name}_minus_{names[0]}"] = component_bootstrap(
                    values, args.bootstrap_draws, 20260816)
            audit["strata"][f"k{k}_{stratum}"] = cell

    # Target-level diagnostics: per-target B1-A0 CI/shape/calibration changes.
    rows = []
    for name in names:
        for k, entries in records[name].items():
            for target, component, sim, pred, truth in entries:
                error = pred - truth
                count = len(truth)
                idx_rows, idx_cols = np.triu_indices(count, 1)
                comparable = np.abs(truth[idx_rows] - truth[idx_cols]) > 0
                signed = (truth[idx_rows] - truth[idx_cols]) * (
                    pred[idx_rows] - pred[idx_cols])
                ci = ((signed > 0).sum() + 0.5 * (signed == 0).sum()) \
                    / comparable.sum()
                rows.append({"arm": name, "k": k, "component": component,
                             "target": target, "ci": float(ci),
                             "calibration_pk": float(error.mean() ** 2),
                             "shape_pk": float(error.var())})
    target_corr = {}
    for k in SUPPORT_SIZES:
        base = {r["target"]: r for r in rows
                if r["arm"] == names[0] and r["k"] == k}
        for name in names[1:]:
            joined = [(r, base[r["target"]]) for r in rows
                      if r["arm"] == name and r["k"] == k
                      and r["target"] in base]
            def corr(field):
                a = np.asarray([r[field] for r, _ in joined])
                b = np.asarray([b[field] for _, b in joined])
                return float(np.corrcoef(a - b, a)[0, 1]) if len(a) > 2 else None
            target_corr[f"k{k}_{name}_vs_{names[0]}"] = {
                "corr_ci_change_with_ci": corr("ci"),
                "corr_ci_change_with_shape_change": corr("shape_pk"),
                "corr_ci_change_with_calibration_change": corr("calibration_pk"),
                "mean_ci_change": float(np.mean([r["ci"] - b["ci"]
                                                 for r, b in joined])),
            }

    payload = {
        "schema": "MetaSieve.StageR9PairAudit.v1",
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "arms": names, "support_sizes": list(SUPPORT_SIZES),
        "strata_definition": {
            "cliff": f"Tanimoto >= {CLIFF_SIM} and |dy| >= {CLIFF_GAP} pK",
            "similar_small_gap": f"Tanimoto >= {CLIFF_SIM} and |dy| < {CLIFF_GAP}",
            "high_sim": f"Tanimoto >= {CLIFF_SIM}",
            "mid_sim": f"{LOW_SIM} <= Tanimoto < {CLIFF_SIM}",
            "low_sim": f"Tanimoto < {LOW_SIM}",
            "large_gap": f"|dy| >= {CLIFF_GAP} pK",
            "small_gap": "|dy| < 0.5 pK",
            "mid_gap": "0.5 <= |dy| < 1.0 pK",
            "top5_pairs": "both ligands in the panel's top-5 affinities",
            "non_cliff": "everything except cliff",
        },
        "strata": audit["strata"],
        "target_correlations": target_corr,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print("strata summary (k=0):")
    for stratum in strata_names:
        cell = audit["strata"][f"k0_{stratum}"]
        line = "  %-18s" % stratum
        for name in names:
            e = cell[name]
            line += " | %s pairs=%d sign=%.3f" % (
                name, e["pairs"], e["sign_accuracy"] or 0.0)
        print(line)


if __name__ == "__main__":
    main()
