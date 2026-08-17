"""Phase 2 analysis: strata, contrasts, component-paired intervals, gates.

Reads only the prediction rows the arms wrote.  No model is loaded and nothing
is trained, so this can be re-run without touching the arms.  The preregistered
thresholds are imported from this module's `GATES` block, which mirrors
`PREREGISTRATION.md` verbatim; the file's SHA-256 is recorded in the output so a
reader can verify that the thresholds were not edited after the fact.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from tools.research.stageS_sar_field.features import ProteinFeatureStore
from tools.research.stageS_sar_field.metrics import (
    Rows, component_paired_bootstrap, metrics_of,
)
from tools.research.stageS_sar_field.pairs import load_data
from tools.research.stageS_sar_field.train import build_banks, pair_id

HERE = Path(__file__).resolve().parent

# Frozen in PREREGISTRATION.md section 7.
GATES = {
    "pearson_gain_threshold": 0.05,
    "shuffled_protein_share_of_gain": 0.5,
    "label_shuffle_pearson_ceiling": 0.10,
    "cliff_sign_floor": 0.50,
    "bootstrap_draws": 2000,
    "novelty_terciles": (0.30136987566947937, 0.5606504082679749),
}
CANDIDATES = ("B_protein", "C_protein_cf")


def load_rows(path: Path) -> Rows:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return Rows(payload["pair_id"], payload["delta_y"], payload["delta_hat"],
                payload["target"], payload["component"])


def pair_novelty(data, banks) -> dict[str, float]:
    """Mean max-Tanimoto-to-fit-ligands of a pair's two ligands."""
    fit, internal = partition_components(data)
    fit_set = set(fit)
    reference: set[str] = set()
    for cell in data.cells:
        if cell["protein_group_40"] in fit_set:
            reference.add(cell["ligand_id"])
    table = data.fingerprints
    keys = sorted(reference)
    matrix = np.stack([table[key].numpy() for key in keys])
    counts = matrix.sum(axis=1)
    cache: dict[str, float] = {}

    def novelty(ligand_id: str) -> float:
        if ligand_id not in cache:
            row = table[ligand_id].numpy()
            intersection = matrix @ row
            union = counts + row.sum() - intersection
            with np.errstate(divide="ignore", invalid="ignore"):
                value = np.where(union > 0,
                                 intersection / np.maximum(union, 1e-12), 0.0)
            cache[ligand_id] = float(value.max()) if value.size else 0.0
        return cache[ligand_id]

    out: dict[str, float] = {}
    for spec in (*banks.internal_same_panel, *banks.internal_cross_panel):
        left = novelty(data.cells[spec.a]["ligand_id"])
        right = novelty(data.cells[spec.b]["ligand_id"])
        out[pair_id(spec)] = 0.5 * (left + right)
    return out


def strata_masks(rows: Rows, specs_by_id: dict, novelty: dict[str, float]
                 ) -> dict[str, np.ndarray]:
    stratum = np.asarray([specs_by_id[key].stratum for key in rows.pair_id])
    values = np.asarray([novelty[key] for key in rows.pair_id])
    low, high = GATES["novelty_terciles"]
    masks = {
        "all": np.ones(len(rows), dtype=bool),
        "cliff": stratum == "cliff",
        "local": stratum == "local",
        "medium": stratum == "medium",
        "distant": stratum == "distant",
        "novelty_low": values < low,
        "novelty_mid": (values >= low) & (values < high),
        "novelty_high": values >= high,
    }
    return {name: mask for name, mask in masks.items() if int(mask.sum()) >= 8}


def contrast(left: Rows, right: Rows, draws: int) -> dict:
    return component_paired_bootstrap(left, right, draws=draws)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=HERE / "runs")
    parser.add_argument("--output", type=Path, default=HERE / "RESULT.json")
    parser.add_argument("--draws", type=int, default=GATES["bootstrap_draws"])
    args = parser.parse_args()

    data = load_data()
    seal = data.seal_record()
    proteins = ProteinFeatureStore(data)
    banks = build_banks(data, proteins, 20260819)
    specs_by_id = {pair_id(spec): spec
                   for spec in (*banks.internal_same_panel,
                                *banks.internal_cross_panel,
                                *banks.fit_unsampled)}
    novelty = pair_novelty(data, banks)

    arms = sorted(path.name for path in args.runs.iterdir() if path.is_dir())
    rows: dict[str, dict[str, Rows]] = {}
    runs: dict[str, dict] = {}
    for arm in arms:
        runs[arm] = json.loads((args.runs / arm / "RUN.json").read_text(
            encoding="utf-8"))
        rows[arm] = {}
        for path in sorted((args.runs / arm).glob("*.rows.json")):
            rows[arm][path.name.replace(".rows.json", "")] = load_rows(path)

    report: dict = {
        "schema": "MetaSieve.StageS.Result.v1",
        "stage": "stageS_sar_field",
        "preregistration_sha256": hashlib.sha256(
            (HERE / "PREREGISTRATION.md").read_bytes()).hexdigest(),
        "gates_frozen": GATES,
        "meta_test": seal,
        "arms": {arm: {"config": runs[arm]["arm"],
                       "parameters": runs[arm]["parameters"],
                       "steps": runs[arm]["steps"],
                       "seed": runs[arm]["seed"],
                       "elapsed_seconds": runs[arm]["elapsed_seconds"],
                       "final_training_loss": runs[arm]["history"][-1]["loss"]}
                 for arm in arms},
        "partition": runs[arms[0]]["partition"],
        "selection": runs[arms[0]]["selection"],
    }

    # -- headline metrics per arm per bank ---------------------------------
    banks_of_interest = sorted({name for arm in arms for name in rows[arm]})
    report["metrics"] = {
        arm: {name: metrics_of(rows[arm][name]) for name in sorted(rows[arm])}
        for arm in arms
    }

    primary = "internal_same_panel_correct"
    wrong = "internal_same_panel_hard_wrong"

    # -- strata -------------------------------------------------------------
    masks = strata_masks(rows[arms[0]][primary], specs_by_id, novelty)
    report["strata_sizes"] = {name: int(mask.sum()) for name, mask in masks.items()}
    report["strata_metrics"] = {
        arm: {name: metrics_of(rows[arm][primary].select(mask))
              for name, mask in masks.items()}
        for arm in arms
    }

    # -- pooled same-panel + cross-panel bank (gate G5) ---------------------
    def pooled(arm: str) -> Rows:
        same = rows[arm][primary]
        cross = rows[arm]["internal_cross_panel_correct"]
        return Rows(same.pair_id + cross.pair_id,
                    np.concatenate([same.delta_y, cross.delta_y]),
                    np.concatenate([same.delta_hat, cross.delta_hat]),
                    same.target + cross.target, same.component + cross.component)

    # -- contrasts ----------------------------------------------------------
    contrasts: dict[str, dict] = {}
    baseline = "A_ligand_only"
    for arm in arms:
        if arm == baseline:
            continue
        contrasts[f"{arm}_vs_{baseline}__same_panel"] = contrast(
            rows[arm][primary], rows[baseline][primary], args.draws)
        contrasts[f"{arm}_vs_{baseline}__pooled_panels"] = contrast(
            pooled(arm), pooled(baseline), args.draws)
        contrasts[f"{arm}_vs_{baseline}__cross_panel"] = contrast(
            rows[arm]["internal_cross_panel_correct"],
            rows[baseline]["internal_cross_panel_correct"], args.draws)
        if wrong in rows[arm]:
            contrasts[f"{arm}__correct_vs_hard_wrong_protein"] = contrast(
                rows[arm][primary], rows[arm][wrong], args.draws)
        for name in ("cliff", "local", "novelty_high"):
            if name in masks:
                contrasts[f"{arm}_vs_{baseline}__{name}"] = contrast(
                    rows[arm][primary].select(masks[name]),
                    rows[baseline][primary].select(masks[name]), args.draws)
    report["contrasts"] = contrasts

    # -- gates --------------------------------------------------------------
    threshold = GATES["pearson_gain_threshold"]
    gates: dict[str, dict] = {}
    for candidate in CANDIDATES:
        if candidate not in rows:
            continue
        same = contrasts[f"{candidate}_vs_{baseline}__same_panel"]
        gain = same["pearson"]["delta"]
        g1 = {"delta": gain, "lo": same["pearson"]["lo"],
              "hi": same["pearson"]["hi"],
              "pass": bool(gain >= threshold and same["pearson"]["lo"] > 0.0)}

        swap = contrasts.get(f"{candidate}__correct_vs_hard_wrong_protein")
        g2 = {"delta": swap["pearson"]["delta"], "lo": swap["pearson"]["lo"],
              "hi": swap["pearson"]["hi"],
              "pass": bool(swap["pearson"]["delta"] >= threshold
                           and swap["pearson"]["lo"] > 0.0)} if swap else {
            "pass": False, "reason": "arm has no protein input"}

        mse = same["equal_component_target_mean_mse"]
        ranking = {name: same[name] for name in ("spearman", "ci", "sign_accuracy")}
        improved = [name for name, value in ranking.items()
                    if value["delta"] > 0 and value["lo"] > 0.0]
        degraded = [name for name, value in ranking.items()
                    if value["delta"] < 0 and value["hi"] < 0.0]
        g3 = {
            "equal_component_target_mean_mse": mse,
            "ranking": ranking,
            "resolved_ranking_improvements": improved,
            "resolved_ranking_degradations": degraded,
            "pass": bool(mse["delta"] < 0 and mse["hi"] < 0.0
                         and improved and not degraded),
        }

        shuffled_gain = contrasts[
            f"D_protein_shuffled_vs_{baseline}__same_panel"]["pearson"]["delta"]
        label_pearson = report["metrics"]["E_label_shuffled"][primary]["pearson"]
        label_contrast = contrast(rows[candidate][primary],
                                  rows["E_label_shuffled"][primary], args.draws)
        g4 = {
            "shuffled_protein_gain": shuffled_gain,
            "candidate_gain": gain,
            "shuffled_protein_pass": bool(
                shuffled_gain <= GATES["shuffled_protein_share_of_gain"] * gain),
            "label_shuffle_pearson": label_pearson,
            "candidate_minus_label_shuffle": label_contrast["pearson"],
            "label_shuffle_pass": bool(
                label_pearson <= GATES["label_shuffle_pearson_ceiling"]
                and label_contrast["pearson"]["delta"] >= threshold
                and label_contrast["pearson"]["lo"] > 0.0),
        }
        g4["pass"] = bool(g4["shuffled_protein_pass"] and g4["label_shuffle_pass"])

        pooled_gain = contrasts[
            f"{candidate}_vs_{baseline}__pooled_panels"]["pearson"]["delta"]
        cross_gain = contrasts[
            f"{candidate}_vs_{baseline}__cross_panel"]["pearson"]["delta"]
        g5 = {
            "same_panel_gain": gain, "pooled_gain": pooled_gain,
            "cross_panel_gain": cross_gain,
            "pass": bool(np.sign(pooled_gain) == np.sign(gain)
                         and gain >= threshold),
        }

        cliff_candidate = report["strata_metrics"][candidate].get("cliff")
        cliff_baseline = report["strata_metrics"][baseline].get("cliff")
        cliff_contrast = contrasts.get(f"{candidate}_vs_{baseline}__cliff")
        if cliff_candidate and cliff_contrast:
            sign_delta = cliff_contrast["sign_accuracy"]
            g6 = {
                "cliff_sign_accuracy": cliff_candidate["sign_accuracy"],
                "cliff_sign_accuracy_baseline": cliff_baseline["sign_accuracy"],
                "cliff_sign_delta": sign_delta,
                "cliff_pearson_gain": cliff_contrast["pearson"]["delta"],
                "pass": bool(
                    cliff_candidate["sign_accuracy"] >= GATES["cliff_sign_floor"]
                    and sign_delta["delta"] >= 0.0
                    and not (sign_delta["hi"] < 0.0)
                    and cliff_contrast["pearson"]["delta"] >= 0.0),
            }
        else:  # pragma: no cover - the stratum is populated on this corpus
            g6 = {"pass": False, "reason": "cliff stratum too small"}

        every = {"G1": g1, "G2": g2, "G3": g3, "G4": g4, "G5": g5, "G6": g6}
        gates[candidate] = {
            **every,
            "all_pass": bool(all(value["pass"] for value in every.values())),
        }
    report["gates"] = gates
    report["verdict"] = {
        "route_passes": bool(any(value["all_pass"] for value in gates.values())),
        "candidates": {name: value["all_pass"] for name, value in gates.items()},
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "verdict": report["verdict"],
        "gates": {arm: {key: value.get("pass")
                        for key, value in gates[arm].items() if key != "all_pass"}
                  for arm in gates},
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
