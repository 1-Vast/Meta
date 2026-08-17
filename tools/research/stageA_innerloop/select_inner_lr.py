"""Choose the inner step size on `meta_train` only, by train-only cross-fitting.

A held-out objective is needed to pick a step size, and `meta_val` is not
available for it. So the objective is built inside `meta_train`: for each
target, adapt on that target's support ligands and score on that target's
*remaining* ligands. Targets are grouped into folds by homology component, so a
fold's score is a genuine unseen-component estimate rather than a memorized one.

The sweep runs on the frozen A0 checkpoint rather than a fresh initialization.
A random model's adaptation curve says nothing about the regime the experiment
will actually train in, and picking a step size there would be picking it for
the wrong problem.

No meta_val is read. No meta_test is constructible.
"""
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

from scripts.qpsmp_data import EpisodeSpec, QPSMPData, stable_seed    # noqa: E402
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.inner_loop import AdaptationConfig  # noqa: E402
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, episode_tensors, predict,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
A0 = (ROOT / "report/meta_fewshot/stageR3R4_level_shape_20260815"
      / "A0_incumbent_seed20260815" / "checkpoint.pt")
OUT = Path(__file__).resolve().parent / "INNER_LR_SELECTION.json"

LR_GRID = (0.003, 0.01, 0.03, 0.1, 0.3, 1.0)
STEP_GRID = (1, 2, 3)
SUPPORT = 5
HELD = 8
FOLDS = 5
FOLD_SEED = 20260818
SELECTION_SEED = 20260817
MAX_TARGETS = 120


def build_tasks(data: QPSMPData) -> list[tuple[str, str, EpisodeSpec]]:
    """One (component, target, spec) per eligible meta_train target."""
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    tasks = []
    for target in sorted(data.tasks["meta_train"]):
        rng = np.random.default_rng(
            stable_seed("inner-lr", SELECTION_SEED, target))
        order = data._unique_ligand_order(data.tasks["meta_train"][target], rng)
        if len(order) < SUPPORT + 2:
            continue
        support = tuple(int(i) for i in order[:SUPPORT])
        held = tuple(int(i) for i in order[SUPPORT:SUPPORT + HELD])
        donor = target
        tasks.append((component_of[target], target,
                      EpisodeSpec("meta_train", component_of[target], target,
                                  support, held, donor)))
    rng = np.random.default_rng(SELECTION_SEED)
    if len(tasks) > MAX_TARGETS:
        pick = rng.choice(len(tasks), size=MAX_TARGETS, replace=False)
        tasks = [tasks[int(i)] for i in sorted(pick)]
    return tasks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    model, kind, seed = load_arm(A0, data, arguments.device)
    model.eval()

    tasks = build_tasks(data)
    components = sorted({component for component, _, _ in tasks})
    order = np.random.default_rng(FOLD_SEED).permutation(len(components))
    fold_of = {components[int(index)]: rank % FOLDS
               for rank, index in enumerate(order)}
    print(f"{len(tasks)} meta_train tasks across {len(components)} components")

    cached = []
    for component, target, spec in tasks:
        episode = compact_episode(normalized_episode(
            data.materialize(spec), label_scale))
        parts = align_atoms(episode_tensors(
            model, episode, arguments.device, torch.float32))
        cached.append((fold_of[component], component, target, parts))

    def score(config: AdaptationConfig) -> list[tuple[int, float]]:
        out = []
        for fold, _, _, parts in cached:
            with torch.no_grad():
                output = predict(model, parts, config)
                error = float(F.mse_loss(output["prediction"],
                                         parts["query_y"]))
            out.append((fold, error * label_scale.scale ** 2))
        return out

    baseline = score(AdaptationConfig(inner_steps=0))
    baseline_mse = float(np.mean([v for _, v in baseline]))
    table = {}
    for steps in STEP_GRID:
        for lr in LR_GRID:
            values = score(AdaptationConfig(inner_steps=steps, inner_lr=lr))
            per_fold = [float(np.mean([v for f, v in values if f == fold]))
                        for fold in range(FOLDS)]
            table[f"steps{steps}_lr{lr:g}"] = {
                "inner_steps": steps, "inner_lr": lr,
                "held_out_mse_pk": float(np.mean([v for _, v in values])),
                "per_fold_mse_pk": per_fold,
                "folds_improved": int(sum(
                    1 for fold in range(FOLDS)
                    if per_fold[fold] < float(np.mean(
                        [v for f, v in baseline if f == fold])))),
            }
            print(f"  steps={steps} lr={lr:<6g} "
                  f"held-out MSE {table[f'steps{steps}_lr{lr:g}']['held_out_mse_pk']:.4f} "
                  f"(baseline {baseline_mse:.4f})")

    best_key = min(table, key=lambda key: table[key]["held_out_mse_pk"])
    best_at_one_step = min(
        (k for k in table if table[k]["inner_steps"] == 1),
        key=lambda key: table[key]["held_out_mse_pk"])
    payload = {
        "schema": "MetaSieve.StageA.InnerLrSelection.v1",
        "date": "2026-08-17",
        "selection_population": "meta_train only, component folds",
        "meta_val_read": False,
        "reference_checkpoint": str(A0.relative_to(ROOT)),
        "reference_kind": kind, "reference_seed": seed,
        "tasks": len(tasks), "components": len(components), "folds": FOLDS,
        "support": SUPPORT, "held_out_ligands": HELD,
        "baseline_no_adaptation_mse_pk": baseline_mse,
        "sweep": table,
        "best_overall": best_key,
        "selected_for_training": table[best_at_one_step],
        "selection_rule": (
            "lowest held-out MSE at the preregistered inner_steps=1; the "
            "multi-step rows are recorded for the evaluation-time sweep and "
            "did not choose the training value"),
        "meta_test": data.seal_record(),
    }
    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nbest overall: {best_key}")
    print(f"selected for training (inner_steps=1): {best_at_one_step}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
