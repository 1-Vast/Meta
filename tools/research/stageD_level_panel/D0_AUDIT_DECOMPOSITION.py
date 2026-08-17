"""Re-audit the Stage C level/shape decomposition.

Questions from the governing task:
1. Is the level/shape split defined per episode, per target, and per draw?
2. Is the calibrated constant fitted on meta_train only?
3. Does 'target level' mix assay history / panel composition into the level term?
4-5. How much of the measured level^2 is within-target panel-sampling variance
   rather than cross-target calibration error?

Method: reuse the Stage C evaluation banks (meta_val, k=0, 16 queries, 2 draws
per target) and the leak-free Stage B T checkpoint. For every episode we
recompute the exact decomposition; then we aggregate per target across its
draws and compare the drawn-panel mean to the canonical target mean over all
the target's unique ligands. The difference between per-episode level^2 and
per-target canonical level^2 is the panel-composition part of 'level'.

No training. meta_test never constructed. meta_val read once.
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
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)
from tools.research.stageB_complementary.arms import (                # noqa: E402
    StageBAdaptation, predict,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
QUERY_SIZE = 16
DRAWS = 2
EVALUATION_SEED = 73101
OUT = Path(__file__).resolve().parent / "D0_AUDIT_DECOMPOSITION.json"


def canonical_target_level(data: QPSMPData, split: str) -> dict[str, float]:
    levels = {}
    for target, indices in data.tasks[split].items():
        seen, values = set(), []
        for index in indices:
            cell = data.cells[int(index)]
            if cell["ligand_id"] in seen:
                continue
            seen.add(cell["ligand_id"])
            values.append(cell["pK"])
        levels[target] = float(np.mean(values))
    return levels


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path,
                        default=ROOT / "report/meta_fewshot/stageB_complementary_20260817/T/checkpoint.pt")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    # Same nested-bank construction as Stage C FEASIBILITY.json:
    # sizes (0,1,2,3,5) with max_support=5, then the k=0 slice.
    specs = data.fixed_nested_episode_banks(
        "meta_val", (0, 1, 2, 3, 5), QUERY_SIZE, DRAWS, EVALUATION_SEED, None)[0]
    banks = tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in specs)

    blob = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    adaptation = StageBAdaptation.from_dict(blob["adaptation"])
    model, _, _ = load_arm(args.checkpoint, data, args.device)
    model.eval()

    scale, mean = label_scale.scale, label_scale.mean
    canonical = canonical_target_level(data, "meta_val")
    rows = []
    for episode in banks:
        spec = episode.spec
        parts = align_atoms(episode_tensors(model, episode, args.device, torch.float32))
        truth = parts["query_y"].squeeze(0).cpu().numpy() * scale + mean
        with torch.no_grad():
            task = encode_parts(model, parts)
            prediction = predict(model, parts, task, adaptation)["prediction"]
        values = prediction.squeeze(0).cpu().numpy() * scale + mean
        level = (values.mean() - truth.mean()) ** 2
        centered = (((values - values.mean()) - (truth - truth.mean())) ** 2).mean()
        rows.append({
            "target": spec.target, "component": spec.component,
            "mse": float(((values - truth) ** 2).mean()),
            "level_squared": float(level), "centered_mse": float(centered),
            "prediction_mean": float(values.mean()), "truth_mean": float(truth.mean()),
            "canonical_target_mean": canonical[spec.target],
        })

    targets = sorted({r["target"] for r in rows})
    by_target = {t: [r for r in rows if r["target"] == t] for t in targets}
    per_target = {}
    for t, rs in by_target.items():
        preds = np.asarray([r["prediction_mean"] for r in rs])
        truths = np.asarray([r["truth_mean"] for r in rs])
        per_target[t] = {
            "component": rs[0]["component"],
            "draw_level_squared_mean": float(np.mean([r["level_squared"] for r in rs])),
            "level_squared_of_target_means": float((preds.mean() - truths.mean()) ** 2),
            "truth_level_error_vs_canonical": float((truths.mean() - rs[0]["canonical_target_mean"]) ** 2),
            "draw_truth_spread": float(np.var(truths)),
            "draw_prediction_spread": float(np.var(preds)),
        }

    ep_level = float(np.mean([r["level_squared"] for r in rows]))
    ep_centered = float(np.mean([r["centered_mse"] for r in rows]))
    ep_mse = float(np.mean([r["mse"] for r in rows]))
    tgt_level = float(np.mean([v["level_squared_of_target_means"] for v in per_target.values()]))
    panel_noise_truth = float(np.mean([v["draw_truth_spread"] for v in per_target.values()]))
    canonical_err = float(np.mean([v["truth_level_error_vs_canonical"] for v in per_target.values()]))

    payload = {
        "schema": "MetaSieve.StageD.AuditDecomposition.v1",
        "date": "2026-08-17",
        "checkpoint": str(args.checkpoint.resolve().relative_to(ROOT)),
        "population": {"split": "meta_val", "targets": len(targets),
                       "episodes": len(rows),
                       "components": len({r["component"] for r in rows})},
        "per_episode": {"mse": ep_mse, "level_squared": ep_level,
                        "centered_mse": ep_centered,
                        "level_share": ep_level / ep_mse},
        "per_target": {
            "level_squared_of_target_means": tgt_level,
            "mean_draw_level_squared": float(np.mean([
                v["draw_level_squared_mean"] for v in per_target.values()])),
            "panel_sampling_variance_of_truth_means": panel_noise_truth,
            "draw_prediction_mean_spread": float(np.mean([
                v["draw_prediction_spread"] for v in per_target.values()])),
            "truth_mean_error_vs_canonical_target_mean": canonical_err,
        },
        "audit_answers": {
            "q1_decomposition_granularity": (
                "Stage C (FEASIBILITY.json) decomposes PER EPISODE (one draw of "
                "a 16-query panel), then aggregates with equal weight per "
                "target and per component. It is therefore NOT the canonical "
                "cross-target level error: the drawn-panel mean differs from "
                "the target's own mean, and that sampling variance is inside "
                "the level term."),
            "q2_calibrated_constant": (
                "The calibrated constant (LEVEL_CEILING.json, "
                "calibrated_constant_REFERENCE) is the meta_val target-level "
                "mean: it READS meta_val labels. It is disclosed as a "
                "REFERENCE, not a legitimate predictor. The meta_train-only "
                "constant is global_mean_meta_train (2.1703), which the ESM "
                "MLP probe (1.6357) and the incumbent (1.7078) both beat."),
            "q3_4_5_panel_mixing": (
                "The level term is a joint property of protein, assay and the "
                "drawn query panel: it equals (drawn-panel prediction mean - "
                "drawn-panel truth mean)^2. The panel-sampling variance of the "
                "truth means and the drawn-vs-canonical target-mean error are "
                "measured below; if they are large, part of 'target level' is "
                "panel composition, not protein."),
        },
        "meta_test": data.seal_record(),
    }
    print(json.dumps({k: v for k, v in payload.items() if k != "audit_answers"},
                     indent=1))
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
