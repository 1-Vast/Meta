"""Stage L2: a fair directional-SAR probe. Supersedes `stage_l_ligand_sar.py`.

Four defects in the first version, each of which flattered the learned arm:

1. **Unbalanced pairs.** Pairs were sampled with independent draws of `i` and
   `j`, so the orientation distribution was arbitrary. A symmetric score such
   as Tanimoto then has an *incidental*, non-zero correlation with the signed
   gap. Balanced construction — every unordered pair contributes both `(i,j)`
   and `(j,i)` — makes the signed target exactly antisymmetric, so any
   symmetric predictor has **identically zero** signed correlation by
   construction rather than by measurement.
2. **A misleading comparison.** Reporting "Tanimoto scores +0.028 on signed
   Δy" invited the reading that it performs poorly. It cannot perform at all:
   symmetry forbids it. Under the balanced construction this is a structural
   zero and is reported as such. Tanimoto is scored where it is meaningful —
   `|Δy|` and applicability — and nowhere else.
3. **An algebraically collapsing head.** `⟨w, U Δe⟩` with `U: D→R` and
   `w ∈ R^R` is exactly `⟨Uᵀw, Δe⟩` — a single linear functional. The
   rank-4/8/16 "selection" was selecting over a reparameterisation of one
   vector. Replaced by an explicit `Linear(D, 1, bias=False)`, which is the
   same hypothesis class stated honestly.
4. **Residualisation fitted on the evaluation split.** The incremental-value
   slope was computed on `meta_val` using `meta_val` labels. Now fitted on
   `meta_train` and applied frozen.

Also removed: the claim that Tanimoto "points in the wrong direction" on
activity cliffs. A symmetric similarity has no direction; the earlier −0.37 was
an artifact of the unbalanced sampling.

Controls are capacity-matched — every arm is one linear functional over a
difference of the same width, so a win cannot come from parameter count.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import stable_seed                          # noqa: E402
from scripts.stageR0_retrieval_falsification import (               # noqa: E402
    component_bootstrap, component_target_mean,
)
from tools.research.a2_readiness_v2 import _frozen                  # noqa: E402

STEPS = 600
LEARNING_RATE = 3e-3
PAIRS_PER_TARGET = 48          # unordered; doubled by the balanced construction
SAMPLING_SEEDS = (20260901, 20260902, 20260903)
FOLDS = 5
WEIGHT_DECAYS = (0.0, 1e-4, 1e-2)     # the only hyperparameter, on meta_train
CLIFF_SIMILARITY = 0.6
CLIFF_GAP_PK = 1.0
MIN_STRATUM_TARGETS = 8
MIN_STRATUM_COMPONENTS = 5


class DirectionalHead(nn.Module):
    """`Δŷ = ⟨w, Δe⟩`. One linear functional, stated without a fake rank."""

    def __init__(self, width: int, freeze: bool = False):
        super().__init__()
        self.weight = nn.Linear(width, 1, bias=False)
        if freeze:
            self.weight.weight.requires_grad_(False)

    def forward(self, delta: torch.Tensor) -> torch.Tensor:
        return self.weight(delta).squeeze(-1)


def balanced_pairs(blocks, seed: int, per_target: int = PAIRS_PER_TARGET):
    """Unique unordered pairs per target, then both orientations of each.

    Enumerates exhaustively when a target has few enough ligands, otherwise
    samples unordered pairs without replacement. Either way each unordered pair
    appears exactly twice, once per orientation, so `Δy` is antisymmetric and
    balanced.
    """
    left, right, owner = [], [], []
    for target, block in enumerate(blocks):
        if len(block) < 2:
            continue
        rng = np.random.default_rng(stable_seed("stage-l2", seed, target) % (2 ** 32))
        rows, columns = np.triu_indices(len(block), 1)
        total = len(rows)
        if total > per_target:
            picks = rng.choice(total, size=per_target, replace=False)
            rows, columns = rows[picks], columns[picks]
        a, b = block[rows], block[columns]
        left.append(np.concatenate([a, b]))
        right.append(np.concatenate([b, a]))
        owner.append(np.full(2 * len(a), target, dtype=np.int64))
    return (np.concatenate(left), np.concatenate(right), np.concatenate(owner))


def tanimoto_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    intersection = (left * right).sum(-1)
    union = left.sum(-1) + right.sum(-1) - intersection
    return intersection / np.maximum(union, 1e-9)


def correlation(prediction: np.ndarray, truth: np.ndarray) -> float:
    p, t = prediction - prediction.mean(), truth - truth.mean()
    denominator = float(np.sqrt((p ** 2).mean()) * np.sqrt((t ** 2).mean()))
    return float((p * t).mean() / denominator) if denominator > 1e-12 else 0.0


def per_target(prediction, truth, owner, components) -> list:
    rows = []
    for target in np.unique(owner):
        select = owner == target
        p, t = prediction[select], truth[select]
        if len(p) < 4 or p.std() < 1e-9 or t.std() < 1e-9:
            continue
        rows.append((int(components[target]), int(target), correlation(p, t)))
    return rows


def summarise(rows) -> dict:
    return {"delta_r": component_target_mean(rows),
            "delta_r_ci": component_bootstrap(
                rows, _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED),
            "targets": len(rows),
            "components": len({component for component, _, _ in rows})}


def fit(width: int, features: np.ndarray, delta: np.ndarray, device: str,
        freeze: bool, seed: int, weight_decay: float) -> DirectionalHead:
    torch.manual_seed(seed)
    head = DirectionalHead(width, freeze).to(device)
    trainable = [p for p in head.parameters() if p.requires_grad]
    if not trainable:
        return head.eval()
    optimiser = torch.optim.Adam(trainable, lr=LEARNING_RATE,
                                 weight_decay=weight_decay)
    x = torch.as_tensor(features, dtype=torch.float32, device=device)
    y = torch.as_tensor(delta, dtype=torch.float32, device=device)
    for _ in range(STEPS):
        optimiser.zero_grad(set_to_none=True)
        (head(x) - y).square().mean().backward()
        optimiser.step()
    return head.eval()


def predict(head: DirectionalHead, features: np.ndarray, device: str) -> np.ndarray:
    with torch.no_grad():
        return head(torch.as_tensor(features, dtype=torch.float32,
                                    device=device)).cpu().numpy()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    def load(name):
        with np.load(arguments.features / f"{name}.npz") as stored:
            data = {key: stored[key] for key in stored.files}
        owner = data["target_index"]
        data["_blocks"] = [np.flatnonzero(owner == i)
                           for i in range(int(owner.max()) + 1)]
        return data

    train, val = load("meta_train"), load("meta_val")
    noise = np.random.default_rng(_frozen.CONTROL_SEED + 7)
    for block in (train, val):
        block["random__feature"] = noise.normal(
            size=block["correct__embed"].shape).astype(np.float32)

    # Difference representations. Every arm is the same hypothesis class — one
    # linear functional over a difference — so capacity is matched by design.
    ARMS = {
        "embed":            "correct__embed",
        "ligand_encoder":   "correct__ligand",
        "morgan_difference": "fingerprint",
        "random_feature":   "random__feature",
    }

    payload: dict = {
        "schema": "MetaSieve.A2ExactProbe.StageL2.v1",
        "split": "meta_val",
        "frozen_design": _frozen.frozen_manifest(),
        "construction": {
            "pairs": "balanced: every unordered pair contributes (i,j) and (j,i)",
            "consequence": ("the signed target is exactly antisymmetric, so any "
                            "symmetric score has identically zero signed "
                            "correlation by construction"),
            "head": "Linear(D, 1, bias=False) — one linear functional, no "
                    "algebraically redundant rank factorisation",
            "hyperparameter": "weight decay only, selected on meta_train "
                              "component folds",
            "residualisation": "slope fitted on meta_train, applied frozen to "
                               "meta_val",
            "sampling_seeds": list(SAMPLING_SEEDS),
        },
        "seeds": {},
    }

    for sampling_seed in SAMPLING_SEEDS:
        train_left, train_right, train_owner = balanced_pairs(
            train["_blocks"], sampling_seed)
        val_left, val_right, val_owner = balanced_pairs(
            val["_blocks"], sampling_seed)
        train_delta = train["y"][train_left] - train["y"][train_right]
        val_delta = val["y"][val_left] - val["y"][val_right]
        components = train["component_of_target"]
        fold_of_pair = np.asarray([int(components[t]) % FOLDS for t in train_owner])

        block: dict = {}
        predictions: dict = {}
        for arm, key in ARMS.items():
            raw_train, raw_val = train[key], val[key]
            spread = raw_train.std(0, keepdims=True) + 1e-6
            x_train = (raw_train[train_left] - raw_train[train_right]) / spread
            x_val = (raw_val[val_left] - raw_val[val_right]) / spread
            width = raw_train.shape[1]

            best, best_score = WEIGHT_DECAYS[0], -np.inf
            for decay in WEIGHT_DECAYS:
                held = []
                for fold in range(FOLDS):
                    inside = fold_of_pair != fold
                    if inside.sum() < 200 or (~inside).sum() < 200:
                        continue
                    head = fit(width, x_train[inside], train_delta[inside],
                               arguments.device, False,
                               _frozen.CONTROL_SEED + fold, decay)
                    held.append(component_target_mean(per_target(
                        predict(head, x_train[~inside], arguments.device),
                        train_delta[~inside], train_owner[~inside], components)))
                value = float(np.mean(held)) if held else -np.inf
                if value > best_score:
                    best, best_score = decay, value

            head = fit(width, x_train, train_delta, arguments.device, False,
                       _frozen.CONTROL_SEED, best)
            predictions[arm] = predict(head, x_val, arguments.device)
            block[arm] = dict(summarise(per_target(
                predictions[arm], val_delta, val_owner,
                val["component_of_target"])),
                weight_decay=best, meta_train_fold_score=best_score,
                trainable_parameters=width)

        # Shuffled labels: refit `embed` on within-target permuted pK.
        shuffle = np.random.default_rng(sampling_seed + 11)
        shuffled = train["y"].copy()
        for rows in train["_blocks"]:
            shuffled[rows] = shuffled[shuffle.permutation(rows)]
        raw_train, raw_val = train["correct__embed"], val["correct__embed"]
        spread = raw_train.std(0, keepdims=True) + 1e-6
        head = fit(raw_train.shape[1],
                   (raw_train[train_left] - raw_train[train_right]) / spread,
                   shuffled[train_left] - shuffled[train_right],
                   arguments.device, False, _frozen.CONTROL_SEED,
                   block["embed"]["weight_decay"])
        predictions["shuffled_labels"] = predict(
            head, (raw_val[val_left] - raw_val[val_right]) / spread,
            arguments.device)
        block["shuffled_labels"] = summarise(per_target(
            predictions["shuffled_labels"], val_delta, val_owner,
            val["component_of_target"]))

        # Tanimoto: symmetric. Reported on the signed target only to
        # demonstrate the structural zero, and on |Δy| where it is meaningful.
        similarity = tanimoto_rows(val["fingerprint"][val_left],
                                   val["fingerprint"][val_right])
        block["tanimoto_signed_structural_zero"] = summarise(per_target(
            -similarity, val_delta, val_owner, val["component_of_target"]))
        block["tanimoto_on_absolute_gap"] = summarise(per_target(
            -similarity, np.abs(val_delta), val_owner,
            val["component_of_target"]))
        block["embed_on_absolute_gap"] = summarise(per_target(
            np.abs(predictions["embed"]), np.abs(val_delta), val_owner,
            val["component_of_target"]))

        # Incremental value, with the slope fitted on meta_train and frozen.
        # Tanimoto is symmetric so it cannot carry signed information; the
        # meaningful incremental question is whether `embed` adds to the
        # protein-blind ligand encoder.
        raw_train_l, raw_val_l = train["correct__ligand"], val["correct__ligand"]
        spread_l = raw_train_l.std(0, keepdims=True) + 1e-6
        ligand_head = fit(
            raw_train_l.shape[1],
            (raw_train_l[train_left] - raw_train_l[train_right]) / spread_l,
            train_delta, arguments.device, False, _frozen.CONTROL_SEED,
            block["ligand_encoder"]["weight_decay"])
        train_ligand = predict(
            ligand_head,
            (raw_train_l[train_left] - raw_train_l[train_right]) / spread_l,
            arguments.device)
        train_embed = predict(
            fit(raw_train.shape[1],
                (raw_train[train_left] - raw_train[train_right]) / spread,
                train_delta, arguments.device, False, _frozen.CONTROL_SEED,
                block["embed"]["weight_decay"]),
            (raw_train[train_left] - raw_train[train_right]) / spread,
            arguments.device)

        def frozen_slope(target: np.ndarray, against: np.ndarray) -> float:
            centred = against - against.mean()
            denominator = float((centred ** 2).sum())
            return (float(((target - target.mean()) * centred).sum() / denominator)
                    if denominator > 1e-12 else 0.0)

        slope_delta = frozen_slope(train_delta, train_ligand)
        slope_embed = frozen_slope(train_embed, train_ligand)
        residual_truth = val_delta - slope_delta * predictions["ligand_encoder"]
        residual_embed = (predictions["embed"]
                          - slope_embed * predictions["ligand_encoder"])
        block["embed_given_ligand_encoder"] = dict(summarise(per_target(
            residual_embed, residual_truth, val_owner,
            val["component_of_target"])),
            slopes_fitted_on="meta_train")
        block["prediction_correlation_embed_vs_ligand_encoder"] = correlation(
            predictions["embed"], predictions["ligand_encoder"])

        # Strata, with adequacy recorded rather than assumed.
        train_fp = train["fingerprint"]

        def max_similarity_to_train(rows: np.ndarray) -> np.ndarray:
            out = np.empty(len(rows), dtype=np.float32)
            chunk_rows = val["fingerprint"][rows]
            for start in range(0, len(rows), 64):
                chunk = chunk_rows[start:start + 64]
                intersection = chunk @ train_fp.T
                union = (chunk.sum(-1)[:, None] + train_fp.sum(-1)[None, :]
                         - intersection)
                out[start:start + 64] = (
                    intersection / np.maximum(union, 1e-9)).max(1)
            return out

        novelty = np.maximum(max_similarity_to_train(val_left),
                             max_similarity_to_train(val_right))
        low, high = np.quantile(novelty, [1 / 3, 2 / 3])
        strata = {
            "cliff_subset": (similarity >= CLIFF_SIMILARITY)
                            & (np.abs(val_delta) >= CLIFF_GAP_PK),
            "low_novelty": novelty >= high,
            "mid_novelty": (novelty >= low) & (novelty < high),
            "high_novelty": novelty < low,
        }
        block["strata"] = {}
        for label, mask in strata.items():
            rows = per_target(predictions["embed"][mask], val_delta[mask],
                              val_owner[mask], val["component_of_target"])
            cell = dict(summarise(rows), pairs=int(mask.sum()))
            cell["adequate"] = (cell["targets"] >= MIN_STRATUM_TARGETS
                                and cell["components"] >= MIN_STRATUM_COMPONENTS)
            cell["status"] = "confirmatory" if cell["adequate"] else "exploratory"
            block["strata"][label] = cell

        # Exact overlap accounting.
        block["overlap"] = {
            "val_ligands": len(set(val["ligand_ids"].tolist())),
            "shared_ligand_ids_with_meta_train": len(
                set(val["ligand_ids"].tolist())
                & set(train["ligand_ids"].tolist())),
            "mean_max_tanimoto_to_meta_train": float(novelty.mean()),
        }
        block["pairs"] = int(len(val_delta))
        payload["seeds"][str(sampling_seed)] = block

    # Across-seed summary of the headline arm.
    payload["across_sampling_seeds"] = {
        arm: {
            "delta_r_mean": float(np.mean([
                payload["seeds"][str(s)][arm]["delta_r"] for s in SAMPLING_SEEDS])),
            "delta_r_sd": float(np.std([
                payload["seeds"][str(s)][arm]["delta_r"] for s in SAMPLING_SEEDS])),
        }
        for arm in ("embed", "ligand_encoder", "morgan_difference",
                    "random_feature", "shuffled_labels")
    }

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")

    first = payload["seeds"][str(SAMPLING_SEEDS[0])]
    print(f"{'arm':<34}{'delta_r':>9}{'[95% CI]':>22}{'tgt':>5}{'cmp':>5}")
    for name, cell in first.items():
        if not isinstance(cell, dict) or "delta_r" not in cell:
            continue
        interval = cell["delta_r_ci"]
        print(f"{name:<34}{cell['delta_r']:>+9.4f}"
              f"   [{interval['lo']:+.4f},{interval['hi']:+.4f}]"
              f"{cell['targets']:>5}{cell['components']:>5}")
    print("\nacross sampling seeds (mean ± sd):")
    for arm, cell in payload["across_sampling_seeds"].items():
        print(f"  {arm:<24}{cell['delta_r_mean']:+.4f} ± {cell['delta_r_sd']:.4f}")
    print("\nstrata (embed):")
    for label, cell in first["strata"].items():
        print(f"  {label:<16}{cell['pairs']:>6} pairs {cell['targets']:>3} targets "
              f"{cell['components']:>3} components  {cell['delta_r']:+.4f}  "
              f"{cell['status']}")
    print(f"\noverlap: {json.dumps(first['overlap'])}")
    print(f"\nwrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
