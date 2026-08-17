"""How well can an unseen protein's mean affinity be predicted?

`feasibility.py` shows the entire path to MSE <= 1.00 runs through target-level
calibration: with a perfect level predictor every k lands below the target, and
the model's within-target ordering is already no better than predicting the
target's own mean. So the research question reduces to one regression problem —
predict `mean_ligands(pK)` for a protein from an unseen homology component.

This script bounds that problem without training the full model:

* `global_mean`     — predict the meta_train grand mean for every target. The
                      no-information baseline, and the number any protein
                      feature must beat.
* `incumbent`       — the trained model's own k=0 target-mean prediction.
* `esm_probe`       — a small SGD probe on the pooled ESM protein embedding,
                      weight decay chosen on meta_train component folds, read
                      once on meta_val. No ridge, no closed form.
* `esm_knn`         — a parameter-free nearest-neighbour level transfer in the
                      same embedding, as a capacity-free reference.
* `oracle`          — the true target mean; error zero by construction, shown
                      only to anchor the scale.

Selection uses meta_train components only. `meta_val` is read once. `meta_test`
is never constructed.
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

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.stageR0_retrieval_falsification import (                 # noqa: E402
    component_bootstrap,
)
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
FOLDS = 5
FOLD_SEED = 20260818
PROBE_STEPS = 600
PROBE_LR = 3e-3
WEIGHT_DECAYS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816
OUT = Path(__file__).resolve().parent / "LEVEL_CEILING.json"


def target_levels(data: QPSMPData, split: str) -> dict[str, float]:
    """Mean pK over each target's unique ligands."""
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


def protein_matrix(data: QPSMPData, targets: list[str]) -> np.ndarray:
    return np.stack([np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
                     for t in targets])


def train_probe(features, values, decay, seed, steps=PROBE_STEPS):
    torch.manual_seed(seed)
    probe = nn.Linear(features.shape[1], 1)
    nn.init.zeros_(probe.weight)
    nn.init.zeros_(probe.bias)
    optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                  weight_decay=decay)
    x = torch.as_tensor(features, dtype=torch.float32)
    y = torch.as_tensor(values, dtype=torch.float32)
    for _ in range(steps):
        loss = nn.functional.mse_loss(probe(x).squeeze(-1), y)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
    return probe


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    train_levels = target_levels(data, "meta_train")
    val_levels = target_levels(data, "meta_val")
    train_targets = sorted(train_levels)
    val_targets = sorted(val_levels)

    grand_mean = float(np.mean([train_levels[t] for t in train_targets]))
    val_values = np.asarray([val_levels[t] for t in val_targets])
    between_sd_train = float(np.std([train_levels[t] for t in train_targets]))
    between_sd_val = float(np.std(val_values))

    print(f"meta_train targets {len(train_targets)}, meta_val {len(val_targets)}")
    print(f"target-level grand mean (meta_train) = {grand_mean:.4f} pK")
    print(f"between-target SD: meta_train {between_sd_train:.4f}, "
          f"meta_val {between_sd_val:.4f} pK")

    train_matrix = protein_matrix(data, train_targets)
    val_matrix = protein_matrix(data, val_targets)
    centre = train_matrix.mean(0, keepdims=True)
    scale = train_matrix.std(0, keepdims=True)
    scale[scale < 1e-8] = 1.0
    train_x = (train_matrix - centre) / scale
    val_x = (val_matrix - centre) / scale
    train_y = np.asarray([train_levels[t] for t in train_targets]) - grand_mean

    # --- weight decay on meta_train component folds only -------------------
    components = sorted({component_of[t] for t in train_targets})
    order = np.random.default_rng(FOLD_SEED).permutation(len(components))
    fold_of = {components[int(i)]: rank % FOLDS for rank, i in enumerate(order)}
    fold_index = np.asarray([fold_of[component_of[t]] for t in train_targets])
    table = {}
    for decay in WEIGHT_DECAYS:
        errors = []
        for fold in range(FOLDS):
            keep = fold_index != fold
            probe = train_probe(train_x[keep], train_y[keep], decay, seed=fold)
            with torch.no_grad():
                held = probe(torch.as_tensor(train_x[~keep],
                                             dtype=torch.float32)).squeeze(-1).numpy()
            errors.append(float(((held - train_y[~keep]) ** 2).mean()))
        table[f"{decay:g}"] = float(np.mean(errors))
        print(f"  decay {decay:<6g} fold MSE {table[f'{decay:g}']:.4f}")
    best_decay = float(min(table, key=lambda key: table[key]))
    print(f"selected weight decay {best_decay:g} (meta_train folds only)")

    probe = train_probe(train_x, train_y, best_decay, seed=0)
    with torch.no_grad():
        probe_prediction = probe(torch.as_tensor(
            val_x, dtype=torch.float32)).squeeze(-1).numpy() + grand_mean

    # --- parameter-free nearest neighbour in the same space ----------------
    normed_train = train_x / np.maximum(
        np.linalg.norm(train_x, axis=1, keepdims=True), 1e-9)
    normed_val = val_x / np.maximum(
        np.linalg.norm(val_x, axis=1, keepdims=True), 1e-9)
    similarity = normed_val @ normed_train.T
    nearest = similarity.argmax(1)
    knn_prediction = np.asarray([train_levels[train_targets[i]] for i in nearest])

    # --- the incumbent's own level -----------------------------------------
    label_scale = training_label_scale(data)
    specs = data.fixed_nested_episode_banks("meta_val", (0,), 16, 2, 73101, None)[0]
    banks = [compact_episode(normalized_episode(data.materialize(s), label_scale))
             for s in specs]
    blob = torch.load(arguments.checkpoint, map_location="cpu", weights_only=False)
    adaptation = StageBAdaptation.from_dict(blob["adaptation"])
    model, _, _ = load_arm(arguments.checkpoint, data, arguments.device)
    model.eval()
    incumbent: dict[str, list[float]] = {}
    for episode in banks:
        parts = align_atoms(episode_tensors(model, episode, arguments.device,
                                            torch.float32))
        with torch.no_grad():
            task = encode_parts(model, parts)
            values = predict(model, parts, task, adaptation)["prediction"]
        pk = values.squeeze(0).cpu().numpy() * label_scale.scale + label_scale.mean
        incumbent.setdefault(episode.spec.target, []).append(float(pk.mean()))
    incumbent_prediction = np.asarray(
        [float(np.mean(incumbent[t])) if t in incumbent else grand_mean
         for t in val_targets])

    # A non-linear probe, so the boundary is not an artifact of linearity.
    torch.manual_seed(1)
    mlp = nn.Sequential(nn.Linear(train_x.shape[1], 64), nn.GELU(),
                        nn.Linear(64, 1))
    optimizer = torch.optim.AdamW(mlp.parameters(), lr=PROBE_LR,
                                  weight_decay=best_decay)
    xt = torch.as_tensor(train_x, dtype=torch.float32)
    yt = torch.as_tensor(train_y, dtype=torch.float32)
    for _ in range(PROBE_STEPS):
        loss = nn.functional.mse_loss(mlp(xt).squeeze(-1), yt)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
    with torch.no_grad():
        mlp_prediction = mlp(torch.as_tensor(
            val_x, dtype=torch.float32)).squeeze(-1).numpy() + grand_mean

    # Sequence length is a crude but genuinely protein-derived covariate; if
    # even it beat ESM, the embedding would be the problem rather than the task.
    lengths_train = np.asarray([[len(data._protein_sequences[t])]
                                for t in train_targets], dtype=np.float32)
    lengths_val = np.asarray([[len(data._protein_sequences[t])]
                              for t in val_targets], dtype=np.float32)
    lc, ls = lengths_train.mean(0, keepdims=True), lengths_train.std(0, keepdims=True)
    ls[ls < 1e-8] = 1.0
    length_probe = train_probe((lengths_train - lc) / ls, train_y, best_decay, 2)
    with torch.no_grad():
        length_prediction = length_probe(torch.as_tensor(
            (lengths_val - lc) / ls, dtype=torch.float32)).squeeze(-1).numpy() + grand_mean

    methods = {
        "global_mean_meta_train": np.full_like(val_values, grand_mean),
        "esm_knn": knn_prediction,
        "esm_probe_linear": probe_prediction,
        "esm_probe_mlp": mlp_prediction,
        "sequence_length_probe": length_prediction,
        "incumbent_model": incumbent_prediction,
        # REFERENCE, not a model: it reads meta_val labels. It is the error a
        # perfectly calibrated but completely uninformative constant would make,
        # and therefore the line every protein feature has to beat before it can
        # be said to carry any target-level information at all.
        "calibrated_constant_REFERENCE": np.full_like(
            val_values, float(val_values.mean())),
        "oracle": val_values,
    }
    payload = {
        "schema": "MetaSieve.StageC.LevelCeiling.v1", "date": "2026-08-17",
        "question": "how well can an unseen protein's mean affinity be predicted?",
        "selection": "weight decay on meta_train component folds; meta_val read once",
        "grand_mean_pk": grand_mean,
        "between_target_sd": {"meta_train": between_sd_train,
                              "meta_val": between_sd_val},
        "weight_decay_folds": table, "selected_weight_decay": best_decay,
        "targets": {"meta_train": len(train_targets), "meta_val": len(val_targets)},
        "methods": {}, "meta_test": data.seal_record(),
    }
    print(f"\n{'method':<30} {'level MSE':>10} {'level RMSE':>11} "
          f"{'vs global':>10} {'vs calib':>9}")
    baseline = float(((methods["global_mean_meta_train"] - val_values) ** 2).mean())
    calibrated = float(((methods["calibrated_constant_REFERENCE"]
                         - val_values) ** 2).mean())
    for name, prediction in methods.items():
        mse = float(((prediction - val_values) ** 2).mean())
        pairs = [(component_of[t], t, float((prediction[i] - val_values[i]) ** 2))
                 for i, t in enumerate(val_targets)]
        payload["methods"][name] = {
            "level_mse": mse, "level_rmse": float(np.sqrt(mse)),
            "relative_to_global_mean": mse / baseline if baseline else float("nan"),
            "relative_to_calibrated_constant": (
                mse / calibrated if calibrated else float("nan")),
            "beats_calibrated_constant": bool(mse < calibrated),
            "component_bootstrap": component_bootstrap(
                pairs, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED),
        }
        print(f"{name:<30} {mse:>10.4f} {np.sqrt(mse):>11.4f} "
              f"{mse / baseline if baseline else float('nan'):>10.3f} "
              f"{mse / calibrated if calibrated else float('nan'):>9.3f}")

    # What each level predictor would give at k=0, holding shape fixed.
    shape = 0.8761   # measured k=0 centered MSE, FEASIBILITY.json
    payload["implied_k0_mse_holding_shape_fixed"] = {
        name: block["level_mse"] + shape
        for name, block in payload["methods"].items()}
    payload["shape_term_used"] = shape
    print(f"\nimplied k=0 MSE (level MSE + measured shape {shape:.4f}):")
    for name, value in payload["implied_k0_mse_holding_shape_fixed"].items():
        flag = "<= 1.00 TARGET MET" if value <= 1.0 else ""
        print(f"  {name:<18} {value:>8.4f}  {flag}")

    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
