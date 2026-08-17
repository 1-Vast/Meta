"""Stage H0: can structure-derived pocket priors predict a cold target level?

No training of the DTA model. For targets with a homologous holo structure
(>=30% identity, >=50% query coverage, local pilot20k corpus), pocket
descriptors (pocket size, volume, amino-acid and element composition, holo
ligand size, mapping identity/coverage) are tested as level predictors with
SGD MLP/linear probes: weight decay selected on meta_train COMPONENT folds,
meta_val read once, plus a shuffled-pocket control. External data lane;
meta_test never constructed.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    normalized_episode, training_label_scale,
)
from scripts.stageR0_retrieval_falsification import component_bootstrap

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
BASE = Path(__file__).resolve().parent
QUERY_SIZE = 16
DRAWS = 2
EVALUATION_SEED = 73101
FOLDS = 5
FOLD_SEED = 20260818
PROBE_STEPS = 600
PROBE_LR = 3e-3
DECAYS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
OUT = BASE / "H0_POCKET_IDENTIFIABILITY.json"


def load_pockets():
    with np.load(BASE / "pocket_descriptors.npz", allow_pickle=False) as store:
        return {str(k): store["vectors"][i].astype(np.float32)
                for i, k in enumerate(store["keys"])}


def main():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    pockets = load_pockets()

    def esm_pooled(target):
        return np.asarray(data.protein_bank.get(target)[0], dtype=np.float32)

    def build(split):
        specs = data.fixed_nested_episode_banks(
            split, (0, 1, 2, 3, 5), QUERY_SIZE, DRAWS, EVALUATION_SEED, None)[0]
        rows = []
        for spec in specs:
            if spec.target not in pockets:
                continue
            episode = normalized_episode(data.materialize(spec), label_scale)
            truth = (episode.query_y.numpy() * label_scale.scale
                     + label_scale.mean)
            rows.append({
                "component": component_of[spec.target], "target": spec.target,
                "pocket": pockets[spec.target],
                "protein": esm_pooled(spec.target),
                "truth_mean": float(truth.mean()),
            })
        return rows

    train_rows = build("meta_train")
    val_rows = build("meta_val")
    print("episodes with pocket descriptors: train", len(train_rows),
          "val", len(val_rows))

    def matrix(rows, key):
        return np.stack([r[key] for r in rows]).astype(np.float32)

    groups = {
        "pocket": (matrix(train_rows, "pocket"), matrix(val_rows, "pocket")),
        "protein_esm150": (matrix(train_rows, "protein"),
                           matrix(val_rows, "protein")),
        "pocket_plus_protein": (
            np.concatenate([matrix(train_rows, "pocket"),
                            matrix(train_rows, "protein")], 1),
            np.concatenate([matrix(val_rows, "pocket"),
                            matrix(val_rows, "protein")], 1)),
    }
    device = "cuda" if torch.cuda.is_available() else "cpu"

    def train_probe(x, y, hidden, decay, seed):
        torch.manual_seed(seed)
        probe = (nn.Sequential(nn.Linear(x.shape[1], hidden), nn.GELU(),
                               nn.Linear(hidden, 1))
                 if hidden else nn.Linear(x.shape[1], 1)).to(device)
        optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                      weight_decay=decay)
        xt = torch.as_tensor(x, dtype=torch.float32, device=device)
        yt = torch.as_tensor(y, dtype=torch.float32, device=device)
        for _ in range(PROBE_STEPS):
            loss = nn.functional.mse_loss(probe(xt).squeeze(-1), yt)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        return probe

    def predict(probe, x):
        with torch.no_grad():
            return probe(torch.as_tensor(x, dtype=torch.float32,
                                         device=device)).squeeze(-1).cpu().numpy()

    grand = float(np.mean([r["truth_mean"] for r in train_rows]))
    val_values = np.asarray([r["truth_mean"] for r in val_rows])
    components = sorted({r["component"] for r in train_rows})
    order = np.random.default_rng(FOLD_SEED).permutation(len(components))
    fold_of = {components[int(order[i])]: i % FOLDS for i in range(len(components))}
    folds = np.asarray([fold_of[r["component"]] for r in train_rows])

    payload = {"schema": "MetaSieve.StageH.PocketIdentifiability.v1",
               "date": "2026-08-17", "grand_mean_pk": grand,
               "methods": {}, "meta_test": data.seal_record()}
    results = {}
    for name, (train_x, val_x) in groups.items():
        train_y = np.asarray([r["truth_mean"] for r in train_rows]) - grand
        table = {}
        for decay in DECAYS:
            errs = []
            for fold in range(FOLDS):
                keep, held = folds != fold, folds == fold
                if not held.any():
                    continue
                probe = train_probe(train_x[keep], train_y[keep], 64, decay,
                                    seed=fold)
                errs.append(float(((predict(probe, train_x[held])
                                    - train_y[held]) ** 2).mean()))
            table[str(decay)] = float(np.mean(errs)) if errs else float("nan")
        best = min(table, key=table.get)
        probe = train_probe(train_x, train_y, 64, float(best), seed=0)
        pred = predict(probe, val_x) + grand
        mse = float(((pred - val_values) ** 2).mean())
        lin = train_probe(train_x, train_y, 0, float(best), seed=0)
        lin_mse = float(((predict(lin, val_x) + grand - val_values) ** 2).mean())
        control_mse = None
        if name == "pocket":
            rng = np.random.default_rng(0)
            shuffled = train_x[rng.permutation(len(train_x))]
            ctrl = train_probe(shuffled, train_y, 64, float(best), seed=0)
            control_mse = float(((predict(ctrl, val_x) + grand
                                  - val_values) ** 2).mean())
        pairs = [(r["component"], r["target"],
                  float((pred[i] - val_values[i]) ** 2))
                 for i, r in enumerate(val_rows)]
        interval = component_bootstrap(pairs, 9999, 20260816)
        results[name] = {"selected_weight_decay": float(best),
                         "fold_table": table, "level_mse": mse,
                         "level_rmse": float(np.sqrt(mse)),
                         "linear_level_mse": lin_mse,
                         "shuffled_pocket_control_mse": control_mse,
                         "component_bootstrap": interval}
        print(name, "decay", best, "mse", round(mse, 4),
              "linear", round(lin_mse, 4), "shuffled",
              round(control_mse, 4) if control_mse else None)
    baseline = float(((np.full_like(val_values, grand) - val_values) ** 2).mean())
    payload["grand_mean_mse_on_covered"] = baseline
    payload["methods"] = results
    payload["verdict"] = {"best_method": min(results, key=lambda k: results[k]["level_mse"]),
                          "best_mse": min(r["level_mse"] for r in results.values())}
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print("baseline", round(baseline, 4), "wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
