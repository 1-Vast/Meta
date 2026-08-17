"""D0: which legal inputs can predict a cold target's level?

The Stage C ceiling tested protein features only (ESM-150M pooled, sequence
length, kNN). This diagnostic tests the two remaining legal families that Stage
C named but never measured:

* panel composition: statistics of the query ligand set itself (molecular
  weight, atom-feature fractions, ring/aromatic content, scaffold diversity,
  Morgan fingerprint means). The query ligands are model INPUTS at every k,
  including k=0, so their set-statistics are legal features. A BindingDB
  target's mean affinity depends on which ligands were tested against it, so
  panel composition is a candidate proxy for assay history.
* assay covariates: endpoint type (Ki vs other), distinct-document count,
  panel_count, replicate_count, parsed from each cell's legal metadata.

Plus a shuffled-panel control: the panel features of meta_train are permuted
across episodes before fitting, so any association the probe finds must come
from genuine panel-level signal rather than capacity.

Protocol: episodes from the fixed nested k=0 bank (16 queries, 2 draws/target),
target = episode truth mean (pK). Probes are SGD-trained MLPs/linear maps with
weight decay selected on meta_train COMPONENT folds; meta_val is read once.
No closed forms. meta_test never constructed.
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
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    normalized_episode, training_label_scale,
)
from scripts.stageR0_retrieval_falsification import component_bootstrap  # noqa: E402

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
QUERY_SIZE = 16
DRAWS = 2
EVALUATION_SEED = 73101
FOLDS = 5
FOLD_SEED = 20260818
PROBE_STEPS = 600
PROBE_LR = 3e-3
WEIGHT_DECAYS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816
OUT = Path(__file__).resolve().parent / "D0_LEVEL_IDENTIFIABILITY.json"

ATOM_FEAT_NAMES = ([f"elem{i}" for i in range(12)]
                   + [f"deg{i}" for i in range(6)]
                   + [f"charge{i}" for i in range(5)]
                   + [f"hyb{i}" for i in range(4)]
                   + ["aromatic", "in_ring"]
                   + [f"chiral{i}" for i in range(3)])


def panel_features(data: QPSMPData, episode, ligand_mw) -> tuple[np.ndarray, dict]:
    """Statistics of the query panel, all computable from model inputs."""
    cells = data.cells
    q_idx = list(episode.spec.query)
    feats = []
    for i in q_idx:
        cell = cells[int(i)]
        ligand = cell["ligand_id"]
        atoms, bonds, mask = data.ligand_bank.get(ligand)
        n = int(mask.sum())
        if n == 0:
            continue
        x = atoms[:n]
        feats.append(x.mean(0))
    feats = np.asarray(feats, dtype=np.float32)          # [Q, 32]
    fp_rows = np.asarray([data.fingerprints[cells[int(i)]["ligand_id"]].numpy()
                          for i in q_idx], dtype=np.float32)
    mw = np.asarray([ligand_mw[cells[int(i)]["ligand_id"]][0]
                     for i in q_idx], dtype=np.float32)
    scaffolds = [ligand_mw[cells[int(i)]["ligand_id"]][1] for i in q_idx]
    endpoints, dois, panel_counts, rep_counts = [], set(), [], []
    for i in q_idx:
        cell = cells[int(i)]
        endpoints.append(str(cell["panel_ids"][0]).split("|")[1]
                          if cell["panel_ids"] else "none")
        for pid in cell["panel_ids"]:
            dois.add(str(pid).split("|")[0])
        panel_counts.append(int(cell["panel_count"]))
        rep_counts.append(int(cell["replicate_count"]))
    stats = {
        "panel_size": float(len(q_idx)),
        "mw_mean": float(mw.mean()), "mw_std": float(mw.std()),
        "unique_scaffolds": float(len({s for s in scaffolds})),
        "unique_scaffold_fraction": float(len(set(scaffolds)) / max(len(scaffolds), 1)),
        "endpoint_ki_fraction": float(np.mean([e == "Ki" for e in endpoints])),
        "distinct_documents": float(len(dois)),
        "panel_count_mean": float(np.mean(panel_counts)),
        "replicate_count_mean": float(np.mean(rep_counts)),
    }
    atom_stats = np.concatenate([feats.mean(0), feats.std(0)])   # [64]
    fp_mean = fp_rows.mean(0)                                     # [1024]
    vector = np.concatenate([
        np.asarray([stats[k] for k in (
            "panel_size", "mw_mean", "mw_std", "unique_scaffolds",
            "unique_scaffold_fraction", "endpoint_ki_fraction",
            "distinct_documents", "panel_count_mean",
            "replicate_count_mean")], dtype=np.float32),
        atom_stats, fp_mean,
    ])
    return vector, stats


def protein_features(data: QPSMPData, episode) -> np.ndarray:
    pooled = data.protein_bank.get(episode.spec.target)[0]
    return np.asarray(pooled, dtype=np.float32)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--esm650", type=Path, default=None,
                        help="optional pooled ESM-650M feature file (.npz)")
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}

    esm650 = None
    if args.esm650 is not None:
        with np.load(args.esm650, allow_pickle=False) as store:
            esm650 = {str(k): store["pooled"][i].astype(np.float32)
                      for i, k in enumerate(store["keys"])}

    ligand_mw = {}
    for row in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
        ligand_mw[str(row["drug_key"])] = (float(row["molecular_weight"]),
                                           str(row.get("scaffold", "")))
    _ = data.fingerprints  # build the Morgan table once, up front

    def build(split):
        specs = data.fixed_nested_episode_banks(
            split, (0, 1, 2, 3, 5), QUERY_SIZE, DRAWS, EVALUATION_SEED, None)[0]
        rows = []
        for spec in specs:
            episode = normalized_episode(data.materialize(spec), label_scale)
            truth = (episode.query_y.numpy() * label_scale.scale
                     + label_scale.mean)
            pf, stats = panel_features(data, episode, ligand_mw)
            pr = protein_features(data, episode)
            p650 = (esm650[spec.target] if esm650 is not None
                    and spec.target in esm650 else None)
            rows.append({
                "component": component_of[spec.target], "target": spec.target,
                "panel": pf, "protein": pr, "protein650": p650,
                "stats": stats, "truth_mean": float(truth.mean()),
                "n": int(len(truth)),
            })
        return rows

    train_rows = build("meta_train")
    val_rows = build("meta_val")
    print(f"meta_train episodes {len(train_rows)}, meta_val {len(val_rows)}")

    def matrix(rows, key, dtype=np.float32):
        return np.stack([r[key] for r in rows]).astype(dtype)

    def standardize(train, apply):
        centre = train.mean(0, keepdims=True)
        scale = train.std(0, keepdims=True)
        scale[scale < 1e-8] = 1.0
        return (apply - centre) / scale

    groups = {
        "protein_esm150": ("protein", None),
        "panel": ("panel", None),
        "protein_plus_panel": ("protein_plus_panel", None),
    }
    for name, (key, _) in groups.items():
        if name == "protein_plus_panel":
            train_x = np.concatenate(
                [matrix(train_rows, "protein"), matrix(train_rows, "panel")], 1)
            val_x = np.concatenate(
                [matrix(val_rows, "protein"), matrix(val_rows, "panel")], 1)
        else:
            train_x, val_x = matrix(train_rows, key), matrix(val_rows, key)
        groups[name] = (train_x, val_x)
    if esm650 is not None and all(r["protein650"] is not None for r in train_rows + val_rows):
        groups["protein_esm650"] = (
            matrix(train_rows, "protein650"), matrix(val_rows, "protein650"))

    probe_device = "cuda" if torch.cuda.is_available() else "cpu"

    def train_probe(x, y, hidden, decay, seed):
        torch.manual_seed(seed)
        probe = (nn.Sequential(nn.Linear(x.shape[1], hidden), nn.GELU(),
                               nn.Linear(hidden, 1))
                 if hidden else nn.Linear(x.shape[1], 1)).to(probe_device)
        optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                      weight_decay=decay)
        xt = torch.as_tensor(x, dtype=torch.float32, device=probe_device)
        yt = torch.as_tensor(y, dtype=torch.float32, device=probe_device)
        for _ in range(PROBE_STEPS):
            loss = nn.functional.mse_loss(probe(xt).squeeze(-1), yt)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        return probe

    def predict(probe, x):
        with torch.no_grad():
            return probe(torch.as_tensor(x, dtype=torch.float32,
                                         device=probe_device)
                         ).squeeze(-1).cpu().numpy()

    grand_mean = float(np.mean([r["truth_mean"] for r in train_rows]))
    val_values = np.asarray([r["truth_mean"] for r in val_rows])

    components = sorted({r["component"] for r in train_rows})
    order = np.random.default_rng(FOLD_SEED).permutation(len(components))
    fold_of = {components[int(i)]: rank % FOLDS for rank, i in enumerate(order)}
    folds = np.asarray([fold_of[r["component"]] for r in train_rows])

    payload = {
        "schema": "MetaSieve.StageD.LevelIdentifiability.v1",
        "date": "2026-08-17",
        "grand_mean_pk": grand_mean,
        "methods": {},
        "meta_test": data.seal_record(),
    }
    results = {}
    for name, (train_x, val_x) in groups.items():
        train_y = np.asarray([r["truth_mean"] for r in train_rows]) - grand_mean
        table = {}
        for decay in WEIGHT_DECAYS:
            errors = []
            for fold in range(FOLDS):
                keep, held = folds != fold, folds == fold
                if not held.any():
                    continue
                probe = train_probe(train_x[keep], train_y[keep], 64, decay,
                                    seed=fold)
                errors.append(float(((predict(probe, train_x[held])
                                      - train_y[held]) ** 2).mean()))
            table[f"{decay:g}"] = float(np.mean(errors)) if errors else float("nan")
        best = min(table, key=lambda k: table[k])
        # MLP probe selected on folds, read once on meta_val
        probe = train_probe(train_x, train_y, 64, float(best), seed=0)
        pred = predict(probe, val_x) + grand_mean
        mse = float(((pred - val_values) ** 2).mean())
        # linear variant
        lin = train_probe(train_x, train_y, 0, float(best), seed=0)
        lin_pred = predict(lin, val_x) + grand_mean
        lin_mse = float(((lin_pred - val_values) ** 2).mean())
        # shuffled-panel control: permute panel rows across train episodes
        if name == "panel":
            rng = np.random.default_rng(0)
            shuffled = train_x[rng.permutation(len(train_x))]
            control = train_probe(shuffled, train_y, 64, float(best), seed=0)
            control_pred = predict(control, val_x) + grand_mean
            control_mse = float(((control_pred - val_values) ** 2).mean())
        else:
            control_mse = None
        pairs = [(r["component"], r["target"],
                  float((pred[i] - val_values[i]) ** 2))
                 for i, r in enumerate(val_rows)]
        interval = component_bootstrap(pairs, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)
        results[name] = {
            "selected_weight_decay": float(best), "fold_table": table,
            "level_mse": mse, "level_rmse": float(np.sqrt(mse)),
            "linear_level_mse": lin_mse,
            "shuffled_panel_control_mse": control_mse,
            "relative_to_grand_mean": mse / float(
                ((np.full_like(val_values, grand_mean) - val_values) ** 2).mean()),
            "component_bootstrap": interval,
        }
        print(f"{name:<24} fold-decay {best:<5} level MSE {mse:8.4f} "
              f"linear {lin_mse:8.4f} shuffled {control_mse if control_mse is not None else float('nan'):8.4f} "
              f"boot {interval['mean']:.4f} [{interval['lo']:.4f},{interval['hi']:.4f}]")

    baseline_mse = float(((np.full_like(val_values, grand_mean)
                           - val_values) ** 2).mean())
    payload["grand_mean_mse"] = baseline_mse
    payload["methods"] = results
    payload["verdict"] = {
        "best_method": min(results, key=lambda k: results[k]["level_mse"]),
        "best_mse": min(r["level_mse"] for r in results.values()),
    }
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print(f"grand-mean baseline level MSE {baseline_mse:.4f}")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
