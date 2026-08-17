"""Stage M0: frozen ChemBERTa ligand-embedding probes (no training of the DTA model)."""
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

from scripts.qpsmp_data import QPSMPData, stable_seed
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    normalized_episode, training_label_scale,
)
from scripts.stageR0_retrieval_falsification import component_bootstrap

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
BASE = Path(__file__).resolve().parent
OUT = BASE / "M0_CHEMBERTA_PROBES.json"
FOLDS = 5


def load_embeddings():
    with np.load(ROOT / "tools/runtime/chemberta_ligand_pooled/embeddings.npz",
                 allow_pickle=False) as store:
        return {str(k): store["pooled"][i].astype(np.float32)
                for i, k in enumerate(store["keys"])}


def panels(data, split):
    from scripts.qpsmp_data import EpisodeSpec
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    out = []
    for target in sorted(data.tasks[split]):
        rng = np.random.default_rng(stable_seed("stageM", 20260818, split, target))
        order = data._unique_ligand_order(data.tasks[split][target], rng)
        if len(order) < 6:
            continue
        cells = tuple(int(i) for i in order[:24])
        out.append((component_of[target], target, EpisodeSpec(
            split, component_of[target], target, (), cells, target)))
    return out


def train_probe(blocks, width, decay, seed):
    torch.manual_seed(seed)
    probe = nn.Linear(width, 1, bias=False)
    nn.init.zeros_(probe.weight)
    optimizer = torch.optim.AdamW(probe.parameters(), lr=3e-3,
                                  weight_decay=decay)
    tensors = [(torch.as_tensor(x, dtype=torch.float32),
                torch.as_tensor(y, dtype=torch.float32)) for x, y in blocks]
    for _ in range(500):
        loss = torch.zeros(())
        for x, y in tensors:
            p = probe(x).squeeze(-1)
            loss = loss + nn.functional.mse_loss(p - p.mean(), y - y.mean())
        (loss / len(tensors)).backward()
        optimizer.step()
        optimizer.zero_grad(set_to_none=True)
    return probe


def within_r(probe, block):
    x, y = block
    with torch.no_grad():
        p = probe(torch.as_tensor(x, dtype=torch.float32)).squeeze(-1).numpy()
    pc, yc = p - p.mean(), y - y.mean()
    d = float(np.sqrt((pc ** 2).sum()) * np.sqrt((yc ** 2).sum()))
    return float((pc * yc).sum() / d) if d > 1e-12 else 0.0


def main():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    embeddings = load_embeddings()
    cells = data.cells

    store = {}
    for split in ("meta_train", "meta_val"):
        rows = []
        for component, target, spec in panels(data, split):
            episode = normalized_episode(data.materialize(spec), label_scale)
            labels = (episode.query_y.numpy() * label_scale.scale
                      + label_scale.mean)
            features = np.stack([embeddings[cells[int(i)]["ligand_id"]]
                                 for i in spec.query]).astype(np.float32)
            rows.append({"component": component, "target": target,
                         "features": features, "labels": labels})
        store[split] = rows
    print("panels: train", len(store["meta_train"]), "val", len(store["meta_val"]))

    width = store["meta_train"][0]["features"].shape[1]
    components = sorted({r["component"] for r in store["meta_train"]})
    order = np.random.default_rng(20260818).permutation(len(components))
    fold_of = {components[int(order[i])]: i % FOLDS for i in range(len(components))}

    train_blocks = [(r["features"], r["labels"]) for r in store["meta_train"]]
    folds = np.asarray([fold_of[r["component"]] for r in store["meta_train"]])
    table = {}
    for decay in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
        scores = []
        for fold in range(FOLDS):
            fit = [b for b, f in zip(train_blocks, folds) if f != fold]
            held = [b for b, f in zip(train_blocks, folds) if f == fold]
            if not fit or not held:
                continue
            probe = train_probe(fit, width, decay, seed=fold)
            scores.extend(within_r(probe, b) for b in held)
        table[str(decay)] = float(np.mean(scores)) if scores else float("nan")
    best = max(table, key=table.get)
    probe = train_probe(train_blocks, width, float(best), seed=0)
    val_blocks = [(r["features"], r["labels"]) for r in store["meta_val"]]
    val_scores = [within_r(probe, b) for b in val_blocks]
    pairs = [(r["component"], r["target"], s)
             for r, s in zip(store["meta_val"], val_scores)]
    q1 = {"fold_table": table, "selected_weight_decay": float(best),
          "meta_val_within_target_r": component_bootstrap(pairs, 9999, 20260816)}

    def level_rows(split):
        out = []
        for spec in data.fixed_nested_episode_banks(
                split, (0, 1, 2, 3, 5), 16, 2, 73101, None)[0]:
            episode = normalized_episode(data.materialize(spec), label_scale)
            truth = (episode.query_y.numpy() * label_scale.scale
                     + label_scale.mean)
            feats = np.stack([embeddings[cells[int(i)]["ligand_id"]]
                              for i in spec.query]).mean(0).astype(np.float32)
            component = cells[int(spec.query[0])]["protein_group_40"]
            out.append({"component": component, "target": spec.target,
                        "features": feats, "truth_mean": float(truth.mean())})
        return out
    lr_train = level_rows("meta_train")
    lr_val = level_rows("meta_val")
    grand = float(np.mean([r["truth_mean"] for r in lr_train]))
    lr_folds = np.asarray([fold_of[r["component"]] for r in lr_train])
    lr_table = {}
    for decay in (1e-4, 1e-3, 1e-2, 1e-1, 1.0):
        errs = []
        for fold in range(FOLDS):
            keep = lr_folds != fold
            held = lr_folds == fold
            if not held.any():
                continue
            fit_blocks = [(lr_train[i]["features"],
                           np.asarray([lr_train[i]["truth_mean"] - grand]))
                          for i in np.flatnonzero(keep)]
            probe = train_probe(fit_blocks, width, decay, seed=fold)
            preds = [probe(torch.as_tensor(lr_train[i]["features"],
                                           dtype=torch.float32)
                          ).squeeze(-1).item()
                     for i in np.flatnonzero(held)]
            truths = [lr_train[i]["truth_mean"] - grand
                      for i in np.flatnonzero(held)]
            errs.append(float(np.mean([(p - t) ** 2
                                       for p, t in zip(preds, truths)])))
        lr_table[str(decay)] = float(np.mean(errs)) if errs else float("nan")
    best_lr = min(lr_table, key=lr_table.get)
    probe_lr = train_probe(
        [(r["features"], np.asarray([r["truth_mean"] - grand]))
         for r in lr_train], width, float(best_lr), seed=0)
    val_pred = np.asarray([probe_lr(torch.as_tensor(r["features"],
                                                    dtype=torch.float32)
                               ).squeeze(-1).item() for r in lr_val]) + grand
    val_values = np.asarray([r["truth_mean"] for r in lr_val])
    q2 = {"fold_table": lr_table, "selected_weight_decay": float(best_lr),
          "level_mse": float(((val_pred - val_values) ** 2).mean()),
          "grand_mean_baseline": float(((np.full_like(val_values, grand)
                                         - val_values) ** 2).mean())}

    payload = {"schema": "MetaSieve.StageM.ChemBERTaProbes.v1",
               "date": "2026-08-18", "q1_within_target_ordering": q1,
               "q2_level": q2, "meta_test": data.seal_record()}
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(json.dumps(payload, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
