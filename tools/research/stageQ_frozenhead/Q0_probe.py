"""Stage Q0: joint frozen-feature level identifiability (journal + panel stats
+ frozen ESM pooled), with no trunk coupling.

Stage L's failure mode was measured: training the level head on k=0 episodes
reshaped the SHARED trunk (the head consumed trunk-derived summary and
ligand encodings), degrading k>=1 ordering. This probe measures the ceiling
of a level head built exclusively from FROZEN features — journal/publisher
bags, handcrafted panel statistics and the frozen ESM-150M pooled vector —
so the head cannot couple to the trunk by construction. If the joint probe
reaches <=1.45 level MSE, a preregistered single-stage training candidate
(Stage Q) follows; otherwise the composition is falsified here.

Selection on meta_train component folds; meta_val read once; meta_test never
constructed.
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
OUT = BASE / "Q0_JOINT_FROZEN_IDENTIFIABILITY.json"


def journal_codes(data, spec):
    cells = data.cells
    from collections import Counter
    counts = Counter()
    for i in spec.query:
        for pid in cells[int(i)]["panel_ids"]:
            body = str(pid).split("|")[0]
            parts = body.split("/")
            if len(parts) >= 2 and parts[0].startswith("doi:"):
                publisher = parts[0].split(".")[-1]
                journal = "".join(c for c in parts[1] if c.isalpha())[:4]
                counts[("pub_" + publisher,)] += 1
                counts[("jnl_" + journal,)] += 1
    return counts


def panel_stats(data, spec):
    cells = data.cells
    ligand_mw = {}
    for row in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
        ligand_mw[str(row["drug_key"])] = float(row["molecular_weight"])
    mw, atom_means = [], []
    for i in spec.query:
        ligand = cells[int(i)]["ligand_id"]
        mw.append(ligand_mw[ligand])
        atoms, _, mask = data.ligand_bank.get(ligand)
        n = int(mask.sum())
        if n:
            atom_means.append(atoms[:n].mean(0))
    mw = np.asarray(mw, dtype=np.float32)
    atom_vec = np.asarray(atom_means, dtype=np.float32).mean(0) if atom_means \
        else np.zeros(32, dtype=np.float32)
    return np.concatenate([
        np.asarray([len(spec.query), mw.mean(), mw.std()], dtype=np.float32),
        atom_vec])


def esm_pooled(data, target):
    return np.asarray(data.protein_bank.get(target)[0], dtype=np.float32)


def main():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}

    def build(split):
        specs = data.fixed_nested_episode_banks(
            split, (0, 1, 2, 3, 5), 16, 2, 73101, None)[0]
        rows = []
        for spec in specs:
            episode = normalized_episode(data.materialize(spec), label_scale)
            truth = (episode.query_y.numpy() * label_scale.scale
                     + label_scale.mean)
            rows.append({"component": component_of[spec.target],
                         "target": spec.target,
                         "journal": journal_codes(data, spec),
                         "panel": panel_stats(data, spec),
                         "protein": esm_pooled(data, spec.target),
                         "truth_mean": float(truth.mean())})
        return rows

    train_rows = build("meta_train")
    val_rows = build("meta_val")
    vocab = sorted({key for r in train_rows + val_rows for key in r["journal"]})
    index = {key: i for i, key in enumerate(vocab)}

    def journal_vec(rows):
        x = np.zeros((len(rows), len(vocab)), dtype=np.float32)
        for j, r in enumerate(rows):
            total = sum(r["journal"].values()) or 1
            for key, value in r["journal"].items():
                x[j, index[key]] = value / total
        return x

    def matrix(rows, key):
        return np.stack([r[key] for r in rows]).astype(np.float32)

    groups = {}
    train_j = journal_vec(train_rows)
    val_j = journal_vec(val_rows)
    groups["journal"] = (train_j, val_j)
    groups["journal_panel"] = (
        np.concatenate([train_j, matrix(train_rows, "panel")], 1),
        np.concatenate([val_j, matrix(val_rows, "panel")], 1))
    groups["joint_all"] = (
        np.concatenate([train_j, matrix(train_rows, "panel"),
                        matrix(train_rows, "protein")], 1),
        np.concatenate([val_j, matrix(val_rows, "panel"),
                        matrix(val_rows, "protein")], 1))

    grand = float(np.mean([r["truth_mean"] for r in train_rows]))
    val_values = np.asarray([r["truth_mean"] for r in val_rows])
    components = sorted({r["component"] for r in train_rows})
    order = np.random.default_rng(20260818).permutation(len(components))
    fold_of = {components[int(order[i])]: i % 5 for i in range(len(components))}
    folds = np.asarray([fold_of[r["component"]] for r in train_rows])

    def fit(x, y, hidden, decay, seed):
        torch.manual_seed(seed)
        probe = (nn.Sequential(nn.Linear(x.shape[1], hidden), nn.GELU(),
                               nn.Linear(hidden, 1))
                 if hidden else nn.Linear(x.shape[1], 1))
        optimizer = torch.optim.AdamW(probe.parameters(), lr=3e-3,
                                      weight_decay=decay)
        xt = torch.as_tensor(x, dtype=torch.float32)
        yt = torch.as_tensor(y, dtype=torch.float32)
        for _ in range(600):
            loss = nn.functional.mse_loss(probe(xt).squeeze(-1), yt)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        return probe

    def predict(probe, x):
        with torch.no_grad():
            return probe(torch.as_tensor(x, dtype=torch.float32)
                         ).squeeze(-1).numpy()

    results = {}
    for name, (train_x, val_x) in groups.items():
        train_y = np.asarray([r["truth_mean"] for r in train_rows]) - grand
        table = {}
        for decay in (1e-3, 1e-2, 1e-1, 1.0):
            errs = []
            for fold in range(5):
                keep, held = folds != fold, folds == fold
                if not held.any():
                    continue
                probe = fit(train_x[keep], train_y[keep], 64, decay, fold)
                errs.append(float(((predict(probe, train_x[held])
                                    - train_y[held]) ** 2).mean()))
            table[str(decay)] = float(np.mean(errs)) if errs else float("nan")
        best = min(table, key=table.get)
        probe = fit(train_x, train_y, 64, float(best), 0)
        pred = predict(probe, val_x) + grand
        mse = float(((pred - val_values) ** 2).mean())
        pairs = [(r["component"], r["target"],
                  float((pred[i] - val_values[i]) ** 2))
                 for i, r in enumerate(val_rows)]
        results[name] = {"fold_table": table, "selected_weight_decay": float(best),
                         "level_mse": mse,
                         "component_bootstrap": component_bootstrap(
                             pairs, 9999, 20260816)}
        print(name, "decay", best, "mse", round(mse, 4), flush=True)

    baseline = float(((np.full_like(val_values, grand) - val_values) ** 2).mean())
    payload = {"schema": "MetaSieve.StageQ.JointFrozenIdentifiability.v1",
               "date": "2026-08-18", "grand_mean_baseline": baseline,
               "methods": results, "meta_test": data.seal_record()}
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(json.dumps(payload, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
