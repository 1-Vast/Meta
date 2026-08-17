"""D0c: does JOURNAL (assay provenance) transfer target level?

The panel_ids carry DOIs like "doi:10.1021/jm070284z|Ki|hash". Document
overlap is zero between splits, but the PUBLISHER and JOURNAL codes repeat
across documents: "1021/jm" = ACS Journal of Medicinal Chemistry,
"jbc" = J. Biol. Chem., etc. A journal is a coarse assay covariate:
medicinal-chemistry campaigns and biochemistry papers test different ligand
populations and report systematically different affinity levels. This
diagnostic tests, with frozen SGD probes on the episode-level target,
whether journal identity carries transferable level signal. Selection on
meta_train component folds; meta_val read once. No training of the DTA
model; meta_test never constructed.
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
OUT = BASE / "D0c_JOURNAL_IDENTIFIABILITY.json"


def journal_features(data, episode):
    """Per-episode journal/publisher bag: fraction of cells from each code."""
    cells = data.cells
    from collections import Counter
    counts = Counter()
    for i in episode.spec.query:
        for pid in cells[int(i)]["panel_ids"]:
            body = str(pid).split("|")[0]
            # "doi:10.1021/jm070284z" -> publisher "1021", journal "jm"
            parts = body.split("/")
            if len(parts) >= 2 and parts[0].startswith("doi:"):
                publisher = parts[0].split(".")[-1]
                journal = "".join(c for c in parts[1] if c.isalpha())[:4]
                counts[("pub_" + publisher,)] += 1
                counts[("jnl_" + journal,)] += 1
    return counts


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
            counts = journal_features(data, episode)
            rows.append({"component": component_of[spec.target],
                         "target": spec.target,
                         "counts": counts, "truth_mean": float(truth.mean())})
        return rows

    train_rows = build("meta_train")
    val_rows = build("meta_val")
    vocab = sorted({key for r in train_rows + val_rows for key in r["counts"]})
    print("journal/publisher vocabulary:", len(vocab))

    def matrix(rows):
        x = np.zeros((len(rows), len(vocab)), dtype=np.float32)
        index = {key: i for i, key in enumerate(vocab)}
        for j, r in enumerate(rows):
            total = sum(r["counts"].values()) or 1
            for key, value in r["counts"].items():
                x[j, index[key]] = value / total
        return x

    train_x, val_x = matrix(train_rows), matrix(val_rows)
    grand = float(np.mean([r["truth_mean"] for r in train_rows]))
    val_values = np.asarray([r["truth_mean"] for r in val_rows])
    components = sorted({r["component"] for r in train_rows})
    order = np.random.default_rng(20260818).permutation(len(components))
    fold_of = {components[int(order[i])]: i % 5 for i in range(len(components))}
    folds = np.asarray([fold_of[r["component"]] for r in train_rows])

    def fit(x, y, decay, seed):
        torch.manual_seed(seed)
        probe = nn.Linear(x.shape[1], 1)
        nn.init.zeros_(probe.weight)
        nn.init.zeros_(probe.bias)
        optimizer = torch.optim.AdamW(probe.parameters(), lr=3e-3,
                                      weight_decay=decay)
        xt = torch.as_tensor(x, dtype=torch.float32)
        yt = torch.as_tensor(y, dtype=torch.float32)
        for _ in range(400):
            loss = nn.functional.mse_loss(probe(xt).squeeze(-1), yt)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        return probe

    def predict(probe, x):
        with torch.no_grad():
            return probe(torch.as_tensor(x, dtype=torch.float32)
                         ).squeeze(-1).numpy()

    train_y = np.asarray([r["truth_mean"] for r in train_rows]) - grand
    table = {}
    for decay in (1e-3, 1e-2, 1e-1, 1.0):
        errs = []
        for fold in range(5):
            keep, held = folds != fold, folds == fold
            if not held.any():
                continue
            probe = fit(train_x[keep], train_y[keep], decay, fold)
            errs.append(float(((predict(probe, train_x[held])
                                - train_y[held]) ** 2).mean()))
        table[str(decay)] = float(np.mean(errs)) if errs else float("nan")
    best = min(table, key=table.get)
    probe = fit(train_x, train_y, float(best), 0)
    pred = predict(probe, val_x) + grand
    mse = float(((pred - val_values) ** 2).mean())
    # shuffled-journal control
    rng = np.random.default_rng(0)
    shuffled = train_x[rng.permutation(len(train_x))]
    ctrl = fit(shuffled, train_y, float(best), 0)
    control_mse = float(((predict(ctrl, val_x) + grand - val_values) ** 2).mean())
    baseline = float(((np.full_like(val_values, grand) - val_values) ** 2).mean())
    pairs = [(r["component"], r["target"],
              float((pred[i] - val_values[i]) ** 2))
             for i, r in enumerate(val_rows)]
    interval = component_bootstrap(pairs, 9999, 20260816)
    payload = {"schema": "MetaSieve.StageD.JournalIdentifiability.v1",
               "date": "2026-08-18", "vocabulary": len(vocab),
               "fold_table": table, "selected_weight_decay": float(best),
               "level_mse": mse, "level_rmse": float(np.sqrt(mse)),
               "shuffled_control_mse": control_mse,
               "grand_mean_baseline_mse": baseline,
               "component_bootstrap": interval,
               "meta_test": data.seal_record()}
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(json.dumps(payload, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
