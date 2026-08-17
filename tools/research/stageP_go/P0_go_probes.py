"""Stage P0: protein function annotation (GO) probes from ProteinKG25.

External knowledge lane: sequence-match the governed targets against the
local ProteinKG25 corpus (protein sequences + GO annotation triples) and
probe whether GO function bags carry target-level signal. No training of
the DTA model; component-fold selection on meta_train; meta_val read once;
meta_test never constructed.
"""
from __future__ import annotations

import hashlib
import json
import zipfile
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
OUT = BASE / "P0_GO_PROBES.json"
ZIP = ROOT / "dataset/raw/protein_knowledge/proteinkg25/ProteinKG25.zip"


def seq_key(sequence):
    return hashlib.sha256(sequence.upper().encode("utf-8")).hexdigest()


def load_kg():
    with zipfile.ZipFile(ZIP) as archive:
        seqs = archive.read("ProteinKG25/protein_seq.txt").decode("utf-8").splitlines()
        go2id = {}
        for line in archive.read("ProteinKG25/go2id.txt").decode("utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 2:
                go2id[parts[0]] = int(parts[1])
        go_type = archive.read("ProteinKG25/go_type.txt").decode("utf-8").splitlines()
        triples = []
        for name in ("protein_go_train_triplet.txt", "protein_go_valid_triplet.txt",
                     "protein_go_test_triplet.txt"):
            for line in archive.read("ProteinKG25/" + name).decode("utf-8").splitlines():
                parts = line.split()
                if len(parts) >= 3:
                    triples.append((int(parts[0]), int(parts[2])))
    # sequence index: sha256(seq) -> set of KG protein ids
    seq_to_ids = {}
    for kg_id, sequence in enumerate(seqs):
        key = seq_key(sequence)
        seq_to_ids.setdefault(key, []).append(kg_id)
    go_of = {}
    for protein, go in triples:
        go_of.setdefault(protein, set()).add(go)
    go_types = {}
    for go_id_str, kg_id in go2id.items():
        if kg_id < len(go_type):
            go_types[kg_id] = go_type[kg_id]
    return seq_to_ids, go_of, go_types, len(go2id)


def main():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    seq_to_ids, go_of, go_types, n_go = load_kg()
    print("KG proteins:", len(seq_to_ids), "GO terms:", n_go)

    targets = {}
    for split in ("meta_train", "meta_val"):
        for target in sorted(data.tasks[split]):
            targets[target] = split

    matched = {}
    for target in targets:
        ids = seq_to_ids.get(target, [])
        if not ids:
            continue
        gos = set()
        for kg_id in ids:
            gos |= go_of.get(kg_id, set())
        if gos:
            matched[target] = gos
    print("targets with GO annotations:", len(matched), "of", len(targets))

    # GO one-hot over the union of matched targets' terms (vocab from train only)
    train_matched = [t for t in matched if targets[t] == "meta_train"]
    vocab = sorted({go for t in train_matched for go in matched[t]})
    index = {go: i for i, go in enumerate(vocab)}
    print("GO vocab (train):", len(vocab))

    def featurize(target):
        x = np.zeros(len(vocab), dtype=np.float32)
        for go in matched.get(target, set()):
            if go in index:
                x[index[go]] = 1.0
        return x

    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}

    def build(split):
        specs = data.fixed_nested_episode_banks(
            split, (0, 1, 2, 3, 5), 16, 2, 73101, None)[0]
        rows = []
        for spec in specs:
            if spec.target not in matched:
                continue
            episode = normalized_episode(data.materialize(spec), label_scale)
            truth = (episode.query_y.numpy() * label_scale.scale
                     + label_scale.mean)
            rows.append({"component": component_of[spec.target],
                         "target": spec.target,
                         "features": featurize(spec.target),
                         "truth_mean": float(truth.mean())})
        return rows

    train_rows = build("meta_train")
    val_rows = build("meta_val")
    print("episodes with GO features: train", len(train_rows), "val", len(val_rows))

    grand = float(np.mean([r["truth_mean"] for r in train_rows]))
    val_values = np.asarray([r["truth_mean"] for r in val_rows])
    train_x = np.stack([r["features"] for r in train_rows])
    val_x = np.stack([r["features"] for r in val_rows])
    train_y = np.asarray([r["truth_mean"] for r in train_rows]) - grand

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
    baseline = float(((np.full_like(val_values, grand) - val_values) ** 2).mean())
    pairs = [(r["component"], r["target"],
              float((pred[i] - val_values[i]) ** 2))
             for i, r in enumerate(val_rows)]
    payload = {
        "schema": "MetaSieve.StageP.GOProbes.v1",
        "date": "2026-08-18",
        "kg_proteins": len(seq_to_ids), "go_terms_kg": n_go,
        "targets_matched": len(matched), "targets_total": len(targets),
        "go_vocab_train": len(vocab),
        "episodes_train": len(train_rows), "episodes_val": len(val_rows),
        "fold_table": table, "selected_weight_decay": float(best),
        "level_mse": mse, "grand_mean_baseline": baseline,
        "component_bootstrap": component_bootstrap(pairs, 9999, 20260816),
        "meta_test": data.seal_record(),
    }
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(json.dumps(payload, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
