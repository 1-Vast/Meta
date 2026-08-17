"""D0: anatomy of target level within meta_train.

Decompose the between-target level variance of meta_train by fitting SGD linear
probes (no closed form) to canonical target level (unique-ligand mean pK) from:

* component identity (one-hot, 258 components)
* document identity (one-hot over DOIs present in the target's cells)
* panel composition (ligand-set statistics over the target's own ligands)
* protein features (ESM-150M pooled + sequence length)
* the joint feature set

Evaluation: 5-fold CV by homology component. Held-out MSE tells how much of the
level variance each covariate family can carry ACROSS components; in-fold vs
out-of-fold contrast shows how much is within-component memorization. This
answers whether 'level' is a protein property, an assay/document property or a
panel property. meta_test never constructed; meta_val not read.
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

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    training_label_scale,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
FOLDS = 5
FOLD_SEED = 20260818
PROBE_STEPS = 400
PROBE_LR = 3e-3
DECAYS = (1e-3, 1e-2, 1e-1, 1.0)
OUT = Path(__file__).resolve().parent / "D0_LEVEL_ANATOMY.json"


def main() -> int:
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    cells = data.cells

    ligand_mw, ligand_scaffold = {}, {}
    for row in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
        ligand_mw[str(row["drug_key"])] = float(row["molecular_weight"])
        ligand_scaffold[str(row["drug_key"])] = str(row.get("scaffold", ""))

    components = sorted(data.components["meta_train"])
    component_index = {c: i for i, c in enumerate(components)}
    documents = sorted({str(pid).split("|")[0]
                        for cell in cells
                        if cell["split"] == "meta_train"
                        for pid in cell["panel_ids"]})
    document_index = {d: i for i, d in enumerate(documents)}

    targets = sorted(data.tasks["meta_train"])
    rows = []
    for target in targets:
        indices = data.tasks["meta_train"][target]
        seen, values = set(), []
        panel_mw, panel_scaffolds, docs = [], set(), set()
        atom_means = []
        for idx in indices:
            cell = cells[int(idx)]
            if cell["ligand_id"] in seen:
                continue
            seen.add(cell["ligand_id"])
            values.append(cell["pK"])
            panel_mw.append(ligand_mw[cell["ligand_id"]])
            panel_scaffolds.add(ligand_scaffold[cell["ligand_id"]])
            for pid in cell["panel_ids"]:
                docs.add(str(pid).split("|")[0])
            atoms, _, mask = data.ligand_bank.get(cell["ligand_id"])
            n = int(mask.sum())
            if n:
                atom_means.append(atoms[:n].mean(0))
        pooled = np.asarray(data.protein_bank.get(target)[0], dtype=np.float32)
        atom_vec = np.asarray(atom_means, dtype=np.float32).mean(0) if atom_means             else np.zeros(32, dtype=np.float32)
        mw = np.asarray(panel_mw, dtype=np.float32)
        panel_vec = np.concatenate([
            np.asarray([len(seen), mw.mean(), mw.std(),
                        float(len(panel_scaffolds)),
                        float(len(panel_scaffolds)) / max(len(seen), 1)],
                       dtype=np.float32),
            atom_vec,
        ])
        comp_onehot = np.zeros(len(components), dtype=np.float32)
        comp_onehot[component_index[cells[int(indices[0])]["protein_group_40"]]] = 1.0
        doc_onehot = np.zeros(len(documents), dtype=np.float32)
        for d in docs:
            if d in document_index:
                doc_onehot[document_index[d]] = 1.0
        length = np.asarray([len(data._protein_sequences[target])], dtype=np.float32)
        rows.append({
            "target": target,
            "component": cells[int(indices[0])]["protein_group_40"],
            "level": float(np.mean(values)),
            "component_onehot": comp_onehot,
            "document_onehot": doc_onehot,
            "panel": panel_vec,
            "protein": np.concatenate([pooled, length]),
        })

    levels = np.asarray([r["level"] for r in rows])
    grand = float(levels.mean())
    y = levels - grand
    total_var = float(np.var(levels))
    n = len(rows)

    comp_order = np.random.default_rng(FOLD_SEED).permutation(len(components))
    fold_of = {components[int(comp_order[rank])]: rank % FOLDS
               for rank in range(len(components))}
    folds = np.asarray([fold_of[r["component"]] for r in rows])

    def fit(x, y_fit, decay, seed=0):
        torch.manual_seed(seed)
        probe = nn.Linear(x.shape[1], 1)
        nn.init.zeros_(probe.weight)
        nn.init.zeros_(probe.bias)
        optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                      weight_decay=decay)
        xt = torch.as_tensor(x, dtype=torch.float32)
        yt = torch.as_tensor(y_fit, dtype=torch.float32)
        for _ in range(PROBE_STEPS):
            loss = nn.functional.mse_loss(probe(xt).squeeze(-1), yt)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
        with torch.no_grad():
            return probe

    def predict(probe, x):
        with torch.no_grad():
            return probe(torch.as_tensor(x, dtype=torch.float32)).squeeze(-1).numpy()

    groups = {
        "component_onehot": np.stack([r["component_onehot"] for r in rows]),
        "document_onehot": np.stack([r["document_onehot"] for r in rows]),
        "panel_composition": np.stack([r["panel"] for r in rows]),
        "protein_esm150_len": np.stack([r["protein"] for r in rows]),
        "joint": np.concatenate([
            np.stack([r["component_onehot"] for r in rows]),
            np.stack([r["document_onehot"] for r in rows]),
            np.stack([r["panel"] for r in rows]),
            np.stack([r["protein"] for r in rows])], 1),
    }

    payload = {
        "schema": "MetaSieve.StageD.LevelAnatomy.v1",
        "date": "2026-08-17",
        "targets": n, "components": len(components),
        "documents": len(documents),
        "between_target_variance_meta_train": total_var,
        "grand_mean_pk": grand,
        "methods": {},
        "meta_test": data.seal_record(),
    }
    print(f"targets {n}, components {len(components)}, docs {len(documents)}, "
          f"between-target var {total_var:.4f}")

    for name, x in groups.items():
        table = {}
        for decay in DECAYS:
            errs = []
            for fold in range(FOLDS):
                keep, held = folds != fold, folds == fold
                if not held.any():
                    continue
                probe = fit(x[keep], y[keep], decay)
                errs.append(float(((predict(probe, x[held])
                                    - y[held]) ** 2).mean()))
            table[f"{decay:g}"] = float(np.mean(errs)) if errs else float("nan")
        best = min(table, key=lambda k: table[k])
        held_mse = table[best]
        # in-fold reference
        probe = fit(x, y, float(best))
        in_mse = float(((predict(probe, x) - y) ** 2).mean())
        payload["methods"][name] = {
            "selected_weight_decay": float(best),
            "fold_table": table,
            "held_out_mse": held_mse,
            "held_out_r2": 1.0 - held_mse / total_var,
            "in_fold_mse": in_mse,
            "in_fold_r2": 1.0 - in_mse / total_var,
            "shared_variance_vs_constant": 1.0 - held_mse / total_var,
        }
        print(f"{name:<22} best decay {best:<5} held MSE {held_mse:8.4f} "
              f"(R2 {1-held_mse/total_var:+.3f})  in-fold R2 {1-in_mse/total_var:+.3f}")

    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
