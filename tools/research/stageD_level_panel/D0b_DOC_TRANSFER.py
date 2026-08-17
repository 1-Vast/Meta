"""D0b: does document (assay) history transfer target level across targets?

Within meta_train only (meta_val is never read here): for targets sharing at
least one document (DOI), predict each target's canonical level by the mean
level of the OTHER targets in the same document (leave-one-target-out). If
this crude transfer beats the meta_train grand-mean constant, the level is
partly a document/assay property that transfers between proteins. This
complements D0_LEVEL_ANATOMY, which measured document one-hot CV.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
OUT = Path(__file__).resolve().parent / "D0b_DOC_TRANSFER.json"


def main() -> int:
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    cells = data.cells
    target_docs = {}
    target_level = {}
    for target, indices in data.tasks["meta_train"].items():
        docs, seen, values = set(), set(), []
        for idx in indices:
            cell = cells[int(idx)]
            for pid in cell["panel_ids"]:
                docs.add(str(pid).split("|")[0])
            if cell["ligand_id"] in seen:
                continue
            seen.add(cell["ligand_id"])
            values.append(cell["pK"])
        target_docs[target] = docs
        target_level[target] = float(np.mean(values))

    doc_targets = {}
    for target, docs in target_docs.items():
        for doc in docs:
            doc_targets.setdefault(doc, []).append(target)
    shared = [t for t, docs in target_docs.items()
              if any(len(doc_targets[d]) > 1 for d in docs)]

    rows = []
    grand = float(np.mean(list(target_level.values())))
    for target in shared:
        donor_levels = []
        for doc in target_docs[target]:
            for other in doc_targets[doc]:
                if other != target:
                    donor_levels.append(target_level[other])
        if not donor_levels:
            continue
        pred = float(np.mean(donor_levels))
        rows.append({"target": target, "level": target_level[target],
                     "pred": pred, "n_donors": len(donor_levels)})

    mse_transfer = float(np.mean([(r["pred"] - r["level"]) ** 2
                                  for r in rows])) if rows else float("nan")
    mse_grand = float(np.mean([(grand - r["level"]) ** 2
                               for r in rows])) if rows else float("nan")
    var = float(np.var([r["level"] for r in rows])) if rows else float("nan")
    payload = {
        "schema": "MetaSieve.StageD.DocTransfer.v1",
        "date": "2026-08-17",
        "targets_total": len(target_level),
        "targets_sharing_document": len(shared),
        "documents": len(doc_targets),
        "doc_transfer_level_mse": mse_transfer,
        "grand_mean_level_mse": mse_grand,
        "between_target_variance": var,
        "doc_transfer_r2": 1.0 - mse_transfer / var,
        "meta_test": data.seal_record(),
    }
    print(json.dumps(payload, indent=1))
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
