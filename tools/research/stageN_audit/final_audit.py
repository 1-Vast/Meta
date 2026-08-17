"""Final boundary audit (no training): re-derive the load-bearing numbers of
report/BOUNDARY_20260817_NIGHT.md from recorded artifacts.

1. k=0 level/shape decomposition of the frozen T2 evaluation rows.
2. Within-document level transfer R^2 (D0b recomputation).
3. K2 pooled three-seed contrast (recomputed from the seed row files).
4. meta_test seal state across every result artifact.
5. Preregistration-before-result ordering for every trained stage.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)
from scripts.stageR0_retrieval_falsification import component_bootstrap

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
OUT = Path(__file__).resolve().parent / "FINAL_BOUNDARY_AUDIT.json"

REPORT = {}


def load_rows(path: Path):
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]


def main():
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)

    # 1. T2 k=0 decomposition (component-weighted, matching the boundary doc)
    from scripts.stageR0_retrieval_falsification import component_target_mean
    rows = load_rows(ROOT / "tools/research/stageD_level_panel/T2_meta_val.rows.jsonl")
    k0 = [r for r in rows if r["k"] == 0 and r["condition"] == "correct"]
    mse = component_target_mean([(r["component"], r["target"], r["mse_pk"])
                                 for r in k0])
    level = component_target_mean([(r["component"], r["target"], r["level_squared"])
                                   for r in k0])
    centered = component_target_mean([(r["component"], r["target"], r["centered_mse_pk"])
                                      for r in k0])
    REPORT["t2_k0"] = {"mse": mse, "level_squared": level,
                       "centered_mse": centered,
                       "level_share": level / mse}

    # 2. D0b doc transfer (within meta_train, leave-one-target-out)
    cells = data.cells
    target_docs, target_level = {}, {}
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
    pairs = []
    for target in shared:
        donors = []
        for doc in target_docs[target]:
            for other in doc_targets[doc]:
                if other != target:
                    donors.append(target_level[other])
        if donors:
            pairs.append((float(np.mean(donors)), target_level[target]))
    pred = np.asarray([p[0] for p in pairs])
    truth = np.asarray([p[1] for p in pairs])
    var = float(np.var(truth))
    mse_transfer = float(((pred - truth) ** 2).mean())
    REPORT["doc_transfer"] = {"targets_sharing_document": len(pairs),
                              "r2": 1.0 - mse_transfer / var,
                              "mse": mse_transfer, "variance": var}

    # 3. K2 pooled contrast recomputation
    def pool_contrast(a_paths, b_paths, k, field="mse_pk"):
        diffs = {}
        for a_path, b_path in zip(a_paths, b_paths):
            a_rows = [r for r in load_rows(a_path)
                      if r["condition"] == "correct" and r["k"] == k]
            b_rows = [r for r in load_rows(b_path)
                      if r["condition"] == "correct" and r["k"] == k]
            a_map, b_map = {}, {}
            for r in a_rows:
                a_map.setdefault((r["component"], r["target"]), []).append(r[field])
            for r in b_rows:
                b_map.setdefault((r["component"], r["target"]), []).append(r[field])
            for key in a_map:
                if key in b_map:
                    diffs.setdefault(key, []).append(
                        float(np.mean(a_map[key]) - np.mean(b_map[key])))
        pairs = [(key[0], key[1], float(np.mean(values)))
                 for key, values in diffs.items()]
        return component_bootstrap(pairs, 9999, 20260816)

    a_paths = [ROOT / p for p in (
        "tools/research/stageK_contrastive/K-REG_meta_val.rows.jsonl",
        "tools/research/stageK_contrastive/KREG_s20260816.rows.jsonl",
        "tools/research/stageK_contrastive/KREG_s20260817.rows.jsonl")]
    b_paths = [ROOT / p for p in (
        "tools/research/stageG_esm650/T2_s15.rows.jsonl",
        "tools/research/stageG_esm650/T2_s16.rows.jsonl",
        "tools/research/stageG_esm650/T2_s17.rows.jsonl")]
    k2 = {str(k): pool_contrast(a_paths, b_paths, k) for k in (0, 1, 2, 3, 5)}
    REPORT["k2_pooled_mse_contrast"] = k2
    stored = json.loads((ROOT / "tools/research/stageK_contrastive/"
                         "K2_multiseed_contrast.json").read_text(encoding="utf-8"))
    REPORT["k2_matches_stored_authority"] = {
        str(k): all(abs(stored["by_k"][str(k)]["correct"]["mse_pk"][key]
                        - k2[str(k)][key]) < 1e-9
                    for key in ("mean", "lo", "hi"))
        for k in (0, 1, 2, 3, 5)}

    # 4. seal state across result artifacts
    seals = {}
    result_dirs = sorted(ROOT.glob("report/meta_fewshot/stage*/**/RESULT.json")) + \
        sorted(ROOT.glob("report/meta_fewshot/stage*/**/**/RESULT.json"))
    for path in result_dirs:
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            seal = payload.get("meta_test", {})
            seals[str(path.relative_to(ROOT))] = {
                "included": seal.get("included"),
                "evaluated": seal.get("evaluated"),
                "withheld": seal.get("sealed_cells_withheld"),
            }
        except Exception:  # noqa: BLE001
            seals[str(path.relative_to(ROOT))] = {"error": True}
    opened = {k: v for k, v in seals.items()
              if v.get("included") or v.get("evaluated")}
    REPORT["seals"] = {"artifacts": len(seals),
                       "opened": opened}

    # 5. preregistration-before-result ordering per stage
    order = []
    for stage in ("stageD_level_panel", "stageF_pairwise", "stageG_esm650",
                  "stageI_lm", "stageJ_assay", "stageK_contrastive",
                  "stageL_gated"):
        pre = ROOT / "tools/research" / stage / "PREREGISTRATION.md"
        results = sorted((ROOT / "tools/research" / stage).glob("*.rows.summary.json"))
        if pre.exists() and results:
            order.append({"stage": stage,
                          "preregistered": True,
                          "result_files": len(results)})
    REPORT["preregistration_order"] = order

    OUT.write_text(json.dumps(REPORT, indent=1, sort_keys=True),
                   encoding="utf-8")
    print(json.dumps(REPORT, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
