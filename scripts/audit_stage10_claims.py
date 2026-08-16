"""Independent verification of the eleven binding audit findings on Stage 9/10.

Every number here is recomputed from the retained artifacts
(`stage10_retrieval_prior_20260815/BLEND_meta_val.rows.jsonl`, the three accepted
`similarity_only` checkpoints, and the governed corpus). Nothing is taken from
the narrative reports. The output is a machine-readable verdict per finding so
that the reports can be corrected against evidence rather than assertion.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)


def paired_component_bootstrap(by_component: dict[str, list[float]],
                               draws: int, seed: int) -> dict:
    """Resample components with replacement; components are the unit of evidence."""
    keys = sorted(by_component)
    if not keys:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"),
                "components": 0}
    per_component = np.asarray([float(np.mean(by_component[key])) for key in keys])
    rng = np.random.default_rng(seed)
    index = rng.integers(0, len(keys), size=(draws, len(keys)))
    samples = per_component[index].mean(-1)
    return {"mean": float(per_component.mean()),
            "lo": float(np.quantile(samples, 0.025)),
            "hi": float(np.quantile(samples, 0.975)),
            "components": len(keys)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--checkpoint", action="append", default=[])
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=20)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = [json.loads(line) for line in
            args.rows.read_text(encoding="utf-8").splitlines() if line.strip()]
    findings: dict[str, dict] = {}

    # ---- 3. population size ------------------------------------------------
    components = sorted({row["component"] for row in rows})
    targets = sorted({row["target"] for row in rows})
    findings["3_population"] = {
        "claim": "Stage C contains 44 targets and 10 protein components",
        "targets": len(targets), "components": len(components),
        "verdict": "CONFIRMED" if (len(targets) == 44 and len(components) == 10)
        else "DIFFERS",
    }

    # ---- 4. how many components actually improve at k=0 --------------------
    def value(arm: str, k: int, field: str) -> dict[tuple[str, str], float]:
        out: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in rows:
            if row["k"] == k and row["arm"] == arm and row.get(field) is not None:
                out[(row["component"], row["target"])].append(float(row[field]))
        return {key: float(np.mean(v)) for key, v in out.items()}

    base = value("blend_w0", 0, "mse_pk")
    treat = value("blend_w0.5", 0, "mse_pk")
    shared = sorted(set(base) & set(treat))
    per_component: dict[str, list[float]] = defaultdict(list)
    for component, target in shared:
        per_component[component].append(base[(component, target)]
                                        - treat[(component, target)])
    improved = {c: float(np.mean(v)) for c, v in per_component.items()}
    findings["4_component_direction"] = {
        "claim": "only 6 of 10 components improve",
        "per_component_mse_reduction": improved,
        "improving": sum(1 for v in improved.values() if v > 0),
        "total": len(improved),
        "verdict": "CONFIRMED" if sum(1 for v in improved.values() if v > 0) == 6
        else "DIFFERS",
    }

    # ---- 5/6. exact ligand overlap with meta_train -------------------------
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK)
    train_ligands = {cell["ligand_id"] for cell in data.cells
                     if cell["split"] == "meta_train"}
    banks = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, args.query_size, 1, args.evaluation_seed, None)
    overlap_cells = 0
    total_cells = 0
    disjoint_targets: set[str] = set()
    overlap_fraction: dict[str, float] = {}
    for spec in banks[0]:
        ids = [data.cells[i]["ligand_id"] for i in spec.query]
        hits = sum(1 for key in ids if key in train_ligands)
        overlap_cells += hits
        total_cells += len(ids)
        overlap_fraction[spec.target] = hits / max(len(ids), 1)
        if hits == 0:
            disjoint_targets.add(spec.target)
    findings["5_exact_ligand_overlap"] = {
        "claim": "305/624 Stage C query cells (48.9%) have the exact ligand in meta_train",
        "query_cells": total_cells, "exact_in_meta_train": overlap_cells,
        "fraction": overlap_cells / max(total_cells, 1),
        "verdict": "CONFIRMED" if (total_cells == 624 and overlap_cells == 305)
        else "DIFFERS",
    }

    disjoint_component: dict[str, list[float]] = defaultdict(list)
    for component, target in shared:
        if target in disjoint_targets:
            disjoint_component[component].append(base[(component, target)]
                                                 - treat[(component, target)])
    findings["6_ligand_disjoint_subset"] = {
        "claim": "on completely ligand-disjoint targets the improvement is small "
                 "and statistically unresolved",
        "disjoint_targets": len(disjoint_targets),
        "evaluated_targets": sum(len(v) for v in disjoint_component.values()),
        "bootstrap": paired_component_bootstrap(
            disjoint_component, args.bootstrap_draws, 20260815),
        "all_targets_bootstrap": paired_component_bootstrap(
            per_component, args.bootstrap_draws, 20260815),
    }
    entry = findings["6_ligand_disjoint_subset"]["bootstrap"]
    findings["6_ligand_disjoint_subset"]["verdict"] = (
        "CONFIRMED" if entry["lo"] <= 0.0 <= entry["hi"] else "DIFFERS")

    # ---- 9. transport beta against each checkpoint's learned scale ---------
    learned = {}
    for path in args.checkpoint:
        payload = torch.load(Path(path), map_location="cpu", weights_only=False)
        state = payload["model_state"]
        learned[Path(path).parent.name] = {
            "similarity_scale": float(state["transport.similarity_scale"]),
            "log_shrinkage": float(state["transport.log_shrinkage"]),
            "arch": payload["config"].get("arch"),
        }
    scales = [v["similarity_scale"] for v in learned.values()]
    findings["9_transport_beta"] = {
        "claim": "w=0 is only numerically close to the checkpoint baseline because "
                 "the script uses beta=8 instead of the learned similarity scale",
        "script_transport_beta": 8.0,
        "checkpoints": learned,
        "max_abs_difference": float(max(abs(s - 8.0) for s in scales)) if scales else None,
        "verdict": "CONFIRMED" if scales and any(s != 8.0 for s in scales) else "DIFFERS",
    }

    # ---- 10. compression of raw pooled ESM cosine similarity ---------------
    train_targets = sorted({cell["target_id"] for cell in data.cells
                            if cell["split"] == "meta_train"})
    val_targets = sorted({cell["target_id"] for cell in data.cells
                          if cell["split"] == "meta_val"})
    pooled = {t: np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
              for t in (*train_targets, *val_targets)}
    train_matrix = np.stack([pooled[t] for t in train_targets])
    val_matrix = np.stack([pooled[t] for t in val_targets])

    def unit(matrix: np.ndarray) -> np.ndarray:
        return matrix / np.maximum(np.linalg.norm(matrix, axis=-1, keepdims=True), 1e-9)

    def describe(cross: np.ndarray) -> dict:
        top = np.sort(cross, axis=-1)[:, -16:]
        return {"mean": float(cross.mean()), "std": float(cross.std()),
                "p01": float(np.quantile(cross, 0.01)),
                "p99": float(np.quantile(cross, 0.99)),
                "max": float(cross.max()),
                "top16_spread": float(np.mean(top.max(-1) - top.min(-1)))}

    raw = describe(unit(val_matrix) @ unit(train_matrix).T)
    center = train_matrix.mean(0, keepdims=True)          # train-only statistic
    centered = describe(unit(val_matrix - center) @ unit(train_matrix - center).T)
    deviation = train_matrix - center
    covariance = deviation.T @ deviation / max(len(train_matrix) - 1, 1)
    values, vectors = np.linalg.eigh(covariance.astype(np.float64))
    whiten = (vectors / np.sqrt(np.maximum(values, 1e-3))) @ vectors.T
    whitened = describe(unit((val_matrix - center) @ whiten.T.astype(np.float32))
                        @ unit((train_matrix - center) @ whiten.T.astype(np.float32)).T)
    findings["10_esm_similarity_compression"] = {
        "claim": "raw pooled ESM cosine similarities are highly compressed, so the "
                 "current experiment does not rule out all protein representations",
        "raw_cosine": raw, "train_centered_cosine": centered,
        "train_whitened_cosine": whitened,
        "verdict": "CONFIRMED" if raw["std"] < 0.5 * centered["std"] else "DIFFERS",
    }

    # ---- findings verified by construction rather than recomputation -------
    findings["1_status"] = {
        "claim": "Stage C's 12.3% is a meta_val development result, not confirmed",
        "verdict": "ACCEPTED_AS_POLICY",
        "basis": "selection and inference used the same population; see finding 2",
    }
    findings["2_selection_on_inference_data"] = {
        "claim": "beta=24, retrieval source and w=0.5 were selected on the same "
                 "meta_val data used for the bootstrap",
        "verdict": "CONFIRMED",
        "basis": "stage10_retrieval_prior.py sweeps three sources x five weights "
                 "on meta_val and BOOTSTRAP_meta_val.json resamples that same "
                 "population; no outer fold separates selection from inference",
    }
    findings["7_transductive_upper_bound"] = {
        "claim": "Stage A's 25.5% composed retrieval uses query-panel means",
        "verdict": "CONFIRMED",
        "basis": "stage9_k0_decomposition.py builds retrieval[shape|level] as "
                 "shape - shape.mean() + level.mean(); both means are over the "
                 "query panel of the episode, so the predictor is transductive",
    }
    findings["8_offline_evaluator"] = {
        "claim": "the retrieval prior is an offline evaluator, not part of the "
                 "model, checkpoint or standard evaluation path",
        "verdict": "CONFIRMED",
        "basis": "the blend exists only in scripts/stage10_retrieval_prior.py; "
                 "no ARCHITECTURES entry, no checkpoint tensor, and "
                 "scripts/evaluate_qpsmp.py does not compute it",
    }
    findings["11_meta_test_status"] = {
        "claim": "the previous meta_test split has been inspected many times and "
                 "may not be used for selection or as a pristine confirmation set",
        "verdict": "ACCEPTED_AS_POLICY",
    }

    payload = {"schema": "MetaSieve.Stage10AuditVerification.v1",
               "rows": str(args.rows), "findings": findings}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    for key in sorted(findings):
        print(f"{key:34s} {findings[key]['verdict']}")
    print(json.dumps({k: v for k, v in findings.items()
                      if k in ("3_population", "4_component_direction",
                               "5_exact_ligand_overlap", "6_ligand_disjoint_subset",
                               "9_transport_beta", "10_esm_similarity_compression")},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
