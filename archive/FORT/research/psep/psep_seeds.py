"""R3: is the ranking-objective effect stable to training randomness?

R2 measured, on 373 homology components with 5-fold cross-fitting, that swapping
an MSE objective for a bounded smoothed-concordance surrogate on within-document
pairs is worth **+0.0108 [+0.0009, +0.0207]** in within-document concordance.  The
interval excludes zero but its lower bound sits below the registered 0.005 MDE,
and the whole thing rests on a single training seed.  Fold-to-fold spread was
~0.016 -- comparable to the effect itself -- so seed variance is not obviously
negligible.

This gate re-runs only the two arms that matter, across seeds.

**The folds are deliberately held fixed.**  `component_fold` keys on the frozen
programme SEED, not on the training seed, so every seed sees exactly the same
component partition and the same test rows.  The per-component paired delta is
therefore comparable across seeds, and the only thing varying is initialisation
and batch order.  That isolates training randomness, which is the question.

Aggregation, in the order that respects the nesting:
  1. per (seed, component): mean paired delta over that component's units
  2. per component: mean across seeds  -> bootstrapped over components
  3. per seed: the component-mean contrast -> gives the between-seed SD

Reported: the seed-averaged contrast with a component bootstrap, the per-seed
contrasts, and the between-seed SD.  A mechanism whose sign flips across seeds is
not a mechanism.

Reads the `discover` role only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits, paired_bootstrap
from research.psep.psep_operator import concordance
from research.psep.psep_representation import (
    DEVICE, FOLDS, MIN_DOC_ROWS, component_fold, train,
)
from research.psep.psep_transfer import MIN_HALF_PAIRS, within_document_pairs

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_seeds_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_seeds_records_2026-08-02.parquet"

SEEDS = (20260802, 20260803, 20260804, 20260805, 20260806)
ARMS = (("mse", "raw"), ("rank", "raw"))
MDE = 0.005


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    from scipy.sparse import load_npz

    rows = pd.read_parquet(substrate_dir / "rows.parquet")
    rows = rows.loc[rows.role == role].reset_index(drop=True)
    bits = load_npz(substrate_dir / "morgan.npz").tocsr()
    descriptors = np.load(substrate_dir / "descriptors.npy")
    index = rows.structure_row.to_numpy()
    scale = np.maximum(descriptors.std(axis=0), 1e-6)
    standardised = ((descriptors - descriptors.mean(axis=0)) / scale)[index]
    dense = np.hstack([np.asarray(bits[index].todense(), dtype=np.float32),
                       standardised.astype(np.float32)])
    design = torch.as_tensor(dense, device=DEVICE)
    del dense

    affinity_np = rows.affinity.to_numpy(dtype=np.float64)
    base_np = rows.base.to_numpy(dtype=np.float64)
    documents = rows.docs.astype(str).to_numpy()
    context_np = pd.factorize(documents)[0]
    affinity = torch.as_tensor(affinity_np, dtype=torch.float32, device=DEVICE)
    context = torch.as_tensor(context_np, dtype=torch.long, device=DEVICE)

    folds = rows.component.map(component_fold).to_numpy()
    splits = [s for s in build_splits(rows) if s.regime == "separated"]
    units = []
    for split in splits:
        left, right = within_document_pairs(documents[split.evaluation], affinity_np[split.evaluation])
        if len(left) < 2 * MIN_HALF_PAIRS:
            continue
        split.pair_left, split.pair_right = left, right
        units.append(split)

    document_groups: dict[int, list[np.ndarray]] = {}
    frame = pd.DataFrame({"doc": context_np, "fold": folds, "row": np.arange(len(rows))})
    for (fold, _), group in frame.groupby(["fold", "doc"], sort=False):
        if len(group) >= MIN_DOC_ROWS:
            document_groups.setdefault(int(fold), []).append(group.row.to_numpy())

    print(f"{len(rows)} rows | {len(units)} units | {rows.component.nunique()} components | "
          f"{len(SEEDS)} seeds x {len(ARMS)} arms x {FOLDS} folds | {DEVICE}", flush=True)

    records: list[dict[str, object]] = []
    for seed in SEEDS:
        for objective, context_mode in ARMS:
            name = f"{objective}_{context_mode}"
            for fold in range(FOLDS):
                test_mask = folds == fold
                validation_fold = (fold + 1) % FOLDS
                train_rows = np.flatnonzero(~test_mask & (folds != validation_fold))
                train_documents = [g for f, groups in document_groups.items()
                                   if f not in (fold, validation_fold) for g in groups]
                val_groups = document_groups.get(validation_fold, [])[:2000]
                model, info = train(objective, context_mode, design, affinity, affinity,
                                    context, train_rows, val_groups, train_documents, seed + fold)
                model.eval()
                with torch.no_grad():
                    for unit in units:
                        if not test_mask[unit.evaluation[0]]:
                            continue
                        batch = torch.as_tensor(unit.evaluation, device=DEVICE)
                        prediction = model(design[batch]).cpu().numpy().astype(np.float64)
                        label = affinity_np[unit.evaluation]
                        records.append({
                            "seed": seed, "arm": name, "fold": fold,
                            "unit": unit.unit, "component": unit.component,
                            "endpoint": unit.endpoint,
                            "ci_within": concordance(prediction, label, unit.pair_left, unit.pair_right),
                            "base_ci_within": concordance(base_np[unit.evaluation], label,
                                                          unit.pair_left, unit.pair_right),
                        })
                print(f"  seed {seed} {name} fold {fold}: val_ci={info['best_val_ci']:.4f}", flush=True)

    table = pd.DataFrame.from_records(records)
    table["gain"] = table.ci_within - table.base_ci_within
    table.to_parquet(records_path, index=False)

    wide = table.pivot_table(index=["seed", "component", "unit"], columns="arm",
                             values="gain", aggfunc="mean").reset_index()
    wide["delta"] = wide["rank_raw"] - wide["mse_raw"]

    per_seed = {}
    for seed, part in wide.groupby("seed"):
        per_seed[int(seed)] = {
            "rank_vs_base": paired_bootstrap(part.rename(columns={"rank_raw": "v"}), "v"),
            "mse_vs_base": paired_bootstrap(part.rename(columns={"mse_raw": "v"}), "v"),
            "rank_minus_mse": paired_bootstrap(part, "delta"),
        }
    seed_contrasts = np.array([per_seed[s]["rank_minus_mse"]["mean"] for s in sorted(per_seed)])

    averaged = wide.groupby(["component", "unit"], as_index=False).agg(
        delta=("delta", "mean"), rank_raw=("rank_raw", "mean"), mse_raw=("mse_raw", "mean"))
    combined = paired_bootstrap(averaged, "delta")
    rank_combined = paired_bootstrap(averaged.rename(columns={"rank_raw": "v"}), "v")
    mse_combined = paired_bootstrap(averaged.rename(columns={"mse_raw": "v"}), "v")

    sign_stable = bool(np.all(seed_contrasts > 0) or np.all(seed_contrasts < 0))
    clears_mde = bool(combined["lower95"] > MDE)
    if not sign_stable:
        verdict = "OBJECTIVE_EFFECT_SIGN_UNSTABLE_ACROSS_SEEDS"
    elif clears_mde:
        verdict = "OBJECTIVE_EFFECT_REPLICATES_AND_CLEARS_MDE"
    else:
        verdict = "OBJECTIVE_EFFECT_REPLICATES_BUT_BELOW_MDE"

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "device": str(DEVICE),
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seeds": list(SEEDS), "arms": [f"{o}_{c}" for o, c in ARMS], "folds": FOLDS,
            "folds_fixed_across_seeds": True,
            "metric": "within-document pair concordance vs the frozen cross-fitted linear base",
            "aggregation": "seed -> component mean -> component bootstrap; between-seed SD reported",
            "mde": MDE,
        },
        "counts": {"units": len(units), "components": int(rows.component.nunique())},
        "per_seed": per_seed,
        "seed_contrast_mean": float(seed_contrasts.mean()),
        "seed_contrast_sd": float(seed_contrasts.std(ddof=1)),
        "seed_contrast_min": float(seed_contrasts.min()),
        "seed_contrast_max": float(seed_contrasts.max()),
        "combined": {
            "rank_vs_base": rank_combined,
            "mse_vs_base": mse_combined,
            "rank_minus_mse": combined,
        },
        "sign_stable_across_seeds": sign_stable,
        "clears_mde": clears_mde,
        "verdict": verdict,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="R3 multi-seed replication")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({"verdict": payload["verdict"],
                      "seed_contrast_mean": payload["seed_contrast_mean"],
                      "seed_contrast_sd": payload["seed_contrast_sd"]}, indent=2))
    print(f"\n{'seed':>10s}{'rank vs base':>22s}{'mse vs base':>22s}{'rank - mse':>22s}")
    for seed, cell in payload["per_seed"].items():
        r, m, d = cell["rank_vs_base"], cell["mse_vs_base"], cell["rank_minus_mse"]
        print(f"{seed:>10}{r['mean']:+.4f}{'':>15s}{m['mean']:+.4f}{'':>15s}"
              f"{d['mean']:+.4f} [{d['lower95']:+.4f},{d['upper95']:+.4f}]")
    c = payload["combined"]
    print(f"\nseed-averaged, component bootstrap (n={c['rank_minus_mse']['components']}):")
    for name, cell in c.items():
        print(f"  {name:<16s}{cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}]")


if __name__ == "__main__":
    main()
