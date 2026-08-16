"""T2: can k support labels *select a source target* even though they cannot
estimate a task vector?

T1 established the shape of the problem.  Leave-one-homology-component-out, on
pair-half B with every selective arm choosing on half A:

    uniform (global multitask floor)      +0.0000
    universal_best (one globally good head) +0.0118
    own head, fitted from ~140 labels     +0.0134
    honest oracle over source heads       +0.1222   <- nine times the own head
    chemistry / protein similarity routing  ~0.000  <- no observable predicts it

So there is a large, genuinely target-specific pool of transferable structure in
the source population, and *no label-free feature locates it*.  That is precisely
the situation in which the few remaining labels should be spent on **selection
rather than estimation**.

The information argument.  Estimating the target's 266-dimensional head from k
labels is hopeless: k=5 spans 5 of ~115 effective dimensions (M2/M3).  But
choosing among N ~ 550 pre-fitted source heads requires only about
`log2(550) ~ 9.1` bits.  Five real-valued measurements can carry that.  The
question is entirely whether the support *fit* of a source head predicts its
*transfer* to the query documents -- and that is what this gate measures, with
the same split and the same firewall.

Selection criterion is computed on the target's own support rows, which are
document- and assay-disjoint from the query rows, and centred, because the
per-context offset carries 68 % of residual variance and within-document
concordance is invariant to it.

Arms
  uniform            mean of all source heads                (floor)
  support_select     argmin centred support MSE over sources (the mechanism)
  support_top5       mean of the five best-fitting sources
  wrong_support      selection driven by *another target's* support   (control)
  random_source      no selection                                     (control)
  ridge_k            the target's own ridge head fitted on the same k (incumbent)
  oracle_split       best source chosen on pair-half A                (ceiling)

Reads the `discover` role only.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits, paired_bootstrap
from research.psep.psep_m0 import rich_basis
from research.psep.psep_m2 import ridge_fit
from research.psep.psep_m4 import CAPACITY, RIDGE
from research.psep.psep_transfer import (
    MIN_HALF_PAIRS, MIN_TRAIN, TOP_M, concordance_matrix, within_document_pairs,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_source_select_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_source_select_records_2026-08-02.parquet"

K_SWEEP = (1, 3, 5, 10, 20)
DRAWS = 12
RIDGE_GRID = (10.0, 100.0, 1000.0)


def centred_support_fit(
    design_support: np.ndarray, residual_support: np.ndarray, matrix: np.ndarray
) -> np.ndarray:
    """Centred MSE of every source head on the support rows.  Lower is better.

    Centring is not cosmetic: the per-context offset is 68 % of residual variance
    and support rows sit in different documents from the query, so an uncentred
    criterion would rank heads by how well they reproduce a nuisance that does
    not transfer and does not affect the metric.
    """

    predicted = design_support @ matrix.T                      # (k, n_sources)
    predicted = predicted - predicted.mean(axis=0, keepdims=True)
    target = (residual_support - residual_support.mean())[:, None]
    return np.mean((predicted - target) ** 2, axis=0)


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    from scipy.stats import spearmanr

    basis, substrate, stats = rich_basis(substrate_dir, role)
    design = np.ascontiguousarray(basis[:, :CAPACITY])
    splits = [s for s in build_splits(substrate.rows) if s.regime == "separated"]
    documents = substrate.rows.docs.astype(str).to_numpy()

    units, heads = [], []
    for split in splits:
        if len(split.train) < MIN_TRAIN:
            continue
        left, right = within_document_pairs(documents[split.evaluation],
                                            substrate.affinity[split.evaluation])
        if len(left) < 2 * MIN_HALF_PAIRS:
            continue
        order = np.random.default_rng(SEED ^ hash(split.unit) % (2**32)).permutation(len(left))
        half = len(order) // 2
        split.pairs_a = (left[order[:half]], right[order[:half]])
        split.pairs_b = (left[order[half:]], right[order[half:]])
        units.append(split)
        heads.append(ridge_fit(design[split.train], substrate.residual[split.train], (RIDGE,))[0])
    matrix = np.stack(heads)
    components = np.asarray([unit.component for unit in units])
    print(f"{len(units)} units / {len(set(components))} components", flush=True)

    records: list[dict[str, object]] = []
    diagnostics: list[dict[str, float]] = []
    for position, unit in enumerate(units):
        left_a, right_a = unit.pairs_a
        left_b, right_b = unit.pairs_b
        label = substrate.affinity[unit.evaluation]
        base = substrate.base[unit.evaluation]
        design_eval = design[unit.evaluation]
        allowed = np.flatnonzero(components != unit.component)
        sources = matrix[allowed]

        predictions = base[:, None] + design_eval @ sources.T
        ci_a = concordance_matrix(predictions, label, left_a, right_a)
        ci_b = concordance_matrix(predictions, label, left_b, right_b)
        reference = float(concordance_matrix(base[:, None], label, left_b, right_b)[0])

        def score_b(head: np.ndarray) -> float:
            return float(concordance_matrix(
                (base + design_eval @ head)[:, None], label, left_b, right_b)[0])

        uniform = score_b(sources.mean(axis=0))
        donor = units[int(np.random.default_rng(
            int(sha256(f"{SEED}:donor:{unit.unit}".encode()).hexdigest()[:8], 16)
        ).choice(np.flatnonzero(components != unit.component)))]

        row_common = {
            "unit": unit.unit, "component": unit.component, "endpoint": unit.endpoint,
            "base": reference, "uniform": uniform,
            "oracle_split": float(ci_b[int(ci_a.argmax())]),
            "n_sources": int(len(allowed)),
        }
        rng = np.random.default_rng(int(sha256(f"{SEED}:sel:{unit.unit}".encode()).hexdigest()[:8], 16))
        for k in K_SWEEP:
            if len(unit.train) < k or len(donor.train) < k:
                continue
            picks = {"support_select": [], "support_top5": [], "wrong_support": [],
                     "random_source": [], "ridge_k": []}
            rank_quality = []
            for _ in range(DRAWS):
                support = rng.choice(unit.train, size=k, replace=False)
                fit = centred_support_fit(design[support], substrate.residual[support], sources)
                order = np.argsort(fit)
                picks["support_select"].append(float(ci_b[order[0]]))
                picks["support_top5"].append(score_b(sources[order[:TOP_M]].mean(axis=0)))
                picks["random_source"].append(float(ci_b[rng.integers(len(allowed))]))

                wrong = rng.choice(donor.train, size=k, replace=False)
                wrong_fit = centred_support_fit(design[wrong], substrate.residual[wrong], sources)
                picks["wrong_support"].append(float(ci_b[int(wrong_fit.argmin())]))

                best = -np.inf
                for ridge in RIDGE_GRID:
                    weight = ridge_fit(design[support], substrate.residual[support], (ridge,))[0]
                    best = max(best, score_b(weight))
                picks["ridge_k"].append(best)

                finite = np.isfinite(fit) & np.isfinite(ci_b)
                if finite.sum() > 10:
                    rank_quality.append(float(spearmanr(-fit[finite], ci_b[finite])[0]))
            row = dict(row_common)
            row["k"] = k
            for name, values in picks.items():
                row[name] = float(np.mean(values)) if values else float("nan")
            row["support_fit_rank_spearman"] = float(np.mean(rank_quality)) if rank_quality else float("nan")
            diagnostics.append({"unit": unit.unit, "component": unit.component, "k": k,
                                "spearman": row["support_fit_rank_spearman"]})
            records.append(row)
        if position % 100 == 0:
            print(f"  {position}/{len(units)}", flush=True)

    frame = pd.DataFrame.from_records(records)
    arms = ("uniform", "support_select", "support_top5", "wrong_support",
            "random_source", "ridge_k", "oracle_split")
    summary: dict[str, object] = {}
    for k in K_SWEEP:
        part = frame.loc[frame.k == k].copy()
        if part.empty:
            continue
        cell: dict[str, object] = {}
        for arm in arms:
            part[f"{arm}__vs_base"] = part[arm] - part["base"]
            part[f"{arm}__vs_uniform"] = part[arm] - part["uniform"]
            cell[arm] = {
                "vs_base": paired_bootstrap(part, f"{arm}__vs_base"),
                "vs_uniform": paired_bootstrap(part, f"{arm}__vs_uniform"),
                "negative_transfer_rate": float((part[f"{arm}__vs_base"] < 0).mean()),
            }
        part["select_minus_wrong"] = part["support_select"] - part["wrong_support"]
        part["select_minus_ridge"] = part["support_select"] - part["ridge_k"]
        cell["support_select_minus_wrong_support"] = paired_bootstrap(part, "select_minus_wrong")
        cell["support_select_minus_ridge_k"] = paired_bootstrap(part, "select_minus_ridge")
        cell["support_fit_rank_spearman"] = paired_bootstrap(part, "support_fit_rank_spearman")
        summary[f"k{k}"] = cell

    def value(k: int, arm: str, field: str = "mean") -> float:
        try:
            return summary[f"k{k}"][arm]["vs_uniform"][field]
        except KeyError:
            return float("nan")

    passes = {
        k: bool(
            value(k, "support_select", "lower95") > 0.005
            and summary[f"k{k}"]["support_select_minus_wrong_support"]["lower95"] > 0.005
            and summary[f"k{k}"]["support_select_minus_ridge_k"]["lower95"] > 0.005
        )
        for k in K_SWEEP if f"k{k}" in summary
    }
    verdict = ("SOURCE_SELECTION_IS_A_LOAD_BEARING_MECHANISM"
               if any(passes.values()) else "SUPPORT_CANNOT_SELECT_THE_SOURCE")

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED, "capacity": CAPACITY, "k_sweep": list(K_SWEEP), "draws": DRAWS,
            "criterion": "centred support MSE of each pre-fitted source head",
            "metric": "within-document pair concordance on pair-half B, component bootstrap",
            "leave_out": "homology component",
            "ridge_k_grant": "own-head ridge granted its best value per episode (optimistic)",
        },
        "counts": {"units": len(units), "components": int(len(set(components)))},
        "summary": summary,
        "gate": passes,
        "verdict": verdict,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="T2 support-driven source selection")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({"verdict": payload["verdict"], "gate": payload["gate"]}, indent=2))
    for key, cell in payload["summary"].items():
        print(f"\n--- {key} ---")
        for arm in ("uniform", "support_select", "support_top5", "wrong_support",
                    "random_source", "ridge_k", "oracle_split"):
            entry = cell[arm]["vs_uniform"]
            print(f"  {arm:<16s}{entry['mean']:+.4f} [{entry['lower95']:+.4f},{entry['upper95']:+.4f}]"
                  f"  neg={cell[arm]['negative_transfer_rate']:.2f}")
        for name in ("support_select_minus_wrong_support", "support_select_minus_ridge_k",
                     "support_fit_rank_spearman"):
            entry = cell[name]
            print(f"  {name:<34s}{entry['mean']:+.4f} [{entry['lower95']:+.4f},{entry['upper95']:+.4f}]")


if __name__ == "__main__":
    main()
