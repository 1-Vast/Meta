"""T3: is the source-transfer headroom a target property or a document lottery?

T1 measured an honest split-sample oracle over source heads at +0.1222 -- nine
times the target's own label-fitted head.  Before treating that as headroom for
any routing mechanism, one thing must be checked, and it is the thing the pair
split cannot check: **pair-halves A and B are drawn from the same query rows in
the same query documents.**  A source head selected on half A is selected partly
for fitting those particular molecules in those particular documents.

Two worlds produce an identical +0.1222:

  (a) *Target structure.*  Some source targets genuinely share this target's
      chemistry-response, so the head that wins does so for every document of
      the target.  Then routing is worth building.
  (b) *Document lottery.*  With ~550 candidate directions, some head aligns with
      whatever idiosyncratic ordering a given document happens to carry.  The
      winner does not replicate on a different document of the same target.
      Then the headroom is unexploitable and no router can reach it.

The discriminator is a **document-level** split of the query rows: choose the
source on documents half 1, score it on documents half 2, both of which belong to
the same target and are both disjoint from the support documents.  If the
per-source concordance ranking replicates across the two halves, we are in world
(a); if the cross-half rank correlation is ~0, we are in world (b).

`rank_correlation_across_documents` is the headline number.  `oracle_doc_split`
is the exploitable ceiling; `oracle_insample_half2` shows what selection bias
alone would have produced.

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
from research.psep.psep_transfer import MIN_TRAIN, concordance_matrix

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_oracle_stability_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_oracle_stability_records_2026-08-02.parquet"

MIN_DOCUMENTS = 4
MIN_PAIRS_PER_HALF = 24
MAX_PAIRS_PER_HALF = 4096


def pairs_within(documents: np.ndarray, label: np.ndarray, mask: np.ndarray, limit: int):
    index = np.flatnonzero(mask)
    if len(index) < 2:
        return None
    left, right = np.triu_indices(len(index), k=1)
    left, right = index[left], index[right]
    keep = (documents[left] == documents[right]) & (label[left] != label[right])
    left, right = left[keep], right[keep]
    if len(left) < MIN_PAIRS_PER_HALF:
        return None
    if len(left) > limit:
        rng = np.random.default_rng(SEED)
        pick = rng.choice(len(left), size=limit, replace=False)
        left, right = left[pick], right[pick]
    return left, right


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
        units.append(split)
        heads.append(ridge_fit(design[split.train], substrate.residual[split.train], (RIDGE,))[0])
    matrix = np.stack(heads)
    components = np.asarray([unit.component for unit in units])

    records: list[dict[str, object]] = []
    for position, unit in enumerate(units):
        label = substrate.affinity[unit.evaluation]
        base = substrate.base[unit.evaluation]
        design_eval = design[unit.evaluation]
        evaluation_docs = documents[unit.evaluation]
        unique = np.unique(evaluation_docs)
        if len(unique) < MIN_DOCUMENTS:
            continue
        rng = np.random.default_rng(int(sha256(f"{SEED}:doc:{unit.unit}".encode()).hexdigest()[:8], 16))
        shuffled = rng.permutation(unique)
        first = set(shuffled[: len(shuffled) // 2].tolist())
        mask_one = np.isin(evaluation_docs, list(first))
        pair_one = pairs_within(evaluation_docs, label, mask_one, MAX_PAIRS_PER_HALF)
        pair_two = pairs_within(evaluation_docs, label, ~mask_one, MAX_PAIRS_PER_HALF)
        if pair_one is None or pair_two is None:
            continue

        allowed = np.flatnonzero(components != unit.component)
        sources = matrix[allowed]
        predictions = base[:, None] + design_eval @ sources.T
        ci_one = concordance_matrix(predictions, label, *pair_one)
        ci_two = concordance_matrix(predictions, label, *pair_two)
        base_two = float(concordance_matrix(base[:, None], label, *pair_two)[0])
        uniform_two = float(concordance_matrix(
            (base + design_eval @ sources.mean(axis=0))[:, None], label, *pair_two)[0])

        finite = np.isfinite(ci_one) & np.isfinite(ci_two)
        records.append({
            "unit": unit.unit, "component": unit.component, "endpoint": unit.endpoint,
            "documents": int(len(unique)),
            "n_sources": int(len(allowed)),
            "base": base_two,
            "uniform": uniform_two,
            "oracle_doc_split": float(ci_two[int(ci_one.argmax())]),
            "oracle_insample_half2": float(ci_two.max()),
            "own": float(concordance_matrix(
                (base + design_eval @ matrix[position])[:, None], label, *pair_two)[0]),
            "rank_correlation_across_documents": (
                float(spearmanr(ci_one[finite], ci_two[finite])[0]) if finite.sum() > 20 else np.nan
            ),
            "top10_overlap": float(len(
                set(np.argsort(-ci_one)[:10].tolist()) & set(np.argsort(-ci_two)[:10].tolist())
            ) / 10.0),
        })
        if position % 100 == 0:
            print(f"  {position}/{len(units)}", flush=True)

    frame = pd.DataFrame.from_records(records)
    summary: dict[str, object] = {}
    for arm in ("uniform", "oracle_doc_split", "oracle_insample_half2", "own"):
        frame[f"{arm}__vs_base"] = frame[arm] - frame["base"]
        frame[f"{arm}__vs_uniform"] = frame[arm] - frame["uniform"]
        summary[arm] = {
            "vs_base": paired_bootstrap(frame, f"{arm}__vs_base"),
            "vs_uniform": paired_bootstrap(frame, f"{arm}__vs_uniform"),
        }
    stability = paired_bootstrap(frame, "rank_correlation_across_documents")
    overlap = paired_bootstrap(frame, "top10_overlap")

    exploitable = summary["oracle_doc_split"]["vs_uniform"]["mean"]
    inflated = summary["oracle_insample_half2"]["vs_uniform"]["mean"]
    stable = bool(stability["lower95"] > 0.10)
    if not stable and exploitable < 0.005:
        verdict = "HEADROOM_IS_A_DOCUMENT_LOTTERY_NOT_TARGET_STRUCTURE"
    elif stable:
        verdict = "SOURCE_QUALITY_IS_A_STABLE_TARGET_PROPERTY"
    else:
        verdict = "PARTIAL_STABILITY_HEADROOM_LARGELY_DOCUMENT_SPECIFIC"

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED, "capacity": CAPACITY,
            "split": "query documents halved; source chosen on half 1, scored on half 2",
            "metric": "within-document pair concordance, component bootstrap",
        },
        "counts": {"units": int(len(frame)), "components": int(frame.component.nunique())},
        "summary": summary,
        "rank_correlation_across_documents": stability,
        "top10_overlap_across_documents": overlap,
        "headline": {
            "oracle_document_split_over_uniform": exploitable,
            "oracle_insample_over_uniform": inflated,
            "document_specific_component": inflated - exploitable,
        },
        "verdict": verdict,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="T3 oracle stability across documents")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({"verdict": payload["verdict"], "headline": payload["headline"],
                      "counts": payload["counts"]}, indent=2))
    for arm, cell in payload["summary"].items():
        entry = cell["vs_uniform"]
        print(f"  {arm:<24s}{entry['mean']:+.4f} [{entry['lower95']:+.4f},{entry['upper95']:+.4f}]")
    for name in ("rank_correlation_across_documents", "top10_overlap_across_documents"):
        cell = payload[name]
        print(f"  {name:<36s}{cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}]")


if __name__ == "__main__":
    main()
