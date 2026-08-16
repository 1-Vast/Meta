"""T1: is there anything for a source router to route to?

The operator gate closed the `support -> adaptation` framing.  The reframed
question is whether the *source task population* carries transferable structure
that helps a scarce recipient target: `{D_i, p_i}_{i=1}^N -> f_t`.

Before building any router this gate bounds the whole mechanism class.  Every
proposal in the family -- task mixtures, conditional expert routing, MoE over
biological tasks, similarity-weighted transfer, learned source selection -- is
some map from target-observable features to a *weighting over source heads*.  So
the family's ceiling is the best achievable weighting and its floor is the
uniform weighting, which is just a global multitask model.

**Selection bias is the trap here and it is handled explicitly.**  A first pass
of this gate reported an oracle ceiling of +0.165 -- ten times the target's own
label-fitted head (+0.0152).  That was not headroom: it was the maximum of ~550
noisy concordance estimates chosen on the very pairs it was scored on.  The
maximum of N noisy estimates is biased upward by roughly `sigma * sqrt(2 ln N)`,
which at N=553 reproduces the entire "effect".

Every selective arm is therefore **split-sample**: the source (or the top-m set)
is chosen on half the within-document pairs and scored on the *other* half.  All
arms, selective or not, are scored on the same half B, so they are comparable.
`oracle_insample` is retained purely to display the bias it removes.

Similarity temperatures are likewise fixed in advance and reported one by one; no
per-unit temperature is selected on evaluation labels.

Arms (all leave-one-homology-component-out, all scored on pair-half B)
  uniform            mean of all source heads        <- global/multitask floor
  random_source      one random source head
  universal_best     head with best mean transfer on *other* components
  protein_top1       most sequence-similar source
  chemistry_top1     most chemistry-similar source
  {protein,chemistry}_softmax_T   fixed-temperature similarity mixtures
  oracle_split       best source chosen on half A     <- honest 1-expert ceiling
  oracle_top5_split  best five chosen on half A
  oracle_insample    best source chosen on half B     <- biased, shown for contrast
  own                the target's own fitted head     <- uses labels; upper bound

Reads the `discover` role only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits, paired_bootstrap
from research.psep.psep_m0 import rich_basis
from research.psep.psep_m2 import ridge_fit
from research.psep.psep_m4 import CAPACITY, RIDGE

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_transfer_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_transfer_records_2026-08-02.parquet"

TEMPERATURES = (0.02, 0.05, 0.1, 0.25)
TOP_M = 5
MIN_TRAIN = 20
MAX_PAIRS = 4096
MIN_HALF_PAIRS = 16


def concordance_matrix(
    predictions: np.ndarray, label: np.ndarray, left: np.ndarray, right: np.ndarray
) -> np.ndarray:
    truth = np.sign(label[left] - label[right])[:, None]
    gap = np.sign(predictions[left] - predictions[right])
    correct = (gap == truth).astype(np.float64) + 0.5 * (gap == 0)
    return correct.mean(axis=0)


def within_document_pairs(documents: np.ndarray, label: np.ndarray, limit: int = MAX_PAIRS):
    left, right = np.triu_indices(len(documents), k=1)
    keep = (documents[left] == documents[right]) & (label[left] != label[right])
    left, right = left[keep], right[keep]
    if len(left) > limit:
        rng = np.random.default_rng(SEED)
        pick = rng.choice(len(left), size=limit, replace=False)
        left, right = left[pick], right[pick]
    return left, right


def similarity_matrices(units, bits: np.ndarray, accessions: list[str]) -> dict[str, np.ndarray]:
    from sklearn.feature_extraction.text import CountVectorizer

    proteins = pd.read_csv(
        ROOT / "dataset" / "public" / "papyrus_05_7" / "raw" / "05.7_combined_set_protein_targets.tsv.xz",
        sep="\t", low_memory=False, usecols=["target_id", "Sequence"],
    )
    lookup: dict[str, str] = {}
    for target_id, sequence in zip(proteins.target_id, proteins.Sequence):
        if isinstance(sequence, str):
            lookup.setdefault(str(target_id).split("_")[0], sequence)
    vectoriser = CountVectorizer(analyzer="char", ngram_range=(4, 4), binary=True, lowercase=False)
    matrix = vectoriser.fit_transform([lookup.get(a, "X") for a in accessions]).astype(np.float32)
    matrix.data[:] = 1.0
    intersection = (matrix @ matrix.T).toarray()
    sizes = np.asarray(matrix.sum(axis=1)).ravel()
    protein = intersection / np.maximum(sizes[:, None] + sizes[None, :] - intersection, 1.0)

    fingerprint = np.stack([bits[unit.train].mean(axis=0) for unit in units])
    fingerprint = fingerprint / np.maximum(np.linalg.norm(fingerprint, axis=1, keepdims=True), 1e-9)
    return {"protein": protein, "chemistry": fingerprint @ fingerprint.T}


def _indicator(chosen: np.ndarray, length: int) -> np.ndarray:
    weight = np.zeros(length)
    weight[chosen] = 1.0 / max(len(chosen), 1)
    return weight


def run(substrate_dir: Path, role: str, output: Path, records_path: Path) -> dict[str, object]:
    started = time.time()
    from scipy.sparse import load_npz
    from scipy.stats import spearmanr

    basis, substrate, stats = rich_basis(substrate_dir, role)
    design = np.ascontiguousarray(basis[:, :CAPACITY])
    splits = [s for s in build_splits(substrate.rows) if s.regime == "separated"]
    documents = substrate.rows.docs.astype(str).to_numpy()
    accession_of = substrate.rows.accession.astype(str).to_numpy()

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
    endpoints = np.asarray([unit.endpoint for unit in units])
    print(f"{len(units)} units / {len(set(components))} components", flush=True)

    bits_all = load_npz(substrate_dir / "morgan.npz").tocsr()
    bits = np.asarray(bits_all[substrate.rows.structure_row.to_numpy()].todense(), dtype=np.float32)
    similarity = similarity_matrices(units, bits, [accession_of[u.evaluation[0]] for u in units])

    transfer_b = np.full((len(units), len(units)), np.nan, dtype=np.float64)
    rng = np.random.default_rng(SEED)
    records: list[dict[str, object]] = []
    routability: list[dict[str, float]] = []

    for position, unit in enumerate(units):
        left_a, right_a = unit.pairs_a
        left_b, right_b = unit.pairs_b
        label = substrate.affinity[unit.evaluation]
        base = substrate.base[unit.evaluation]
        design_eval = design[unit.evaluation]
        allowed = np.flatnonzero(components != unit.component)

        source_predictions = base[:, None] + design_eval @ matrix[allowed].T
        ci_a = concordance_matrix(source_predictions, label, left_a, right_a)
        ci_b = concordance_matrix(source_predictions, label, left_b, right_b)
        transfer_b[position, allowed] = ci_b

        def score_b(weight: np.ndarray) -> float:
            head = weight @ matrix[allowed]
            return float(concordance_matrix(
                (base + design_eval @ head)[:, None], label, left_b, right_b)[0])

        uniform_weight = np.full(len(allowed), 1.0 / len(allowed))
        top5_a = np.argsort(-ci_a)[:TOP_M]
        row: dict[str, object] = {
            "unit": unit.unit, "component": unit.component, "endpoint": unit.endpoint,
            "n_sources": int(len(allowed)),
            "n_pairs_b": int(len(left_b)),
            "base": float(concordance_matrix(base[:, None], label, left_b, right_b)[0]),
            "uniform": score_b(uniform_weight),
            "random_source": float(ci_b[rng.integers(len(allowed))]),
            "oracle_split": float(ci_b[int(ci_a.argmax())]),
            "oracle_top5_split": score_b(_indicator(top5_a, len(allowed))),
            "oracle_insample": float(ci_b.max()),
            "own": float(concordance_matrix(
                (base + design_eval @ matrix[position])[:, None], label, left_b, right_b)[0]),
        }
        for name in ("protein", "chemistry"):
            values = similarity[name][position, allowed]
            row[f"{name}_top1"] = float(ci_b[int(values.argmax())])
            row[f"{name}_top5"] = score_b(_indicator(np.argsort(-values)[:TOP_M], len(allowed)))
            for temperature in TEMPERATURES:
                weight = np.exp((values - values.max()) / temperature)
                row[f"{name}_softmax_t{temperature:g}"] = score_b(weight / weight.sum())
            finite = np.isfinite(ci_b) & np.isfinite(values)
            if finite.sum() > 10:
                routability.append({
                    "unit": unit.unit, "component": unit.component, "feature": name,
                    "spearman": float(spearmanr(values[finite], ci_b[finite])[0]),
                })
        same = endpoints[allowed] == unit.endpoint
        if same.any() and (~same).any():
            row["endpoint_match_delta"] = float(ci_b[same].mean() - ci_b[~same].mean())
        records.append(row)
        if position % 100 == 0:
            print(f"  {position}/{len(units)}", flush=True)

    frame = pd.DataFrame.from_records(records)

    universal = []
    for position, unit in enumerate(units):
        others = components != unit.component
        mean_transfer = np.nanmean(transfer_b[others], axis=0)
        mean_transfer[components == unit.component] = -np.inf
        universal.append(transfer_b[position, int(np.nanargmax(mean_transfer))])
    frame["universal_best"] = universal

    arms = ["uniform", "random_source", "universal_best", "protein_top1", "protein_top5",
            "chemistry_top1", "chemistry_top5", "oracle_split", "oracle_top5_split",
            "oracle_insample", "own"]
    arms += [f"{name}_softmax_t{t:g}" for name in ("protein", "chemistry") for t in TEMPERATURES]

    summary: dict[str, object] = {"versus_base": {}, "versus_uniform": {}, "negative_transfer_rate": {}}
    for arm in arms:
        frame[f"{arm}__vs_base"] = frame[arm] - frame["base"]
        frame[f"{arm}__vs_uniform"] = frame[arm] - frame["uniform"]
        summary["versus_base"][arm] = paired_bootstrap(frame, f"{arm}__vs_base")
        summary["versus_uniform"][arm] = paired_bootstrap(frame, f"{arm}__vs_uniform")
        summary["negative_transfer_rate"][arm] = float((frame[f"{arm}__vs_base"] < 0).mean())

    routability_frame = pd.DataFrame.from_records(routability)
    routing_signal = {
        feature: paired_bootstrap(part, "spearman")
        for feature, part in routability_frame.groupby("feature")
    } if not routability_frame.empty else {}

    honest_ceiling = summary["versus_uniform"]["oracle_split"]["mean"]
    biased_ceiling = summary["versus_uniform"]["oracle_insample"]["mean"]
    universal_gain = summary["versus_uniform"]["universal_best"]["mean"]
    router_arms = [a for a in arms if "softmax" in a or a.endswith("_top1") or a.endswith("_top5")]
    best_router = max(summary["versus_uniform"][a]["mean"] for a in router_arms)
    best_router_arm = max(router_arms, key=lambda a: summary["versus_uniform"][a]["mean"])

    target_specific = honest_ceiling - universal_gain
    if honest_ceiling <= 0.005:
        verdict = "NO_SOURCE_ROUTING_HEADROOM_EVEN_WITH_HINDSIGHT"
    elif target_specific <= 0.005:
        verdict = "SOURCE_HEADROOM_IS_NOT_TARGET_SPECIFIC"
    elif best_router > 0.005:
        verdict = "ROUTABLE_TARGET_SPECIFIC_SOURCE_STRUCTURE"
    else:
        verdict = "HEADROOM_EXISTS_BUT_NO_OBSERVABLE_PREDICTS_THE_SOURCE"

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "substrate": stats,
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False},
        "protocol": {
            "seed": SEED, "capacity": CAPACITY, "ridge": RIDGE,
            "leave_out": "homology component",
            "metric": "within-document pair concordance on pair-half B, component bootstrap",
            "selection": "every selective arm chooses on pair-half A and is scored on half B",
            "temperatures": list(TEMPERATURES), "top_m": TOP_M,
        },
        "counts": {"units": len(units), "components": int(len(set(components)))},
        "summary": summary,
        "routing_signal_spearman": routing_signal,
        "headline": {
            "honest_oracle_ceiling_over_uniform": honest_ceiling,
            "insample_oracle_over_uniform_BIASED": biased_ceiling,
            "selection_bias_removed": biased_ceiling - honest_ceiling,
            "universal_best_over_uniform": universal_gain,
            "target_specific_headroom": target_specific,
            "best_observable_router_over_uniform": best_router,
            "best_observable_router_arm": best_router_arm,
        },
        "verdict": verdict,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="T1 source-population transfer bound")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records)
    print(json.dumps({"verdict": payload["verdict"], "headline": payload["headline"]}, indent=2))
    print(f"\n{'arm':<26s}{'vs base':>26s}{'vs uniform':>26s}{'neg':>7s}")
    for arm, cell in payload["summary"]["versus_base"].items():
        other = payload["summary"]["versus_uniform"][arm]
        rate = payload["summary"]["negative_transfer_rate"][arm]
        print(f"{arm:<26s}{cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}]"
              f"{other['mean']:+.4f} [{other['lower95']:+.4f},{other['upper95']:+.4f}]{rate:>7.2f}")
    print("\nrouting signal (Spearman, similarity vs realised transfer on half B):")
    for feature, cell in payload["routing_signal_spearman"].items():
        print(f"  {feature:<12s}{cell['mean']:+.4f} [{cell['lower95']:+.4f},{cell['upper95']:+.4f}]")


if __name__ == "__main__":
    main()
