"""Adjudicate crossed-interaction existence from frozen X1A-R rows."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import t


ROOT = Path(__file__).resolve().parents[2]
X1AR_ROOT = ROOT / "report" / "crossed_interaction" / "x1ar_direct_dd"
X1AR_GATE = X1AR_ROOT / "gate.json"
X1AR_ROWS = X1AR_ROOT / "dd_rows.jsonl"
PREREG = ROOT / "research" / "crossed_interaction" / "PREREG_X1B_INTERACTION_EXISTENCE.md"
OUTPUT = ROOT / "report" / "crossed_interaction" / "x1b_interaction"
MIN_RATIO = 0.5


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def cluster_t_inference(values_by_cluster: dict[str, list[float]]) -> dict:
    keys = sorted(values_by_cluster)
    if len(keys) < 2:
        raise ValueError("at least two dependency components are required")
    sizes = np.asarray([len(values_by_cluster[key]) for key in keys], dtype=np.float64)
    means = np.asarray([np.mean(values_by_cluster[key]) for key in keys], dtype=np.float64)
    total = float(sizes.sum())
    weights = sizes / total
    estimate = float(weights @ means)
    contributions = weights * (means - estimate)
    standard_error = float(math.sqrt(len(keys) / (len(keys) - 1) * (contributions @ contributions)))
    critical = float(t.ppf(0.95, len(keys) - 1))
    lcb95 = estimate - critical * standard_error
    loo = []
    for index, key in enumerate(keys):
        keep = np.arange(len(keys)) != index
        loo.append({
            "excluded_component": key,
            "estimate": float((sizes[keep] @ means[keep]) / sizes[keep].sum()),
        })
    return {
        "components": len(keys),
        "rectangles": int(total),
        "estimate_tau2": estimate,
        "cluster_robust_standard_error": standard_error,
        "one_sided_t_critical": critical,
        "lcb95_tau2": lcb95,
        "leave_one_component_out": loo,
    }


def exact_rademacher_pvalue(values_by_cluster: dict[str, list[float]]) -> float:
    keys = sorted(values_by_cluster)
    if len(keys) > 20:
        raise ValueError("exact sign enumeration is limited to 20 components")
    total = sum(len(values_by_cluster[key]) for key in keys)
    contributions = np.asarray([sum(values_by_cluster[key]) / total for key in keys],
                               dtype=np.float64)
    observed = float(contributions.sum())
    exceed = 0
    draws = 0
    for signs in itertools.product((-1.0, 1.0), repeat=len(keys)):
        draws += 1
        if float(np.asarray(signs) @ contributions) >= observed - 1e-15:
            exceed += 1
    return exceed / draws


def adjudicate_endpoint(rows: list[dict], x1ar: dict, noise: dict, endpoint: str) -> dict:
    by_cluster: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        if row["endpoint"] == endpoint:
            by_cluster[row["dependency_cluster"]].append(float(row["Z"]))
    inference = cluster_t_inference(by_cluster)
    ratio = (math.sqrt(max(inference["estimate_tau2"], 0.0)) /
             float(noise["sigma_ucb95"]))
    gates = {
        "B1_contract": bool(by_cluster) and all(math.isfinite(value)
                                                 for values in by_cluster.values()
                                                 for value in values),
        "B2_dependence": bool(x1ar["pass"]),
        "B3_positive_lcb": inference["lcb95_tau2"] > 0.0,
        "B4_design_margin": ratio >= MIN_RATIO,
    }
    if not gates["B1_contract"] or not gates["B2_dependence"]:
        verdict = "X1B_DATA_OR_NOISE_CONTRACT_INVALID"
    elif not gates["B3_positive_lcb"]:
        verdict = "X1B_INTERACTION_NOT_DETECTED"
    elif not gates["B4_design_margin"]:
        verdict = "X1B_INTERACTION_PRESENT_BELOW_DESIGN_MARGIN"
    else:
        verdict = "X1B_REAL_CROSSED_INTERACTION_IDENTIFIED"
    result = {
        **inference,
        "endpoint": endpoint,
        "sigma_rep_ucb95": float(noise["sigma_ucb95"]),
        "interaction_rms_over_sigma_rep_ucb95": ratio,
        "minimum_ratio": MIN_RATIO,
        "gates": gates,
        "terminal_verdict": verdict,
    }
    if endpoint == "Kd":
        result["exact_rademacher_one_sided_pvalue"] = exact_rademacher_pvalue(by_cluster)
    return result


def run(output: Path = OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"no-clobber output exists: {output}")
    x1ar_gate = json.loads(X1AR_GATE.read_text(encoding="utf-8"))
    rows = read_jsonl(X1AR_ROWS)
    if sha256_file(X1AR_ROWS) != x1ar_gate["hashes"]["dd_rows"]:
        raise RuntimeError("X1A-R row hash mismatch")
    passing = set(x1ar_gate["passing_endpoints"])
    endpoints = {}
    for endpoint in ("Ki", "Kd"):
        if endpoint not in passing:
            endpoints[endpoint] = {
                "endpoint": endpoint,
                "terminal_verdict": "NOT_RUN_X1A_R_PRECONDITION_FAILED",
            }
            continue
        endpoints[endpoint] = adjudicate_endpoint(
            rows, x1ar_gate["endpoints"][endpoint], x1ar_gate["noise"][endpoint], endpoint)
    identified = [endpoint for endpoint, result in endpoints.items()
                  if result["terminal_verdict"] == "X1B_REAL_CROSSED_INTERACTION_IDENTIFIED"]
    terminal = ("X1B_REAL_CROSSED_INTERACTION_IDENTIFIED_" + "_AND_".join(identified).upper()
                if identified else "X1B_NO_ENDPOINT_IDENTIFIED")
    output.mkdir(parents=True)
    report = {
        "stage": "E-AFF-X1B_CROSSED_INTERACTION_EXISTENCE",
        "terminal_verdict": terminal,
        "identified_endpoints": identified,
        "endpoints": endpoints,
        "hashes": {
            "x1ar_gate": sha256_file(X1AR_GATE),
            "x1ar_rows": sha256_file(X1AR_ROWS),
            "prereg": sha256_file(PREREG),
            "runner": sha256_file(Path(__file__)),
        },
        "label_reads": "no additional reads; consumes frozen X1A-R rows",
        "training_performed": False,
        "gpu_used": False,
        "x2_preregistration_authorized_endpoints": identified,
    }
    (output / "gate.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                                      encoding="utf-8")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    try:
        result = run(args.output)
    except Exception as exc:
        print(json.dumps({"terminal_verdict": "X1B_DATA_OR_NOISE_CONTRACT_INVALID",
                          "error": f"{type(exc).__name__}: {exc}"}, indent=2), file=sys.stderr)
        return 2
    print(json.dumps({"terminal_verdict": result["terminal_verdict"],
                      "identified_endpoints": result["identified_endpoints"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
