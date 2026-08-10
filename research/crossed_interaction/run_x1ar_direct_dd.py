"""Execute X1A-R direct-DD dependence on the frozen exact-assay design."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "research" / "crossed_interaction"
ASSAY_ROOT = HERE / "artifacts" / "x1a_r_assays"
ASSAY_ROWS = ASSAY_ROOT / "assay_selection.jsonl"
ASSAY_MANIFEST = ASSAY_ROOT / "manifest.json"
ASSAY_SHA256 = "e293fa2bf689bca63825c96f0503491b9c470208dea0947068e95b1f44d02fb4"
PREREG = HERE / "PREREG_X1A_R_DIRECT_DD_DEPENDENCE.md"
X1B_PREREG = HERE / "PREREG_X1B_INTERACTION_EXISTENCE.md"
CHEMBL_DB = (ROOT / "dataset" / "raw" / "source_affinity" / "chembl37_sqlite_v1" /
             "extracted" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db")
OUTPUT = ROOT / "report" / "crossed_interaction" / "x1ar_direct_dd"
RHO_STAR = {"Ki": 0.0915, "Kd": 0.0164}
REQUIRED_N = 245
MAX_SHARE = 0.25
LR_ONE_SIDED_95 = float(chi2.ppf(0.90, 1))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _profile_loglik(groups: list[np.ndarray], rho: float) -> float:
    if not 0 <= rho < 1 or len(groups) < 2:
        return -math.inf
    one_minus = 1.0 - rho
    denom = np.array([one_minus + rho * len(group) for group in groups])
    sums = np.array([group.sum() for group in groups])
    ns = np.array([len(group) for group in groups], dtype=np.float64)
    weights = ns / denom
    mu = float(np.sum(sums / denom) / np.sum(weights))
    qform = 0.0
    logdet = 0.0
    total = 0
    for group, group_denom in zip(groups, denom):
        residual = group - mu
        residual_sum = float(residual.sum())
        qform += (float(residual @ residual) / one_minus -
                  rho * residual_sum ** 2 / (one_minus * group_denom))
        logdet += (len(group) - 1) * math.log(one_minus) + math.log(group_denom)
        total += len(group)
    if qform <= 0 or total <= 1:
        return -math.inf
    sigma2 = qform / total
    return -0.5 * (total * (math.log(2 * math.pi) + 1.0 + math.log(sigma2)) + logdet)


def profile_icc(groups: list[np.ndarray]) -> dict:
    groups = [np.asarray(group, dtype=np.float64) for group in groups if len(group)]
    result = minimize_scalar(lambda rho: -_profile_loglik(groups, rho),
                             bounds=(0.0, 0.999), method="bounded",
                             options={"xatol": 1e-10})
    rho_hat = float(result.x)
    ll_max = -float(result.fun)
    target = ll_max - 0.5 * LR_ONE_SIDED_95

    def equation(rho):
        return _profile_loglik(groups, rho) - target

    if equation(0.999) >= 0:
        ucb = 0.999
    else:
        ucb = float(brentq(equation, max(rho_hat, 1e-12), 0.999))
    return {"rho_mle": rho_hat, "rho_ucb95": ucb,
            "profile_loglik_max": ll_max, "groups": len(groups),
            "rows": sum(map(len, groups))}


def load_values(activity_ids: set[int]) -> tuple[dict[int, float], dict]:
    connection = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("CREATE TEMP TABLE wanted(activity_id INTEGER PRIMARY KEY)")
    connection.executemany("INSERT INTO wanted VALUES (?)", ((value,) for value in activity_ids))
    cursor = connection.execute(
        "SELECT a.activity_id,a.standard_relation,a.standard_type,a.pchembl_value "
        "FROM activities a JOIN wanted w ON w.activity_id=a.activity_id")
    values = {}
    relations = Counter()
    types = Counter()
    rows = 0
    for activity_id, relation, standard_type, value in cursor:
        rows += 1
        relations[str(relation)] += 1
        types[str(standard_type)] += 1
        if relation == "=" and value is not None:
            values[int(activity_id)] = float(value)
    connection.close()
    return values, {"rows_returned": rows, "usable_rows": len(values),
                    "relations": dict(relations), "types": dict(types)}


def build_dd_rows(selections: list[dict], values: dict[int, float]) -> tuple[list[dict], dict]:
    cell_values = []
    for row in selections:
        arrays = []
        for identifiers in row["cell_activity_ids"]:
            observed = np.asarray([values[value] for value in identifiers if value in values],
                                  dtype=np.float64)
            if len(observed) != len(identifiers) or not len(observed):
                raise RuntimeError(f"missing selected activity for {row['rectangle_id']}")
            arrays.append(observed)
            cell_values.append((row["endpoint"], observed))

    noise = {}
    for endpoint in ("Ki", "Kd"):
        arrays = [array for family, array in cell_values if family == endpoint and len(array) >= 2]
        degrees = sum(len(array) - 1 for array in arrays)
        if degrees <= 0:
            raise RuntimeError(f"no exact-assay replicate degrees of freedom for {endpoint}")
        sse = sum(float(((array - array.mean()) ** 2).sum()) for array in arrays)
        sigma2 = sse / degrees
        sigma2_ucb = degrees * sigma2 / float(chi2.ppf(0.05, degrees))
        noise[endpoint] = {"replicate_cells": len(arrays), "degrees_of_freedom": degrees,
                           "sigma2": sigma2, "sigma2_ucb95": sigma2_ucb,
                           "sigma_ucb95": math.sqrt(sigma2_ucb)}

    output = []
    cursor = 0
    for row in selections:
        arrays = [cell_values[cursor + offset][1] for offset in range(4)]
        cursor += 4
        means = [float(array.mean()) for array in arrays]
        dd = means[0] - means[1] - means[2] + means[3]
        d_value = dd / 2.0
        sigma2_ucb = noise[row["endpoint"]]["sigma2_ucb95"]
        v_d_ucb = sigma2_ucb * sum(1.0 / len(array) for array in arrays) / 4.0
        output.append({
            "rectangle_id": row["rectangle_id"], "endpoint": row["endpoint"],
            "dependency_cluster": row["dependency_cluster"],
            "D": d_value, "v_D_ucb": v_d_ucb, "Z": d_value ** 2 - v_d_ucb,
            "replicate_counts": [len(array) for array in arrays],
            "all_four_cells_replicated": all(len(array) >= 2 for array in arrays),
        })
    return output, noise


def endpoint_dependence(rows: list[dict], endpoint: str) -> dict:
    endpoint_rows = [row for row in rows if row["endpoint"] == endpoint]
    by_cluster = defaultdict(list)
    for row in endpoint_rows:
        by_cluster[row["dependency_cluster"]].append(row["Z"])
    groups = [np.asarray(by_cluster[key], dtype=np.float64) for key in sorted(by_cluster)]
    profile = profile_icc(groups)
    loo_ucbs = []
    if endpoint == "Kd":
        for index in range(len(groups)):
            loo_ucbs.append(profile_icc(groups[:index] + groups[index + 1:])["rho_ucb95"])
    conservative_ucb = max([profile["rho_ucb95"], *loo_ucbs])
    sizes = np.asarray([len(group) for group in groups], dtype=np.float64)
    total = int(sizes.sum())
    mean_influence = float((sizes @ sizes) / total)
    design_effect = 1.0 + (mean_influence - 1.0) * conservative_ucb
    effective_n = total / design_effect
    largest_share = float(sizes.max() / total)
    passed = (conservative_ucb < RHO_STAR[endpoint] and effective_n >= REQUIRED_N
              and largest_share <= MAX_SHARE)
    return {
        "profile": profile,
        "small_g_leave_one_cluster_out_ucb95": loo_ucbs,
        "conservative_rho_ucb95": conservative_ucb,
        "rho_star": RHO_STAR[endpoint],
        "rectangles": total, "clusters": len(groups),
        "largest_cluster_share": largest_share,
        "mean_cluster_influence": mean_influence,
        "design_effect_at_conservative_ucb": design_effect,
        "effective_n_at_conservative_ucb": effective_n,
        "replicate_supported_rectangles": sum(row["all_four_cells_replicated"]
                                               for row in endpoint_rows),
        "pass": passed,
    }


def run(output: Path = OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"no-clobber output exists: {output}")
    if sha256_file(ASSAY_ROWS) != ASSAY_SHA256:
        raise RuntimeError("exact-assay selection hash mismatch")
    selections = [row for row in read_jsonl(ASSAY_ROWS)
                  if row["selected_at_frozen_cap"] and row["eligible_primary"]]
    wanted = {int(value) for row in selections for ids in row["cell_activity_ids"] for value in ids}
    values, read_audit = load_values(wanted)
    dd_rows, noise = build_dd_rows(selections, values)
    endpoints = {endpoint: endpoint_dependence(dd_rows, endpoint) for endpoint in ("Ki", "Kd")}
    passing = [endpoint for endpoint, result in endpoints.items() if result["pass"]]
    verdict = ("X1A_R_DEPENDENCE_PRECONDITION_PASSED_" + "_AND_".join(passing).upper()
               if passing else "X1A_R_DEPENDENCE_PRECONDITION_FAILED")

    output.mkdir(parents=True)
    dd_path = output / "dd_rows.jsonl"
    with dd_path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in sorted(dd_rows, key=lambda item: (item["endpoint"], item["rectangle_id"])):
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    report = {
        "stage": "E-AFF-X1A-R_DIRECT_DD_DEPENDENCE",
        "terminal_verdict": verdict,
        "passing_endpoints": passing,
        "endpoints": endpoints,
        "noise": noise,
        "label_reads": {"chembl37": read_audit, "bindingdb": 0, "davis": 0,
                        "kiba": 0, "pdbbind": 0, "recipient": 0},
        "hashes": {"assay_selection": ASSAY_SHA256,
                   "dd_rows": sha256_file(dd_path),
                   "prereg": sha256_file(PREREG),
                   "x1b_prereg": sha256_file(X1B_PREREG),
                   "runner": sha256_file(Path(__file__))},
        "training_performed": False,
        "gpu_used": False,
        "x1b_execution_authorized_endpoints": passing,
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
        print(json.dumps({"terminal_verdict": "X1A_R_CONTRACT_INVALID",
                          "error": f"{type(exc).__name__}: {exc}"}, indent=2), file=sys.stderr)
        return 2
    print(json.dumps({"terminal_verdict": result["terminal_verdict"],
                      "passing_endpoints": result["passing_endpoints"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
