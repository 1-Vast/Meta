"""X1B crossed protein-by-ligand interaction existence audit.

Registered by
  research/crossed_interaction/PREREG_X1B_INTERACTION_EXISTENCE.md
  (sha256 fbcbc9c5..., commit pending at write time) committed BEFORE any DD
  value was computed.

Audit only: trains nothing, adds no module, introduces no parameter. Tests
interaction VARIANCE against replicate noise, never whether mean(DD) differs
from zero. Ki and Kd stay completely separate.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "crossed_interaction"
sys.path.insert(0, str(HERE))

from x1a_icc_audit import (  # noqa: E402
    CHEMBL_DB, ENDPOINTS, OUT, RECOVERED, X1ContractError, adjust_additive,
    git_head, read_jsonl, sha_file, write_json,
)

PREREG = HERE / "PREREG_X1B_INTERACTION_EXISTENCE.md"
PREREG_SHA = "fbcbc9c5d6ca959545d14fc4897bcce81614fe504f97d26fd262e6916414d595"

# ---- frozen contract (prereg sections 3, 4, 5, 7) --------------------------
VAR_REPLICATE = {"Ki": 0.38189950748311347, "Kd": 3.458710327272734}
X0B_UNITS = {"Ki": 11168, "Kd": 1041}
EFFECT_FLOOR = 1.0
GATE_I_REAL = 0.30
GATE_INR = 0.50
SEED_BOOT = 20260903
SEED_NULL = 20260904
N_BOOT = 10_000
N_NULL = 200


def component_bootstrap(per_cluster: dict, seed: int = SEED_BOOT,
                        n_boot: int = N_BOOT) -> dict:
    """Resample dependency clusters. Rectangles are never resampled as IID."""
    keys = sorted(per_cluster)
    if len(keys) < 2:
        return {"mean": float("nan"), "lcb95_one_sided": float("nan"),
                "clusters": len(keys)}
    values = np.array([per_cluster[k] for k in keys], dtype=np.float64)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), (n_boot, len(keys)))
    draws = values[idx].mean(1)
    return {"mean": float(values.mean()),
            "lcb95_one_sided": float(np.percentile(draws, 5)),
            "ci95": [float(np.percentile(draws, 2.5)),
                     float(np.percentile(draws, 97.5))],
            "clusters": len(keys)}


def pack_rectangles(cells_by_panel: dict) -> list:
    """Deterministic greedy cell-disjoint packing, X0-B's own rule."""
    units = []
    for panel in sorted(cells_by_panel):
        grid = cells_by_panel[panel]
        targets = sorted({t for t, _l in grid})
        ligands = sorted({l for _t, l in grid})
        used = set()
        for i in range(len(targets)):
            for j in range(i + 1, len(targets)):
                t1, t2 = targets[i], targets[j]
                for a in range(len(ligands)):
                    for b in range(a + 1, len(ligands)):
                        la, lb = ligands[a], ligands[b]
                        corners = ((t1, la), (t1, lb), (t2, la), (t2, lb))
                        if any(c not in grid for c in corners):
                            continue
                        if any(c in used for c in corners):
                            continue
                        used.update(corners)
                        units.append((panel, corners))
    return units


def evaluate(units, grid_value, grid_noise, cluster_of_panel, endpoint) -> dict:
    """DD, its noise budget and the derived statistics, aggregated by cluster."""
    by_cluster_dd2, by_cluster_noise, by_cluster_prev = (defaultdict(list),
                                                         defaultdict(list),
                                                         defaultdict(list))
    by_pair = defaultdict(list)
    for panel, corners in units:
        (t1, la), (t1b, lb), (t2, la2), (t2b, lb2) = corners
        dd = (grid_value[(t1, la)] - grid_value[(t1b, lb)]
              - grid_value[(t2, la2)] + grid_value[(t2b, lb2)])
        noise = sum(grid_noise[c] for c in corners)
        cluster = cluster_of_panel[panel]
        by_cluster_dd2[cluster].append(dd * dd)
        by_cluster_noise[cluster].append(noise)
        by_cluster_prev[cluster].append(1.0 if abs(dd) >= EFFECT_FLOOR else 0.0)
        by_pair[(panel, t1, t2)].append(dd)

    clusters = sorted(by_cluster_dd2)
    dd2 = {c: float(np.mean(by_cluster_dd2[c])) for c in clusters}
    noise = {c: float(np.mean(by_cluster_noise[c])) for c in clusters}
    prevalence = {c: float(np.mean(by_cluster_prev[c])) for c in clusters}
    i_real = {c: float(np.sqrt(max(0.0, dd2[c] - noise[c]))) for c in clusters}
    inr = {c: (i_real[c] / np.sqrt(noise[c]) if noise[c] > 0 else 0.0)
           for c in clusters}

    consistent, pairs_used = [], 0
    for values in by_pair.values():
        if len(values) < 4:
            continue
        arr = np.asarray(values)
        sign = np.sign(arr.mean())
        if sign == 0:
            continue
        pairs_used += 1
        consistent.append(float(np.mean(np.sign(arr) == sign)))

    return {
        "units": len(units),
        "mean_dd2": float(np.mean([dd2[c] for c in clusters])),
        "mean_noise": float(np.mean([noise[c] for c in clusters])),
        "i_real": float(np.mean([i_real[c] for c in clusters])),
        "inr": float(np.mean([inr[c] for c in clusters])),
        "prevalence_above_floor": float(np.mean([prevalence[c] for c in clusters])),
        "rank1_sign_consistency": (float(np.mean(consistent)) if consistent
                                   else float("nan")),
        "rank1_target_pairs": pairs_used,
        "per_cluster": {"i_real": i_real, "inr": inr, "prevalence": prevalence},
    }


def run(n_boot: int = N_BOOT, n_null: int = N_NULL) -> dict:
    started = time.time()
    if sha_file(PREREG) != PREREG_SHA:
        raise X1ContractError("X1B preregistration hash mismatch")
    x1a = json.loads((OUT / "X1A_ICC_AUDIT.json").read_text(encoding="utf-8"))
    if x1a["TERMINAL_VERDICT"] != "X1_ICC_PRECONDITION_PASSED":
        raise X1ContractError("X1A did not authorize X1B")

    cells = list(read_jsonl(RECOVERED / "eaff__x0_v1_cells.jsonl"))
    cluster_of_panel = {}
    for component in read_jsonl(RECOVERED / "eaff__x0_v1_dependency_components.jsonl"):
        if component.get("stratum") != "all_rectangles":
            continue
        for panel in component["panels"]:
            cluster_of_panel[(component["endpoint_family"], panel)] = \
                component["dependency_component_id"]

    import sqlite3
    wanted = set()
    for cell in cells:
        wanted.update(int(a) for a in cell["activity_ids"])
    connection = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    connection.execute("CREATE TEMP TABLE wanted(activity_id INTEGER PRIMARY KEY)")
    connection.executemany("INSERT OR IGNORE INTO wanted VALUES (?)",
                           ((a,) for a in wanted))
    values = {}
    for activity_id, relation, pchembl in connection.execute(
            "SELECT a.activity_id, a.standard_relation, a.pchembl_value "
            "FROM activities a JOIN wanted w ON w.activity_id = a.activity_id"):
        if relation == "=" and pchembl is not None:
            values[int(activity_id)] = float(pchembl)
    connection.close()

    per_endpoint = {}
    for endpoint in ENDPOINTS:
        var_rep = VAR_REPLICATE[endpoint]
        cells_by_panel = defaultdict(dict)
        flat = []
        for cell in cells:
            if cell["endpoint_family"] != endpoint:
                continue
            if (endpoint, cell["panel_id"]) not in cluster_of_panel:
                continue
            observed = [values[int(a)] for a in cell["activity_ids"]
                        if int(a) in values]
            if not observed:
                continue
            key = (cell["protein_sequence_sha256"], cell["ligand_connectivity_key"])
            cells_by_panel[cell["panel_id"]][key] = {
                "mean": float(np.mean(observed)), "n": len(observed)}
            flat.append({"panel": cell["panel_id"], "cluster": cluster_of_panel[
                (endpoint, cell["panel_id"])], "cell": key,
                "target": key[0], "ligand": key[1],
                "values": observed, "value_mean": float(np.mean(observed))})

        units = pack_rectangles({p: set(g) for p, g in cells_by_panel.items()})
        units = [(p, c) for p, c in units
                 if (endpoint, p) in cluster_of_panel]
        if not units:
            per_endpoint[endpoint] = {"status": "NO_UNITS"}
            continue

        grid_value = {}
        grid_noise = {}
        for panel, grid in cells_by_panel.items():
            for key, entry in grid.items():
                grid_value[key] = entry["mean"]
                grid_noise[key] = var_rep / entry["n"]
        panel_of_cluster = {p: cluster_of_panel[(endpoint, p)]
                            for p in cells_by_panel}
        observed_stats = evaluate(units, grid_value, grid_noise,
                                  panel_of_cluster, endpoint)

        # ---- additive null: identical design, values regenerated
        adjusted = adjust_additive(flat, within_panel=False)
        fitted = {}
        residual_by_cell = defaultdict(list)
        for row in adjusted:
            residual_by_cell[row["cell"]].append(row["r"])
        for row in flat:
            r = float(np.mean(residual_by_cell[row["cell"]]))
            fitted[row["cell"]] = row["value_mean"] - r
        rng = np.random.default_rng(SEED_NULL)
        null_i_real, null_inr, null_prev, null_consistency = [], [], [], []
        counts = {key: cells_by_panel[p][key]["n"]
                  for p, g in cells_by_panel.items() for key in g}
        for _ in range(n_null):
            null_value = {k: fitted.get(k, 0.0) +
                          rng.normal(0.0, np.sqrt(var_rep / counts[k]))
                          for k in grid_value}
            stats = evaluate(units, null_value, grid_noise, panel_of_cluster,
                             endpoint)
            null_i_real.append(stats["i_real"])
            null_inr.append(stats["inr"])
            null_prev.append(stats["prevalence_above_floor"])
            null_consistency.append(stats["rank1_sign_consistency"])

        null_i_real = np.asarray(null_i_real)
        null_p95 = float(np.percentile(null_i_real, 95))
        boot_i = component_bootstrap(observed_stats["per_cluster"]["i_real"],
                                     n_boot=n_boot)
        boot_inr = component_bootstrap(observed_stats["per_cluster"]["inr"],
                                       n_boot=n_boot)

        gate_a = bool(observed_stats["i_real"] >= GATE_I_REAL and
                      boot_i["lcb95_one_sided"] > 0)
        gate_b = bool(observed_stats["inr"] >= GATE_INR and
                      boot_inr["lcb95_one_sided"] > 0)
        gate_c = bool(observed_stats["i_real"] > null_p95)

        per_endpoint[endpoint] = {
            "status": "EVALUATED",
            "reconstructed_units": observed_stats["units"],
            "x0b_reported_units": X0B_UNITS[endpoint],
            "unit_reconstruction_ratio":
                observed_stats["units"] / X0B_UNITS[endpoint],
            "clusters": len(observed_stats["per_cluster"]["i_real"]),
            "mean_dd_squared": observed_stats["mean_dd2"],
            "mean_noise_budget": observed_stats["mean_noise"],
            "i_real": observed_stats["i_real"],
            "i_real_bootstrap": boot_i,
            "inr": observed_stats["inr"],
            "inr_bootstrap": boot_inr,
            "prevalence_above_1_log_unit": observed_stats["prevalence_above_floor"],
            "rank1_sign_consistency": observed_stats["rank1_sign_consistency"],
            "rank1_target_pairs": observed_stats["rank1_target_pairs"],
            "additive_null": {
                "replicates": int(null_i_real.size),
                "i_real_mean": float(null_i_real.mean()),
                "i_real_p95": null_p95,
                "inr_mean": float(np.mean(null_inr)),
                "prevalence_mean": float(np.mean(null_prev)),
                "rank1_sign_consistency_mean": float(np.nanmean(null_consistency)),
            },
            "gates": {
                "X1B_A_i_real": {"observed": observed_stats["i_real"],
                                 "required_at_least": GATE_I_REAL,
                                 "lcb95": boot_i["lcb95_one_sided"],
                                 "pass": gate_a},
                "X1B_B_inr": {"observed": observed_stats["inr"],
                              "required_at_least": GATE_INR,
                              "lcb95": boot_inr["lcb95_one_sided"],
                              "pass": gate_b},
                "X1B_C_above_additive_null": {
                    "observed": observed_stats["i_real"],
                    "null_p95": null_p95, "pass": gate_c},
            },
            "endpoint_pass": bool(gate_a and gate_b and gate_c),
        }
        print(f"[{endpoint}] units={observed_stats['units']} "
              f"(X0-B {X0B_UNITS[endpoint]}) I_real={observed_stats['i_real']:.4f} "
              f"INR={observed_stats['inr']:.4f} null_p95={null_p95:.4f} "
              f"pass={per_endpoint[endpoint]['endpoint_pass']}", flush=True)

    evaluated = [e for e in ENDPOINTS
                 if per_endpoint.get(e, {}).get("status") == "EVALUATED"]
    passing = [e for e in evaluated if per_endpoint[e]["endpoint_pass"]]
    verdict = ("REAL_CROSSED_AFFINITY_INTERACTION_IDENTIFIED" if passing
               else "REAL_CROSSED_AFFINITY_INTERACTION_NOT_IDENTIFIED")

    result = {
        "schema": "MetaSieve.CrossedInteraction.X1B.Gate.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "trains_nothing": True, "trainable_parameters_introduced": 0,
        "mean_dd_not_tested": "opposing selectivity effects cancel; only "
                              "interaction variance is tested",
        "endpoints": per_endpoint,
        "endpoints_passing": passing,
        "TERMINAL_VERDICT": verdict,
        "meta_learning_consequence": (
            "a real interaction makes the README target section a_t identifiable"
            if passing else
            "with interaction at or below noise the support set identifies only a "
            "per-target scalar offset, so the target section a_t is not "
            "identifiable from this source and X2 must not be trained"),
        "x2_authorized": bool(passing),
        "elapsed_seconds": round(time.time() - started, 1),
    }
    write_json(OUT / "X1B_INTERACTION_AUDIT.json", result)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    parser.add_argument("--n-null", type=int, default=N_NULL)
    args = parser.parse_args(argv)
    try:
        result = run(n_boot=args.n_boot, n_null=args.n_null)
        print(json.dumps({"TERMINAL_VERDICT": result["TERMINAL_VERDICT"],
                          "endpoints_passing": result["endpoints_passing"]},
                         indent=2), flush=True)
        return 0
    except Exception as exc:
        failure = {"schema": "MetaSieve.CrossedInteraction.X1B.FailClosed.v1",
                   "error_type": type(exc).__name__, "error": str(exc),
                   "TERMINAL_VERDICT": "X1B_CONTRACT_FAIL_CLOSED"}
        write_json(OUT / "X1B_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
