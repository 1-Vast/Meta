"""X1A crossed-interaction ICC precondition and data-contract audit.

Registered by
  research/crossed_interaction/PREREG_X1_ICC_AND_DATA_CONTRACT.md
  (sha256 67e3c651..., commit 008c82a) committed BEFORE any ChEMBL37 value
  was read.

Trains nothing. Reads only activity_id, standard_relation, standard_type and
pchembl_value for activity ids already enumerated by the label-blind X0 census.
Ki and Kd are analysed completely separately and are never merged. BindingDB,
DAVIS, KIBA, PDBbind and recipient labels are not opened.
"""
from __future__ import annotations

import argparse
import json
import sqlite3
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "crossed_interaction"
RECOVERED = HERE / "recovered"
OUT = ROOT / "report" / "crossed_interaction"

PREREG = HERE / "PREREG_X1_ICC_AND_DATA_CONTRACT.md"
PREREG_SHA = "67e3c651de8d3f932934d76fb955f4554ded668551200ca439253d0548549bbc"
CHEMBL_DB = (ROOT / "dataset" / "raw" / "source_affinity" / "chembl37_sqlite_v1" /
             "extracted" / "chembl_37" / "chembl_37_sqlite" / "chembl_37.db")

CELLS_SHA = "898df88235401a2be2341ae1ab222e6c5903202796c8312d8e9091cf76741562"
COMPONENTS_SHA = "8970d059bcbf4dd1ca0f3b7ec9e5ab3f594740ec37bcb5c173fdb66ffb8a7779"
PANELS_SHA = "f378cdd610205fc02850e84e5530830d2255fcf539f4ecca81518eb19bf036bc"

# ---- frozen contract (prereg sections 3, 5, 6) -----------------------------
ENDPOINTS = ("Ki", "Kd")
RHO_STAR = {"Ki": 0.0915, "Kd": 0.0164}
X0B_UNITS = {"Ki": 11168, "Kd": 1041}
X0B_CLUSTERS = {"Ki": 36, "Kd": 12}
X0B_CAP_AT_RHO_STAR = {"Ki": 32, "Kd": 125}
REQUIRED_EFFECTIVE_N = 245
MAX_CLUSTER_SHARE = 0.25
SEED_BOOT = 20260903
N_BOOT = 10_000


class X1ContractError(RuntimeError):
    pass


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str),
                    encoding="utf-8")


def sha_file(path: Path) -> str:
    import hashlib
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


def read_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


# --------------------------------------------------------------- variance
def one_way_components(groups: list[np.ndarray]) -> tuple[float, float, float]:
    """Unbalanced one-way random-effects moment estimator (Henderson III).

    Returns (between variance, within variance, effective group size). Negative
    between-variance is truncated at zero, as registered.
    """
    groups = [g for g in groups if g.size >= 1]
    k = len(groups)
    sizes = np.array([g.size for g in groups], dtype=np.float64)
    total = sizes.sum()
    if k < 2 or total <= k:
        return 0.0, float(np.var(np.concatenate(groups), ddof=1)) if total > 1 else 0.0, 0.0
    means = np.array([g.mean() for g in groups])
    grand = float(np.concatenate(groups).mean())
    ss_within = float(sum(((g - g.mean()) ** 2).sum() for g in groups))
    ss_between = float((sizes * (means - grand) ** 2).sum())
    ms_within = ss_within / (total - k)
    ms_between = ss_between / (k - 1)
    n0 = (total - (sizes ** 2).sum() / total) / (k - 1)
    between = (ms_between - ms_within) / n0 if n0 > 0 else 0.0
    return max(0.0, float(between)), float(ms_within), float(n0)


def nested_rho(records: list[dict]) -> dict:
    """Nested cluster / panel / cell / replicate decomposition of the
    additively-adjusted residual, per prereg section 5."""
    by_cell = defaultdict(list)
    for row in records:
        by_cell[(row["cluster"], row["panel"], row["cell"])].append(row["r"])

    replicate_groups = [np.asarray(v, dtype=np.float64)
                        for v in by_cell.values() if len(v) >= 2]
    if replicate_groups:
        _between_cell, var_e, _n0 = one_way_components(replicate_groups)
    else:
        var_e = 0.0

    cell_mean = {key: float(np.mean(v)) for key, v in by_cell.items()}
    by_panel = defaultdict(list)
    for (cluster, panel, _cell), value in cell_mean.items():
        by_panel[(cluster, panel)].append(value)
    var_w, _resid_panel, _ = one_way_components(
        [np.asarray(v) for v in by_panel.values() if len(v) >= 2])

    panel_mean = {key: float(np.mean(v)) for key, v in by_panel.items()}
    by_cluster = defaultdict(list)
    for (cluster, _panel), value in panel_mean.items():
        by_cluster[cluster].append(value)
    var_v, _resid_cluster, _ = one_way_components(
        [np.asarray(v) for v in by_cluster.values() if len(v) >= 2])

    cluster_mean = np.asarray(sorted(float(np.mean(v)) for v in by_cluster.values()))
    var_u = float(np.var(cluster_mean, ddof=1)) if cluster_mean.size > 1 else 0.0

    total = var_u + var_v + var_w + var_e
    return {"var_cluster": var_u, "var_panel": var_v, "var_cell": var_w,
            "var_replicate": var_e, "total": total,
            "rho": (var_u / total) if total > 0 else float("nan"),
            "cells": len(by_cell), "panels": len(by_panel),
            "clusters": len(by_cluster),
            "replicate_supported_cells": len(replicate_groups)}


def adjust_additive(cells: list[dict]) -> list[dict]:
    """Remove additive target and ligand effects within each panel by one
    least-squares two-way fit. DD cancels exactly these effects, so the residual
    is the part of the measurement that DD inherits."""
    out = []
    by_panel = defaultdict(list)
    for cell in cells:
        by_panel[cell["panel"]].append(cell)
    for panel, rows in by_panel.items():
        targets = sorted({r["target"] for r in rows})
        ligands = sorted({r["ligand"] for r in rows})
        if len(targets) < 2 or len(ligands) < 2:
            continue
        t_index = {t: i for i, t in enumerate(targets)}
        l_index = {l: i for i, l in enumerate(ligands)}
        design = np.zeros((len(rows), 1 + len(targets) + len(ligands)))
        y = np.zeros(len(rows))
        for i, row in enumerate(rows):
            design[i, 0] = 1.0
            design[i, 1 + t_index[row["target"]]] = 1.0
            design[i, 1 + len(targets) + l_index[row["ligand"]]] = 1.0
            y[i] = row["value_mean"]
        coef, *_ = np.linalg.lstsq(design, y, rcond=None)
        fitted = design @ coef
        for i, row in enumerate(rows):
            for value in row["values"]:
                out.append({"cluster": row["cluster"], "panel": panel,
                            "cell": row["cell"], "r": float(value - fitted[i])})
    return out


def cluster_bootstrap_ucb(records: list[dict], seed: int = SEED_BOOT,
                          n_boot: int = N_BOOT) -> dict:
    """One-sided 95% upper bound by resampling dependency clusters, the frozen
    inference unit. Rectangles, cells and measurements are never resampled."""
    by_cluster = defaultdict(list)
    for row in records:
        by_cluster[row["cluster"]].append(row)
    clusters = sorted(by_cluster)
    if len(clusters) < 2:
        return {"ucb95": float("nan"), "draws": 0, "clusters": len(clusters)}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        picked = rng.integers(0, len(clusters), len(clusters))
        resampled = []
        for offset, index in enumerate(picked):
            source = clusters[index]
            for row in by_cluster[source]:
                resampled.append({**row, "cluster": f"{source}#{offset}"})
        value = nested_rho(resampled)["rho"]
        if np.isfinite(value):
            draws.append(value)
    if not draws:
        return {"ucb95": float("nan"), "draws": 0, "clusters": len(clusters)}
    array = np.asarray(draws)
    return {"ucb95": float(np.percentile(array, 95)),
            "lcb95": float(np.percentile(array, 5)),
            "median": float(np.median(array)),
            "draws": int(array.size), "clusters": len(clusters)}


# --------------------------------------------------------------- run
def run(n_boot: int = N_BOOT) -> dict:
    started = time.time()
    if sha_file(PREREG) != PREREG_SHA:
        raise X1ContractError("X1 preregistration hash mismatch")
    files = {
        "cells": (RECOVERED / "eaff__x0_v1_cells.jsonl", CELLS_SHA),
        "dependency_components": (RECOVERED / "eaff__x0_v1_dependency_components.jsonl",
                                  COMPONENTS_SHA),
        "panels": (RECOVERED / "eaff__x0_v1_panels.jsonl", PANELS_SHA),
    }
    recovered_hashes = {}
    for name, (path, expected) in files.items():
        actual = sha_file(path)
        recovered_hashes[name] = actual
        if actual != expected:
            raise X1ContractError(f"recovered {name} hash mismatch")
    if not CHEMBL_DB.is_file():
        raise X1ContractError("pinned ChEMBL37 archive is missing")

    cells = list(read_jsonl(files["cells"][0]))
    cluster_of_panel = {}
    for component in read_jsonl(files["dependency_components"][0]):
        if component.get("stratum") != "all_rectangles":
            continue
        for panel in component["panels"]:
            cluster_of_panel[(component["endpoint_family"], panel)] = \
                component["dependency_component_id"]

    wanted = set()
    for cell in cells:
        wanted.update(int(a) for a in cell["activity_ids"])
    print(f"cells={len(cells)} activity ids requested={len(wanted)}", flush=True)

    # ---- the first and only read of affinity values, per prereg section 4
    connection = sqlite3.connect(f"file:{CHEMBL_DB}?mode=ro", uri=True)
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("CREATE TEMP TABLE wanted(activity_id INTEGER PRIMARY KEY)")
    connection.executemany("INSERT OR IGNORE INTO wanted VALUES (?)",
                           ((a,) for a in wanted))
    cursor = connection.execute(
        "SELECT a.activity_id, a.standard_relation, a.standard_type, a.pchembl_value "
        "FROM activities a JOIN wanted w ON w.activity_id = a.activity_id")
    values, relation_census, type_census = {}, Counter(), Counter()
    rows_returned = 0
    for activity_id, relation, standard_type, pchembl in cursor:
        rows_returned += 1
        relation_census[relation if relation is not None else "NULL"] += 1
        type_census[standard_type if standard_type is not None else "NULL"] += 1
        if relation != "=" or pchembl is None:
            continue
        values[int(activity_id)] = float(pchembl)
    connection.close()
    print(f"activities returned={rows_returned} usable pchembl '='={len(values)}",
          flush=True)

    per_endpoint = {}
    for endpoint in ENDPOINTS:
        usable = []
        dropped = Counter()
        for cell in cells:
            if cell["endpoint_family"] != endpoint:
                continue
            cluster = cluster_of_panel.get((endpoint, cell["panel_id"]))
            if cluster is None:
                dropped["no_dependency_cluster"] += 1
                continue
            observed = [values[int(a)] for a in cell["activity_ids"]
                        if int(a) in values]
            if not observed:
                dropped["no_usable_value"] += 1
                continue
            usable.append({
                "cluster": cluster, "panel": cell["panel_id"],
                "cell": (cell["panel_id"], cell["protein_sequence_sha256"],
                         cell["ligand_connectivity_key"]),
                "target": cell["protein_sequence_sha256"],
                "ligand": cell["ligand_connectivity_key"],
                "values": observed, "value_mean": float(np.mean(observed)),
                "replicates": len(observed),
            })
        if not usable:
            per_endpoint[endpoint] = {"status": "NO_USABLE_CELLS",
                                      "dropped": dict(dropped)}
            continue

        adjusted = adjust_additive(usable)
        if not adjusted:
            per_endpoint[endpoint] = {"status": "NO_CROSSED_PANEL",
                                      "dropped": dict(dropped)}
            continue
        components = nested_rho(adjusted)
        bootstrap = cluster_bootstrap_ucb(adjusted, n_boot=n_boot)

        rho_point = components["rho"]
        rho_star = RHO_STAR[endpoint]
        cap = X0B_CAP_AT_RHO_STAR[endpoint]
        sizes = Counter(row["cluster"] for row in adjusted)
        capped = {k: min(v, cap) for k, v in sizes.items()}
        capped_total = sum(capped.values())
        largest_capped_share = (max(capped.values()) / capped_total
                                if capped_total else float("nan"))
        largest_uncapped_share = (max(sizes.values()) / sum(sizes.values())
                                  if sizes else float("nan"))
        mean_influence = capped_total / max(len(capped), 1)
        design_effect = 1.0 + (mean_influence - 1.0) * rho_point
        n_eff = capped_total / design_effect if design_effect > 0 else 0.0

        gate_icc = bool(np.isfinite(bootstrap["ucb95"]) and
                        bootstrap["ucb95"] < rho_star)
        gate_domination = bool(np.isfinite(largest_capped_share) and
                               largest_capped_share <= MAX_CLUSTER_SHARE)
        gate_power = bool(n_eff >= REQUIRED_EFFECTIVE_N)
        per_endpoint[endpoint] = {
            "status": "EVALUATED",
            "cells_with_values": len(usable),
            "measurements": sum(len(r["values"]) for r in usable),
            "replicate_supported_cells": sum(1 for r in usable if r["replicates"] >= 2),
            "adjusted_measurements": len(adjusted),
            "dropped": dict(dropped),
            "variance_components": components,
            "rho_point": rho_point,
            "rho_bootstrap": bootstrap,
            "rho_star": rho_star,
            "x0b_units": X0B_UNITS[endpoint],
            "x0b_clusters": X0B_CLUSTERS[endpoint],
            "cap": cap,
            "clusters_present": len(sizes),
            "largest_cluster_share_uncapped": largest_uncapped_share,
            "largest_cluster_share_capped": largest_capped_share,
            "mean_capped_cluster_influence": mean_influence,
            "design_effect": design_effect,
            "effective_n": n_eff,
            "gates": {
                "icc_upper_bound_below_rho_star": {
                    "observed_ucb95": bootstrap["ucb95"],
                    "required_below": rho_star, "pass": gate_icc},
                "no_cluster_dominates": {
                    "observed_capped_share": largest_capped_share,
                    "required_at_most": MAX_CLUSTER_SHARE, "pass": gate_domination},
                "effective_power": {
                    "observed_effective_n": n_eff,
                    "required_at_least": REQUIRED_EFFECTIVE_N, "pass": gate_power},
            },
            "endpoint_pass": bool(gate_icc and gate_domination and gate_power),
        }
        print(f"[{endpoint}] rho={rho_point:.4f} ucb95={bootstrap['ucb95']:.4f} "
              f"(rho*={rho_star}) capped_share={largest_capped_share:.4f} "
              f"n_eff={n_eff:.1f} pass={per_endpoint[endpoint]['endpoint_pass']}",
              flush=True)

    evaluated = [e for e in ENDPOINTS
                 if per_endpoint[e].get("status") == "EVALUATED"]
    passing = [e for e in evaluated if per_endpoint[e]["endpoint_pass"]]
    if not evaluated:
        verdict = "X1_DATA_CONTRACT_INVALID"
    elif not passing:
        verdict = "X1_INTERACTION_UNDERDETERMINED"
    else:
        verdict = "X1_ICC_PRECONDITION_PASSED"

    result = {
        "schema": "MetaSieve.CrossedInteraction.X1A.Gate.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "recovered_x0_hashes": recovered_hashes,
        "chembl37_sha256_pinned":
            "4be13df3b68e25dcd0bff44bf094033b5aebe98f415acdc8c1cdf380e0c15142",
        "trains_nothing": True,
        "trainable_parameters_introduced": 0,
        "label_reads": {
            "chembl37_activity_rows_returned": rows_returned,
            "chembl37_usable_pchembl_equals": len(values),
            "fields_read": ["activity_id", "standard_relation", "standard_type",
                            "pchembl_value"],
            "bindingdb": 0, "davis": 0, "kiba": 0, "pdbbind": 0, "recipient": 0,
        },
        "relation_census": dict(relation_census),
        "standard_type_census": dict(type_census.most_common(12)),
        "endpoints": per_endpoint,
        "endpoints_passing": passing,
        "TERMINAL_VERDICT": verdict,
        "x1b_authorized_for": passing,
        "authorized_next_action": (
            f"preregister X1B for {passing}" if passing else
            "none; stop at the earliest failed precondition"),
        "elapsed_seconds": round(time.time() - started, 1),
    }
    write_json(OUT / "X1A_ICC_AUDIT.json", result)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    args = parser.parse_args(argv)
    try:
        result = run(n_boot=args.n_boot)
        print(json.dumps({"TERMINAL_VERDICT": result["TERMINAL_VERDICT"],
                          "endpoints_passing": result["endpoints_passing"]},
                         indent=2), flush=True)
        return 0
    except Exception as exc:
        failure = {"schema": "MetaSieve.CrossedInteraction.X1A.FailClosed.v1",
                   "error_type": type(exc).__name__, "error": str(exc),
                   "TERMINAL_VERDICT": "X1_DATA_CONTRACT_INVALID"}
        write_json(OUT / "X1A_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
