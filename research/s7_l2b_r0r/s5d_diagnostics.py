"""S5D ligand-direction collapse and symmetric-difference conditional estimand.

Registered by
  research/s7_l2b_r0r/PREREG_PHASE2B_S5D_ESTIMAND_AND_COLLAPSE_DIAGNOSTICS.md
  (sha256 1ea639e0..., commit ba2390d) committed BEFORE this file existed.

Trains nothing. Reuses the frozen S4R views and checkpoints byte-for-byte and
adds no capacity, no representation and no arm that was not registered. The
S4R verdict `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED` is not reopened by
this stage; only the estimand and its metric are interrogated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "s7_l2b_r0r"
sys.path.insert(0, str(HERE))

from p2b_residue_residual import (  # noqa: E402
    SEED_BOOT, ap_exact, aggregate, component_bootstrap, project_np, sha_file,
)
import s4r_run as S4  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S5D_ESTIMAND_AND_COLLAPSE_DIAGNOSTICS.md"
PREREG_SHA = "1ea639e09a43c93e80e9969e09a45ca1bc7bb517fc040438e43b14db8d5386bb"

OUT = ROOT / "report" / "s7_l2b_r0r"
EXEC = S4.EXEC

# ---- frozen S5D contract (preregistration sections 4, 6) -------------------
MIN_PAIRS_FOR_RHO = 3
RHO_COLLAPSE_MEDIAN = 0.80
RHO_EXCESS_OVER_DATA = 0.10
E1_MARGIN = 0.05
E2_MARGIN = 0.03
E3_MARGIN = 0.03
ARMS = ("candidate", "baseline41", "foreign", "permuted")


class S5DContractError(RuntimeError):
    pass


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


def top_principal_energy_fraction(vectors: np.ndarray) -> float:
    """Share of mean-centred variance on the first principal direction of a set
    of unit vectors. 1.0 means every vector points the same way."""
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    keep = norms[:, 0] > 0
    if keep.sum() < MIN_PAIRS_FOR_RHO:
        return float("nan")
    unit = vectors[keep] / norms[keep]
    centred = unit - unit.mean(axis=0)
    singular = np.linalg.svd(centred, compute_uv=False)
    power = float((singular ** 2).sum())
    if power <= 0:
        return 1.0
    return float(singular[0] ** 2 / power)


def conditional_metrics(score: np.ndarray, gain_set, loss_set, length: int):
    """AP restricted to the residues that changed, plus its constant-score
    chance. Pocket membership is constant inside this comparison, so it cancels
    exactly. Eligibility matches the registered `ap_symdiff_conditional`."""
    gain = np.zeros(length, dtype=np.int8)
    loss = np.zeros(length, dtype=np.int8)
    for residue in gain_set:
        gain[residue] = 1
    for residue in loss_set:
        loss[residue] = 1
    idx = np.flatnonzero((gain + loss) > 0)
    if idx.size < 2:
        return None
    labels = gain[idx]
    positives = int(labels.sum())
    if positives == 0 or positives == idx.size:
        return None
    return {
        "ap_cond": ap_exact(np.asarray(score, dtype=np.float64)[idx], labels),
        "chance_cond": ap_exact(np.zeros(idx.size, dtype=np.float64), labels),
        "changed": int(idx.size),
        "gain_fraction": positives / float(idx.size),
    }


def run() -> dict:
    S4.require_absent([OUT / "PHASE2B_S5D_GATE.json",
                       OUT / "PHASE2B_S5D_REPORT.md"])
    if sha_file(PREREG) != PREREG_SHA:
        raise S5DContractError("S5D preregistration hash mismatch")
    s4r = json.loads((OUT / "PHASE2B_S4R_GATE.json").read_text())
    if s4r["TERMINAL_VERDICT"] != "REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED":
        raise S5DContractError("S5D expects the recorded S4R terminal verdict")

    ctx = S4.RuntimeContext("heldoutA")
    S4.validate_runtime_inputs(ctx, "heldoutA")
    if len(S4.LABEL_READ_LOG) != 1 or not S4.LABEL_READ_LOG[0].endswith(
            "heldoutA_residue_masks.json.gz"):
        raise S5DContractError(f"label-view firewall failed: {S4.LABEL_READ_LOG}")

    weights = {name: S4.load_head(name).W.detach().double().numpy()
               for name in ("candidate", "baseline41", "permuted")}
    ligand_table = {"candidate": ctx.gvec_graph, "baseline41": ctx.gvec,
                    "foreign": ctx.gvec_graph, "permuted": ctx.gvec_graph}
    arm_weight = {"candidate": weights["candidate"],
                  "baseline41": weights["baseline41"],
                  "foreign": weights["candidate"],
                  "permuted": weights["permuted"]}

    by_sk = defaultdict(list)
    for row in ctx.pairs:
        by_sk[row["sk"]].append(row)

    rho_rows = []
    conditional = {name: {} for name in ARMS}
    conditional_chance, pair_meta = {}, {}
    foreign_cosines, changed_sizes, gain_fractions = [], [], []

    for sk in sorted(by_sk):
        rows = by_sk[sk]
        q = ctx.Qs[sk]
        h = np.asarray(ctx.h(sk), dtype=np.float64)
        z = {name: h @ arm_weight[name] for name in ARMS}
        fields = {"candidate": [], "baseline41": [], "dg": []}
        for row in rows:
            a, b = row["a"], row["b"]
            pid = f"{a}|{b}"
            gk_a = ctx.records[a]["graph_key"]
            gk_b = ctx.records[b]["graph_key"]
            fk_a, fk_b = ctx.control["foreign_pair_map"][pid]
            length = int(ctx.records[a]["n_res"])
            gain, loss = S4.labels_for(ctx, a, b)

            scores = {}
            for name in ARMS:
                table = ligand_table[name]
                if name == "foreign":
                    delta = table[fk_a] - table[fk_b]
                else:
                    delta = table[gk_a] - table[gk_b]
                scores[name] = project_np(q, z[name] @ delta)

            fields["candidate"].append(scores["candidate"])
            fields["baseline41"].append(scores["baseline41"])
            fields["dg"].append(ctx.gvec_graph[gk_a] - ctx.gvec_graph[gk_b])
            left, right = scores["candidate"], scores["foreign"]
            denominator = np.linalg.norm(left) * np.linalg.norm(right)
            if denominator > 0:
                foreign_cosines.append(float(left @ right / denominator))

            reference = conditional_metrics(scores["candidate"], gain, loss, length)
            if reference is None:
                continue
            conditional_chance[pid] = reference["chance_cond"]
            pair_meta[pid] = sk
            changed_sizes.append(reference["changed"])
            gain_fractions.append(reference["gain_fraction"])
            for name in ARMS:
                metrics = conditional_metrics(scores[name], gain, loss, length)
                if metrics is None:
                    raise S5DContractError(
                        f"conditional eligibility differs for arm {name}")
                conditional[name][pid] = metrics["ap_cond"]

        if len(rows) >= MIN_PAIRS_FOR_RHO:
            rho_rows.append({
                "construct": sk,
                "component": ctx.construct_component[sk],
                "pairs": len(rows),
                "rho_dg": top_principal_energy_fraction(np.stack(fields["dg"])),
                "rho_graph": top_principal_energy_fraction(
                    np.stack(fields["candidate"])),
                "rho_base": top_principal_energy_fraction(
                    np.stack(fields["baseline41"])),
            })

    if not conditional["candidate"]:
        raise S5DContractError("empty conditional panel")
    reference_keys = set(conditional["candidate"])
    for name in ARMS:
        if set(conditional[name]) != reference_keys:
            raise S5DContractError(f"conditional mask differs for arm {name}")
    common_mask_sha = hashlib.sha256("\n".join(
        f"{pid}|{pair_meta[pid]}|{ctx.construct_component[pair_meta[pid]]}"
        for pid in sorted(reference_keys)).encode()).hexdigest()

    component, macro = {}, {}
    for name in ARMS:
        component[name], macro[name] = aggregate(
            conditional[name], pair_meta, ctx.construct_component)
    component["chance"], macro["chance"] = aggregate(
        conditional_chance, pair_meta, ctx.construct_component)

    def contrast(left, right, margin=None):
        result = component_bootstrap(component[left], component[right],
                                     n_boot=10000, seed=SEED_BOOT)
        if margin is None:
            result["gating"] = False
            return result
        result["margin"] = margin
        result["gating"] = True
        result["pass"] = bool(result["delta"] >= margin and
                              result["lcb95_one_sided"] > 0)
        return result

    gates = {
        "E1_vs_conditional_chance": contrast("candidate", "chance", E1_MARGIN),
        "E2_vs_foreign_ligand_pair": contrast("candidate", "foreign", E2_MARGIN),
        "E3_vs_trained_permuted_learner": contrast("candidate", "permuted", E3_MARGIN),
    }
    non_gating = {
        "E4_candidate_minus_baseline41": contrast("candidate", "baseline41"),
        "E5_baseline41_minus_chance": contrast("baseline41", "chance"),
    }

    rho_graph = np.array([r["rho_graph"] for r in rho_rows], dtype=np.float64)
    rho_dg = np.array([r["rho_dg"] for r in rho_rows], dtype=np.float64)
    rho_base = np.array([r["rho_base"] for r in rho_rows], dtype=np.float64)
    finite = np.isfinite(rho_graph) & np.isfinite(rho_dg) & np.isfinite(rho_base)
    d1 = {
        "constructs_with_at_least_3_pairs": int(finite.sum()),
        "median_rho_dg": float(np.median(rho_dg[finite])),
        "median_rho_graph": float(np.median(rho_graph[finite])),
        "median_rho_baseline41": float(np.median(rho_base[finite])),
        "median_rho_excess_over_data": float(
            np.median(rho_graph[finite] - rho_dg[finite])),
        "fraction_constructs_rho_graph_above_0p90": float(
            np.mean(rho_graph[finite] > 0.90)),
        "median_true_vs_foreign_field_cosine": float(np.median(foreign_cosines)),
        "pairs_with_foreign_cosine": len(foreign_cosines),
        "collapse_confirmed": None,
        "rule": (f"median rho_graph >= {RHO_COLLAPSE_MEDIAN} and "
                 f"median rho_graph >= median rho_dg + {RHO_EXCESS_OVER_DATA}"),
    }
    d1["collapse_confirmed"] = bool(
        d1["median_rho_graph"] >= RHO_COLLAPSE_MEDIAN and
        d1["median_rho_graph"] >= d1["median_rho_dg"] + RHO_EXCESS_OVER_DATA)

    if not d1["collapse_confirmed"]:
        verdict = "LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED"
    elif not all(gate["pass"] for gate in gates.values()):
        verdict = "POSE_FREE_LIGAND_CONDITIONED_DIRECTION_ABSENT_UNDER_CONDITIONAL_ESTIMAND"
    else:
        verdict = "CONDITIONAL_ESTIMAND_RECOVERS_LIGAND_SPECIFIC_DIRECTION_IN_DEVELOPMENT"
    gates_all_pass = all(gate["pass"] for gate in gates.values())

    result = {
        "schema": "MetaSieve.S7L2B.P2B.S5D.Gate.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "trains_nothing": True,
        "trainable_parameters_introduced": 0,
        "reuses_s4r_checkpoints": sorted(weights),
        "conditional_panel": {
            "pairs": len(reference_keys),
            "components": len(component["candidate"]),
            "common_mask_sha256": common_mask_sha,
            "median_changed_residues": float(np.median(changed_sizes)),
            "median_gain_fraction": float(np.median(gain_fractions)),
        },
        "primary_panel_pairs_before_conditional_eligibility": len(ctx.pairs),
        "D1_ligand_steering_collapse": d1,
        "D1_per_construct": rho_rows,
        "D2_macro_ap_conditional": {name: macro[name]
                                    for name in (*ARMS, "chance")},
        "gates": gates,
        "non_gating_contrasts": non_gating,
        "gates_all_pass": gates_all_pass,
        "s4r_verdict_unchanged": "REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED",
        "heldoutA_status": "ALREADY_CONSUMED_BY_S3R_AND_S4R_DEVELOPMENT_EVIDENCE_ONLY",
        "heldoutB_status": "NOT_CREATED_AND_NOT_READ",
        "R6_status": "NOT_RUN",
        "affinity_value_reads": 0,
        "label_view_paths_opened": list(S4.LABEL_READ_LOG),
        "TERMINAL_VERDICT": verdict,
        "claims_not_made": ["exact residue-atom coupling", "physical energy",
                            "affinity", "selectivity", "few-shot section",
                            "biological z", "validated end-to-end DTA",
                            "reopening of the ligand representation route"],
        "authorized_next_action": (
            "preregister a conditional-estimand confirmation on a panel that is "
            "not heldout-A"
            if verdict == "CONDITIONAL_ESTIMAND_RECOVERS_LIGAND_SPECIFIC_DIRECTION_IN_DEVELOPMENT"
            else "none; stop at the earliest failed boundary"),
    }
    S4.write_json(OUT / "PHASE2B_S5D_GATE.json", result)
    write_report(result)
    return result


def write_report(result) -> None:
    d1 = result["D1_ligand_steering_collapse"]
    lines = [
        "# Phase 2B S5D estimand and collapse diagnostics", "",
        f"Terminal verdict: `{result['TERMINAL_VERDICT']}`", "",
        "Trains nothing; reuses the frozen S4R checkpoints. The S4R verdict",
        f"`{result['s4r_verdict_unchanged']}` is unchanged by this stage.", "",
        "## D1 ligand-steering collapse", "",
        f"{d1['constructs_with_at_least_3_pairs']} heldout-A constructs have at "
        "least three eligible pairs.", "",
        "| quantity | median |", "|---|---:|",
        f"| `rho_dg`, unit ligand differences (data-side upper bound) | "
        f"{d1['median_rho_dg']:.4f} |",
        f"| `rho_graph`, candidate residue fields | {d1['median_rho_graph']:.4f} |",
        f"| `rho_base`, baseline41 residue fields | "
        f"{d1['median_rho_baseline41']:.4f} |",
        f"| `rho_graph - rho_dg` | {d1['median_rho_excess_over_data']:.4f} |",
        f"| true-vs-foreign field cosine | "
        f"{d1['median_true_vs_foreign_field_cosine']:.4f} |",
        "",
        f"Rule: {d1['rule']}. Collapse confirmed: "
        f"`{d1['collapse_confirmed']}`.", "",
        "## D2 symmetric-difference conditional estimand", "",
        f"{result['conditional_panel']['pairs']} eligible pairs across "
        f"{result['conditional_panel']['components']} closure components, from "
        f"{result['primary_panel_pairs_before_conditional_eligibility']} primary "
        "pairs. Median changed residues per pair "
        f"{result['conditional_panel']['median_changed_residues']:.1f}, median "
        f"gain fraction {result['conditional_panel']['median_gain_fraction']:.4f}.",
        "", "| arm | component-macro AP_cond |", "|---|---:|",
    ]
    lines += [f"| {name} | {value:.6f} |"
              for name, value in result["D2_macro_ap_conditional"].items()]
    lines += ["", "| Gate | delta | LCB95 | margin | PASS |",
              "|---|---:|---:|---:|:---:|"]
    lines += [f"| {name} | {gate['delta']:.6f} | {gate['lcb95_one_sided']:.6f} | "
              f"{gate['margin']:.2f} | {gate['pass']} |"
              for name, gate in result["gates"].items()]
    lines += ["", "| non-gating contrast | delta | LCB95 |", "|---|---:|---:|"]
    lines += [f"| {name} | {row['delta']:.6f} | {row['lcb95_one_sided']:.6f} |"
              for name, row in result["non_gating_contrasts"].items()]
    lines += ["",
              "Heldout-A was already consumed by S3R and S4R, so every number here",
              "is development evidence and none of it confirms anything. Heldout-B",
              "was neither created nor read, R6 was not opened, no affinity value",
              "was read and the frozen law operator was not modified.", ""]
    (OUT / "PHASE2B_S5D_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", nargs="?", default="run", choices=("run",))
    parser.parse_args(argv)
    started = time.time()
    try:
        result = run()
        print(json.dumps({"TERMINAL_VERDICT": result["TERMINAL_VERDICT"],
                          "D1": result["D1_ligand_steering_collapse"],
                          "D2": result["D2_macro_ap_conditional"],
                          "gates": result["gates"],
                          "non_gating": result["non_gating_contrasts"],
                          "elapsed_seconds": round(time.time() - started, 3)},
                         indent=2, default=str), flush=True)
        return 0
    except Exception as exc:
        failure = {
            "schema": "MetaSieve.S7L2B.P2B.S5D.FailClosed.v1",
            "created_utc": "2026-08-10",
            "error_type": type(exc).__name__, "error": str(exc),
            "TERMINAL_VERDICT": "S5D_CONTRACT_FAIL_CLOSED",
            "affinity_value_reads": 0,
        }
        S4.write_json(OUT / "PHASE2B_S5D_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
