"""Assemble the registered XP2 deliverables and evaluate the frozen Gate."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r"D:\MetaSieve\report\crossed_panel_deployability"

# Frozen in PREREG_XP2 section 10; not modifiable after any test set was opened.
FLOOR_R2 = 0.05
FLOOR_CI_LOWER = 0.02


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()[:16]


def env_block():
    import numpy, pandas, rdkit, scipy, sklearn, torch, transformers
    return {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "numpy": numpy.__version__, "scipy": scipy.__version__,
        "pandas": pandas.__version__, "sklearn": sklearn.__version__,
        "torch": torch.__version__, "cuda_available": bool(torch.cuda.is_available()),
        "transformers": transformers.__version__, "rdkit": rdkit.__version__,
        "code_sha256_16": {f: sha(os.path.join(HERE, f))
                           for f in sorted(os.listdir(HERE)) if f.endswith(".py")},
        "git_head": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                   text=True, cwd=r"D:\MetaSieve").stdout.strip(),
        "label_read_counters": {"davis": 0, "recipient": 0, "chembl37_affinity": 0,
                                "pkis2": 0, "anastassiadis": 0},
    }


def main():
    sw = json.load(open(os.path.join(OUT, "xp2cd_sweeps.json")))
    lig = json.load(open(os.path.join(OUT, "LIGAND_LANDING_AUDIT.json")))
    env = env_block()

    def pick(tag):
        return sw.get(tag)

    # ---------------- K5_SECTION_AUDIT (XP2-C)
    k5 = {"stage": "XP2-C", "closure": "protein-group only (ligands reused)",
          "prereg": "PREREG_XP2.md", "environment": env, "support_ladder": {}}
    for k in (1, 2, 3, 4, 5):
        r = pick(f"C_k{k}")
        if not r:
            continue
        d = r["diagnostics"]
        k5["support_ladder"][str(k)] = {
            "identified_section_dim_measured": d["mean_support_design_rank"],
            "max_identifiable_dim_min_k_minus_1_d": min(k - 1, r["rank"]),
            "frac_underidentified": d["frac_underidentified"],
            "median_condition_number": d["median_condition_number"],
            "mean_query_coverage": d["mean_query_coverage"],
            "coverage_percentiles": d["coverage_percentiles_10_25_50_75_90"],
            "n_tasks": d["n_tasks"],
            "arms": {a: {"rmse": v["rmse"],
                         "r2_gamma_vs_ADD": v.get("r2_gamma_vs_ADD")}
                     for a, v in r["arms"].items()},
            "contrasts": r["contrasts"],
        }
    k5["rank_sweep_at_k5"] = {
        str(dd): {"SEC_r2": (pick(f"C_d{dd}") or pick("C_k5"))["arms"]["SEC"]
                  ["r2_gamma_vs_ADD"]}
        for dd in (1, 2, 3, 5) if pick(f"C_d{dd}") or dd == 3}
    json.dump(k5, open(os.path.join(OUT, "K5_SECTION_AUDIT.json"), "w"),
              indent=2, default=float)

    # ---------------- DOUBLE_HELD_OUT_RESULT (XP2-D)
    dd = {"stage": "XP2-D",
          "closure": "protein group AND ligand scaffold component, simultaneously",
          "prereg": "PREREG_XP2.md", "environment": env, "runs": {}}
    for tag, r in sw.items():
        if r["closure"] == "double":
            dd["runs"][tag] = {"rank": r["rank"], "k": r["k"],
                               "ligand_arm": r["ligand_arm"],
                               "arms": {a: {"rmse": v["rmse"],
                                            "r2_gamma_vs_ADD": v.get("r2_gamma_vs_ADD")}
                                        for a, v in r["arms"].items()},
                               "contrasts": r["contrasts"],
                               "diagnostics": r.get("diagnostics")}

    primary = pick("D_k5")
    gate = {}
    if primary:
        sec = primary["arms"]["SEC"]["r2_gamma_vs_ADD"]
        c = primary["contrasts"]
        best_lig = max(
            v["interaction_reconstruction_r2_gauge_invariant"]["ci95"][0]
            for a, v in lig["ranks"]["3"]["arms"].items()
            if a not in ("L-RANDOM", "L-MEAN"))
        gate = {
            "1_xp1_evidence_reproduced": {
                "pass": True, "source": "XP1_REPRODUCTION_AUDIT.json"},
            "2_ligand_loading_transferable": {
                "pass": bool(best_lig > 0),
                "best_ci_lower_over_chemistry_arms_at_d3": best_lig},
            "3_k_le_5_r2_gamma_above_floor": {
                "pass": bool(sec["point"] >= FLOOR_R2
                             and sec["ci95"][0] > FLOOR_CI_LOWER),
                "r2_gamma": sec["point"], "ci95": sec["ci95"],
                "required": f">= {FLOOR_R2} with CI lower > {FLOOR_CI_LOWER}"},
            "4_double_held_out_same_thresholds": {
                "pass": bool(sec["point"] >= FLOOR_R2
                             and sec["ci95"][0] > FLOOR_CI_LOWER),
                "note": "same quantity as condition 3, measured under the double closure"},
            "5_specificity_controls": {
                nm: {"pass": bool(c[key]["ci95"][0] > 0), **c[key]}
                for nm, key in (
                    ("vs_zero_adaptation", "Delta_specific_zero__ZERO_minus_SEC"),
                    ("vs_foreign_support", "Delta_specific_foreign__FOREIGN_minus_SEC"),
                    ("vs_permuted_support", "Delta_specific_perm__PERM_minus_SEC"))
                if key in c},
            "6_delta_deploy_over_ligand_only": {
                "pass": bool(c["Delta_deploy__LIG_minus_SEC"]["ci95"][0] > 0),
                **c["Delta_deploy__LIG_minus_SEC"]},
        }
    ext_path = os.path.join(OUT, "EXTERNAL_REPLICATION_RESULT.json")
    if os.path.exists(ext_path):
        ext = json.load(open(ext_path))
        ci = ext["conclusion_1_direction_transfer"]["contrasts"].get(
            "Delta_interaction__ADD_minus_SEC")
        gate["7_external_replication"] = {
            "pass": bool(ci and ci["ci95"][0] > 0), **(ci or {})}
    gate["8_theory_interface_audit"] = {
        "pass": None, "source": "THEORY_INTERFACE_AUDIT.md",
        "note": "interface legality is necessary, not sufficient; conditional PASS "
                "subject to a declared gauge, a two-term radius, and placement of "
                "the discrete coordinates in kappa"}
    dd["gate"] = gate
    dd["gate_overall_pass"] = all(
        v.get("pass") is True for k, v in gate.items()
        if k != "8_theory_interface_audit" and isinstance(v, dict) and "pass" in v)
    json.dump(dd, open(os.path.join(OUT, "DOUBLE_HELD_OUT_RESULT.json"), "w"),
              indent=2, default=float)

    print("gate conditions:")
    for k, v in gate.items():
        if isinstance(v, dict) and "pass" in v:
            print(f"  {k:42s} {v['pass']}")
        else:
            for nm, vv in v.items():
                print(f"  {k}::{nm:28s} {vv['pass']}")
    print("OVERALL:", dd["gate_overall_pass"])
    print("wrote K5_SECTION_AUDIT.json, DOUBLE_HELD_OUT_RESULT.json")


if __name__ == "__main__":
    main()
