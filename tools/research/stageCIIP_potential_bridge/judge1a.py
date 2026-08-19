"""Stage CIIP-1A gate adjudication (frozen thresholds, prereg 31d3eeaf...).

Reads RESULT_SCREENING.json or RESULT_FULL.json; computes the frozen
gates: pair-mean Spearman >= 0.30, dead-zone sign accuracy >= 0.65,
bootstrap 2.5% lower bounds of (unified_local - family_shuffle) >= 0.05
and (unified_local - ligand_only) >= 0.05 (cluster = parent, 2000 draws,
SHA-256 keyed). Free pairwise gap > 0.10 -> expression-insufficiency
note. Screening mode: reports the same quantities for seed 1 only and
the structural-readiness statement; it does not decide PASS (full
3-seed median does).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import train1a as T  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", default=str(HERE / "RESULT_SCREENING.json"))
    ap.add_argument("--mode", choices=["screening", "full"], default="screening")
    args = ap.parse_args()
    res = json.loads(Path(args.result).read_text(encoding="utf-8"))
    arms = res["arms"]
    seeds = sorted({int(s) for a in arms.values() for s in a})
    ul_rows = [r for s in seeds for r in arms["unified_local"][str(s)]["test_rows"]]
    lig_rows = [r for s in seeds for r in arms["ligand_only"][str(s)]["test_rows"]]
    fs_rows = [r for s in seeds for r in arms["family_shuffle"][str(s)]["test_rows"]]
    fp_rows = [r for s in seeds for r in arms["free_pairwise"][str(s)]["test_rows"]]
    sp = float(np.nanmean([r["spearman"] for r in ul_rows]))
    sa = float(np.nanmean([r["sign_acc"] for r in ul_rows]))
    gap_fs = T.bootstrap_gap(ul_rows, fs_rows, "spearman", "screen_v_family")
    gap_lig = T.bootstrap_gap(ul_rows, lig_rows, "spearman", "screen_v_lig")
    gap_fp = float(np.nanmean([r["spearman"] for r in fp_rows]) - sp)
    out = {
        "schema": "MetaSieve.StageCIIP1A.Gate.v1",
        "mode": args.mode, "seeds": seeds,
        "unified_local": {"spearman": sp, "sign_acc": sa,
                          "pearson": float(np.nanmean(
                              [r["pearson"] for r in ul_rows])),
                          "mse": float(np.mean([r["mse"] for r in ul_rows])),
                          "scale_median": float(np.nanmedian(
                              [r["scale"] for r in ul_rows]))},
        "controls": {
            "ligand_only_spearman": float(np.nanmean(
                [r["spearman"] for r in lig_rows])),
            "family_shuffle_spearman": float(np.nanmean(
                [r["spearman"] for r in fs_rows])),
            "free_pairwise_spearman": float(np.nanmean(
                [r["spearman"] for r in fp_rows])),
        },
        "gates": {
            "spearman_ge_0.30": sp >= 0.30,
            "sign_acc_ge_0.65": sa >= 0.65,
            "gap_vs_family_lo2.5_ge_0.05": gap_fs >= 0.05,
            "gap_vs_ligand_lo2.5_ge_0.05": gap_lig >= 0.05,
        },
        "bootstrap_lo2.5": {"vs_family_shuffle": gap_fs,
                            "vs_ligand_only": gap_lig},
        "free_pairwise_expression_note": (
            "free_pairwise exceeds unified by >0.10 -> pairwise signal may "
            "exist but the integrable potential expression is insufficient"
            if gap_fp > 0.10 else "no expression-insufficiency flag"),
        "gap_free_pairwise_vs_unified": gap_fp,
    }
    if args.mode == "full":
        out["PASS"] = all(out["gates"].values())
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
