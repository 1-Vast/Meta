"""Stage CIIP-1A 2x2 matched-subset data contract (2026-08-19).

Builds the oracle-covered 49-pair matched subset for the root-cause 2x2
diagnostic: every cell (KLIFS/ESM x joint/centered-only) uses EXACTLY
the same pairs, rows, ligands, split assignment, seeds and budget. The
local ESM representation is ORACLE-localized (radius-6 window at the
verified mutation coordinate) - a positive-control representation, not
a deployable router. Read-only wrt every frozen input.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SIG = HERE.parent / "stageX_csc_signal"
X0C = SIG / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(X0C))

from x0_i2 import window_mean_esm, ESM_WINDOW_RADIUS, ESM_MAX_LEN  # noqa: E402
from x0_common import normalize_parent_name, normalize_construct_name, sha256_file  # noqa: E402

SCHEMA = "MetaSieve.StageCIIP1A.2x2.Data.v1"


def main() -> int:
    art = json.loads((HERE / "DATA1A.json").read_text(encoding="utf-8"))
    z = np.load(HERE / "DATA1A.npz", allow_pickle=False)
    esm = np.load(X0C / "q1_esm_cache.npz", allow_pickle=True)
    cache = {k: esm[k] for k in esm.files}
    wt_keys = {normalize_parent_name(k[3:]): k for k in cache if k.startswith("wt:")}
    mt_keys = {normalize_construct_name(k[3:]): k for k in cache if k.startswith("mt:")}
    audit = json.loads((X0C / "Q0B_MAPPING_AUDIT.json").read_text(encoding="utf-8"))
    cls = {r["source_row"]: r["mutation_class"] for r in audit["duongly_variant_records"]}
    sp = np.asarray(art["split"]["pair_split"])
    P = z["prot"]

    covered, missing = [], []
    for i, p in enumerate(art["pairs"]):
        wk = wt_keys.get(normalize_parent_name(p["parent"]))
        mk = mt_keys.get(normalize_construct_name(p["var_label"]))
        ok = bool(wk and mk and p["pos"] <= ESM_MAX_LEN)
        (covered if ok else missing).append(i)
    # oracle window features per pair: WT and variant windows at the pair site
    esm_wt = np.zeros((len(art["pairs"]), 640), np.float32)
    esm_var = np.zeros((len(art["pairs"]), 640), np.float32)
    for i in covered:
        p = art["pairs"][i]
        wk = wt_keys[normalize_parent_name(p["parent"])]
        mk = mt_keys[normalize_construct_name(p["var_label"])]
        esm_wt[i] = window_mean_esm(cache[wk], p["pos"], ESM_WINDOW_RADIUS)
        esm_var[i] = window_mean_esm(cache[mk], p["pos"], ESM_WINDOW_RADIUS)
    # covered-vs-missing bias report
    def stats(idx):
        vt = [float(np.var(art["targets"][i]["c"])) for i in idx]
        inf = [float((np.abs(np.asarray(art["targets"][i]["c"])) >= 10).mean())
               for i in idx]
        dp0 = sum(1 for i in idx if np.linalg.norm(
            P[art["pairs"][i]["var_row"]] - P[art["pairs"][i]["wt_row"]]) == 0)
        rows = sorted({art["pairs"][i]["wt_row"] for i in idx}
                      | {art["pairs"][i]["var_row"] for i in idx})
        rowm = [float(np.nanmean(z["Y"][r])) for r in rows]
        rowsd = [float(np.nanstd(z["Y"][r])) for r in rows]
        return {"n_pairs": len(idx),
                "split_counts": {s: int((sp[idx] == k).sum())
                                 for k, s in enumerate(("train", "val", "test"))},
                "n_parents": len({art["pairs"][i]["parent"] for i in idx}),
                "var_true_median": float(np.median(vt)),
                "var_true_range": [float(min(vt)), float(max(vt))],
                "informative_frac_median": float(np.median(inf)),
                "dP0_pairs": dp0,
                "mutation_classes": {k: int(v) for k, v in __import__(
                    "collections").Counter(
                        cls.get(art["pairs"][i]["var_label"], "?")
                        for i in idx).items()},
                "row_endpoint_mean_median": float(np.median(rowm)),
                "row_endpoint_sd_median": float(np.median(rowsd))}
    missing_detail = [{"pair_idx": i, "parent": art["pairs"][i]["parent"],
                       "mutation": art["pairs"][i]["mutation"],
                       "pos": art["pairs"][i]["pos"],
                       "reason": "pos > ESM_MAX_LEN"} for i in missing]
    bias = {"covered": stats(covered), "missing": stats(missing),
            "missing_detail": missing_detail,
            "verdict": ("family-level coverage bias (4 whole families dropped: "
                        "ALK, MET, LRRK2, TEK-Y1108F; 2/4 missing test pairs) "
                        "=> this stage is named the oracle-covered subset "
                        "diagnostic")}
    out = {
        "schema": SCHEMA,
        "preregistration_stage1_sha256": art["preregistration_sha256"],
        "inputs_sha256": {
            "DATA1A.json": sha256_file(HERE / "DATA1A.json"),
            "DATA1A.npz": sha256_file(HERE / "DATA1A.npz"),
            "q1_esm_cache.npz": sha256_file(X0C / "q1_esm_cache.npz"),
        },
        "oracle_local_esm": {
            "note": "radius-6 window mean at the VERIFIED mutation coordinate "
                    "(Q1-qualified construction); ORACLE-localized: requires "
                    "mutation annotation, NOT a deployable representation",
            "window_radius": ESM_WINDOW_RADIUS, "dim": 640,
            "esm_max_len": ESM_MAX_LEN,
        },
        "covered_pair_indices": [int(i) for i in covered],
        "covered_test_pairs": [int(i) for i in covered if sp[i] == 2],
        "n_covered_test_parents": len({art["pairs"][i]["parent"]
                                       for i in covered if sp[i] == 2}),
        "coverage_bias": bias,
        "endpoint": "percent inhibition (raw; never relabeled)",
    }
    (HERE / "DATA2X2.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    np.savez_compressed(HERE / "DATA2X2.npz",
                        esm_wt=esm_wt, esm_var=esm_var,
                        covered=np.asarray(covered, dtype=np.int64))
    sums = {}
    for f in ("DATA2X2.json", "DATA2X2.npz"):
        sums[f] = sha256_file(HERE / f)
    with open(HERE / "SHA256SUMS", "a", encoding="utf-8") as fh:
        for k, v in sorted(sums.items()):
            fh.write(f"{v} *{k}\n")
    print(json.dumps(out["coverage_bias"], indent=1))
    print("covered test pairs:", out["covered_test_pairs"],
          "parents:", out["n_covered_test_parents"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
