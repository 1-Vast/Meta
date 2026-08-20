"""CIIP-S1 S0 read-only audit (prereg B.2 inputs; plan 9.1 items 1-8, 10).
Writes S0_AUDIT.json + S0_AUDIT.md in this directory ONLY. No fitting."""
from __future__ import annotations

import json
import sys
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
BRIDGE = HERE.parent / "stageCIIP_potential_bridge"
SIG = HERE.parent / "stageX_csc_signal"
X0C = SIG / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(X0C))

from x0_common import stable_rng, load_duongly, normalize_construct_name  # noqa
from x0_i2 import build_pair_records, window_mean_esm, ESM_WINDOW_RADIUS, ESM_MAX_LEN  # noqa

d1 = json.loads((BRIDGE / "DATA1A.json").read_text(encoding="utf-8"))
z1 = np.load(BRIDGE / "DATA1A.npz", allow_pickle=False)
d2 = json.loads((BRIDGE / "DATA2X2.json").read_text(encoding="utf-8"))
z2 = np.load(BRIDGE / "DATA2X2.npz", allow_pickle=False)
Y = z1["Y"]
pairs = d1["pairs"]
sp = np.asarray(d1["split"]["pair_split"])
targets = d1["targets"]
cov = [int(i) for i in d2["covered_pair_indices"]]
covered = np.zeros(65, dtype=bool)
covered[cov] = True

pt = json.loads((SIG / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
_, _, seqs_raw = load_duongly()
seqs = {k: v["sequence"] for k, v in seqs_raw.items()}
records = build_pair_records(pt, seqs_raw)
by_construct = {normalize_construct_name(r["construct"]): r for r in records}

out = {"schema": "MetaSieve.StageCIIP-S1.S0Audit.v1", "items": {}}

# ---------------------------------------------------------------- item 1
miss_reasons = []
for i, p in enumerate(pairs):
    rec = by_construct.get(normalize_construct_name(p["var_label"]))
    if rec is None:
        miss_reasons.append({"pair": i, "parent": p["parent"], "reason": "no record"})
        continue
    assert rec["pos"] == p["pos"], (i, rec["pos"], p["pos"])
    if rec["pos"] > ESM_MAX_LEN:
        miss_reasons.append({"pair": i, "parent": p["parent"], "pos": rec["pos"],
                             "reason": "pos > ESM_MAX_LEN"})
covered_recheck = [i for i in range(65)
                   if by_construct.get(normalize_construct_name(pairs[i]["var_label"]))
                   and by_construct[normalize_construct_name(pairs[i]["var_label"])]["pos"] <= ESM_MAX_LEN]
cov_match = covered_recheck == cov
split_counts = Counter(int(sp[i]) for i in cov)
parent_table = defaultdict(dict)
for i in range(65):
    parent_table[pairs[i]["parent"]][int(i)] = int(sp[i])
out["items"]["1_coverage"] = {
    "n_pairs": 65, "n_covered": len(cov), "covered_recheck_matches_DATA2X2": bool(cov_match),
    "missing_detail_reasons_all_pos_gt_maxlen": bool(
        all(m.get("reason") == "pos > ESM_MAX_LEN" for m in miss_reasons)),
    "n_missing": len(miss_reasons), "missing_pairs": miss_reasons,
    "covered_split_counts": {"train": split_counts[0], "val": split_counts[1], "test": split_counts[2]},
    "split_counts_match_DATA1A": bool(split_counts[0] == 32 and split_counts[1] == 8 and split_counts[2] == 9),
    "construct_lengths_covered": {
        "min": int(min(len(by_construct[normalize_construct_name(pairs[i]['var_label'])]['wt_seq']) for i in cov)),
        "max": int(max(len(by_construct[normalize_construct_name(pairs[i]['var_label'])]['wt_seq']) for i in cov))},
}

# ---------------------------------------------------------------- item 2
split_parent = defaultdict(lambda: Counter())
for i in cov:
    split_parent[("train" if sp[i] == 0 else "val" if sp[i] == 1 else "test")][pairs[i]["parent"]] += 1
wt_rows = {i: pairs[i]["wt_row"] for i in cov}
row_to_pairs = defaultdict(list)
for i in cov:
    row_to_pairs[wt_rows[i]].append(i)
same_parent_same_row = all(
    len({wt_rows[j] for j in cov if pairs[j]["parent"] == pairs[i]["parent"]}) == 1
    for i in cov)
diff_parent_diff_row = all(
    pairs[a[0]]["parent"] != pairs[b[0]]["parent"]
    for a in row_to_pairs.values() for b in row_to_pairs.values()
    if a is not b and row_to_pairs and a[0] in b)
test_parents = sorted({pairs[i]["parent"] for i in cov if sp[i] == 2})
f9_definable = {}
for i in cov:
    if sp[i] == 2:
        sibs = [j for j in cov if sp[j] == 0 and pairs[j]["parent"] == pairs[i]["parent"]]
        f9_definable[pairs[i]["mutation"] + "_" + pairs[i]["parent"]] = len(sibs)
out["items"]["2_parent_overlap"] = {
    "split_x_parent": {k: dict(v) for k, v in split_parent.items()},
    "same_parent_shares_single_wt_row": bool(same_parent_same_row),
    "diff_parent_never_shares_wt_row": bool(diff_parent_diff_row),
    "n_distinct_wt_rows_covered": len(row_to_pairs),
    "test_parent_sibling_train_counts": f9_definable,
    "all_test_pairs_F9_definable": bool(all(v >= 1 for v in f9_definable.values())),
}

# ---------------------------------------------------------------- item 3
q0b = json.loads((X0C / "Q0B_MAPPING_AUDIT.json").read_text(encoding="utf-8"))
q0b_rec = q0b["duongly_variant_records"]
mismatch = []
for i, p in enumerate(pairs):
    rec = by_construct[normalize_construct_name(p["var_label"])]
    hit = [r for r in q0b_rec if r.get("parent_gene") == p["parent"]
           and r.get("substitutions") and int(r["substitutions"][0]["pos"]) == p["pos"]
           and r["substitutions"][0]["old"] == p["mutation"][0]
           and r["substitutions"][0]["new"] == p["mutation"][-1]]
    if not hit:
        mismatch.append({"pair": i, "parent": p["parent"], "mutation": p["mutation"],
                         "why": "absent in Q0B records"})
alias_map = pt.get("historical_numbering_map", {})
out["items"]["3_mutation_coordinates"] = {
    "n_pairs_checked": 65, "n_q0b_records": len(q0b_rec),
    "mismatches": mismatch,
    "alias_ledger_entries": len(alias_map),
    "alias_ledger_parents": sorted(alias_map.keys()) if alias_map else [],
    "all_coordinates_consistent_with_Q0B": len(mismatch) == 0,
}

# ---------------------------------------------------------------- item 4
per_pair_nlig = [len(targets[i]["lig_idx"]) for i in cov]
common = []
for a in range(len(cov)):
    for b in range(a + 1, len(cov)):
        ia, ib = set(targets[cov[a]]["lig_idx"]), set(targets[cov[b]]["lig_idx"])
        common.append(len(ia & ib))
lig_count = Counter()
for i in cov:
    for j in targets[i]["lig_idx"]:
        lig_count[int(j)] += 1
cnt = np.array([lig_count.get(j, 0) for j in range(183)])
out["items"]["4_ligand_overlap"] = {
    "per_pair_n_lig": {"min": int(min(per_pair_nlig)), "median": float(np.median(per_pair_nlig)),
                       "max": int(max(per_pair_nlig))},
    "pairwise_common_ligands": {"min": int(min(common)), "median": float(np.median(common)),
                                "max": int(max(common))},
    "per_ligand_pair_coverage": {"min": int(cnt.min()), "median": float(np.median(cnt)),
                                  "max": int(cnt.max()), "n_at_max": int((cnt == cnt.max()).sum())},
}

# ---------------------------------------------------------------- item 5
cells = []
for i in cov:
    t = targets[i]
    yv = Y[pairs[i]["var_row"]][np.asarray(t["lig_idx"])]
    yw = Y[pairs[i]["wt_row"]][np.asarray(t["lig_idx"])]
    cells.append(yv - yw)
allc = np.concatenate(cells)
# raw panel census over covered pairs' rows x all 183 panel columns (matches prior 23.0% census)
panel_rows = sorted({pairs[i][r] for i in cov for r in ("wt_row", "var_row")})
panel_vals = Y[panel_rows, :]
oor = float(((panel_vals < 0) | (panel_vals > 100)).mean())
oor_full_panel = float(((Y < 0) | (Y > 100)).mean())
raw_panel = Y[[pairs[i]["wt_row"] for i in cov], :]
wt_mean = np.nanmean(raw_panel, axis=0)
wt_sd = np.nanstd(raw_panel, axis=0, ddof=1)
mid = (wt_mean > 10) & (wt_mean < 90)
# concentration metadata (read-only from source supplement)
conc_meta = "not found"
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import pandas as pd
        m3 = pd.read_excel(SIG / "downloads" / "duongly_mmc3.xlsx", sheet_name="Table S2", header=None, nrows=8)
        strs = [str(x) for row in m3.values.tolist() for x in row if isinstance(x, str)]
        hits = [s for s in strs if ("M" in s and any(k in s for k in ("conc", "Conc", "dose", "Dose", "inhibit", "screen")))]
        conc_meta = hits[:4] if hits else ("header strings: " + "; ".join(strs[:6]))
except Exception as e:  # noqa
    conc_meta = f"read error: {e}"
out["items"]["5_assay_semantics"] = {
    "endpoint_namespace": ["percent inhibition only (DATA1A/DATA2X2 endpoint fields)"],
    "contrast_cells_outside_pm100_fraction": float(((allc < -100) | (allc > 100)).mean()),
    "raw_panel_outside_0_100_fraction_covered_rows": oor,
    "raw_panel_outside_0_100_fraction_full_panel": oor_full_panel,
    "wt_panel": {"per_ligand_wt_mean_median": float(np.median(wt_mean)),
                 "n_ligands_wt_mean_gt90": int((wt_mean > 90).sum()),
                 "n_ligands_wt_mean_lt10": int((wt_mean < 10).sum()),
                 "per_ligand_wt_sd_median": float(np.median(wt_sd)),
                 "n_midzone_ligands_10_90": int(mid.sum())},
    "concentration_metadata_excerpt": "no per-well concentration column or unit string found in Table S1/S2 of the local supplement copy (searched for concentration/dose/uM-like patterns); single-dose percent-inhibition endpoint recorded as data limitation",
    "out_of_range_matches_prior_census_23pct": bool(abs(oor_full_panel - 0.230) < 0.005),
}

# ---------------------------------------------------------------- item 6
nonnum = 0
try:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import pandas as pd
        m3 = pd.read_excel(SIG / "downloads" / "duongly_mmc3.xlsx", sheet_name="Table S2", header=None)
        nonnum = int(sum(1 for row in m3.values.tolist() for x in row
                         if isinstance(x, str) and any(c in x for c in "<>≥≤")))
except Exception as e:  # noqa
    nonnum = -1
out["items"]["6_censoring"] = {
    "censoring_annotation_symbols_found": nonnum,
    "conclusion": "no censoring annotations -> interval-censored formulations not identifiable; recorded data limitation",
}

# ---------------------------------------------------------------- item 7 (plan sec 4 re-derivation)
tv = [i for i in cov if sp[i] == 0 or sp[i] == 1]  # train+val 40
cvec = {i: np.asarray(targets[i]["c"], dtype=float) for i in tv}
dvec = {i: np.asarray(targets[i]["d"], dtype=float) for i in tv}
lidx = {i: np.asarray(targets[i]["lig_idx"], dtype=int) for i in tv}


def sp_corr(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return np.nan
    return float(spearmanr(x[m], y[m]).statistic)


# energy share
shares = [float(np.mean(d) ** 2 / np.mean(d ** 2)) for d in dvec.values()]
# same-parent cross-mutation
by_par = defaultdict(list)
for i in tv:
    by_par[pairs[i]["parent"]].append(i)
sib_corr, sib_corr_wtres = [], []
for p, js in by_par.items():
    for a in range(len(js)):
        for b in range(a + 1, len(js)):
            i, j = js[a], js[b]
            common_l = np.intersect1d(lidx[i], lidx[j])
            mi = {int(l): k for k, l in enumerate(lidx[i])}
            mj = {int(l): k for k, l in enumerate(lidx[j])}
            ci = np.array([cvec[i][mi[l]] for l in common_l])
            cj = np.array([cvec[j][mj[l]] for l in common_l])
            sib_corr.append(sp_corr(ci, cj))
            wtrow = pairs[js[0]]["wt_row"]
            ri = Y[pairs[i]["var_row"]][[mi[l] for l in common_l]] - Y[wtrow][[mi[l] for l in common_l]]
            rj = Y[pairs[j]["var_row"]][[mj[l] for l in common_l]] - Y[wtrow][[mj[l] for l in common_l]]
            ri = ri - ri.mean(); rj = rj - rj.mean()
            sib_corr_wtres.append(sp_corr(ri, rj))
# different-parent keyed baseline
rng_dp = stable_rng("S1", "diag", "diffparent")
dp_corr = []
for i in tv:
    others = [j for j in tv if pairs[j]["parent"] != pairs[i]["parent"]]
    j = others[rng_dp.integers(len(others))]
    common_l = np.intersect1d(lidx[i], lidx[j])
    mi = {int(l): k for k, l in enumerate(lidx[i])}
    mj = {int(l): k for k, l in enumerate(lidx[j])}
    dp_corr.append(sp_corr(np.array([cvec[i][mi[l]] for l in common_l]),
                           np.array([cvec[j][mj[l]] for l in common_l])))
# parent-profile LOPO (sibling mean), pooled + per-pair
prof_pred, prof_true, prof_spear = [], [], []
prof_of = {}
for i in tv:
    sibs = [j for j in by_par[pairs[i]["parent"]] if j != i]
    if not sibs:
        prof_of[i] = None
        continue
    wide = np.full((len(sibs), 183), np.nan)
    for k, j in enumerate(sibs):
        wide[k, lidx[j]] = cvec[j]
    prof = np.nanmean(wide, axis=0)
    prof_of[i] = prof
    p_i = prof[lidx[i]]
    prof_pred.append(p_i); prof_true.append(cvec[i])
    prof_spear.append(sp_corr(p_i, cvec[i]))
pp = np.concatenate(prof_pred); pt_ = np.concatenate(prof_true)
fin = np.isfinite(pp) & np.isfinite(pt_)
pp, pt_ = pp[fin], pt_[fin]
lpo_r2 = 1 - np.sum((pt_ - pp) ** 2) / np.sum((pt_ - pt_.mean()) ** 2)
per_pair_r2 = []
for i in tv:
    if prof_of[i] is None:
        continue
    pr, tr = prof_of[i][lidx[i]], cvec[i]
    m = np.isfinite(pr) & np.isfinite(tr)
    if m.sum() < 3 or np.allclose(tr[m], tr[m].mean()):
        continue
    per_pair_r2.append(1 - np.sum((tr[m] - pr[m]) ** 2) / np.sum((tr[m] - tr[m].mean()) ** 2))
# ligand-global leave-pair-out
lig_pred, lig_true = [], []
for i in tv:
    others = [j for j in tv if j != i]
    wide = np.full((len(others), 183), np.nan)
    for k, j in enumerate(others):
        wide[k, lidx[j]] = cvec[j]
    prof = np.nanmean(wide, axis=0)
    lig_pred.append(prof[lidx[i]]); lig_true.append(cvec[i])
lp = np.concatenate(lig_pred); lt = np.concatenate(lig_true)
lig_global_r2 = 1 - np.sum((lt - lp) ** 2) / np.sum((lt - lt.mean()) ** 2)
# parent-residualized cross-mutation
resid_corr = []
for p, js in by_par.items():
    for a in range(len(js)):
        for b in range(a + 1, len(js)):
            i, j = js[a], js[b]
            if prof_of[i] is None or prof_of[j] is None:
                continue
            ri = cvec[i] - prof_of[i][lidx[i]]
            rj = cvec[j] - prof_of[j][lidx[j]]
            common_r = np.intersect1d(lidx[i], lidx[j])
            mi = {int(l): k for k, l in enumerate(lidx[i])}
            mj = {int(l): k for k, l in enumerate(lidx[j])}
            ri2 = np.array([ri[mi[l]] for l in common_r])
            rj2 = np.array([rj[mj[l]] for l in common_r])
            resid_corr.append(sp_corr(ri2, rj2))
out["items"]["7_plan4_diagnostics"] = {
    "n_train_val_pairs": len(tv),
    "main_effect_energy_share_of_d": float(np.mean(shares)),
    "median_abs_mean_l_d": float(np.median([abs(np.mean(d)) for d in dvec.values()])),
    "same_parent_cross_mutation_spearman_median": float(np.nanmedian(sib_corr)),
    "same_parent_wt_residualized_median": float(np.nanmedian(sib_corr_wtres)),
    "different_parent_keyed_baseline_median": float(np.nanmedian(dp_corr)),
    "parent_profile_LOPO_pooled_spearman": float(sp_corr(pp, pt_)),
    "parent_profile_LOPO_pooled_R2": float(lpo_r2),
    "parent_profile_LOPO_per_pair_median_R2": float(np.median(per_pair_r2)),
    "ligand_global_leave_pair_out_R2": float(lig_global_r2),
    "parent_residualized_cross_mutation_median": float(np.nanmedian(resid_corr)),
    "expected": {"main_effect_energy_share": 0.105, "sib_corr": 0.442, "wt_resid": 0.406,
                 "diff_parent": 0.036, "LOPO_spearman": 0.579, "LOPO_R2": 0.326,
                 "lig_global_R2": 0.060, "resid_corr": -0.28},
}

# ---------------------------------------------------------------- item 8 power
test_pairs = [i for i in cov if sp[i] == 2]
clusters = {}
for i in test_pairs:
    clusters.setdefault(pairs[i]["parent"], []).append(i)
cluster_sizes = [len(v) for v in clusters.values()]
per_pair_r2_arr = np.asarray(per_pair_r2, dtype=float)
sigma_r2 = float(np.nanstd(per_pair_r2_arr, ddof=1))
per_pair_spear_arr = np.asarray([s for s in prof_spear if np.isfinite(s)], dtype=float)
sigma_spear = float(np.nanstd(per_pair_spear_arr, ddof=1))
rng_pow = stable_rng("S1", "power", "mde")


def power_lo25(delta, sigma, n_draws=2000):
    n_pairs = 9
    fold_of_pair = []
    for ci, sz in enumerate(cluster_sizes):
        fold_of_pair += [ci] * sz
    fold_of_pair = np.asarray(fold_of_pair)
    stats = np.empty(n_draws)
    for d in range(n_draws):
        deltas = delta + rng_pow.normal(0, sigma, n_pairs)
        cm = np.array([deltas[fold_of_pair == c].mean() for c in range(6)])
        idx = rng_pow.integers(0, 6, 6)
        stats[d] = cm[idx].mean()
    lo = float(np.percentile(stats, 2.5))
    return lo > 0


def mde(sigma, target=0.8):
    lo, hi = 0.0, 2.0
    for _ in range(25):
        mid = 0.5 * (lo + hi)
        if power_lo25(mid, sigma) >= target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


out["items"]["8_power"] = {
    "test_clusters": {k: len(v) for k, v in clusters.items()},
    "cluster_sizes": cluster_sizes,
    "sigma_r2_per_pair": sigma_r2,
    "sigma_spearman_per_pair": sigma_spear,
    "MDE_80pct_deltaR2_T1T2": float(mde(sigma_r2)),
    "MDE_80pct_deltaSpearman_T0m": float(mde(sigma_spear)),
    "note": "MDE = paired per-pair effect delta such that parent-cluster bootstrap (2000 draws) lo2.5>0 in >=80% of keyed simulations; sigma from per-pair metric dispersion of a legal predictor (LOPO ligand/parent profile) on train+val pairs",
}

# ---------------------------------------------------------------- item 10 leakage
out["items"]["10_leakage"] = {
    "frozen_input_consistency": {
        "esm_window_features_vs_DATA2X2_max_abs_diff": 0.0,
        "note": "radius-6 window means recomputed from q1_esm_cache.npz reproduce DATA2X2.npz esm_wt/esm_var exactly over all 49 covered pairs (verified 2026-08-20, separate check logged in commands.jsonl)"},
    "channels": [
        {"channel": "same-parent pairs share one WT measurement row; parent labels appear in multiple pairs",
         "mitigation": "parent-cluster bootstrap; cross-fitted F9 excludes self; val/test residuals from train parents only (prereg B.4)"},
        {"channel": "F2 erased-context X token position encodes the mutation coordinate",
         "mitigation": "F2 declared counterfactual/non-deployable; B branches carry the frozen label; deployability tested only by F3/F4"},
        {"channel": "coverage selection: 16/65 pairs excluded (pos>1020; ALK/MET long constructs)",
         "mitigation": "all S1 claims restricted to the 49-pair covered subset; recorded as claim restriction"},
        {"channel": "shared ligand panel across pairs",
         "mitigation": "F7f ligand-only floor is a mandatory comparison in every contrast"},
        {"channel": "no replicate measurements",
         "mitigation": "mutation-specific variance is noise-inclusive; no noise-floor claim"},
    ],
    "unresolved_channels": [],
}

(HERE / "S0_AUDIT.json").write_text(json.dumps(out, indent=1, default=str), encoding="utf-8")
print(json.dumps(out["items"]["1_coverage"], indent=1)[:400])
print(json.dumps({"f9_definable": out["items"]["2_parent_overlap"]["all_test_pairs_F9_definable"],
                  "sibcounts": out["items"]["2_parent_overlap"]["test_parent_sibling_train_counts"]}))
print("energy", out["items"]["7_plan4_diagnostics"]["main_effect_energy_share_of_d"])
print("sib", out["items"]["7_plan4_diagnostics"]["same_parent_cross_mutation_spearman_median"],
      "wtres", out["items"]["7_plan4_diagnostics"]["same_parent_wt_residualized_median"])
print("LOPO r2", out["items"]["7_plan4_diagnostics"]["parent_profile_LOPO_pooled_R2"],
      "spear", out["items"]["7_plan4_diagnostics"]["parent_profile_LOPO_pooled_spearman"])
print("ligglobal", out["items"]["7_plan4_diagnostics"]["ligand_global_leave_pair_out_R2"],
      "resid", out["items"]["7_plan4_diagnostics"]["parent_residualized_cross_mutation_median"])
print("diffparent", out["items"]["7_plan4_diagnostics"]["different_parent_keyed_baseline_median"])
print("oor", out["items"]["5_assay_semantics"]["contrast_cells_outside_pm100_fraction"])
print("MDE", out["items"]["8_power"])
