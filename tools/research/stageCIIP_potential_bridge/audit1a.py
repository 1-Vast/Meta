"""Stage CIIP-1A READ-ONLY collapse audit (2026-08-19).

Diagnoses why the unified potential collapsed (3/13 nonconstant test
pairs, pair-mean Spearman -0.0409) WITHOUT retraining: output collapse
table, potential variance sources, gradient competition at
initialization, KLIFS-vs-local-ESM representation comparison, and the
true meaning of the free-pairwise diagnostic. Writes ONLY the audit
artifacts; every frozen input is opened read-only and pinned by SHA.
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

import torch  # noqa: E402

import potential as POT  # noqa: E402
import train1a as T  # noqa: E402
from x0_i2 import window_mean_esm, ESM_WINDOW_RADIUS, ESM_MAX_LEN  # noqa: E402
from x0_common import normalize_parent_name, normalize_construct_name  # noqa: E402

DEAD_ZONE = 10.0
SCHEMA = "MetaSieve.StageCIIP1A.CollapseAudit.v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_frozen():
    art = json.loads((HERE / "DATA1A.json").read_text(encoding="utf-8"))
    res = json.loads((HERE / "RESULT_SCREENING.json").read_text(encoding="utf-8"))
    z = np.load(HERE / "DATA1A.npz", allow_pickle=False)
    q1 = json.loads((X0C / "Q1_SELECTIVITY.json").read_text(encoding="utf-8"))
    esm = np.load(X0C / "q1_esm_cache.npz", allow_pickle=True)
    return art, res, z, q1, esm


def collapse_table(art, res):
    te = [i for i, s in enumerate(art["split"]["pair_split"]) if s == 2]
    table = {}
    for arm in ("unified_local", "free_pairwise", "ligand_only"):
        rows = res["arms"][arm]["1"]["test_rows"]
        by_pair = {r["pair"]: r for r in rows}
        recs = []
        for i in te:
            r = by_pair[i]
            c = np.asarray(art["targets"][i]["c"])
            var_true = float(c.var())
            mse = r["mse"]
            scale = r["scale"]
            nonconst = scale == scale and abs(scale) > 0  # var(pred) > 0
            recs.append({
                "pair": i, "parent": r["parent"], "mutation": r["mutation"],
                "n_lig": r["n"],
                "var_true": var_true,
                "mse": mse,
                "r2": float(1 - mse / var_true) if var_true > 0 else None,
                "ols_slope": scale if nonconst else None,
                "sign_acc": r["sign_acc"],
                "spearman": r["spearman"] if r["spearman"] == r["spearman"] else None,
                "pearson": r["pearson"] if r["pearson"] == r["pearson"] else None,
                "N_target_informative": int((np.abs(c) >= DEAD_ZONE).sum()),
                "prediction_nonconstant": bool(nonconst),
                "rank_evaluable": bool(nonconst),
            })
        table[arm] = recs
    counts = {}
    for arm, recs in table.items():
        nonc = sum(1 for r in recs if r["prediction_nonconstant"])
        counts[arm] = {
            "N_total": len(recs),
            "N_target_informative_median": float(np.median(
                [r["N_target_informative"] for r in recs])),
            "N_prediction_nonconstant": nonc,
            "N_rank_evaluable": nonc,
            "finite_pairs/total_pairs": f"{nonc}/{len(recs)}",
        }
    return table, counts


def erank(singular_values):
    sv = np.asarray(singular_values, dtype=np.float64)
    if sv.sum() == 0:
        return 0.0
    p = sv / sv.sum()
    return float(1.0 / (p @ p))


def klifs_diff_matrix(art, z):
    P = z["prot"]
    D = np.stack([P[p["var_row"]] - P[p["wt_row"]] for p in art["pairs"]])
    norms = np.linalg.norm(D, axis=1)
    zero = norms == 0
    split = np.asarray(art["split"]["pair_split"])
    zero_idx = [int(i) for i in np.where(zero)[0]]
    u, s, _ = np.linalg.svd(D, full_matrices=False)
    out = {"shape": list(D.shape), "norm_min": float(norms.min()),
           "norm_max": float(norms.max()),
           "n_zero_difference_pairs": int(zero.sum()),
           "zero_difference_pairs": zero_idx,
           "zero_difference_by_split": {
               "train": int((zero & (split == 0)).sum()),
               "val": int((zero & (split == 1)).sum()),
               "test": int((zero & (split == 2)).sum())},
           "nonzero_norms_are_sqrt2_only": bool(
               np.allclose(norms[~zero], np.sqrt(2.0))),
           "singular_values": [float(x) for x in s[:12]],
           "effective_rank": erank(s), "rank": int((s > 1e-6).sum())}
    return out


def esm_diff_matrix(art, esm_cache):
    """Local-ESM mutation delta at the verified site (radius-6 window),
    same construction as Q1 pair_centered_local_esm (read-only)."""
    cache = {k: esm_cache[k] for k in esm_cache.files}
    wt_keys = {}
    mt_keys = {}
    for k in cache:
        if k.startswith("wt:"):
            wt_keys[normalize_parent_name(k[3:])] = k
        elif k.startswith("mt:"):
            mt_keys[normalize_construct_name(k[3:])] = k
    rows, meta = [], []
    for p in art["pairs"]:
        wk = wt_keys.get(normalize_parent_name(p["parent"]))
        mk = mt_keys.get(normalize_construct_name(p["var_label"]))
        if wk and mk and p["pos"] <= ESM_MAX_LEN:
            dw = window_mean_esm(cache[wk], p["pos"], ESM_WINDOW_RADIUS)
            dm = window_mean_esm(cache[mk], p["pos"], ESM_WINDOW_RADIUS)
            rows.append((dm - dw).ravel())
            meta.append({"pair": p["parent"] + " " + p["mutation"],
                         "norm": float(np.linalg.norm(dm - dw))})
        else:
            meta.append({"pair": p["parent"] + " " + p["mutation"],
                         "norm": None})
    if rows:
        D = np.stack(rows)
        u, s, _ = np.linalg.svd(D, full_matrices=False)
        norms = np.asarray([m["norm"] for m in meta if m["norm"] is not None])
        return {"n_pairs_in_cache": len(rows), "n_pairs_total": len(art["pairs"]),
                "shape": list(D.shape),
                "norm_min": float(norms.min()), "norm_median": float(
                    np.median(norms)), "norm_max": float(norms.max()),
                "singular_values": [float(x) for x in s[:12]],
                "effective_rank": erank(s), "rank": int((s > 1e-6).sum()),
                "per_pair": meta}
    return {"n_pairs_in_cache": 0}


def q1_evidence(q1):
    out = {}
    task_a = q1.get("results", {}).get("task_A", {})
    for name in ("pair_centered_local_esm", "klifs_pocket",
                 "global_pooled_esm", "edit_descriptor"):
        v = task_a.get(name, "NOT_IN_Q1")
        if isinstance(v, dict):
            out[name] = {
                "selectivity": v.get("selectivity"),
                "bootstrap_ci_lo": (v.get("bootstrap") or {}).get("ci_lo"),
                "bootstrap_ci_hi": (v.get("bootstrap") or {}).get("ci_hi"),
                "n_parents": v.get("n_parents"),
                "note": v.get("note"),
            }
        else:
            out[name] = v
    out["q1_pass"] = q1.get("q1_pass")
    out["frozen_pass_rule"] = q1.get("frozen_pass_rule")
    return out


def init_model_analysis(art, z, seed=1):
    torch.manual_seed(seed)
    m = POT.UnifiedPotential()
    P = torch.from_numpy(z["prot"]).float()
    L = torch.from_numpy(z["lig"]).float()
    pairs = art["pairs"]
    with torch.no_grad():
        # Cov(psi) over all 183 ligands at init
        psi = m.psi(torch.relu(m.l_enc(L))).numpy()  # (183, 8)
        cov = np.cov(psi, rowvar=False)
        s_cov = np.linalg.svd(cov, compute_uv=False)
        alpha_wt = m.alpha(torch.relu(m.p_enc(P[[p["wt_row"] for p in pairs]]))).numpy()
        alpha_v = m.alpha(torch.relu(m.p_enc(P[[p["var_row"] for p in pairs]]))).numpy()
        da = alpha_v - alpha_wt  # (65, 8)
        da_norms = np.linalg.norm(da, axis=1)
        pred_var = np.einsum("ik,kl,il->i", da, cov, da)  # da_i^T cov da_i
        # potential output scale at init (sample of cells)
        rows = np.arange(97)
        ligs = np.arange(183)
        rr = np.repeat(rows, min(8, 183))
        ll = np.tile(ligs[:min(8, 183)], 97)
        s_out = m.potential(P[rr], L[ll]).numpy()
        f_out = m(P[rr], L[ll]).numpy()
        return {
            "cov_psi_trace": float(np.trace(cov)),
            "cov_psi_singular_values": [float(x) for x in s_cov],
            "cov_psi_effective_rank": erank(s_cov),
            "delta_alpha_norm_min": float(da_norms.min()),
            "delta_alpha_norm_median": float(np.median(da_norms)),
            "delta_alpha_norm_max": float(da_norms.max()),
            "n_delta_alpha_near_zero": int((da_norms < 1e-6).sum()),
            "pred_var_per_pair_min/median/max": [
                float(np.percentile(pred_var, 0)), float(np.median(pred_var)),
                float(np.percentile(pred_var, 100))],
            "potential_output_var_at_init": float(s_out.var()),
            "f_output_var_at_init": float(f_out.var()),
            "note": "at INITIALIZATION (seed-1 init; the training trajectory "
                    "was not persisted, so no 10/25/50/75/100% checkpoints exist)",
        }


def gradient_audit(art, z, seed=1):
    """g_abs vs g_ctr on the interaction params at initialization,
    reproducing the trainer's epoch-0 first batch (same rng streams)."""
    torch.manual_seed(seed)
    m = POT.UnifiedPotential()
    pairs = art["pairs"]
    targets = art["targets"]
    tr_idx = [i for i, s in enumerate(art["split"]["pair_split"]) if s == 0]
    va_te_rows = ({pairs[i]["wt_row"] for i, s in enumerate(
        art["split"]["pair_split"]) if s != 0}
        | {pairs[i]["var_row"] for i, s in enumerate(
            art["split"]["pair_split"]) if s != 0})
    abs_rows = sorted({r for i in tr_idx
                       for r in (pairs[i]["wt_row"], pairs[i]["var_row"])}
                      - va_te_rows)
    P = torch.from_numpy(z["prot"]).float()
    L = torch.from_numpy(z["lig"]).float()
    Y = torch.from_numpy(z["Y"]).float()
    cells = [(i, j) for i in tr_idx for j in range(len(targets[i]["c"]))]
    rng = T.stable_rng("stageCIIP1A", "order", "unified_local", seed, "epoch", 0)
    order = np.asarray(cells)[rng.permutation(len(cells))][:T.BATCH]
    # contrast loss with gradient
    hats = {}
    for i in {i for i, _ in order}:
        t = targets[i]
        Lm = L[t["lig_idx"]]
        Pw = P[pairs[i]["wt_row"]:pairs[i]["wt_row"] + 1].expand(len(t["lig_idx"]), -1)
        Pv = P[pairs[i]["var_row"]:pairs[i]["var_row"] + 1].expand(len(t["lig_idx"]), -1)
        hats[i] = m.centered_mutation_effect(Pw, Pv, Lm)
    hi = torch.stack([hats[i][j] for i, j in order])
    ci = torch.tensor([targets[i]["c"][j] for i, j in order], dtype=torch.float32)
    loss_ctr = ((hi - ci) ** 2).mean()
    loss_ctr.backward()
    g_ctr = {n: p.grad.detach().clone() for n, p in m.named_parameters() if p.grad is not None}
    m.zero_grad()
    # absolute loss (same frozen row mask + rng stream as trainer epoch 0)
    rng_a = T.stable_rng("stageCIIP1A", "abs", "unified_local", seed, "epoch", 0)
    ra = np.asarray(abs_rows)[rng_a.permutation(len(abs_rows))[:T.BATCH]]
    la = rng_a.integers(0, Y.shape[1], size=len(ra))
    fin = np.isfinite(Y[ra, la])
    ra, la = ra[fin], la[fin]
    f = m(P[ra], L[la])
    loss_abs = ((f - Y[ra, la]) ** 2).mean()
    loss_abs.backward()
    g_abs = {n: p.grad.detach().clone() for n, p in m.named_parameters() if p.grad is not None}

    def norm(gdict, prefix):
        return float(np.sqrt(sum(float((g ** 2).sum()) for n, g in gdict.items()
                                  if n.startswith(prefix))))
    s_params = ("alpha.", "psi.")
    g_ctr_s = sum((g_ctr[n] ** 2).sum() for n in g_ctr if n.startswith(s_params))
    g_abs_s = sum((g_abs[n] ** 2).sum() for n in g_abs if n.startswith(s_params))
    g_ctr_s = float(g_ctr_s.sqrt()); g_abs_s = float(g_abs_s.sqrt())
    cos = float((sum((g_ctr[n] * g_abs[n]).sum() for n in g_ctr if n.startswith(s_params))
                 / (g_ctr_s * g_abs_s + 1e-12)))
    return {
        "loss_ctr_value": float(loss_ctr.detach()),
        "loss_abs_value": float(loss_abs.detach()),
        "R_g": g_abs_s / (g_ctr_s + 1e-12),
        "C_g": cos,
        "per_branch": {
            "g_ctr": {"s": g_ctr_s, "b_P": norm(g_ctr, "b_P."),
                      "b_L": norm(g_ctr, "b_L."),
                      "enc": norm(g_ctr, "p_enc.") + norm(g_ctr, "l_enc.")},
            "g_abs": {"s": g_abs_s, "b_P": norm(g_abs, "b_P."),
                      "b_L": norm(g_abs, "b_L."),
                      "enc": norm(g_abs, "p_enc.") + norm(g_abs, "l_enc.")},
        },
        "note": "initialization only (no persisted trajectory)",
    }


def free_vs_unified(art, res, table):
    te = [i for i, s in enumerate(art["split"]["pair_split"]) if s == 2]
    fu = {r["pair"]: r for r in table["free_pairwise"]}
    ul = {r["pair"]: r for r in table["unified_local"]}
    joint = [i for i in te if fu[i]["rank_evaluable"] and ul[i]["rank_evaluable"]]
    by_parent = {}
    for i in joint:
        by_parent.setdefault(art["pairs"][i]["parent"], []).append(
            fu[i]["spearman"] - ul[i]["spearman"])
    parent_med = {p: float(np.median(v)) for p, v in by_parent.items()}
    rows_u = table["unified_local"]
    rows_f = table["free_pairwise"]
    # bootstrap on spearman over parents (both arms, cluster=parent)
    parents = sorted({r["parent"] for r in rows_u})
    rng = T.stable_rng("stageCIIP1A", "boot", "audit", "spearman", T.BOOT_SEED)
    vals_u = {p: [] for p in parents}
    vals_f = {p: [] for p in parents}
    for r in rows_u:
        if r["spearman"] is not None:
            vals_u[r["parent"]].append(r["spearman"])
    for r in rows_f:
        if r["spearman"] is not None:
            vals_f[r["parent"]].append(r["spearman"])
    lo = []
    for _ in range(T.BOOT_DRAWS):
        idx = rng.integers(len(parents), size=len(parents))
        su = [x for pi in idx for x in vals_u[parents[pi]]]
        sf = [x for pi in idx for x in vals_f[parents[pi]]]
        if su and sf:
            lo.append(float(np.mean(sf) - np.mean(su)))
    return {
        "jointly_evaluable_pairs": joint,
        "n_jointly_evaluable": len(joint),
        "per_joint_pair_delta": {int(i): float(fu[i]["spearman"] - ul[i]["spearman"])
                                 for i in joint},
        "parent_level_median_spearman_unified": {p: float(np.median(
            [r["spearman"] for r in rows_u if r["parent"] == p and r["spearman"] is not None]))
            for p in parents if any(r["parent"] == p and r["spearman"] is not None for r in rows_u)},
        "parent_level_median_spearman_free": {p: float(np.median(
            [r["spearman"] for r in rows_f if r["parent"] == p and r["spearman"] is not None]))
            for p in parents if any(r["parent"] == p and r["spearman"] is not None for r in rows_f)},
        "bootstrap_2.5_lo_free_minus_unified_spearman": float(np.percentile(lo, 2.5)) if lo else None,
        "bootstrap_50_free_minus_unified_spearman": float(np.percentile(lo, 50)) if lo else None,
        "parent_level_median_delta": parent_med,
        "driven_by_few_pairs": len(joint) <= 3,
        "centered_mse_improved": bool(any(
            fu[i]["mse"] < ul[i]["mse"] - 1e-9 for i in te)),
        "scale_near_1": {int(i): fu[i]["ols_slope"] for i in joint},
    }


def main() -> int:
    art, res, z, q1, esm = load_frozen()
    table, counts = collapse_table(art, res)
    klifs = klifs_diff_matrix(art, z)
    esmrep = esm_diff_matrix(art, esm)
    q1ev = q1_evidence(q1)
    init_an = init_model_analysis(art, z)
    grad = gradient_audit(art, z)
    fvu = free_vs_unified(art, res, table)
    # structural identity: constant prediction <-> zero KLIFS difference
    te_idx = [i for i, s in enumerate(art["split"]["pair_split"]) if s == 2]
    klifs_zero_test = {i for i in te_idx if np.linalg.norm(
        z["prot"][art["pairs"][i]["var_row"]] - z["prot"][art["pairs"][i]["wt_row"]]) == 0}
    const_test = {r["pair"] for r in table["unified_local"]
                  if not r["prediction_nonconstant"]}
    structural = {
        "test_pairs_with_zero_input_difference": sorted(klifs_zero_test),
        "test_pairs_with_constant_prediction": sorted(const_test),
        "constant_prediction_set_equals_zero_input_set": klifs_zero_test == const_test,
        "note": "for zero-input-difference pairs the potential contrast is "
                "IDENTICALLY 0 by antisymmetry (s(P,L)-s(P,L)=0), so constant "
                "prediction is forced by the representation, not by training",
    }
    out = {
        "schema": SCHEMA,
        "structural_collapse_identity": structural,
        "frozen_inputs_sha256": {
            "RESULT_SCREENING.json": sha(HERE / "RESULT_SCREENING.json"),
            "DATA1A.json": sha(HERE / "DATA1A.json"),
            "DATA1A.npz": sha(HERE / "DATA1A.npz"),
            "q1_esm_cache.npz": sha(X0C / "q1_esm_cache.npz"),
            "Q1_SELECTIVITY.json": sha(X0C / "Q1_SELECTIVITY.json"),
            "PREREGISTRATION_STAGE1_CIIP1A.md": sha(HERE / "PREREGISTRATION_STAGE1_CIIP1A.md"),
        },
        "collapse": {"per_arm": table, "counts": counts},
        "variance_sources": {
            "klifs_diff": klifs,
            "local_esm_diff": esmrep,
            "q1_evidence": q1ev,
            "init_model": init_an,
        },
        "gradient_audit": grad,
        "free_vs_unified": fvu,
    }
    path = HERE / "STAGE1_COLLAPSE_AUDIT.json"
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("collapse", "gradient_audit")}, indent=1)[:6000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
