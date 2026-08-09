"""Phase 2A / Phase 2 — teacher ligand-conditionality ceiling. Label-side only.

Registered by PREREG_S7_L2B_PHASE2A.md (sha 4e01401d...), section 6, plus
computational amendment 01.

The controlling idea: the NOISE FLOOR is measurable from the data. Two crystal
structures of the same exact construct with the SAME ligand differ only by
experimental variation. That replicate Jaccard is the correct comparator for
alternative-ligand Jaccard. An arbitrary foreign ligand has no teacher-level
counterpart at all (amendment A5), which is exactly why the wrong-ligand model
control cannot settle this question.

No model is loaded. No affinity source is opened.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import component_macro, paired_bootstrap  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
PREREG_SHA = "4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e"

DJ_MIN = 0.05
SEED_BOOT_T1 = 20260819
SEED_PERM = 20260820
N_PERM = 200
MIN_PAIRS_FOR_RHO = 5
MEANINGFUL_J = 0.5
MEANINGFUL_SYMDIFF = 3


# ------------------------------------------------------------ exact tied AP
def expected_ap(scores: np.ndarray, labels: np.ndarray) -> float | None:
    """Exact expectation of AP under a uniformly random ordering within each
    tied score block (amendment A1)."""
    P = float(labels.sum())
    if P == 0 or P == labels.size:
        return None
    order = np.argsort(-scores, kind="stable")
    s = scores[order]
    y = labels[order].astype(np.float64)
    bounds = np.flatnonzero(np.r_[True, s[1:] != s[:-1], True])
    total = 0.0
    a = 0
    b = 0.0
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        n = hi - lo
        k = float(y[lo:hi].sum())
        if k > 0:
            j = np.arange(1, n + 1, dtype=np.float64)
            grow = 0.0 if n == 1 else (k - 1.0) * (j - 1.0) / (n - 1.0)
            total += k * float(np.mean((b + 1.0 + grow) / (a + j)))
        a += n
        b += k
    return float(total / P)


def mc_ap(scores, labels, rng, reps=200):
    P = labels.sum()
    if P == 0:
        return None
    vals = []
    for _ in range(reps):
        order = np.lexsort((rng.random(scores.size), -scores))
        yy = labels[order].astype(np.float64)
        tp = np.cumsum(yy)
        prec = tp / np.arange(1, yy.size + 1, dtype=np.float64)
        vals.append(float((prec * yy).sum() / P))
    return float(np.mean(vals))


def selftest_expected_ap():
    rng = np.random.default_rng(7)
    worst = 0.0
    cases = []
    for _ in range(12):
        L = int(rng.integers(30, 120))
        src = (rng.random(L) < 0.2).astype(np.float64)          # binary tied source
        tgt = (rng.random(L) < 0.15).astype(np.int8)
        if tgt.sum() == 0:
            continue
        e = expected_ap(src, tgt)
        m = mc_ap(src, tgt, np.random.default_rng(11), reps=400)
        worst = max(worst, abs(e - m))
        cases.append({"L": L, "closed_form": e, "monte_carlo": m})
    return {"max_abs_diff_vs_monte_carlo": worst,
            "monte_carlo_reps": 400, "cases": len(cases),
            "pass": bool(worst < 0.02)}


def spearman(x, y):
    n = len(x)
    if n < 3:
        return None
    rx = np.argsort(np.argsort(np.asarray(x, float), kind="stable"), kind="stable").astype(float)
    ry = np.argsort(np.argsort(np.asarray(y, float), kind="stable"), kind="stable").astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    d = float(np.sqrt((rx * rx).sum() * (ry * ry).sum()))
    return None if d < 1e-12 else float((rx * ry).sum() / d)


def main():
    st = selftest_expected_ap()
    print("expected-AP self-test:", json.dumps({k: st[k] for k in
                                                ("max_abs_diff_vs_monte_carlo", "pass")}),
          flush=True)
    if not st["pass"]:
        (OUT / "PHASE2A_TEACHER_CONDITIONALITY.json").write_text(
            json.dumps({"verdict": "AP_ESTIMATOR_CONTRACT_FAIL", "selftest": st}, indent=2),
            encoding="utf-8")
        return 1

    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    mgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    kept, _q, _c, _f = build()
    comp_of = protein_components(kept)
    _tr, _ha, held_A, _hb = make_split(kept, comp_of)
    heldA = {r["source_key"] for r in held_A}

    # fingerprints keyed by CCD, from the same molecule that defined graph_key
    import pickle
    md = [pickle.load((MONN / "mol_dict").open("rb"), encoding="bytes"),
          pickle.load((MONN / "independent_dataset_mol_dict").open("rb"), encoding="bytes")]
    fps, natoms = {}, {}
    for r in kept:
        ccd = r["ligand_ccd"]
        if ccd in fps:
            continue
        mol = None
        for d in md:
            m = d.get(ccd.encode("ascii", "ignore"))
            if m is not None:
                mol = m
                break
        if mol is None:
            continue
        w = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(w)
            fps[ccd] = mgen.GetFingerprint(w)
            natoms[ccd] = w.GetNumAtoms()
        except Exception:
            pass
    print(f"fingerprints: {len(fps)} ligands", flush=True)

    by_construct = defaultdict(list)
    for r in kept:
        by_construct[r["seq_key"]].append(r)

    masks = {r["source_key"]: frozenset(i for i, _j in r["edges"]) for r in kept}

    # ---------------------------------------------------------------- pairs
    rep_J = defaultdict(list)       # component -> [jaccard]
    alt_J = defaultdict(list)
    rep_J_A = defaultdict(list)
    alt_J_A = defaultdict(list)
    gain, loss = [], []
    meaningful, meaningful_n = 0, 0
    per_construct_rho_in = defaultdict(lambda: ([], [], []))   # dJ, dT, ccds
    n_pairs_rep = n_pairs_alt = 0
    tan_alt = []

    for sk, recs in by_construct.items():
        if len(recs) < 2:
            continue
        for a, b in combinations(recs, 2):
            ra, rb = masks[a["source_key"]], masks[b["source_key"]]
            if not ra or not rb:
                continue
            inter = len(ra & rb)
            uni = len(ra | rb)
            J = inter / uni if uni else 1.0
            comp = comp_of[a["source_key"]]
            both_A = a["source_key"] in heldA and b["source_key"] in heldA
            if a["graph_key"] == b["graph_key"]:
                if a["pdb_id"] == b["pdb_id"]:
                    continue
                rep_J[comp].append(J)
                if both_A:
                    rep_J_A[comp].append(J)
                n_pairs_rep += 1
            else:
                sa, sb = a["scaffold"], b["scaffold"]
                if not (sa and sb and sa != sb):
                    continue
                alt_J[comp].append(J)
                if both_A:
                    alt_J_A[comp].append(J)
                n_pairs_alt += 1
                gain.append(len(rb - ra) / max(len(ra), 1))
                loss.append(len(ra - rb) / max(len(ra), 1))
                meaningful_n += 1
                if J <= MEANINGFUL_J and len(ra ^ rb) >= MEANINGFUL_SYMDIFF:
                    meaningful += 1
                fa, fb = fps.get(a["ligand_ccd"]), fps.get(b["ligand_ccd"])
                if fa is not None and fb is not None:
                    t = DataStructs.TanimotoSimilarity(fa, fb)
                    tan_alt.append(t)
                    dj, dt, cc = per_construct_rho_in[sk]
                    dj.append(1.0 - J)
                    dt.append(1.0 - t)
                    cc.append((a["ligand_ccd"], b["ligand_ccd"]))
        if len(by_construct) > 0 and len(rep_J) + len(alt_J) and False:
            pass
    print(f"pairs: replicate={n_pairs_rep} scaffold-distinct={n_pairs_alt}", flush=True)

    comp_rep = {c: float(np.mean(v)) for c, v in rep_J.items() if v}
    comp_alt = {c: float(np.mean(v)) for c, v in alt_J.items() if v}
    comp_rep_A = {c: float(np.mean(v)) for c, v in rep_J_A.items() if v}
    comp_alt_A = {c: float(np.mean(v)) for c, v in alt_J_A.items() if v}

    t1 = paired_bootstrap(comp_rep, comp_alt, n_boot=10000, seed=SEED_BOOT_T1)
    t1["pass"] = bool(t1["delta"] >= DJ_MIN and t1["lcb95_one_sided"] > 0)
    t1_A = (paired_bootstrap(comp_rep_A, comp_alt_A, n_boot=10000, seed=SEED_BOOT_T1)
            if len(set(comp_rep_A) & set(comp_alt_A)) >= 3 else None)

    # unpaired fallback
    rng = np.random.default_rng(SEED_BOOT_T1)
    ra_ = np.array(list(comp_rep.values()))
    aa_ = np.array(list(comp_alt.values()))
    bs = [ra_[rng.integers(0, ra_.size, ra_.size)].mean()
          - aa_[rng.integers(0, aa_.size, aa_.size)].mean() for _ in range(10000)]
    unpaired = {"delta": float(ra_.mean() - aa_.mean()),
                "lcb95_one_sided": float(np.percentile(bs, 5)),
                "units_replicate": int(ra_.size), "units_alternative": int(aa_.size)}

    # ------------------------------------------------------- T2 AP retention
    ap_self, ap_rep, ap_alt, ap_loo, ap_prev = {}, {}, {}, {}, {}
    for sk, recs in by_construct.items():
        if len(recs) < 2:
            continue
        L = recs[0]["n_res"]
        vec = {}
        for r in recs:
            v = np.zeros(L, dtype=np.float64)
            idx = np.fromiter(masks[r["source_key"]], dtype=np.int64) if masks[r["source_key"]] else None
            if idx is not None:
                v[idx] = 1.0
            vec[r["source_key"]] = v
        for t in recs:
            kt = t["source_key"]
            yt = vec[kt].astype(np.int8)
            if yt.sum() == 0:
                continue
            others = [r for r in recs if r["source_key"] != kt]
            if not others:
                continue
            ap_self[kt] = expected_ap(vec[kt], yt)
            ap_prev[kt] = float(yt.sum() / L)
            reps = [vec[r["source_key"]] for r in others
                    if r["graph_key"] == t["graph_key"] and r["pdb_id"] != t["pdb_id"]]
            alts = [vec[r["source_key"]] for r in others
                    if r["graph_key"] != t["graph_key"] and r["scaffold"] and t["scaffold"]
                    and r["scaffold"] != t["scaffold"]]
            if reps:
                vals = [expected_ap(s, yt) for s in reps]
                ap_rep[kt] = float(np.mean([v for v in vals if v is not None]))
            if alts:
                vals = [expected_ap(s, yt) for s in alts]
                ap_alt[kt] = float(np.mean([v for v in vals if v is not None]))
                loo = np.mean(alts, axis=0)
                ap_loo[kt] = expected_ap(loo, yt)

    def macro(d):
        return component_macro({k: v for k, v in d.items() if v is not None}, comp_of)

    cm_self, m_self = macro(ap_self)
    cm_rep, m_rep = macro(ap_rep)
    cm_alt, m_alt = macro(ap_alt)
    cm_loo, m_loo = macro(ap_loo)
    cm_prev, m_prev = macro(ap_prev)
    t2_rep_minus_alt = paired_bootstrap(cm_rep, cm_alt, n_boot=10000, seed=SEED_BOOT_T1)
    t2_alt_minus_prev = paired_bootstrap(cm_alt, cm_prev, n_boot=10000, seed=SEED_BOOT_T1)
    t2_loo_minus_alt = paired_bootstrap(cm_loo, cm_alt, n_boot=10000, seed=SEED_BOOT_T1)

    # -------------------------------------------------- T4 between-ligand var
    var_multi, var_rep_only = [], []
    for sk, recs in by_construct.items():
        gks = {r["graph_key"] for r in recs}
        if len(recs) < 2:
            continue
        L = recs[0]["n_res"]
        Y = np.zeros((len(recs), L), dtype=np.float64)
        for i, r in enumerate(recs):
            idx = np.fromiter(masks[r["source_key"]], dtype=np.int64) if masks[r["source_key"]] else None
            if idx is not None:
                Y[i, idx] = 1.0
        ybar = Y.mean(0)
        num = float((ybar * (1 - ybar)).sum())
        den = float(ybar.sum())
        if den <= 0:
            continue
        (var_multi if len(gks) >= 2 else var_rep_only).append(num / den)

    # ------------------------------------------- T5/T6 scaffold-distance rho
    rhos, rho_comp = {}, defaultdict(list)
    excluded_small = 0
    perm_rng = np.random.default_rng(SEED_PERM)
    perm_z, perm_p_terms = [], []
    for sk, (dj, dt, cc) in per_construct_rho_in.items():
        if len(dj) < MIN_PAIRS_FOR_RHO:
            excluded_small += 1
            continue
        r = spearman(dt, dj)
        if r is None:
            excluded_small += 1
            continue
        rhos[sk] = r
        comp = comp_of[by_construct[sk][0]["source_key"]]
        rho_comp[comp].append(r)
        null = []
        dtv = np.asarray(dt, float)
        for _ in range(N_PERM):
            rn = spearman(perm_rng.permutation(dtv), dj)
            if rn is not None:
                null.append(rn)
        if len(null) >= N_PERM // 2:
            null = np.asarray(null)
            perm_z.append((r - null.mean()) / max(null.std(), 1e-9))
            perm_p_terms.append(float((null >= r).mean()))

    comp_rho = {c: float(np.mean(v)) for c, v in rho_comp.items() if v}
    macro_rho = float(np.mean(list(comp_rho.values()))) if comp_rho else float("nan")
    # component-level bootstrap of the macro rho
    if comp_rho:
        vals = np.array(list(comp_rho.values()))
        rng2 = np.random.default_rng(SEED_BOOT_T1)
        bs2 = vals[rng2.integers(0, vals.size, (10000, vals.size))].mean(1)
        rho_lcb = float(np.percentile(bs2, 5))
    else:
        rho_lcb = float("nan")
    # global permutation p-value: fraction of constructs whose observed rho does
    # not exceed its own null, combined by the median per-construct p
    median_perm_p = float(np.median(perm_p_terms)) if perm_p_terms else float("nan")
    median_perm_z = float(np.median(perm_z)) if perm_z else float("nan")
    frac_pos_z = float(np.mean(np.asarray(perm_z) > 0)) if perm_z else float("nan")
    # Fisher combination over independent constructs is not valid here (constructs
    # inside one component are dependent), so the registered decision uses the
    # component-macro rho lower bound plus the sign consistency of the z values.
    t6_pass = bool(rho_lcb > 0 and median_perm_p <= 0.05)

    res = {
        "schema": "MetaSieve.S7L2B.P2A.TeacherConditionality.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "amendment": "PREREG_S7_L2B_PHASE2A_AMENDMENT_01.md",
        "repo_commit": "623602e76b7d4f445af069014782278163183d59",
        "scope": "label-side only; no model loaded; no affinity source opened",
        "expected_ap_selftest": st,
        "pairs": {"replicate": n_pairs_rep, "scaffold_distinct": n_pairs_alt,
                  "components_replicate": len(comp_rep),
                  "components_alternative": len(comp_alt),
                  "components_paired": len(set(comp_rep) & set(comp_alt))},
        "T1_jaccard": {
            "definition": "component-macro mean Jaccard of residue masks",
            "replicate_mean": float(np.mean(list(comp_rep.values()))),
            "alternative_ligand_mean": float(np.mean(list(comp_alt.values()))),
            "dJ_paired": t1,
            "dJ_unpaired_fallback": unpaired,
            "heldout_A_only": t1_A,
            "threshold_dJ_min": DJ_MIN,
            "pass": t1["pass"],
        },
        "T2_ap_retention": {
            "self_mask_degenerate": m_self,
            "replicate_mask": m_rep,
            "alternative_ligand_mask": m_alt,
            "leave_one_out_protein_marginal": m_loo,
            "prevalence": m_prev,
            "replicate_minus_alternative": t2_rep_minus_alt,
            "alternative_minus_prevalence": t2_alt_minus_prev,
            "loo_marginal_minus_alternative": t2_loo_minus_alt,
            "corruption_floor_note": "no teacher-level foreign-ligand control exists; "
                                     "residue indices are comparable only within one "
                                     "construct (amendment A5). Prevalence is the floor.",
        },
        "T3_gain_loss": {
            "mean_gain_rate": float(np.mean(gain)) if gain else None,
            "mean_loss_rate": float(np.mean(loss)) if loss else None,
            "median_gain_rate": float(np.median(gain)) if gain else None,
            "median_loss_rate": float(np.median(loss)) if loss else None,
            "n_pairs": len(gain),
        },
        "T4_between_ligand_variance": {
            "multi_ligand_constructs_mean": float(np.mean(var_multi)) if var_multi else None,
            "multi_ligand_constructs_median": float(np.median(var_multi)) if var_multi else None,
            "replicate_only_constructs_mean": float(np.mean(var_rep_only)) if var_rep_only else None,
            "replicate_only_constructs_median": float(np.median(var_rep_only)) if var_rep_only else None,
            "n_multi": len(var_multi), "n_replicate_only": len(var_rep_only),
        },
        "T5_scaffold_distance_sensitivity": {
            "statistic": "Spearman rho between (1 - Jaccard) and (1 - Tanimoto), "
                         "within construct, component-macro",
            "constructs_used": len(rhos),
            "constructs_excluded_fewer_than_5_pairs": excluded_small,
            "component_macro_rho": macro_rho,
            "component_bootstrap_lcb95": rho_lcb,
            "mean_tanimoto_alternative_pairs": float(np.mean(tan_alt)) if tan_alt else None,
        },
        "T6_ligand_permutation_control": {
            "permutations": N_PERM, "seed": SEED_PERM,
            "median_per_construct_permutation_p": median_perm_p,
            "median_per_construct_z": median_perm_z,
            "fraction_constructs_with_positive_z": frac_pos_z,
            "pass": t6_pass,
        },
        "T7_meaningful_change": {
            "definition": f"Jaccard <= {MEANINGFUL_J} and symmetric difference >= "
                          f"{MEANINGFUL_SYMDIFF} residues",
            "fraction_of_scaffold_distinct_pairs": (meaningful / meaningful_n) if meaningful_n else None,
            "n": meaningful_n,
        },
        "teacher_verdict": {
            "rule": "ligand-conditioned iff T1 and T6 both pass",
            "T1_pass": t1["pass"], "T6_pass": t6_pass,
            "ligand_conditioned": bool(t1["pass"] and t6_pass),
        },
    }
    (OUT / "PHASE2A_TEACHER_CONDITIONALITY.json").write_text(
        json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps({k: res[k] for k in
                      ("pairs", "T1_jaccard", "T2_ap_retention", "T3_gain_loss",
                       "T4_between_ligand_variance", "T5_scaffold_distance_sensitivity",
                       "T6_ligand_permutation_control", "T7_meaningful_change",
                       "teacher_verdict")}, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
