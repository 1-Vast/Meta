"""P1 stratification + P2 screening analysis on the frozen bank.

Recomputes the two cheap arms through the same code path as
p_baselines.py (determinism asserted against the SHA-pinned artifact),
then:
- P1 strata: support-query mean Tanimoto bands, Bemis-Murcko scaffold
  novelty vs p_train, activity-cliff cells (|d pKi| >= 2.0 to any support
  neighbour).
- P2 screening (frozen threshold: active = pKi >= 6.0): EF@1/5/10%,
  BEDROC(alpha=20), PR-AUC over pooled query cells per split and arm,
  k in {0,5,10,20,40} (k=0 is the degenerate global-mean ranker, reported
  as such).
Artifacts: P1_STRATA.json, P2_SCREENING.json (+ manifests).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import DataStructs

import p_baselines as PB

HERE = Path(__file__).resolve().parent
OUT = HERE / "artifacts"
SCHEMA_S = "MetaSieve.StageP.P1Strata.v1"
SCHEMA_C = "MetaSieve.StageP.P2Screening.v1"
K_LIST = (0, 1, 2, 3, 5, 10, 20, 40)
CLIFF_DELTA = 2.0
ACTIVE_THRESHOLD = 6.0
BANDS = {"low": (-1.0, 0.30), "mid": (0.30, 0.60), "high": (0.60, 2.0)}


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def ef(yhat, y, active, frac):
    n = len(y)
    top = int(np.ceil(frac * n))
    if top <= 0:
        return float("nan")
    hits = active[np.argsort(-yhat)[:top]].sum()
    total_hits = active.sum()
    if total_hits == 0:
        return float("nan")
    return float(hits / total_hits / frac)


def bedroc(yhat, active, alpha=20.0):
    n = len(yhat)
    order = np.argsort(-yhat)
    rank = np.empty(n, dtype=np.float64)
    rank[order] = np.arange(1, n + 1)
    ri = np.sum(1.0 / np.power(rank[active], alpha))
    ra = active.sum()
    ri_max = sum(1.0 / np.power(float(i), alpha) for i in range(1, int(ra) + 1))
    ri_min = sum(1.0 / np.power(float(n - ra + i), alpha)
                 for i in range(1, int(ra) + 1))
    if ri_max == ri_min:
        return float("nan")
    return float((ri - ri_min) / (ri_max - ri_min))


def pr_auc(yhat, active):
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(active, yhat))


def predict_record(rec, pki, split_art, lig, gmean, fp_cache):
    """Same predictions as p_baselines (verified against artifact)."""
    k = rec["k"]
    q_ids = rec["query_cell_ids"]
    q_y = np.array([pki[c] for c in q_ids], dtype=np.float64)
    q_ligs = [split_art["cell_split"][c]["ligand_id"] for c in q_ids]
    if k == 0:
        yh_lo = np.full(len(q_ids), gmean)
        yh_t = np.full(len(q_ids), gmean)
        sims_all = None
    else:
        sup_y = np.array([pki[c] for c in rec["support_cell_ids"]], dtype=np.float64)
        yh_lo = np.full(len(q_ids), sup_y.mean())
        sup_ligs = [split_art["cell_split"][c]["ligand_id"]
                    for c in rec["support_cell_ids"]]
        yh_t = []
        sims_all = []
        for ql in q_ligs:
            qf = fp_cache.setdefault(ql, PB.ecfp(lig[ql]))
            sims = []
            for sl in sup_ligs:
                sf = fp_cache.setdefault(sl, PB.ecfp(lig[sl]))
                sims.append(0.0 if qf is None or sf is None else
                            DataStructs.TanimotoSimilarity(qf, sf))
            sims = np.asarray(sims)
            top = np.argsort(-sims)[:3]
            w = sims[top]
            yh_t.append(float(sup_y.mean()) if w.sum() <= 0 else
                        float((w * sup_y[top]).sum() / w.sum()))
            sims_all.append(float(sims.mean()))
        yh_t = np.asarray(yh_t)
        sims_all = np.asarray(sims_all)
    return q_y, yh_lo, yh_t, sims_all


def main() -> int:
    bank = json.loads((OUT / "P_BANK.json").read_text(encoding="utf-8"))
    split_art = json.loads((OUT / "P_SPLIT.json").read_text(encoding="utf-8"))
    pki = PB.load_labels()
    lig = PB.load_ligands()
    fp_cache = {}
    train_ids = [cid for cid, rec in split_art["cell_split"].items()
                 if rec["split"] == "p_train"]
    gmean = float(np.mean([pki[c] for c in train_ids]))
    # scaffolds from the corpus ligands.jsonl (Bemis-Murcko, label-free)
    scaffolds = {}
    with open(PB.CORPUS / "ligands.jsonl", encoding="utf-8") as f:
        for line in f:
            d = json.loads(line)
            scaffolds[d["drug_key"]] = d.get("scaffold")
    train_scaffolds = {scaffolds[split_art["cell_split"][cid]["ligand_id"]]
                       for cid in train_ids
                       if split_art["cell_split"][cid]["ligand_id"] in scaffolds}

    strata = {arm: {str(k): [] for k in K_LIST} for arm in ("ligand_only", "tanimoto")}
    screening = {arm: {sp: {str(k): None for k in K_LIST}
                       for sp in ("p_val", "p_test")}
                 for arm in ("ligand_only", "tanimoto")}
    pool = {arm: {sp: {str(k): {"yhat": [], "y": []} for k in K_LIST}
                  for sp in ("p_val", "p_test")}
            for arm in ("ligand_only", "tanimoto")}
    for rec in bank["records"]:
        q_y, yh_lo, yh_t, sims = predict_record(rec, pki, split_art, lig, gmean,
                                                fp_cache)
        k = str(rec["k"])
        sp = rec["split"]
        for arm, yh in (("ligand_only", yh_lo), ("tanimoto", yh_t)):
            pool[arm][sp][k]["yhat"].extend(yh.tolist())
            pool[arm][sp][k]["y"].extend(q_y.tolist())
            mse = float(np.mean((yh - q_y) ** 2))
            strata[arm][k].append({
                "split": sp, "k": rec["k"], "mse": mse,
                "sim_mean": float(sims.mean()) if sims is not None else None,
                "scaffold_novel_frac": float(np.mean(
                    [scaffolds.get(split_art["cell_split"][c]["ligand_id"])
                     not in train_scaffolds for c in rec["query_cell_ids"]]))
                if rec["k"] > 0 else None,
                "cliff_frac": float(np.mean(
                    [max(abs(pki[c] - pki[s]) for s in rec["support_cell_ids"])
                     >= CLIFF_DELTA for c in rec["query_cell_ids"]]))
                if rec["k"] > 0 else None,
            })
    # determinism check against the pinned baseline artifact
    bl = json.loads((OUT / "P1_BASELINES.json").read_text(encoding="utf-8"))
    for arm in ("ligand_only", "tanimoto"):
        for k in K_LIST:
            key = f"p_test:k{k}"
            a = bl["arms"][arm].get(key)
            if a is None:
                continue
            rows = [r for r in strata[arm][str(k)] if r["split"] == "p_test"]
            mse = float(np.mean([r["mse"] for r in rows]))
            assert abs(mse - a["mse"]) < 1e-9, (arm, k, mse, a["mse"])
    print("determinism vs P1_BASELINES.json: OK")

    out_s = {"schema": SCHEMA_S,
             "bank_sha256": sha256_file(OUT / "P_BANK.json"),
             "baselines_sha256": sha256_file(OUT / "P1_BASELINES.json"),
             "bands": BANDS, "cliff_delta": CLIFF_DELTA, "strata": {}}
    for arm in ("ligand_only", "tanimoto"):
        out_s["strata"][arm] = {}
        for sp in ("p_val", "p_test"):
            for k in K_LIST:
                rows = [r for r in strata[arm][str(k)] if r["split"] == sp
                        and r["sim_mean"] is not None]
                if not rows:
                    continue
                agg = {}
                for band, (lo, hi) in BANDS.items():
                    sel = [r for r in rows if lo <= r["sim_mean"] < hi]
                    agg[band] = {
                        "n": len(sel),
                        "mse": float(np.mean([r["mse"] for r in sel])) if sel else None,
                        "scaffold_novel_frac": float(np.mean(
                            [r["scaffold_novel_frac"] for r in sel])) if sel else None,
                        "cliff_frac": float(np.mean(
                            [r["cliff_frac"] for r in sel])) if sel else None,
                    }
                for key, fn in (("cliff", lambda r: r["cliff_frac"] >= 0.5),
                                ("non_cliff", lambda r: r["cliff_frac"] < 0.5)):
                    sel = [r for r in rows if fn(r)]
                    agg[key] = {
                        "n": len(sel),
                        "mse": float(np.mean([r["mse"] for r in sel])) if sel else None,
                    }
                for key, fn in (("novel", lambda r: r["scaffold_novel_frac"] >= 0.5),
                                ("known", lambda r: r["scaffold_novel_frac"] < 0.5)):
                    sel = [r for r in rows if fn(r)]
                    agg[key] = {
                        "n": len(sel),
                        "mse": float(np.mean([r["mse"] for r in sel])) if sel else None,
                    }
                out_s["strata"][arm][f"{sp}:k{k}"] = agg

    out_c = {"schema": SCHEMA_C,
             "bank_sha256": sha256_file(OUT / "P_BANK.json"),
             "active_threshold_pki": ACTIVE_THRESHOLD,
             "screening": {}}
    for arm in ("ligand_only", "tanimoto"):
        out_c["screening"][arm] = {}
        for sp in ("p_val", "p_test"):
            for k in K_LIST:
                yh = np.asarray(pool[arm][sp][str(k)]["yhat"])
                y = np.asarray(pool[arm][sp][str(k)]["y"])
                if len(yh) == 0:
                    continue
                active = y >= ACTIVE_THRESHOLD
                out_c["screening"][arm][f"{sp}:k{k}"] = {
                    "n": int(len(yh)),
                    "active_frac": float(active.mean()),
                    "ef_1pct": ef(yh, y, active, 0.01),
                    "ef_5pct": ef(yh, y, active, 0.05),
                    "ef_10pct": ef(yh, y, active, 0.10),
                    "bedroc_a20": bedroc(yh, active, 20.0),
                    "pr_auc": pr_auc(yh, active),
                    "degenerate_ranker": bool(k == 0),
                }
    for name, obj, schema in (("P1_STRATA.json", out_s, SCHEMA_S),
                              ("P2_SCREENING.json", out_c, SCHEMA_C)):
        text = json.dumps(obj, indent=1, sort_keys=True)
        path = OUT / name
        path.write_text(text, encoding="utf-8")
        art_sha = sha256_file(path)
        (OUT / (name + ".manifest.json")).write_text(
            json.dumps({"schema": schema + ".Manifest", "file": name,
                        "sha256": art_sha}, indent=1), encoding="utf-8")
        print("wrote", name, "sha", art_sha)

    # console table (p_test)
    print("=== P1 strata MSE (p_test, tanimoto) ===")
    for k in K_LIST:
        key = f"p_test:k{k}"
        if key not in out_s["strata"]["tanimoto"]:
            continue
        agg = out_s["strata"]["tanimoto"][key]
        print(f"k={k:<2}", " ".join(
            f"{b}={agg[b]['mse']:.3f}(n={agg[b]['n']})" if agg[b]['mse'] is not None
            else f"{b}=-(n=0)" for b in ("low", "mid", "high")),
            f"cliff={agg['cliff']['mse']:.3f}(n={agg['cliff']['n']})"
            if agg['cliff']['mse'] is not None else "cliff=-",
            f"noncliff={agg['non_cliff']['mse']:.3f}(n={agg['non_cliff']['n']})"
            if agg['non_cliff']['mse'] is not None else "noncliff=-")
    print("=== P2 screening (p_test, tanimoto) ===")
    for k in K_LIST:
        a = out_c["screening"]["tanimoto"].get(f"p_test:k{k}")
        if a is None:
            continue
        print(f"k={k:<2} n={a['n']} act={a['active_frac']:.3f} "
              f"EF1={a['ef_1pct']:.3f} EF5={a['ef_5pct']:.3f} "
              f"EF10={a['ef_10pct']:.3f} BEDROC={a['bedroc_a20']:.3f} "
              f"PR-AUC={a['pr_auc']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
