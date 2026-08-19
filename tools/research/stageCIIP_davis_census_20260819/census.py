"""Stage CIIP-0b Davis census (read-only; prereg 2dd8b708...)."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DL = HERE.parent / "stageX_csc_signal" / "downloads"
OUT = HERE
PREREG_SHA = "2dd8b70829c2fdfef68f9201e0764d0cc10d7d9809934a5ac195b2b922f8c3fc"


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    meta = pd.read_excel(DL / "davis_MOESM3.xls", "SuppTable1-050511")
    kd = pd.read_excel(DL / "davis_MOESM5.xls", "Sheet1")
    # align rows by order (both 442 rows, same ordering asserted below)
    assert len(meta) == len(kd) == 442
    assert (meta["Accession Number"].astype(str).tolist()
            == kd["Accession Number"].astype(str).tolist())
    mut = meta["Mutant"].astype(str).str.upper().str.strip()
    gene = meta["Entrez Gene Symbol"].astype(str)
    kdmat = kd.iloc[:, 3:].apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float64)
    n_lig = kdmat.shape[1]

    is_wt = (mut == "NO").to_numpy()
    is_mut = (mut == "YES").to_numpy()
    quarantined = ~is_wt & ~is_mut

    # duplicates
    dup_rows = {}
    for i in range(442):
        key = (gene[i], mut[i])
        dup_rows.setdefault(key, []).append(i)
    dup_keys = {k: v for k, v in dup_rows.items() if len(v) > 1}

    # fusion heuristic: kinase names containing fusion markers
    kname = meta["Kinase"].astype(str).str.upper()
    fusion_markers = ["-ABL", "BCR-", "FUSION", "TEL-", "NPM-", "EML4-", "/"]
    fusion_rows = kname.apply(lambda s: any(m in s for m in fusion_markers))

    pairs = []  # (wt_i, mut_j, gene)
    for g in sorted(set(gene)):
        wts = [i for i in range(442) if gene[i] == g and is_wt[i]]
        mus = [i for i in range(442) if gene[i] == g and is_mut[i]]
        for w in wts:
            for m_ in mus:
                common = np.isfinite(kdmat[w]) & np.isfinite(kdmat[m_])
                if int(common.sum()) >= 1:
                    pairs.append((w, m_, g, int(common.sum())))

    pair_stats = []
    for w, m_, g, nc in pairs:
        fin = np.isfinite(kdmat[w]) & np.isfinite(kdmat[m_])
        if fin.sum() >= 3:
            lr = np.log10(kdmat[w][fin] / kdmat[m_][fin])
            var = float(np.var(lr, ddof=1))
            n_eff = float(np.var(lr, ddof=1) / (fin.sum() ** 0.5)) if fin.sum() > 0 else np.nan
        else:
            var, n_eff = np.nan, np.nan
        pair_stats.append({
            "gene": g, "wt_row": int(w), "mut_row": int(m_),
            "common_ligands": int(nc),
            "logratio_var": var, "logratio_se": n_eff,
        })

    per_parent = {}
    for ps in pair_stats:
        per_parent.setdefault(ps["gene"], []).append(ps)
    parents_multi = {g: v for g, v in per_parent.items() if len(v) >= 2}

    na_frac = float(np.isnan(kdmat).mean())
    cap_9900_frac = float((kdmat == 9900.0).sum() / kdmat.size)

    out = {
        "schema": "MetaSieve.StageCIIP0b.DavisCensus.v1",
        "preregistration_sha256": PREREG_SHA,
        "input_files": {
            "davis_MOESM3.xls_sha256": sha256_file(DL / "davis_MOESM3.xls"),
            "davis_MOESM5.xls_sha256": sha256_file(DL / "davis_MOESM5.xls"),
        },
        "items": {
            "1_usable_wt_variant_pairs": int(len(pairs)),
            "2_identical_ligands_per_pair": {
                "median": float(np.median([p["common_ligands"] for p in pair_stats])),
                "p10": float(np.percentile([p["common_ligands"] for p in pair_stats], 10)),
                "min": int(min(p["common_ligands"] for p in pair_stats)),
                "max": int(max(p["common_ligands"] for p in pair_stats)),
            },
            "3_parent_mutation_fusion_counts": {
                "genes": int(gene.nunique()),
                "kinase_rows": 442,
                "wt_rows": int(is_wt.sum()),
                "mutant_rows": int(is_mut.sum()),
                "quarantined_rows": int(quarantined.sum()),
                "fusion_flagged_rows": int(fusion_rows.sum()),
                "duplicate_keys": {str(k): [int(i) for i in v] for k, v in dup_keys.items()},
            },
            "4_condition_completeness": {
                "state": "PARTIAL",
                "note": "single-assay competition-binding Kd panel; fixed protocol "
                        "(one condition class); measured-cell fraction = "
                        f"{1 - na_frac:.3f}; no per-cell ATP/construct table "
                        "in these SI files",
            },
            "5_duplicate_and_saturation_fraction": {
                "na_fraction": na_frac,
                "cap_9900_fraction": cap_9900_frac,
                "note": "9900 nM is the largest observed value; whether NA means "
                        "'not tested' or '>10 uM' is not distinguishable in this "
                        "file -> recorded UNKNOWN",
            },
            "6_endpoint": {
                "quantity": "Kd [nM], competition binding assay",
                "direction": "larger = weaker binding",
                "note": "never relabeled pK/Ki/DTA",
            },
            "7_parent_connectivity": {
                "parents_with_pairs": int(len(per_parent)),
                "pairs_per_parent_median": float(np.median([len(v) for v in per_parent.values()])),
            },
            "8_heldout_parent_folds": {
                "parents_with_2plus_mutants": int(len(parents_multi)),
                "note": "each such parent supports leave-one-parent-out folds",
            },
            "9_ligand_coverage_variance": {
                "pairs_with_3plus_common_ligands": int(
                    sum(1 for p in pair_stats if p["common_ligands"] >= 3)),
                "logratio_var_median": float(np.nanmedian(
                    [p["logratio_var"] for p in pair_stats if p["common_ligands"] >= 3])),
                "logratio_var_p10_p90": [
                    float(np.nanpercentile([p["logratio_var"] for p in pair_stats
                                            if p["common_ligands"] >= 3], 10)),
                    float(np.nanpercentile([p["logratio_var"] for p in pair_stats
                                            if p["common_ligands"] >= 3], 90))],
            },
            "10_data_availability": {
                "local": True,
                "provenance": "Nature Biotechnology 2011 (Davis et al.) supplementary "
                              "tables, kept under tools/research/stageX_csc_signal/downloads/",
                "urls": ["https://www.nature.com/articles/nbt.1990"],
            },
        },
        "pair_details": pair_stats,
        "stop_rule_check": {
            "pairs_ge_20": len(pairs) >= 20,
            "median_common_ligands_ge_5": float(np.median(
                [p["common_ligands"] for p in pair_stats])) >= 5,
            "heldout_parents_ge_10": len(parents_multi) >= 10,
        },
    }
    (OUT / "CENSUS.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps(out["items"], indent=1)[:3000])
    print("stop_rule_check:", out["stop_rule_check"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
