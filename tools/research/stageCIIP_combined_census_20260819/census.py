"""Stage CIIP-0c combined census (read-only; prereg 6952ed1a...).

Panels audited: Anastassiadis 2011 MOESM23 (% remaining activity),
Duong-Ly mmc2/mmc3 (% inhibition). Duong-Ly row-name parsing admits
C-KIT/C-MET/C-SRC/P38a/PDGFRa/TIE2 (hyphen + Greek). Multi-mutant tags
(slash-separated) are excluded from the single-mutant estimand and
counted separately. Davis numbers re-quoted from stageCIIP_davis_census.
"""
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DL = HERE.parent / "stageX_csc_signal" / "downloads"
PREREG_SHA = "6952ed1aef7dedda92c1be45a9c10192d5640f6f208733bc4c31e993b0fec2fd"

NAME_RE = re.compile(r"^([A-Z0-9α]+(?:/MAPK[0-9]+)?(?:/TEK)?)(?:\(([^)]+)\))?$",
                     re.IGNORECASE)


def sha256_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def parse_duongly():
    m3 = pd.read_excel(DL / "duongly_mmc3.xlsx", "Table S2", header=0)
    names = m3.iloc[:, 0].astype(str).tolist()
    mat = m3.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float64)
    kin = names[1:]
    M = mat[1:]
    wt = {}
    muts = {}
    multi = {}
    unparsed = []
    for i, k in enumerate(kin):
        s = (k.strip().replace("C-KIT", "CKIT").replace("C-MET", "CMET")
             .replace("C-SRC", "CSRC").replace("P38α", "P38A")
             .replace("PDGFRα", "PDGFRA").replace("P38A/MAPK14", "P38A")
             .replace("TIE2/TEK", "TIE2"))
        m = re.match(r"^([A-Z0-9]+)(?:\(([^)]+)\))?$", s)
        if not m:
            unparsed.append(k)
            continue
        g, mu = m.group(1), m.group(2)
        if mu is None or mu.upper() in ("WT", "WILD TYPE", "WILD-TYPE"):
            wt[g] = i
        elif "/" in mu:
            multi.setdefault(g, []).append((i, mu))
        else:
            muts.setdefault(g, []).append((i, mu))
    pairs = []
    for g, w_i in wt.items():
        for m_i, mu in muts.get(g, []):
            fin = np.isfinite(M[w_i]) & np.isfinite(M[m_i])
            if int(fin.sum()) >= 1:
                var = (float(np.nanvar(M[w_i][fin] - M[m_i][fin], ddof=1))
                       if fin.sum() >= 2 else None)
                pairs.append({"gene": g, "mutation": mu,
                              "common_ligands": int(fin.sum()),
                              "effect_var": var})
    return {
        "kinase_rows": len(kin),
        "ligand_columns": int(M.shape[1]),
        "wt_rows": len(wt),
        "single_mutant_rows": sum(len(v) for v in muts.values()),
        "multi_mutant_rows": sum(len(v) for v in multi.values()),
        "multi_mutant_tags": {g: [m for _, m in v] for g, v in multi.items()},
        "unparsed_rows": unparsed,
        "value_range": [float(np.nanmin(M)), float(np.nanmax(M))],
        "na_fraction": float(np.isnan(M).mean()),
        "pairs": pairs,
        "endpoint": "% inhibition (0-100 nominal; observed beyond bounds)",
        "direction": "larger = stronger inhibition",
    }


def anastassiadis():
    xl = pd.ExcelFile(DL / "anastassiadis_MOESM23.xls")
    df = xl.parse("Sheet1", header=None)
    names = df.iloc[:, 0].astype(str).tolist()
    mat = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(
        dtype=np.float64)
    return {
        "kinase_rows": len(names),
        "ligand_columns": int(mat.shape[1]),
        "mutant_rows": int(sum("(" in n for n in names)),
        "wt_rows": int(sum("(" not in n for n in names)),
        "na_fraction": float(np.isnan(mat).mean()),
        "value_range": [float(np.nanmin(mat)), float(np.nanmax(mat))],
        "endpoint": "% remaining activity (larger = weaker inhibition)",
        "pairs_within_panel": 0,
        "role": "WT-heavy screen; no within-panel WT-variant pairs; "
                "cross-endpoint replication reference only",
    }


def main() -> int:
    dl = parse_duongly()
    an = anastassiadis()
    davis = json.loads((HERE.parent / "stageCIIP_davis_census_20260819"
                        / "CENSUS.json").read_text(encoding="utf-8"))
    pairs = dl["pairs"]
    per_parent = {}
    for p in pairs:
        per_parent.setdefault(p["gene"], []).append(p)
    parents_2plus = {g: len(v) for g, v in per_parent.items() if len(v) >= 2}
    common = [p["common_ligands"] for p in pairs]
    out = {
        "schema": "MetaSieve.StageCIIP0c.CombinedCensus.v1",
        "preregistration_sha256": PREREG_SHA,
        "inputs_sha256": {
            "duongly_mmc3.xlsx": sha256_file(DL / "duongly_mmc3.xlsx"),
            "anastassiadis_MOESM23.xls": sha256_file(DL / "anastassiadis_MOESM23.xls"),
        },
        "duongly": {
            "rows": dl["kinase_rows"], "ligands": dl["ligand_columns"],
            "wt_rows": dl["wt_rows"],
            "single_mutant_rows": dl["single_mutant_rows"],
            "multi_mutant_rows": dl["multi_mutant_rows"],
            "multi_mutant_tags": dl["multi_mutant_tags"],
            "unparsed_rows": dl["unparsed_rows"],
            "na_fraction": dl["na_fraction"],
            "value_range": dl["value_range"],
            "endpoint": dl["endpoint"], "direction": dl["direction"],
            "usable_pairs": len(pairs),
            "common_ligands_median": float(np.median(common)) if common else None,
            "common_ligands_min": int(min(common)) if common else None,
            "common_ligands_max": int(max(common)) if common else None,
            "parents_with_pairs": len(per_parent),
            "parents_with_2plus_pairs": int(len(parents_2plus)),
            "per_parent_pairs": {g: len(v) for g, v in
                                 sorted(per_parent.items())},
            "effect_var_median_pairs_2plus_ligands": float(np.nanmedian(
                [p["effect_var"] for p in pairs if p["common_ligands"] >= 2]))
            if any(p["common_ligands"] >= 2 for p in pairs) else None,
            "single_mutant_pairs_ge_20": len(pairs) >= 20,
            "median_common_ge_5": bool(common) and np.median(common) >= 5,
            "heldout_parents_ge_10": len(parents_2plus) >= 10,
            "ciiP1A_admissible": len(pairs) >= 20 and bool(common)
                                  and np.median(common) >= 5,
            "ciiP1B_admissible": (len(pairs) >= 20 and bool(common)
                                   and np.median(common) >= 5
                                   and len(parents_2plus) >= 10),
        },
        "anastassiadis": an,
        "davis_requoted": {
            "pairs": davis["items"]["1_usable_wt_variant_pairs"],
            "heldout_parents": davis["items"]["8_heldout_parent_folds"],
            "verdict": "INSUFFICIENT ALONE",
        },
        "admissibility_classification": {
            "davis": "same-platform (Kd binding): CIIP-1A yes, CIIP-1B no (7 parents)",
            "duongly": ("same-platform (% inhibition): CIIP-1A "
                        + ("yes" if (len(pairs) >= 20 and bool(common)
                                     and np.median(common) >= 5) else "no")
                        + ", CIIP-1B "
                        + ("yes" if len(parents_2plus) >= 10 else "no")
                        + f" ({len(parents_2plus)} parents)"),
            "anastassiadis": "cross-endpoint replication reference only",
            "cross_panel_pooling": "FORBIDDEN (Kd vs % inhibition are different endpoints)",
        },
        "verdict": None,
    }
    # verdict per frozen stop rules
    if out["duongly"]["ciiP1B_admissible"]:
        out["verdict"] = "SUFFICIENT: Duong-Ly alone supports CIIP-1A and CIIP-1B"
    elif out["duongly"]["ciiP1A_admissible"]:
        out["verdict"] = ("PARTIAL: Duong-Ly supports CIIP-1A; CIIP-1B "
                          "UNRESOLVED (held-out parents < 10; no admissible "
                          "pooling with Davis)")
    else:
        out["verdict"] = "UNRESOLVED/DATA-BLOCKED"
    def _ser(o):
        if isinstance(o, (np.floating, np.integer)):
            return float(o) if isinstance(o, np.floating) else int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return str(o)
    (HERE / "CENSUS.json").write_text(
        json.dumps(out, indent=1, default=_ser), encoding="utf-8")
    print("duongly: pairs", len(pairs), "median ligands", np.median(common),
          "parents2plus", len(parents_2plus), "per-parent", out["duongly"]["per_parent_pairs"])
    print("unparsed:", dl["unparsed_rows"])
    print("multi tags:", dl["multi_mutant_tags"])
    print("VERDICT:", out["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
