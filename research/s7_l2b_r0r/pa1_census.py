"""Phase 2A / Phase 1 — data identifiability census. Label-blind power.

Registered by research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md
(SHA-256 4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e), s.5.

Answers: at the EXACT construct level (identical UniProt sequence, so residue
indices are comparable), how many groups carry at least two scaffold-distinct
ligands, how deep are they, how independent are they under protein closure, and
is there enough of that structure to detect the preregistered minimum meaningful
effect dJ_min = 0.05 in Jaccard units?

Opens no affinity source.
"""
from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
RCSB = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "rcsb" / "rcsb_normalized.json"
PREREG_SHA = "4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e"

DJ_MIN = 0.05
SIGMA_GRID = (0.10, 0.15, 0.20, 0.25, 0.30)
ALPHA = 0.05
D1_MIN_COMPONENTS = 30
D2_MIN_PAIRS = 100


def power_two_sided(delta, sigma, n, alpha=ALPHA):
    """Normal-approximation power for a mean-difference test over n independent
    inference units. Label-blind: uses only n and an ASSUMED sigma."""
    from math import erf, sqrt
    if n < 2 or sigma <= 0:
        return 0.0
    zc = 1.959963984540054 if abs(alpha - 0.05) < 1e-12 else None
    if zc is None:
        raise ValueError("only alpha=0.05 is registered")
    lam = delta / (sigma / sqrt(n))
    cdf = lambda x: 0.5 * (1.0 + erf(x / sqrt(2.0)))
    return float((1.0 - cdf(zc - lam)) + cdf(-zc - lam))


def required_n(delta, sigma, target=0.80, alpha=ALPHA):
    for n in range(2, 100001):
        if power_two_sided(delta, sigma, n, alpha) >= target:
            return n
    return None


def main():
    kept, _q, _c, _f = build()
    comp_of = protein_components(kept)
    train, _held_all, held_A, _held_B = make_split(kept, comp_of)
    heldA_keys = {r["source_key"] for r in held_A}

    rcsb = json.loads(RCSB.read_text(encoding="utf-8"))["records"]

    by_construct = defaultdict(list)
    for r in kept:
        by_construct[r["seq_key"]].append(r)

    proteins = {r["uniprot_id"] for r in kept}
    print(f"records={len(kept)} constructs={len(by_construct)} uniprot={len(proteins)}",
          flush=True)

    # ------------------------------------------------ per-construct structure
    groups = {}
    for sk, recs in by_construct.items():
        ccds = {r["ligand_ccd"] for r in recs}
        gks = {r["graph_key"] for r in recs}
        scafs = {r["scaffold"] for r in recs if r["scaffold"]}
        comps = {comp_of[r["source_key"]] for r in recs}
        years = []
        for r in recs:
            md = rcsb.get(r["pdb_id"].upper())
            if md and md.get("release_date"):
                years.append(int(md["release_date"][:4]))
        groups[sk] = {
            "records": len(recs),
            "uniprot": sorted({r["uniprot_id"] for r in recs}),
            "n_ccd": len(ccds), "n_graph_key": len(gks), "n_scaffold": len(scafs),
            "components": sorted(comps),
            "n_in_heldA": sum(1 for r in recs if r["source_key"] in heldA_keys),
            "release_years": years,
            "all_have_mask": all(len(r["edges"]) > 0 for r in recs),
        }

    # ------------------------------------------------------------ pair census
    pair_counts = Counter()
    scaffold_distinct_pairs = 0
    alt_pairs = 0
    replicate_pairs = 0
    comps_with_scaffold_distinct = set()
    comps_with_replicate = set()
    constructs_with_scaffold_distinct = []
    constructs_with_replicate = []
    heldA_scaffold_distinct_pairs = 0
    heldA_replicate_pairs = 0
    per_construct_pairs = {}

    for sk, recs in by_construct.items():
        if len(recs) < 2:
            continue
        n_rep = n_alt = n_sd = 0
        n_rep_A = n_sd_A = 0
        for a, b in combinations(recs, 2):
            same_graph = a["graph_key"] == b["graph_key"]
            both_A = (a["source_key"] in heldA_keys) and (b["source_key"] in heldA_keys)
            if same_graph:
                if a["pdb_id"] != b["pdb_id"]:
                    n_rep += 1
                    if both_A:
                        n_rep_A += 1
                    pair_counts["replicate"] += 1
                else:
                    pair_counts["same_pdb_same_graph"] += 1
            else:
                n_alt += 1
                pair_counts["alternative_ligand"] += 1
                sa, sb = a["scaffold"], b["scaffold"]
                if sa and sb and sa != sb:
                    n_sd += 1
                    if both_A:
                        n_sd_A += 1
                    pair_counts["scaffold_distinct"] += 1
        replicate_pairs += n_rep
        alt_pairs += n_alt
        scaffold_distinct_pairs += n_sd
        heldA_scaffold_distinct_pairs += n_sd_A
        heldA_replicate_pairs += n_rep_A
        per_construct_pairs[sk] = {"replicate": n_rep, "alternative": n_alt,
                                   "scaffold_distinct": n_sd}
        cs = {comp_of[r["source_key"]] for r in recs}
        if n_sd > 0:
            comps_with_scaffold_distinct |= cs
            constructs_with_scaffold_distinct.append(sk)
        if n_rep > 0:
            comps_with_replicate |= cs
            constructs_with_replicate.append(sk)

    # ------------------------------------------------------- closure topology
    comp_sizes = Counter(comp_of[r["source_key"]] for r in kept)
    largest = max(comp_sizes.values())
    both = comps_with_scaffold_distinct & comps_with_replicate

    # depth distribution over constructs that carry scaffold-distinct pairs
    depth_lig = [groups[sk]["n_graph_key"] for sk in constructs_with_scaffold_distinct]
    depth_scf = [groups[sk]["n_scaffold"] for sk in constructs_with_scaffold_distinct]

    def dist(v):
        if not v:
            return {}
        v = np.asarray(v, dtype=float)
        return {"n": int(v.size), "min": float(v.min()), "p25": float(np.percentile(v, 25)),
                "median": float(np.median(v)), "p75": float(np.percentile(v, 75)),
                "max": float(v.max()), "mean": float(v.mean())}

    all_years = [y for g in groups.values() for y in g["release_years"]]
    yr_sd = [y for sk in constructs_with_scaffold_distinct for y in groups[sk]["release_years"]]

    # -------------------------------------------------------------- power
    n_comp_paired = len(both)
    n_comp_sd = len(comps_with_scaffold_distinct)
    power = {}
    for s in SIGMA_GRID:
        power[f"sigma_{s:.2f}"] = {
            "achieved_power_at_paired_components": power_two_sided(DJ_MIN, s, n_comp_paired),
            "achieved_power_at_scaffold_distinct_components": power_two_sided(DJ_MIN, s, n_comp_sd),
            "components_required_for_80pct": required_n(DJ_MIN, s),
        }

    d1 = n_comp_sd >= D1_MIN_COMPONENTS
    d2 = scaffold_distinct_pairs >= D2_MIN_PAIRS
    d3 = all(g["all_have_mask"] for g in groups.values())

    res = {
        "schema": "MetaSieve.S7L2B.P2A.DataIdentifiabilityCensus.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "repo_commit": "623602e76b7d4f445af069014782278163183d59",
        "scope": "label-side census over the FULL admitted corpus; held-out A "
                 "subsets reported separately",
        "corpus": {"records": len(kept), "exact_constructs": len(by_construct),
                   "uniprot_ids": len(proteins), "train": len(train),
                   "heldout_A": len(held_A)},
        "closure": {
            "components_total": len(comp_sizes),
            "largest_component_records": largest,
            "largest_component_fraction": float(largest / len(kept)),
            "components_with_scaffold_distinct_pair": n_comp_sd,
            "components_with_replicate_pair": len(comps_with_replicate),
            "components_with_both_pair_types": n_comp_paired,
        },
        "constructs": {
            "with_at_least_two_graph_keys": sum(1 for g in groups.values() if g["n_graph_key"] >= 2),
            "with_at_least_two_scaffolds": sum(1 for g in groups.values() if g["n_scaffold"] >= 2),
            "with_scaffold_distinct_pair": len(constructs_with_scaffold_distinct),
            "with_replicate_pair": len(constructs_with_replicate),
            "singleton": sum(1 for g in groups.values() if g["records"] == 1),
            "ligand_depth_distribution_in_scaffold_distinct_constructs": dist(depth_lig),
            "scaffold_depth_distribution_in_scaffold_distinct_constructs": dist(depth_scf),
        },
        "pairs": {
            "replicate_same_construct_same_ligand_diff_pdb": replicate_pairs,
            "alternative_ligand_same_construct": alt_pairs,
            "scaffold_distinct_same_construct": scaffold_distinct_pairs,
            "same_pdb_same_graph_key": pair_counts["same_pdb_same_graph"],
            "heldout_A_only_replicate": heldA_replicate_pairs,
            "heldout_A_only_scaffold_distinct": heldA_scaffold_distinct_pairs,
            "overlap_levels": "exact ligand = ligand_ccd; connectivity key = "
                              "graph_key (canonical non-isomeric SMILES hash); "
                              "scaffold = Bemis-Murcko",
        },
        "publication_time": {
            "records_with_release_year": len(all_years),
            "release_year_distribution_all": dist(all_years),
            "release_year_distribution_scaffold_distinct_constructs": dist(yr_sd),
            "release_year_histogram_decade": dict(sorted(
                Counter((y // 5) * 5 for y in all_years).items())),
        },
        "residue_mask_availability": {
            "constructs_all_records_masked": sum(1 for g in groups.values() if g["all_have_mask"]),
            "constructs_total": len(groups),
            "all_masked": d3,
        },
        "within_construct_guarantee": "every comparison in Phase 2 is within one "
                                      "seq_key, so residue indices are identical by "
                                      "construction; cross-construct residue "
                                      "comparison is never performed",
        "label_blind_power": {
            "minimum_meaningful_effect_dJ": DJ_MIN,
            "alpha_two_sided": ALPHA,
            "assumed_sigma_grid": list(SIGMA_GRID),
            "note": "no observed label variance enters this computation",
            "results": power,
        },
        "sufficiency": {
            "D1_components_ge_30": {"value": n_comp_sd, "threshold": D1_MIN_COMPONENTS,
                                    "pass": bool(d1)},
            "D2_scaffold_distinct_pairs_ge_100": {"value": scaffold_distinct_pairs,
                                                  "threshold": D2_MIN_PAIRS, "pass": bool(d2)},
            "D3_all_records_masked": {"pass": bool(d3)},
            "verdict": ("DATA_IDENTIFIABLE" if (d1 and d2 and d3)
                        else "PHASE2A_DATA_NOT_IDENTIFIABLE"),
        },
    }
    (OUT / "PHASE2A_DATA_IDENTIFIABILITY_CENSUS.json").write_text(
        json.dumps(res, indent=2), encoding="utf-8")
    # construct-level table for Phase 2 to consume, keeps the pair definition frozen
    (OUT / "PHASE2A_CONSTRUCT_GROUPS.json").write_text(
        json.dumps({"per_construct_pairs": per_construct_pairs,
                    "constructs_with_scaffold_distinct": constructs_with_scaffold_distinct,
                    "constructs_with_replicate": constructs_with_replicate}, indent=1),
        encoding="utf-8")
    print(json.dumps({k: res[k] for k in
                      ("closure", "constructs", "pairs", "sufficiency")}, indent=2))
    print(json.dumps(res["label_blind_power"]["results"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
