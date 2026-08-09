"""Phase 2A / Phase 5 — label-semantics audit.

Registered by PREREG_S7_L2B_PHASE2A.md (sha 4e01401d...), section 9, plus
computational amendment 03 (the dense-distance comparator).

Answers, using only already-available evidence:
  - what interaction types make up the teacher, and how much of it is INDIRECT
    (water- or metal-mediated) rather than a direct residue-atom contact;
  - how reproducible the teacher is (replicate floor, carried from Phase 2);
  - whether a locally computable dense-distance teacher agrees with it;
  - how the answer depends on the distance threshold;
  - the registered sensitivity: does the T1 teacher-conditionality conclusion
    survive removal of indirect edges?

Nothing is trained. No affinity source is opened.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import component_macro, paired_bootstrap  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
MMCIF_DIRS = [ROOT / "dataset/raw/open_structures/pilot20k/mmcif",
              ROOT / "dataset/raw/open_structures/pilot15k/mmcif",
              ROOT / "dataset/raw/ssl_b2_independent/mmcif"]
PREREG_SHA = "4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e"

INDIRECT_TYPES = {"Water Bridges"}
MEDIATED_SECONDARY = {"Metal Complexes"}
THRESHOLDS = (4.0, 4.5, 5.0)
MIN_MATCHED = 50
MIN_IDENTITY = 0.95
SEED_BOOT = 20260819
DJ_MIN = 0.05


def map_chain(st, seq, gemmi):
    """Exhaustive integer-offset scan (amendment C2). Returns (offset, chain,
    matched, identity) or None."""
    seq_arr = np.frombuffer(seq.encode(), dtype="S1")
    L = seq_arr.size
    best = None
    for ch in st[0]:
        idx, aa = [], []
        for res in ch:
            if res.label_seq is None:
                continue
            t = gemmi.find_tabulated_residue(res.name)
            if t is None or not t.is_amino_acid():
                continue
            idx.append(res.label_seq)
            aa.append(t.one_letter_code.upper())
        if len(idx) < MIN_MATCHED:
            continue
        ii = np.asarray(idx)
        av = np.frombuffer("".join(aa).encode(), dtype="S1")
        lo, hi = 1 - int(ii.max()), L - int(ii.min())
        for off in range(lo, hi + 1):
            k = ii - 1 + off
            m = (k >= 0) & (k < L)
            nm = int(m.sum())
            if nm < MIN_MATCHED:
                continue
            hit = int((seq_arr[k[m]] == av[m]).sum())
            if hit / nm < MIN_IDENTITY:
                continue
            if best is None or hit > best[2] or (hit == best[2] and ch.name < best[1]):
                best = (off, ch.name, hit, hit / nm)
    return best


def main():
    t0 = time.time()
    kept, quarantine, contract, _f = build()
    comp_of = protein_components(kept)
    _tr, _ha, held_A, _hb = make_split(kept, comp_of)

    # ------------------------------------------------- interaction-type census
    type_edges = Counter()
    type_complexes = Counter()
    total_typed = 0
    per_complex_types = []
    for r in kept:
        ts = set()
        for e in r.get("positive_typed_edges", []):
            type_edges[e[2]] += 1
            ts.add(e[2])
            total_typed += 1
        for t in ts:
            type_complexes[t] += 1
        per_complex_types.append(len(ts))
    indirect_edges = sum(v for k, v in type_edges.items() if k in INDIRECT_TYPES)
    mediated_edges = sum(v for k, v in type_edges.items() if k in MEDIATED_SECONDARY)

    # -------------------------------------------------------------- coverage
    pos = np.array([len(r["edges"]) for r in kept], dtype=float)
    cells = np.array([r["n_res"] * r["n_atoms"] for r in kept], dtype=float)
    res_hit = np.array([len({i for i, _ in r["edges"]}) for r in kept], dtype=float)
    atom_hit = np.array([len({j for _, j in r["edges"]}) for r in kept], dtype=float)
    nres = np.array([r["n_res"] for r in kept], dtype=float)
    natom = np.array([r["n_atoms"] for r in kept], dtype=float)

    # -------------- registered sensitivity: T1 with indirect edges removed
    def masks_from(rec, drop_types, source):
        """source='validated' uses the contract-validated binary edges; 'typed'
        uses positive_typed_edges. The sensitivity must compare typed-vs-typed,
        never validated-vs-typed, or the edge source would confound the effect
        of removing indirect interactions."""
        if source == "validated":
            return frozenset(i for i, _j in rec["edges"])
        keep = set()
        for e in rec.get("positive_typed_edges", []):
            if e[2] not in drop_types:
                keep.add(int(e[0]))
        return frozenset(keep)

    def t1_for(drop, source="typed"):
        by_c = defaultdict(list)
        for r in kept:
            by_c[r["seq_key"]].append(r)
        rep, alt = defaultdict(list), defaultdict(list)
        for sk, recs in by_c.items():
            if len(recs) < 2:
                continue
            mk = {r["source_key"]: masks_from(r, drop, source) for r in recs}
            for a, b in combinations(recs, 2):
                ra, rb = mk[a["source_key"]], mk[b["source_key"]]
                if not ra or not rb:
                    continue
                J = len(ra & rb) / len(ra | rb)
                c = comp_of[a["source_key"]]
                if a["graph_key"] == b["graph_key"]:
                    if a["pdb_id"] != b["pdb_id"]:
                        rep[c].append(J)
                else:
                    sa, sb = a["scaffold"], b["scaffold"]
                    if sa and sb and sa != sb:
                        alt[c].append(J)
        cr = {c: float(np.mean(v)) for c, v in rep.items() if v}
        ca = {c: float(np.mean(v)) for c, v in alt.items() if v}
        bs = paired_bootstrap(cr, ca, n_boot=10000, seed=SEED_BOOT)
        bs["replicate_mean"] = float(np.mean(list(cr.values())))
        bs["alternative_mean"] = float(np.mean(list(ca.values())))
        bs["pass"] = bool(bs["delta"] >= DJ_MIN and bs["lcb95_one_sided"] > 0)
        return bs

    t1_validated = t1_for(frozenset(), source="validated")   # cross-check only
    t1_all = t1_for(frozenset(), source="typed")             # sensitivity baseline
    t1_direct = t1_for(INDIRECT_TYPES, source="typed")
    t1_direct_nometal = t1_for(INDIRECT_TYPES | MEDIATED_SECONDARY, source="typed")
    print(f"T1 sensitivity done {time.time()-t0:.0f}s", flush=True)

    # ----------------------------------------- dense-distance comparator
    import gemmi
    local = {}
    for d in MMCIF_DIRS:
        if d.is_dir():
            for f in os.listdir(d):
                local.setdefault(f.split(".")[0].lower(), d / f)
    cand = [r for r in kept if r["pdb_id"].lower() in local]
    print(f"mmCIF comparator candidates: {len(cand)}", flush=True)

    geo = {t: {"jac": {}, "prec": {}, "rec": {}} for t in THRESHOLDS}
    beyond5 = []
    identities = []
    fail = Counter()
    used = 0
    for n_i, r in enumerate(cand):
        try:
            st = gemmi.read_structure(str(local[r["pdb_id"].lower()]))
            st.setup_entities()
            st.remove_alternative_conformations()
            st.remove_hydrogens()
        except Exception:
            fail["parse"] += 1
            continue
        mp = map_chain(st, r["uniprot_sequence"], gemmi)
        if mp is None:
            fail["no_mapped_chain"] += 1
            continue
        off, chname, matched, ident = mp
        identities.append(ident)
        chain = None
        for ch in st[0]:
            if ch.name == chname:
                chain = ch
                break
        if chain is None:
            fail["chain_lost"] += 1
            continue
        prot = []
        for res in chain:
            if res.label_seq is None:
                continue
            t = gemmi.find_tabulated_residue(res.name)
            if t is None or not t.is_amino_acid():
                continue
            k = res.label_seq - 1 + off
            if not (0 <= k < r["n_res"]):
                continue
            for at in res:
                if at.element == gemmi.Element("H"):
                    continue
                prot.append((k, at.pos.x, at.pos.y, at.pos.z))
        if len(prot) < 50:
            fail["too_few_protein_atoms"] += 1
            continue
        pk = np.array([p[0] for p in prot])
        pc = np.array([[p[1], p[2], p[3]] for p in prot])

        copies = []
        for ch in st[0]:
            for res in ch:
                if res.name.upper() != r["ligand_ccd"].upper():
                    continue
                co = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in res
                               if a.element != gemmi.Element("H")])
                if co.size:
                    copies.append(co)
        if not copies:
            fail["ligand_copy_absent"] += 1
            continue
        best_copy, best_n = None, -1
        for co in copies:
            d = np.linalg.norm(pc[:, None, :] - co[None, :, :], axis=2)
            n = int((d.min(1) <= 4.5).sum())
            if n > best_n:
                best_n, best_copy = n, co
        d = np.linalg.norm(pc[:, None, :] - best_copy[None, :, :], axis=2).min(1)
        plip = set(i for i, _j in r["edges"])
        key = r["source_key"]
        for t in THRESHOLDS:
            g = set(pk[d <= t].tolist())
            if not g and not plip:
                continue
            inter = len(plip & g)
            geo[t]["jac"][key] = inter / max(len(plip | g), 1)
            geo[t]["prec"][key] = inter / max(len(plip), 1)
            geo[t]["rec"][key] = inter / max(len(g), 1)
        g5 = set(pk[d <= 5.0].tolist())
        beyond5.append(len(plip - g5) / max(len(plip), 1))
        used += 1
        if (n_i + 1) % 250 == 0:
            print(f"  geo {n_i+1}/{len(cand)} used={used} {time.time()-t0:.0f}s", flush=True)

    def mac(d):
        return component_macro(d, comp_of)[1]

    geo_out = {}
    for t in THRESHOLDS:
        geo_out[f"{t:.1f}A"] = {
            "complexes": len(geo[t]["jac"]),
            "jaccard_component_macro": mac(geo[t]["jac"]) if geo[t]["jac"] else None,
            "plip_covered_by_geometry_component_macro": mac(geo[t]["prec"]) if geo[t]["prec"] else None,
            "geometry_covered_by_plip_component_macro": mac(geo[t]["rec"]) if geo[t]["rec"] else None,
        }

    ambiguity_indirect = indirect_edges / max(total_typed, 1)
    t1_reverses = bool(t1_all["pass"] and not t1_direct["pass"])
    beyond5_mean = float(np.mean(beyond5)) if beyond5 else None
    disagreement_20 = bool(beyond5_mean is not None and beyond5_mean >= 0.20)
    ambiguous = bool((ambiguity_indirect >= 0.20 and t1_reverses) or disagreement_20)

    res = {
        "schema": "MetaSieve.S7L2B.P2A.LabelSemantics.v1",
        "created_utc": "2026-08-10",
        "preregistration_sha256": PREREG_SHA,
        "amendments": ["PREREG_S7_L2B_PHASE2A_AMENDMENT_01.md",
                       "PREREG_S7_L2B_PHASE2A_AMENDMENT_02.md",
                       "PREREG_S7_L2B_PHASE2A_AMENDMENT_03.md"],
        "repo_commit": "623602e76b7d4f445af069014782278163183d59",
        "interaction_type_census": {
            "total_typed_edges": total_typed,
            "edges_by_type": dict(type_edges.most_common()),
            "complexes_by_type": dict(type_complexes.most_common()),
            "indirect_water_mediated_edge_fraction": ambiguity_indirect,
            "metal_mediated_edge_fraction": mediated_edges / max(total_typed, 1),
            "mean_types_per_complex": float(np.mean(per_complex_types)),
        },
        "coverage": {
            "complexes": len(kept),
            "mean_positives_per_complex": float(pos.mean()),
            "median_positives_per_complex": float(np.median(pos)),
            "overall_positive_density": float(pos.sum() / cells.sum()),
            "mean_fraction_of_residues_touched": float((res_hit / nres).mean()),
            "mean_fraction_of_atoms_touched": float((atom_hit / natom).mean()),
        },
        "mapping_failures_from_contract": {
            "atom_contract_census": contract,
            "records_quarantined": len(quarantine),
            "quarantine_reasons": dict(Counter(q["reason"] for q in quarantine)),
        },
        "registered_sensitivity_T1_without_indirect_edges": {
            "note": "baseline and sensitivity are both built from "
                    "positive_typed_edges so the edge source cannot confound the "
                    "effect of removing indirect interactions; the "
                    "contract-validated-edge value is reported only as a cross-check",
            "validated_edges_cross_check": t1_validated,
            "all_edges": t1_all,
            "water_bridges_removed": t1_direct,
            "water_and_metal_removed": t1_direct_nometal,
            "conclusion_reverses": t1_reverses,
        },
        "dense_distance_comparator": {
            "source": "local mmCIF acquired under earlier governed stages",
            "candidates_with_local_coordinates": len(cand),
            "complexes_used": used,
            "mapping_failures": dict(fail),
            "sequence_identity_median": float(np.median(identities)) if identities else None,
            "agreement_by_threshold": geo_out,
            "plip_positives_beyond_5A_mean_fraction": beyond5_mean,
            "reading": "PLIP requires geometric AND chemical criteria, so R_plip "
                       "being a strict subset of R_geom is expected and is NOT "
                       "evidence of missing positives. Only PLIP positives with no "
                       "heavy atom within 5.0 A would indicate a mapping or label "
                       "defect.",
        },
        "second_frozen_interaction_tool": "ABSENT — no second interaction-annotation "
                                          "tool output exists locally for these "
                                          "complexes; that specific comparison remains "
                                          "UNRESOLVED",
        "ambiguity_verdict": {
            "criterion": "indirect fraction >= 0.20 AND T1 reverses; OR a second "
                         "teacher disagrees on >= 0.20 of edges",
            "indirect_fraction": ambiguity_indirect,
            "T1_reverses": t1_reverses,
            "plip_beyond_5A_disagreement_ge_20pct": disagreement_20,
            "LABEL_SEMANTICS_AMBIGUOUS": ambiguous,
        },
        "elapsed_sec": round(time.time() - t0, 1),
    }
    (OUT / "PHASE2A_LABEL_SEMANTICS.json").write_text(json.dumps(res, indent=2),
                                                      encoding="utf-8")
    print(json.dumps({k: res[k] for k in
                      ("interaction_type_census", "coverage",
                       "registered_sensitivity_T1_without_indirect_edges",
                       "dense_distance_comparator", "ambiguity_verdict")},
                     indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
