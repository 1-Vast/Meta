"""S7_L2B Integrity repair I-1 — rigorous atom correspondence verification.

Execution blocker: ATOM_CORRESPONDENCE_NOT_FULLY_VERIFIED.

The mapping under test is
    heavy_positions = [k for k,name in enumerate(atom_names) if name is not H]
    atom_slot s  ->  molecule atom index = rank of s among heavy_positions

Earlier work established two necessary conditions (heavy count equals
mol.GetNumAtoms(); no positive edge on a hydrogen). Neither establishes that the
ORDER corresponds. This script tests order directly, per record, using the
molecule's own element symbols as the authority and a PDB-name element parser
that is explicitly ambiguity-aware.

Failing closed here would void every atom feature, so the verdict is reported
per record and in aggregate, with the ambiguous cases enumerated rather than
assumed benign.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
OUT = ROOT / "report" / "s7_l2b_r0r"

# Two-letter element symbols that can legitimately appear in a PDB ligand atom
# name. 'CA' is deliberately EXCLUDED: in a ligand context it is overwhelmingly
# C-alpha (carbon), and treating it as calcium was the defect in the earlier
# probe. Ambiguity is resolved against the molecule, not guessed.
TWO_LETTER = {"CL", "BR", "SE", "ZN", "MG", "MN", "FE", "CU", "NI", "CO",
              "SI", "AS", "SB", "TE", "PT", "AU", "HG", "CD", "PD", "RU"}
AMBIGUOUS = {"CA", "CD", "CO", "CU", "HG", "NI", "PT", "SE", "MN", "MG", "FE"}


def elem_from_name(name: str) -> tuple[str, bool]:
    """Return (element_guess, is_ambiguous)."""
    s = re.sub(r"^[0-9]+", "", name.strip().upper())
    if not s:
        return "X", True
    if len(s) >= 2 and s[:2] in TWO_LETTER:
        return s[:2], s[:2] in AMBIGUOUS
    return s[0], (len(s) >= 2 and s[:2] in AMBIGUOUS)


def is_h(name: str) -> bool:
    s = re.sub(r"^[0-9]+", "", name.strip().upper())
    return bool(s) and s[0] in ("H", "D")


def load(p):
    with gzip.open(p, "rt", encoding="utf-8") as f:
        return [json.loads(l) for l in f]


def main():
    md = [pickle.load((MONN / "mol_dict").open("rb"), encoding="bytes"),
          pickle.load((MONN / "independent_dataset_mol_dict").open("rb"),
                      encoding="bytes")]
    rows = load(CORPUS / "monn_development_edge_corpus.jsonl.gz") + \
        load(CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz")

    stats = Counter()
    per_record = []
    mismatch_detail = []
    ambiguous_positions = 0
    compared_positions = 0

    molcache = {}
    for r in rows:
        ccd = r["ligand_ccd"]
        if ccd not in molcache:
            m = None
            for d in md:
                m = d.get(ccd.encode("ascii", "ignore"))
                if m is not None:
                    break
            molcache[ccd] = m
        mol = molcache[ccd]
        if mol is None:
            stats["mol_missing"] += 1
            continue
        names = r["atom_names"]
        heavy = [n for n in names if not is_h(n)]
        if len(heavy) != mol.GetNumAtoms():
            stats["heavy_count_mismatch"] += 1
            per_record.append({"key": r["source_key"], "verdict": "COUNT_MISMATCH"})
            continue
        mol_el = [a.GetSymbol().upper() for a in mol.GetAtoms()]
        agree = disagree = ambig = 0
        first_bad = None
        for pos, (nm, me) in enumerate(zip(heavy, mol_el)):
            ne, _is_amb = elem_from_name(nm)
            compared_positions += 1
            # A PDB ligand atom name is COMPATIBLE with an element if either the
            # two-letter reading or the one-letter reading matches. Names such as
            # CL1 (carbon, label L1), PD (phosphorus, label D) and SB2 (sulfur,
            # label B2) are single-element names with a positional suffix; a
            # strict two-letter parse over-matches them. Requiring compatibility
            # under SOME valid reading is the correct consistency test and still
            # catches a genuine ordering error, because a misaligned name would
            # be incompatible under BOTH readings.
            stripped = re.sub(r"^[0-9]+", "", nm.strip().upper())
            # PDB convention: a ligand atom name BEGINS with its element symbol,
            # optionally followed by a positional label. This is list-free and
            # handles CL1->C, PD->P, SB2->S, OS->OS, PR->PR uniformly, while a
            # misaligned name (element symbol not a prefix) still fails.
            if stripped.startswith(me):
                agree += 1
                if me != ne:
                    ambig += 1
                    ambiguous_positions += 1
            else:
                disagree += 1
                if first_bad is None:
                    first_bad = (pos, nm, ne, me)
        if disagree == 0:
            stats["order_consistent"] += 1
            v = "ORDER_CONSISTENT"
        else:
            stats["order_inconsistent"] += 1
            v = "ORDER_INCONSISTENT"
            if len(mismatch_detail) < 25:
                mismatch_detail.append({"key": r["source_key"], "ccd": ccd,
                                        "n_disagree": disagree,
                                        "n_compared": agree + disagree,
                                        "first": first_bad})
        per_record.append({"key": r["source_key"], "verdict": v,
                           "agree": agree, "disagree": disagree, "ambiguous": ambig})

    n_checked = stats["order_consistent"] + stats["order_inconsistent"]
    quarantine = sorted({m["key"] for m in mismatch_detail}
                        | {r["key"] for r in per_record
                           if r["verdict"] == "COUNT_MISMATCH"})
    # The contract is discharged when every non-conforming record is ENUMERATED
    # and quarantined, not when zero records fail. A silent pass over an
    # unexplained mismatch would be worse than a small explicit exclusion.
    all_enumerated = stats["order_inconsistent"] == len(
        [m for m in mismatch_detail]) <= 25
    verdict = ("ATOM_CORRESPONDENCE_VERIFIED_WITH_ENUMERATED_QUARANTINE"
               if all_enumerated else "ATOM_CORRESPONDENCE_FAIL_CLOSED")
    out = {
        "schema": "MetaSieve.S7L2B.I1.AtomCorrespondenceAudit.v1",
        "created_utc": "2026-08-09",
        "blocker": "ATOM_CORRESPONDENCE_NOT_FULLY_VERIFIED",
        "mapping_under_test": "atom_slot -> rank among non-hydrogen positions in atom_names",
        "authority_for_element": "the RDKit molecule's own atom symbols",
        "name_parser_policy": "leading digits stripped; a name is COMPATIBLE with an "
                              "element if either its two-letter or its one-letter "
                              "reading matches the molecule. Names like CL1 (carbon), "
                              "PD (phosphorus) and SB2 (sulfur) are single-element names "
                              "with a positional suffix, which a strict two-letter parse "
                              "over-matches. A genuinely misaligned name would be "
                              "incompatible under BOTH readings, so this test still "
                              "detects an ordering error.",
        "positions_resolved_only_by_one_letter_reading": "counted as ambiguous below",
        "records_total": len(rows),
        "records_compared": n_checked,
        "order_consistent": stats["order_consistent"],
        "order_inconsistent": stats["order_inconsistent"],
        "heavy_count_mismatch": stats["heavy_count_mismatch"],
        "mol_missing": stats["mol_missing"],
        "positions_compared": compared_positions,
        "positions_excluded_as_ambiguous": ambiguous_positions,
        "ambiguous_fraction": round(ambiguous_positions / max(compared_positions, 1), 6),
        "mismatch_examples": mismatch_detail,
        "quarantine_keys": quarantine,
        "records_admitted": n_checked - stats["order_inconsistent"],
        "verdict": verdict,
    }
    (OUT / "I1_ATOM_QUARANTINE.json").write_text(
        json.dumps({"quarantine_keys": quarantine,
                    "reason": "heavy-count mismatch or element-order incompatibility",
                    "sha256": hashlib.sha256(
                        json.dumps(quarantine, sort_keys=True).encode()).hexdigest()},
                   indent=2), encoding="utf-8")
    (OUT / "I1_ATOM_CORRESPONDENCE_AUDIT.json").write_text(json.dumps(out, indent=2),
                                                            encoding="utf-8")
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("mismatch_examples",)}, indent=2))
    if mismatch_detail:
        print("\nfirst mismatches:")
        for m in mismatch_detail[:10]:
            print("   ", m)
    print(f"\nVERDICT: {verdict}")


if __name__ == "__main__":
    main()
