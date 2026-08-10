"""Registered mapping-equivalence check for rule M4 (prereg section 4.2).

M4 replaces P1B's parasail alignment with `sequence_index = label_seq_id - 1`.
The preregistration requires this to be verified against the parasail mapping on
already-exposed structures and to fail closed otherwise. Exposed structures are
used deliberately: the untouched corpus must not be consumed by an engineering
check.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "correspondence_router"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))

from c0_corpus import (AMENDMENT, AMENDMENT_SHA, ATOM_TAGS, OUT,  # noqa: E402
                       PREREG, PREREG_SHA, RAW, SLOTS, C0ContractError,
                       _canonical_rows, git_head, read_jsonl, sequence_mapping,
                       sha_file, write_json)

GOVERNED = (ROOT / "dataset" / "processed" / "open_structures" /
            "pilot20k_holo_governed_v2" / "complexes.jsonl")
MMCIF = RAW / "mmcif"


def parasail_mapping(rows, sequence):
    """The exact P1B path, imported in spirit from scripts/build_holo_complex_index."""
    import gemmi
    import parasail
    names = {}
    for row in rows:
        names.setdefault(row["label_seq_id"], row["label_comp_id"])
    label_seq_ids = sorted(names, key=int)
    structure_sequence = ""
    for label in label_seq_ids:
        code = gemmi.find_tabulated_residue(names[label]).one_letter_code
        structure_sequence += code.upper() if code and code != " " else "X"
    alignment = parasail.nw_trace_striped_16(
        structure_sequence, sequence, 10, 1, parasail.blosum62).traceback
    structure_index = sequence_index = 0
    mapped = {}
    for structure_code, sequence_code in zip(alignment.query, alignment.ref):
        if structure_code != "-" and sequence_code != "-":
            mapped[label_seq_ids[structure_index]] = sequence_index
        if structure_code != "-":
            structure_index += 1
        if sequence_code != "-":
            sequence_index += 1
    return mapped


def run(sample: int = 40) -> dict:
    import gemmi
    if sha_file(PREREG) != PREREG_SHA or sha_file(AMENDMENT) != AMENDMENT_SHA:
        raise C0ContractError("C0/C1 preregistration or amendment hash mismatch")
    checked = agree = slot_agree = 0
    disagreements, skipped = [], 0
    for record in read_jsonl(GOVERNED):
        if checked >= sample:
            break
        pdb_id = record["pdb_id"].lower()
        path = MMCIF / f"{pdb_id}.cif.gz"
        if not path.is_file():
            skipped += 1
            continue
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            block = gemmi.cif.read_string(handle.read()).sole_block()
        table = block.find("_atom_site.", list(ATOM_TAGS))
        atoms = [{tag: str(row[i]) for i, tag in enumerate(ATOM_TAGS)} for row in table]
        rows = [a for a in atoms if a["group_PDB"] == "ATOM"
                and a["label_asym_id"] == record["protein_asym_id"]
                and a["label_seq_id"] not in {"", ".", "?"}
                and a["type_symbol"].upper() != "H"]
        rows = _canonical_rows(rows, ("label_asym_id", "label_seq_id", "label_atom_id"))
        if not rows:
            skipped += 1
            continue
        sequence = record["sequence"]
        length = len(sequence)
        reference = parasail_mapping(rows, sequence)
        checked += 1
        labels, indices, _coverage = sequence_mapping(rows, sequence)
        candidate = dict(zip(labels, indices))
        if candidate == reference:
            agree += 1
        else:
            disagreements.append({"pdb_id": pdb_id, "differing": sum(
                1 for k in reference if reference[k] != candidate.get(k))})
        reference_slots = {k: min(SLOTS - 1, v * SLOTS // length)
                           for k, v in reference.items()}
        candidate_slots = {k: min(SLOTS - 1, v * SLOTS // length)
                           for k, v in candidate.items()}
        if reference_slots == candidate_slots:
            slot_agree += 1

    result = {
        "schema": "MetaSieve.Correspondence.C0.MappingEquivalence.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "panel": "already-exposed pilot20k governed complexes; the untouched "
                 "corpus is deliberately not consumed by this check",
        "structures_checked": checked, "structures_skipped": skipped,
        "sequence_index_agreement": agree,
        "slot_assignment_agreement": slot_agree,
        "disagreements": disagreements[:20],
        "pass": bool(checked > 0 and slot_agree == checked),
    }
    write_json(OUT / "C0_MAPPING_EQUIVALENCE.json", result)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", type=int, default=40)
    args = parser.parse_args(argv)
    result = run(sample=args.sample)
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
