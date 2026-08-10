"""C0 exposure registry, untouched corpus construction and deconvolution census.

Registered by
  research/correspondence_router/PREREG_C0_C1_CORRESPONDENCE_INFORMATION_AUDIT.md
  (sha256 007f8439..., commit f844679) committed BEFORE this file existed.

Audit only. Trains nothing. Opens no affinity field from any source. BioLiP2 is
used solely as a biological-relevance filter; PLINDER is not used at all.
Heldout-A is permanently consumed and is never referenced.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "correspondence_router"
sys.path.insert(0, str(HERE))

PREREG = HERE / "PREREG_C0_C1_CORRESPONDENCE_INFORMATION_AUDIT.md"
PREREG_SHA = "007f8439609078649cf7751b588716492f59c93bc27e2dad997b11afd7172c1e"

OUT = ROOT / "report" / "correspondence_router"
EXEC = ROOT / "dataset" / "processed" / "correspondence_router"
RAW = ROOT / "dataset" / "raw" / "open_structures" / "pilot20k"
MMCIF = RAW / "mmcif"
BIOLIP = RAW / "biolip2" / "BioLiP.txt.gz"
OPEN_PROC = ROOT / "dataset" / "processed" / "open_structures"
MONN_CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
SSL_ROWS = ROOT / "dataset" / "processed" / "ssl_b2" / "parsed_rows_1476.json"

# ---- frozen contract (prereg sections 2, 4, 5) -----------------------------
SLOTS = 128
CONTACT_THRESHOLD = 6.0
SENSITIVITY_THRESHOLD = 4.5
MIN_LIGAND_ATOMS, MAX_LIGAND_ATOMS = 6, 80
MIN_SEQUENCE, MAX_SEQUENCE = 150, 1200
MAX_RESOLUTION = 2.5
MIN_MAPPING_COVERAGE = 0.90

ADDITIVE_BLACKLIST = frozenset("""HOH GOL EDO PEG PG4 PGE SO4 PO4 CL NA K MG CA
ZN MN FE CU NI CD HG ACT ACY DMS MPD TRS EPE IMD BME NAG BMA MAN FUC GAL GLC
XYL IOD BR NO3 CO3 FMT OXL TLA CIT MES CAC SCN AZI UNL UNX""".split())

ATOM_TAGS = ("group_PDB", "label_asym_id", "label_seq_id", "label_comp_id",
             "label_atom_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z",
             "label_alt_id", "occupancy", "pdbx_PDB_model_num")


class C0ContractError(RuntimeError):
    pass


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str),
                    encoding="utf-8")


def write_jsonl(path: Path, rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def read_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha_json(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


# --------------------------------------------------------------- exposure
def exposure_registry() -> dict:
    """Every PDB id consumed by any prior stage, in any role (prereg s.3)."""
    sources: dict[str, set[str]] = {}
    for directory in sorted(OPEN_PROC.glob("*")):
        path = directory / "complexes.jsonl"
        if path.is_file():
            sources[f"open_structures/{directory.name}"] = {
                row["pdb_id"].lower() for row in read_jsonl(path) if row.get("pdb_id")}
    monn = set()
    for path in sorted(MONN_CORPUS.glob("*.jsonl.gz")):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                if row.get("pdb_id"):
                    monn.add(row["pdb_id"].lower())
    sources["monn_s7_b5_s3r_s4r_s5d"] = monn
    if SSL_ROWS.is_file():
        payload = json.loads(SSL_ROWS.read_text(encoding="utf-8"))
        rows = payload if isinstance(payload, list) else payload.get("rows", [])
        ids = set()
        for row in rows:
            if not isinstance(row, dict):
                continue
            for key in ("pdb_id", "pdb", "entry_id", "system_id"):
                if row.get(key):
                    ids.add(str(row[key]).lower()[:4])
                    break
        sources["ssl_b2_independent_structural"] = ids
    union = set().union(*sources.values()) if sources else set()
    return {"per_source_counts": {k: len(v) for k, v in sorted(sources.items())},
            "union_count": len(union), "union": union}


def biolip_ligands(candidate_ids: set[str]) -> dict[str, set[str]]:
    """PDB id -> set of biologically relevant ligand comp ids. Annotation only:
    no numeric BioLiP2 field, and no affinity column, is read."""
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(BIOLIP, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            pdb_id = parts[0].strip().lower()
            if pdb_id not in candidate_ids:
                continue
            comp = parts[4].strip().upper()
            if comp and comp not in ADDITIVE_BLACKLIST:
                out[pdb_id].add(comp)
    return dict(out)


# --------------------------------------------------------------- parsing
def _canonical_rows(rows, key_fields):
    groups: dict[tuple, list] = {}
    for row in rows:
        if row["pdbx_PDB_model_num"] not in {"", ".", "?", "1"}:
            continue
        groups.setdefault(tuple(row[f] for f in key_fields), []).append(row)
    selected = []
    for values in groups.values():
        values.sort(key=lambda row: (
            row["label_alt_id"] not in {"", ".", "?"},
            -float(row["occupancy"] if row["occupancy"] not in {"", ".", "?"} else 0),
            row["label_alt_id"]))
        selected.append(values[0])
    return selected


def parse_entry(pdb_id: str, allowed_comps: set[str]) -> tuple[list, list]:
    """Return (systems, exclusions) for one raw mmCIF entry. Deterministic."""
    import gemmi
    path = MMCIF / f"{pdb_id}.cif.gz"
    exclusions = []
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            block = gemmi.cif.read_string(handle.read()).sole_block()
    except Exception as exc:
        return [], [{"pdb_id": pdb_id, "reason": "unparsable_mmcif",
                     "detail": type(exc).__name__}]

    method = " ".join(block.find_values("_exptl.method")).upper()
    resolution = None
    for tag in ("_refine.ls_d_res_high", "_reflns.d_resolution_high"):
        for value in block.find_values(tag):
            try:
                resolution = float(gemmi.cif.as_string(value))
                break
            except (ValueError, TypeError):
                continue
        if resolution is not None:
            break
    if "X-RAY" not in method:
        return [], [{"pdb_id": pdb_id, "reason": "not_xray", "detail": method[:40]}]
    if resolution is None or resolution > MAX_RESOLUTION:
        return [], [{"pdb_id": pdb_id, "reason": "resolution_filter",
                     "detail": resolution}]

    deposited = ""
    for value in block.find_values(
            "_pdbx_database_status.recvd_initial_deposition_date"):
        deposited = gemmi.cif.as_string(value)
        break

    entity_sequence = {}
    for eid, seq in zip(block.find_values("_entity_poly.entity_id"),
                        block.find_values("_entity_poly.pdbx_seq_one_letter_code_can")):
        entity_sequence[gemmi.cif.as_string(eid).strip()] = "".join(
            gemmi.cif.as_string(seq).split())
    asym_entity = {gemmi.cif.as_string(a).strip(): gemmi.cif.as_string(e).strip()
                   for a, e in zip(block.find_values("_struct_asym.id"),
                                   block.find_values("_struct_asym.entity_id"))}

    covalent = set()
    tags = ["_struct_conn.conn_type_id", "_struct_conn.ptnr1_label_asym_id",
            "_struct_conn.ptnr2_label_asym_id"]
    try:
        for row in block.find(tags):
            if str(row[0]).lower().startswith("coval"):
                covalent.add(str(row[1]))
                covalent.add(str(row[2]))
    except Exception:
        pass

    table = block.find("_atom_site.", list(ATOM_TAGS))
    atoms = [{tag: str(row[i]) for i, tag in enumerate(ATOM_TAGS)} for row in table]
    if not atoms:
        return [], [{"pdb_id": pdb_id, "reason": "no_atom_site"}]

    protein_atoms = [a for a in atoms if a["group_PDB"] == "ATOM"
                     and a["label_seq_id"] not in {"", ".", "?"}
                     and a["type_symbol"].upper() != "H"]
    hetatms = [a for a in atoms if a["group_PDB"] == "HETATM"
               and a["type_symbol"].upper() != "H"
               and a["label_comp_id"].upper() in allowed_comps]
    if not protein_atoms or not hetatms:
        return [], [{"pdb_id": pdb_id, "reason": "no_protein_or_relevant_ligand"}]

    protein_by_asym = defaultdict(list)
    for atom in _canonical_rows(protein_atoms,
                                ("label_asym_id", "label_seq_id", "label_atom_id")):
        protein_by_asym[atom["label_asym_id"]].append(atom)
    ligand_by_key = defaultdict(list)
    for atom in _canonical_rows(hetatms, ("label_asym_id", "label_seq_id",
                                          "label_atom_id")):
        ligand_by_key[(atom["label_asym_id"], atom["label_seq_id"],
                       atom["label_comp_id"].upper())].append(atom)

    systems = []
    for (lig_asym, lig_seq, comp), lig_atoms in sorted(ligand_by_key.items()):
        if not (MIN_LIGAND_ATOMS <= len(lig_atoms) <= MAX_LIGAND_ATOMS):
            exclusions.append({"pdb_id": pdb_id, "ligand": comp,
                               "reason": "ligand_atom_count",
                               "detail": len(lig_atoms)})
            continue
        if lig_asym in covalent:
            exclusions.append({"pdb_id": pdb_id, "ligand": comp,
                               "reason": "covalent_link"})
            continue
        lig_xyz = np.asarray([[float(a["Cartn_x"]), float(a["Cartn_y"]),
                               float(a["Cartn_z"])] for a in lig_atoms],
                             dtype=np.float64)
        for prot_asym in sorted(protein_by_asym):
            sequence = entity_sequence.get(asym_entity.get(prot_asym, ""), "")
            length = len(sequence)
            if not (MIN_SEQUENCE <= length <= MAX_SEQUENCE):
                continue
            residues = defaultdict(list)
            for atom in protein_by_asym[prot_asym]:
                try:
                    index = int(atom["label_seq_id"]) - 1
                except ValueError:
                    continue
                if 0 <= index < length:
                    residues[index].append([float(atom["Cartn_x"]),
                                            float(atom["Cartn_y"]),
                                            float(atom["Cartn_z"])])
            if not residues:
                continue
            coverage = len(residues) / float(length)
            if coverage < MIN_MAPPING_COVERAGE:
                exclusions.append({"pdb_id": pdb_id, "ligand": comp,
                                   "protein_asym": prot_asym,
                                   "reason": "mapping_coverage", "detail": coverage})
                continue
            indices = sorted(residues)
            distances = np.full((len(lig_xyz), len(indices)), np.inf)
            for column, index in enumerate(indices):
                block_xyz = np.asarray(residues[index], dtype=np.float64)
                distances[:, column] = np.sqrt(
                    ((lig_xyz[:, None, :] - block_xyz[None, :, :]) ** 2).sum(-1)).min(-1)
            if not (distances <= CONTACT_THRESHOLD).any():
                continue
            systems.append({
                "pdb_id": pdb_id, "protein_asym_id": prot_asym,
                "ligand_asym_id": lig_asym, "ligand_label_seq_id": lig_seq,
                "ligand_comp_id": comp, "sequence": sequence,
                "sequence_length": length, "resolution_angstrom": resolution,
                "deposition_date": deposited,
                "ligand_heavy_atoms": int(len(lig_xyz)),
                "resolved_residues": len(indices),
                "mapping_coverage": coverage,
                "residue_indices": indices,
                "distances": distances.astype(np.float32),
            })
    return systems, exclusions


def _worker(payload):
    pdb_id, comps = payload
    try:
        return parse_entry(pdb_id, comps)
    except Exception as exc:                    # fail-closed, never silent
        return [], [{"pdb_id": pdb_id, "reason": "worker_exception",
                     "detail": f"{type(exc).__name__}: {exc}"[:200]}]


# --------------------------------------------------------------- census
def census_system(system) -> dict:
    """Deconvolution-unit census for one system (prereg section 5)."""
    distances = system["distances"]
    indices = np.asarray(system["residue_indices"])
    length = system["sequence_length"]
    slot_of = np.minimum(SLOTS - 1, indices * SLOTS // length)
    contact = distances <= CONTACT_THRESHOLD

    per_slot = Counter(slot_of.tolist())
    candidate_slots = sorted(s for s, n in per_slot.items() if n >= 2)
    units = positives = multi = 0
    checkerboards = 0
    for slot in candidate_slots:
        columns = np.flatnonzero(slot_of == slot)
        sub = contact[:, columns]
        counts = sub.sum(1)
        units += sub.shape[0]
        positives += int((counts >= 1).sum())
        multi += int((counts >= 2).sum())
        # valid 2x2: rows i!=j, columns r!=r' with a diagonal contact pattern
        rows = np.flatnonzero((counts >= 1) & (counts < sub.shape[1]))
        for a_pos in range(len(rows)):
            for b_pos in range(a_pos + 1, len(rows)):
                left, right = sub[rows[a_pos]], sub[rows[b_pos]]
                only_left = np.flatnonzero(left & ~right)
                only_right = np.flatnonzero(right & ~left)
                checkerboards += int(len(only_left) * len(only_right))
    return {
        "candidate_slots": len(candidate_slots),
        "slots_used": len(per_slot),
        "max_residues_per_slot": max(per_slot.values()) if per_slot else 0,
        "units": units, "positive_units": positives, "multi_contact_units": multi,
        "valid_checkerboards": checkerboards,
        "exact_contacts": int(contact.sum()),
        "contacted_residues": int(contact.any(0).sum()),
        "contact_density": float(contact.mean()),
    }


def run(limit: int = 0, workers: int = 8) -> dict:
    if sha_file(PREREG) != PREREG_SHA:
        raise C0ContractError("C0/C1 preregistration hash mismatch")
    started = time.time()
    registry = exposure_registry()
    local = sorted({p.name.split(".")[0].lower() for p in MMCIF.glob("*.cif.gz")})
    untouched = sorted(set(local) - registry["union"])
    relevant = biolip_ligands(set(untouched))
    candidates = sorted(pid for pid in untouched if relevant.get(pid))
    if limit:
        candidates = candidates[:limit]
    print(f"exposed={registry['union_count']} local={len(local)} "
          f"untouched={len(untouched)} biolip_relevant={len(candidates)}", flush=True)

    systems, exclusions = [], []
    payloads = [(pid, relevant[pid]) for pid in candidates]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for done, (built, dropped) in enumerate(
                pool.map(_worker, payloads, chunksize=8), start=1):
            systems.extend(built)
            exclusions.extend(dropped)
            if done % 250 == 0:
                print(f"  parsed {done}/{len(payloads)} entries, "
                      f"{len(systems)} systems", flush=True)

    if not systems:
        raise C0ContractError("no admissible systems in the untouched corpus")

    for system in systems:
        system["census"] = census_system(system)
    systems.sort(key=lambda s: (s["pdb_id"], s["protein_asym_id"],
                                s["ligand_asym_id"], s["ligand_label_seq_id"]))
    for index, system in enumerate(systems):
        system["system_id"] = (f"{system['pdb_id']}_{system['protein_asym_id']}"
                               f"_{system['ligand_asym_id']}"
                               f"_{system['ligand_label_seq_id']}")
        system["row"] = index

    EXEC.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        EXEC / "c0_geometry.npz",
        **{f"d_{s['row']}": s["distances"] for s in systems},
        **{f"i_{s['row']}": np.asarray(s["residue_indices"], dtype=np.int32)
           for s in systems})
    metadata = [{k: v for k, v in s.items()
                 if k not in {"distances", "residue_indices"}} for s in systems]
    write_jsonl(EXEC / "c0_systems.jsonl", metadata)
    write_jsonl(EXEC / "c0_exclusions.jsonl", exclusions)

    totals = Counter()
    for system in systems:
        for key, value in system["census"].items():
            if isinstance(value, (int, float)) and key != "contact_density":
                totals[key] += value
    positive_units = totals["positive_units"]
    result = {
        "schema": "MetaSieve.Correspondence.C0.Census.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "exposure_registry": {
            "per_source_counts": registry["per_source_counts"],
            "union_exposed_pdb_ids": registry["union_count"],
            "union_sha256": sha_json(sorted(registry["union"])),
        },
        "corpus": {
            "local_mmcif_entries": len(local),
            "untouched_entries": len(untouched),
            "biolip_relevant_untouched_entries": len(candidates),
            "parsed_entries": len(payloads),
            "admissible_systems": len(systems),
            "distinct_pdb_entries": len({s["pdb_id"] for s in systems}),
            "distinct_sequences": len({s["sequence"] for s in systems}),
            "distinct_ligand_comp_ids": len({s["ligand_comp_id"] for s in systems}),
            "exclusions": len(exclusions),
            "exclusion_reasons": dict(Counter(e["reason"] for e in exclusions)),
        },
        "deconvolution_census": {
            "units": totals["units"],
            "positive_units": positive_units,
            "multi_contact_units": totals["multi_contact_units"],
            "multi_contact_rate": (totals["multi_contact_units"] /
                                   max(positive_units, 1)),
            "valid_checkerboards": totals["valid_checkerboards"],
            "candidate_slots": totals["candidate_slots"],
            "exact_contacts": totals["exact_contacts"],
            "median_residues_per_candidate_slot": float(np.median(
                [s["census"]["max_residues_per_slot"] for s in systems])),
            "systems_with_zero_candidate_slots": sum(
                1 for s in systems if s["census"]["candidate_slots"] == 0),
        },
        "affinity_value_reads": 0,
        "plinder_used": False,
        "biolip2_role": "ANNOTATION_ONLY_biological_relevance_filter",
        "heldoutA_referenced": False,
        "trainable_parameters_introduced": 0,
        "elapsed_seconds": round(time.time() - started, 1),
    }
    write_json(OUT / "C0_CORPUS_AND_CENSUS.json", result)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args(argv)
    try:
        result = run(limit=args.limit, workers=args.workers)
        print(json.dumps({k: v for k, v in result.items()
                          if k not in {"exposure_registry"}}, indent=2, default=str),
              flush=True)
        return 0
    except Exception as exc:
        failure = {"schema": "MetaSieve.Correspondence.C0.FailClosed.v1",
                   "error_type": type(exc).__name__, "error": str(exc),
                   "TERMINAL_VERDICT": "CORRESPONDENCE_DATA_OR_CLOSURE_NOT_IDENTIFIABLE",
                   "affinity_value_reads": 0}
        write_json(OUT / "C0_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
