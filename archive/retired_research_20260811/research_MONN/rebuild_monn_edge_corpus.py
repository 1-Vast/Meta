#!/usr/bin/env python3
"""Rebuild a label-blind residue-atom edge corpus from fixed MONN pickles.

Inputs
------
An authorized non-commercial checkout of https://github.com/lishuya17/MONN at
commit f2b62ccf49c18a9502aa0eb0d582c6e0735ef200.

Outputs
-------
Two deterministic gzip JSONL files containing only complexes with at least one
mapped positive edge, plus a summary JSON with source and output hashes.

This program deliberately never opens either MONN TSV file, so affinity values
cannot enter the localization-only R0 corpus. Python pickle is unsafe for
untrusted inputs; use only the hash-verified official source files.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import pickle
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


COMMIT = "f2b62ccf49c18a9502aa0eb0d582c6e0735ef200"
EXPECTED = {
    "development": {
        "path": "data/out7_final_pairwise_interaction_dict",
        "sha256": "9e7d1128a79139cb3a43d077ba5d19cce6376ddc9cf35d65db925e7f5e7e9d82",
        "raw_records": 12987,
        "mapped_records": 12738,
        "binary_edges": 195798,
        "typed_edges": 202766,
    },
    "additional_pdb": {
        "path": "data/independent_dataset_interaction_dict",
        "sha256": "377b83080190e56a5ceea09101b73234596fbf069acc4be872806c92d4d68598",
        "raw_records": 1853,
        "mapped_records": 1851,
        "binary_edges": 9832,
        "typed_edges": 9832,
    },
}
EVENT_RE = re.compile(r"^(?P<kind>.+)_(?P<number>[0-9]+)$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_trusted_pickle(path: Path) -> dict[str, dict[str, Any]]:
    with path.open("rb") as stream:
        try:
            value = pickle.load(stream)
        except UnicodeDecodeError:
            stream.seek(0)
            value = pickle.load(stream, encoding="latin1")
    if not isinstance(value, dict):
        raise TypeError(f"{path} did not contain a dictionary")
    return value


def deterministic_gzip_text(path: Path) -> io.TextIOWrapper:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")


def sorted_unique(items: Iterable[tuple[Any, ...]]) -> list[tuple[Any, ...]]:
    return sorted(set(items))


def rebuild_record(source_key: str, value: dict[str, Any]) -> tuple[dict[str, Any], int]:
    required = {
        "atom_bond_type",
        "atom_idx",
        "atom_name",
        "ligand",
        "residue_bond_type",
        "uniprot_id",
        "uniprot_seq",
    }
    missing = required.difference(value)
    if missing:
        raise KeyError(f"{source_key}: missing fields {sorted(missing)}")

    atom_names = list(value["atom_name"])
    atom_indices = [int(item) for item in value["atom_idx"]]
    if len(atom_names) != len(atom_indices):
        raise ValueError(f"{source_key}: atom_name/atom_idx length mismatch")
    if len(set(atom_names)) != len(atom_names):
        raise ValueError(f"{source_key}: ligand atom names are not unique")
    atom_slot = {name: index for index, name in enumerate(atom_names)}

    atoms_by_event: dict[str, list[str]] = defaultdict(list)
    residues_by_event: dict[str, list[int]] = defaultdict(list)
    for atom_name, event_id in value["atom_bond_type"]:
        atoms_by_event[str(event_id)].append(str(atom_name))
    for residue_index, event_id in value["residue_bond_type"]:
        residues_by_event[str(event_id)].append(int(residue_index))

    missing_atom_references = 0
    event_edges: list[tuple[int, int, str, str]] = []
    for event_id in sorted(atoms_by_event.keys() & residues_by_event.keys()):
        match = EVENT_RE.match(event_id)
        if match is None:
            raise ValueError(f"{source_key}: invalid full event id {event_id!r}")
        interaction_type = match.group("kind")
        for atom_name in atoms_by_event[event_id]:
            if atom_name not in atom_slot:
                missing_atom_references += 1
                continue
            slot = atom_slot[atom_name]
            for residue_index in residues_by_event[event_id]:
                event_edges.append((residue_index, slot, interaction_type, event_id))

    event_edges = sorted_unique(event_edges)
    binary_edges = sorted_unique((residue, atom) for residue, atom, _, _ in event_edges)
    typed_edges = sorted_unique((residue, atom, kind) for residue, atom, kind, _ in event_edges)

    record = {
        "source_key": source_key,
        "pdb_id": source_key.split("_", 1)[0].lower(),
        "ligand_ccd": str(value["ligand"]),
        "uniprot_id": str(value["uniprot_id"]),
        "uniprot_sequence": str(value["uniprot_seq"]),
        "atom_names": atom_names,
        "source_atom_indices": atom_indices,
        "positive_binary_edges": [list(item) for item in binary_edges],
        "positive_typed_edges": [list(item) for item in typed_edges],
        "positive_event_edges": [list(item) for item in event_edges],
    }
    return record, missing_atom_references


def rebuild_split(
    split: str,
    source_path: Path,
    output_path: Path,
) -> dict[str, Any]:
    source = load_trusted_pickle(source_path)
    mapped_records = 0
    binary_edges = 0
    typed_edges = 0
    event_edges = 0
    missing_atom_references = 0

    with deterministic_gzip_text(output_path) as stream:
        for source_key in sorted(source):
            record, missing = rebuild_record(str(source_key), source[source_key])
            missing_atom_references += missing
            if not record["positive_binary_edges"]:
                continue
            mapped_records += 1
            binary_edges += len(record["positive_binary_edges"])
            typed_edges += len(record["positive_typed_edges"])
            event_edges += len(record["positive_event_edges"])
            stream.write(json.dumps(record, sort_keys=True, separators=(",", ":")))
            stream.write("\n")

    return {
        "split": split,
        "source_path": str(source_path),
        "source_sha256": sha256(source_path),
        "output_path": str(output_path),
        "output_sha256": sha256(output_path),
        "raw_records": len(source),
        "mapped_records": mapped_records,
        "binary_edges": binary_edges,
        "typed_edges": typed_edges,
        "event_edges": event_edges,
        "missing_atom_references": missing_atom_references,
        "affinity_tables_opened": 0,
    }


def check_expected(split: str, result: dict[str, Any], strict_hashes: bool) -> None:
    expected = EXPECTED[split]
    if strict_hashes and result["source_sha256"] != expected["sha256"]:
        raise AssertionError(
            f"{split}: source SHA-256 mismatch: {result['source_sha256']} != {expected['sha256']}"
        )
    for key in ("raw_records", "mapped_records", "binary_edges", "typed_edges"):
        if result[key] != expected[key]:
            raise AssertionError(
                f"{split}: {key} mismatch: {result[key]} != {expected[key]}"
            )
    if result["missing_atom_references"] != 0:
        raise AssertionError(f"{split}: missing atom references were observed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--monn-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--strict-hashes", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.monn_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for split, spec in EXPECTED.items():
        source_path = root / spec["path"]
        if not source_path.is_file():
            raise FileNotFoundError(source_path)
        output_path = output_dir / f"monn_{split}_edge_corpus.jsonl.gz"
        result = rebuild_split(split, source_path, output_path)
        check_expected(split, result, args.strict_hashes)
        results.append(result)

    summary = {
        "schema_version": "1.0",
        "source_repository": "https://github.com/lishuya17/MONN",
        "required_commit": COMMIT,
        "edge_join": "complete PLIP event identifier",
        "affinity_tables_opened": 0,
        "results": results,
    }
    summary_path = output_dir / "raw_corpus_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

