"""Create deterministic protein-homology splits for P1B structural learning."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from scripts.data_contract import read_jsonl, write_jsonl
from scripts.govern_structure_homology import IDENTITY_THRESHOLD, _local_identity, _sequence_id
from scripts.structure_sources.rcsb import sha256_file


SPLIT_FRACTIONS = {"train": 0.8, "val": 0.1, "test": 0.1}


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in sorted(sequences.items()):
            handle.write(f">{identifier}\n{sequence}\n")


def split_structure_corpus(complexes_path: str | Path, output_dir: str | Path, *,
                           mmseqs: str | Path) -> dict:
    records = read_jsonl(complexes_path)
    if not records:
        raise ValueError("structural split requires canonical holo records")
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"structural split output already exists: {output}")
    output.mkdir(parents=True)
    work = output / "work"
    work.mkdir()
    sequences = {_sequence_id(record["sequence"]): record["sequence"] for record in records}
    fasta, prefix, temporary = work / "sequences.fasta", work / "cluster", work / "tmp"
    _write_fasta(fasta, sequences)
    executable = Path(mmseqs).resolve()
    command = [str(executable), "easy-cluster", str(fasta), str(prefix), str(temporary),
               "--min-seq-id", str(IDENTITY_THRESHOLD), "-c", "0.8", "--cov-mode", "0",
               "--threads", "4"]
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(f"MMseqs2 clustering failed ({process.returncode}): {process.stderr[-2000:]}")
    cluster_path = prefix.with_name(prefix.name + "_cluster.tsv")
    if not cluster_path.is_file():
        raise FileNotFoundError(f"MMseqs2 did not emit cluster TSV: {cluster_path}")
    groups: dict[str, set[str]] = {}
    confirmations, rejected_edges = 0, 0
    with cluster_path.open(encoding="utf-8") as handle:
        for line in handle:
            representative, member = line.rstrip("\n").split("\t")[:2]
            if representative != member:
                confirmations += 1
                if _local_identity(sequences[representative], sequences[member]) < IDENTITY_THRESHOLD:
                    rejected_edges += 1
                    groups.setdefault(member, set()).add(member)
                    continue
            groups.setdefault(representative, set()).add(member)
    assigned_members = set().union(*groups.values()) if groups else set()
    for identifier in set(sequences) - assigned_members:
        groups[identifier] = {identifier}

    target = {name: len(sequences) * fraction for name, fraction in SPLIT_FRACTIONS.items()}
    counts = Counter()
    group_split: dict[str, str] = {}
    for representative, members in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
        split = min(SPLIT_FRACTIONS, key=lambda name: (
            (counts[name] + len(members)) / target[name], name))
        group_split[representative] = split
        counts[split] += len(members)
    member_assignment = {member: (representative, group_split[representative])
                         for representative, members in groups.items() for member in members}
    output_records = []
    for record in records:
        sequence_id = _sequence_id(record["sequence"])
        group_id, split = member_assignment[sequence_id]
        output_records.append(record | {"homology_group_id": group_id,
                                        "source_split": split})
    write_jsonl(output / "complexes.jsonl", output_records)
    output_counts = Counter(record["source_split"] for record in output_records)
    version = subprocess.run([str(executable), "version"], capture_output=True,
                             text=True, check=True).stdout.strip()
    manifest = {
        "schema": "MetaSieve.StructureHomologySplit.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "identity_threshold": IDENTITY_THRESHOLD,
        "coverage_threshold": 0.8,
        "coverage_mode": 0,
        "confirmation": "parasail_sw_blosum62_gap_open10_extend1_matches_per_local_alignment_length",
        "records": len(records), "unique_sequences": len(sequences),
        "homology_groups": len(groups), "confirmed_representative_edges": confirmations,
        "rejected_representative_edges": rejected_edges,
        "split_fractions": SPLIT_FRACTIONS, "split_records": dict(output_counts),
        "split_groups": dict(Counter(group_split.values())),
        "input": {"path": str(Path(complexes_path).resolve()),
                  "sha256": sha256_file(complexes_path)},
        "output_sha256": sha256_file(output / "complexes.jsonl"),
        "mmseqs": {"path": str(executable), "version": version,
                   "sha256": sha256_file(executable), "command": command},
        "gate_status": "PASS" if all(output_counts[name] for name in SPLIT_FRACTIONS) else "FAIL",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("complexes")
    parser.add_argument("output")
    parser.add_argument("--mmseqs", default="tools/mmseqs2/mmseqs/bin/mmseqs.exe")
    args = parser.parse_args()
    result = split_structure_corpus(args.complexes, args.output, mmseqs=args.mmseqs)
    print(json.dumps(result, indent=2))
    return 0 if result["gate_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
