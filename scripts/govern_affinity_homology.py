"""Label-free 40% homology governance for EnergyPilot source targets."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.source_affinity.common import sha256_file, write_canonical_json


THRESHOLD = 0.40
CANDIDATE_THRESHOLD = 0.30


def read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def read_fasta(path: Path) -> dict[str, str]:
    records, identifier, chunks = {}, None, []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if identifier is not None:
                records[identifier] = "".join(chunks)
            identifier, chunks = line[1:].split()[0], []
        elif identifier is not None:
            chunks.append(line.strip())
    if identifier is not None:
        records[identifier] = "".join(chunks)
    return records


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in sorted(sequences.items()):
            handle.write(f">{identifier}\n{sequence}\n")


def local_identity(pair: tuple[str, str]) -> tuple[tuple[str, str], float]:
    import parasail

    left, right = pair
    result = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
    return pair, result.matches / result.length if result.length else 0.0


class Components:
    def __init__(self, nodes):
        self.parent = {node: node for node in nodes}

    def find(self, node):
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)


def candidate_pairs(
    mmseqs: Path, query: Path, target: Path, output: Path, temporary: Path, threads: int,
) -> set[tuple[str, str]]:
    command = [
        str(mmseqs), "easy-search", str(query), str(target), str(output), str(temporary),
        "--min-seq-id", str(CANDIDATE_THRESHOLD), "-c", "0.0", "--max-seqs", "100000",
        "-s", "7.5", "--format-output", "query,target", "--threads", str(threads),
    ]
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(f"MMseqs2 failed ({process.returncode}): {process.stderr[-2000:]}")
    pairs = set()
    if output.is_file():
        for line in output.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) >= 2 and fields[0] != fields[1]:
                pairs.add((fields[0], fields[1]))
    return pairs


def govern(
    rows_path: Path, governance_path: Path, protected_fasta: Path,
    output: Path, mmseqs: Path, threads: int = 8,
) -> dict:
    if output.exists():
        raise FileExistsError(f"governance output already exists: {output}")
    output.mkdir(parents=True)
    work = output / "work"
    work.mkdir()
    sequences = {}
    for row in read_jsonl(rows_path):
        key, sequence = row["protein_sequence_sha256"], row["protein_sequence"]
        if hashlib.sha256(sequence.encode("ascii")).hexdigest() != key:
            raise ValueError(f"source sequence hash mismatch: {key}")
        if sequences.setdefault(key, sequence) != sequence:
            raise ValueError(f"conflicting source sequence: {key}")

    governance = list(read_jsonl(governance_path))
    protected_hashes = {
        row["target_key"] for row in governance if row.get("split") in {"metaval", "recipient"}
    }
    fasta = read_fasta(protected_fasta)
    protected = {}
    for key in protected_hashes:
        identifier = "s_" + key[:24]
        sequence = fasta.get(identifier)
        if sequence is None or hashlib.sha256(sequence.encode("ascii")).hexdigest() != key:
            raise ValueError(f"protected sequence unavailable or invalid: {key}")
        protected[key] = sequence

    source_ids = {"s_" + key[:24]: value for key, value in sequences.items()}
    source_hash = {"s_" + key[:24]: key for key in sequences}
    protected_ids = {"p_" + key[:24]: value for key, value in protected.items()}
    protected_hash = {"p_" + key[:24]: key for key in protected}
    source_fasta, protected_path = work / "source.fasta", work / "protected.fasta"
    write_fasta(source_fasta, source_ids)
    write_fasta(protected_path, protected_ids)

    executable = mmseqs.resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"MMseqs2 executable not found: {executable}")
    source_candidates = candidate_pairs(
        executable, source_fasta, source_fasta, work / "source_candidates.tsv",
        work / "source_tmp", threads,
    )
    protected_candidates = candidate_pairs(
        executable, source_fasta, protected_path, work / "protected_candidates.tsv",
        work / "protected_tmp", threads,
    )

    components = Components(sequences)
    source_jobs = [
        (sequences[source_hash[left]], sequences[source_hash[right]])
        for left, right in sorted(source_candidates) if left < right
    ]
    source_keys = [
        (source_hash[left], source_hash[right])
        for left, right in sorted(source_candidates) if left < right
    ]
    with ThreadPoolExecutor(max_workers=threads) as pool:
        source_results = list(pool.map(local_identity, source_jobs))
    confirmed_source = 0
    for keys, (_, identity) in zip(source_keys, source_results):
        if identity >= THRESHOLD:
            components.union(*keys)
            confirmed_source += 1

    protected_jobs = [
        (sequences[source_hash[left]], protected[protected_hash[right]])
        for left, right in sorted(protected_candidates)
    ]
    protected_keys = [
        (source_hash[left], protected_hash[right]) for left, right in sorted(protected_candidates)
    ]
    excluded = {}
    with ThreadPoolExecutor(max_workers=threads) as pool:
        protected_results = list(pool.map(local_identity, protected_jobs))
    for keys, (_, identity) in zip(protected_keys, protected_results):
        if identity >= THRESHOLD and identity > excluded.get(keys[0], {}).get("max_identity", 0):
            excluded[keys[0]] = {"max_identity": identity, "protected_sequence_sha256": keys[1]}

    assignments_path = output / "homology_assignments.jsonl"
    with assignments_path.open("w", encoding="utf-8", newline="\n") as handle:
        for key in sorted(sequences):
            row = {
                "protein_sequence_sha256": key,
                "homology_component_id": components.find(key),
                "excluded_by_davis_protected_homology": key in excluded,
                "protected_match": excluded.get(key),
            }
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    version = subprocess.run(
        [str(executable), "version"], capture_output=True, text=True, check=True,
    ).stdout.strip()
    manifest = {
        "schema": "MetaSieve.AffinityHomologyGovernance.v1",
        "stage": "P1R2B-D1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "identity_threshold": THRESHOLD,
        "candidate_identity_threshold": CANDIDATE_THRESHOLD,
        "identity_definition": "parasail_sw_blosum62_gap_open10_extend1_matches_per_local_alignment_length",
        "source_sequences": len(sequences),
        "protected_sequences": len(protected),
        "source_candidate_pairs": len(source_candidates),
        "source_confirmed_edges": confirmed_source,
        "homology_components": len({components.find(key) for key in sequences}),
        "protected_candidate_pairs": len(protected_candidates),
        "excluded_source_sequences": len(excluded),
        "assignments_sha256": sha256_file(assignments_path),
        "inputs": {
            "canonical_rows_sha256": sha256_file(rows_path),
            "davis_governance_sha256": sha256_file(governance_path),
            "protected_fasta_sha256": sha256_file(protected_fasta),
        },
        "mmseqs": {"path": str(executable), "version": version, "sha256": sha256_file(executable)},
        "recipient_labels_read": False,
        "training_authorized": False,
    }
    write_canonical_json(output / "governance_manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument(
        "--davis-governance", type=Path,
        default=ROOT / "dataset/sealed/DAVIS_mechanism_v2/governance.jsonl",
    )
    parser.add_argument(
        "--protected-fasta", type=Path,
        default=ROOT / "dataset/processed/open_structures/pilot20k_governance_v2/work/benchmark.fasta",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mmseqs", type=Path,
                        default=ROOT / "tools/runtime/mmseqs2/mmseqs/bin/mmseqs.exe")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    result = govern(
        args.rows, args.davis_governance, args.protected_fasta,
        args.output, args.mmseqs, args.threads,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
