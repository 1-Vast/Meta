"""Close R0-C protein dependencies before any coordinate download."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from research.correspondence_router.prepare_r0c_candidates import (
    P1B_CANDIDATE_LIMIT,
)
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.biolip import pilot_candidates, regular_ligand_ids
from scripts.structure_sources.rcsb import sha256_file


IDENTITY_THRESHOLD = 0.40
COVERAGE_THRESHOLD = 0.80
MMSEQS_CANDIDATE_IDENTITY = 0.30
DOWNLOAD_LIMIT = 512
MIN_SCORABLE_COMPONENTS = 120
REPRESENTATIVE_NAMESPACE = "R0C-DOWNLOAD-REP-v1"


def _sequence_id(sequence: str) -> str:
    return "s_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()[:24]


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in sorted(sequences.items()):
            handle.write(f">{identifier}\n{sequence}\n")


def _homologous(left: str, right: str) -> bool:
    import parasail

    result = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
    if not result.length:
        return False
    identity = result.matches / result.length
    left_coverage = result.length / len(left)
    right_coverage = result.length / len(right)
    return (
        identity >= IDENTITY_THRESHOLD
        and left_coverage >= COVERAGE_THRESHOLD
        and right_coverage >= COVERAGE_THRESHOLD
    )


class _UnionFind:
    def __init__(self, values: set[str]):
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            value, self.parent[value] = self.parent[value], root
        return root

    def union(self, left: str, right: str) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root == right_root:
            return
        low, high = sorted((left_root, right_root))
        self.parent[high] = low


def select_fresh_components(
    records: list[dict],
    exposure_hits: set[str],
    internal_edges: set[tuple[str, str]],
    *,
    limit: int = DOWNLOAD_LIMIT,
) -> tuple[list[dict], dict]:
    if limit < MIN_SCORABLE_COMPONENTS:
        raise ValueError("download limit is below the frozen scorable-component target")
    by_sequence = {_sequence_id(str(row["sequence"])): row for row in records}
    if len(by_sequence) != len(records):
        raise ValueError("R0-C candidate sequences must be unique")
    union = _UnionFind(set(by_sequence))
    for left, right in sorted(internal_edges):
        if left not in by_sequence or right not in by_sequence:
            raise ValueError("internal homology edge references an unknown candidate")
        union.union(left, right)
    components: dict[str, list[str]] = {}
    for identifier in by_sequence:
        components.setdefault(union.find(identifier), []).append(identifier)
    contaminated_roots = {union.find(identifier) for identifier in exposure_hits}
    clean_components = {
        root: members
        for root, members in components.items()
        if root not in contaminated_roots
    }

    representatives = []
    for members in clean_components.values():
        representative = min(
            members,
            key=lambda identifier: hashlib.sha256(
                f"{REPRESENTATIVE_NAMESPACE}|{by_sequence[identifier]['source_entry_id']}".encode(
                    "utf-8"
                )
            ).hexdigest(),
        )
        component_id = min(members)
        representatives.append(
            {
                **by_sequence[representative],
                "r0c_component_id": component_id,
                "r0_split": "heldout_b",
            }
        )
    representatives.sort(
        key=lambda row: hashlib.sha256(
            f"{REPRESENTATIVE_NAMESPACE}|PANEL|{row['source_entry_id']}".encode("utf-8")
        ).hexdigest()
    )
    if len(representatives) < limit:
        raise RuntimeError(
            f"only {len(representatives)} fresh components remain, below {limit}"
        )
    selected = representatives[:limit]
    audit = {
        "candidate_records": len(records),
        "internal_components": len(components),
        "exposure_hit_sequences": len(exposure_hits),
        "exposure_contaminated_components": len(contaminated_roots),
        "fresh_components": len(clean_components),
        "selected_records": len(selected),
        "selected_components": len(selected),
        "largest_selected_component_share": 1.0 / len(selected),
        "minimum_scorable_components_after_quarantine": MIN_SCORABLE_COMPONENTS,
    }
    return sorted(selected, key=lambda row: row["source_entry_id"]), audit


def _run_search(
    executable: Path,
    query: Path,
    target: Path,
    output: Path,
    temporary: Path,
) -> list[str]:
    command = [
        str(executable),
        "easy-search",
        str(query),
        str(target),
        str(output),
        str(temporary),
        "--min-seq-id",
        str(MMSEQS_CANDIDATE_IDENTITY),
        "-c",
        str(COVERAGE_THRESHOLD),
        "--cov-mode",
        "0",
        "--format-output",
        "query,target,fident,qcov,tcov,evalue",
        "--threads",
        "4",
    ]
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(f"MMseqs2 failed ({process.returncode}): {process.stderr[-2000:]}")
    return command


def _candidate_pairs(path: Path) -> set[tuple[str, str]]:
    pairs = set()
    if not path.is_file():
        return pairs
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if len(columns) >= 2:
                pairs.add((columns[0], columns[1]))
    return pairs


def prepare(
    candidates_path: str | Path,
    annotation_path: str | Path,
    ligand_summary_path: str | Path,
    output_dir: str | Path,
    *,
    mmseqs: str | Path,
    limit: int = DOWNLOAD_LIMIT,
) -> dict:
    candidates_file = Path(candidates_path)
    annotation = Path(annotation_path)
    ligand_summary = Path(ligand_summary_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"R0-C download panel output already exists: {output}")
    output.mkdir(parents=True)
    work = output / "work"
    work.mkdir()

    records = read_jsonl(candidates_file)
    candidate_sequences = {
        _sequence_id(str(row["sequence"])): str(row["sequence"]) for row in records
    }
    allowed = regular_ligand_ids(ligand_summary)
    exposure_entries = pilot_candidates(
        annotation, limit=P1B_CANDIDATE_LIMIT, allowed_ligands=allowed
    )
    exposure_sequences = {
        _sequence_id(entry.sequence): entry.sequence for entry in exposure_entries
    }
    candidate_fasta = work / "candidates.fasta"
    exposure_fasta = work / "p1b_exposure.fasta"
    _write_fasta(candidate_fasta, candidate_sequences)
    _write_fasta(exposure_fasta, exposure_sequences)

    executable = Path(mmseqs).resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"MMseqs2 executable not found: {executable}")
    exposure_hits_path = work / "exposure_hits.tsv"
    internal_hits_path = work / "internal_hits.tsv"
    exposure_command = _run_search(
        executable,
        candidate_fasta,
        exposure_fasta,
        exposure_hits_path,
        work / "tmp_exposure",
    )
    internal_command = _run_search(
        executable,
        candidate_fasta,
        candidate_fasta,
        internal_hits_path,
        work / "tmp_internal",
    )

    exposure_hits = set()
    confirmed_exposure_edges = 0
    for query, target in sorted(_candidate_pairs(exposure_hits_path)):
        if _homologous(candidate_sequences[query], exposure_sequences[target]):
            exposure_hits.add(query)
            confirmed_exposure_edges += 1
    internal_edges = set()
    for left, right in sorted(_candidate_pairs(internal_hits_path)):
        if left == right:
            continue
        edge = tuple(sorted((left, right)))
        if edge in internal_edges:
            continue
        if _homologous(candidate_sequences[left], candidate_sequences[right]):
            internal_edges.add(edge)

    selected, audit = select_fresh_components(
        records, exposure_hits, internal_edges, limit=limit
    )
    write_jsonl(output / "download_panel.jsonl", selected)
    version = subprocess.run(
        [str(executable), "version"], capture_output=True, text=True, check=True
    ).stdout.strip()
    manifest = {
        "schema": "MetaSieve.R0C.PredownloadPanel.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "identity_threshold": IDENTITY_THRESHOLD,
        "coverage_threshold": COVERAGE_THRESHOLD,
        "coverage_rule": "parasail local-alignment length divided by both full lengths",
        "candidate_identity_threshold": MMSEQS_CANDIDATE_IDENTITY,
        "representative_namespace": REPRESENTATIVE_NAMESPACE,
        **audit,
        "confirmed_exposure_edges": confirmed_exposure_edges,
        "confirmed_internal_edges": len(internal_edges),
        "affinity_value_reads": 0,
        "coordinate_value_reads": 0,
        "inputs": {
            "candidates": {
                "path": str(candidates_file.resolve()),
                "sha256": sha256_file(candidates_file),
            },
            "annotation": {
                "path": str(annotation.resolve()),
                "sha256": sha256_file(annotation),
            },
            "ligand_summary": {
                "path": str(ligand_summary.resolve()),
                "sha256": sha256_file(ligand_summary),
            },
        },
        "download_panel_sha256": sha256_file(output / "download_panel.jsonl"),
        "mmseqs": {
            "path": str(executable),
            "version": version,
            "sha256": sha256_file(executable),
            "exposure_command": exposure_command,
            "internal_command": internal_command,
        },
    }
    (output / "panel_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("candidates")
    parser.add_argument("annotation")
    parser.add_argument("ligand_summary")
    parser.add_argument("output")
    parser.add_argument("--mmseqs", default="tools/mmseqs2/mmseqs/bin/mmseqs.exe")
    parser.add_argument("--limit", type=int, default=DOWNLOAD_LIMIT)
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(
                args.candidates,
                args.annotation,
                args.ligand_summary,
                args.output,
                mmseqs=args.mmseqs,
                limit=args.limit,
            ),
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
