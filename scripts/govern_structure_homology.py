"""Label-free MMseqs2/parasail governance for the structural teacher corpus."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from scripts.build_holo_complex_index import record_set_sha256
from scripts.data_contract import read_jsonl
from scripts.structure_sources.rcsb import sha256_file


IDENTITY_THRESHOLD = 0.40
MMSEQS_CANDIDATE_IDENTITY = 0.30


def _sequence_id(sequence: str) -> str:
    return "s_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()[:24]


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for identifier, sequence in sorted(sequences.items()):
            handle.write(f">{identifier}\n{sequence}\n")


def _benchmark_data(davis_rows: Path, kiba_tab: Path) -> tuple[dict[str, str], set[str]]:
    sequences: dict[str, str] = {}
    smiles: set[str] = set()
    for row in read_jsonl(davis_rows):
        if row.get("split") not in {"val", "test", "metaval", "recipient"}:
            continue
        sequence = str(row["sequence"])
        sequences.setdefault(_sequence_id(sequence), sequence)
        smiles.add(str(row["smiles"]))
    with kiba_tab.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequence = row["X2"].strip('"')
            sequences.setdefault(_sequence_id(sequence), sequence)
            smiles.add(row["X1"].strip('"'))
    return sequences, smiles


def _local_identity(left: str, right: str) -> float:
    import parasail
    result = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
    if not result.length:
        return 0.0
    return result.matches / result.length


def _canonical_ligands(smiles_values: set[str]) -> tuple[set[str], set[str]]:
    exact, scaffolds = set(), set()
    for value in smiles_values:
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            continue
        molecule = Chem.RemoveHs(molecule)
        exact.add(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True))
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        scaffolds.add(Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False))
    return exact, scaffolds


def govern_structure_homology(complexes_path: str | Path, davis_rows: str | Path,
                              kiba_tab: str | Path, output_dir: str | Path, *,
                              mmseqs: str | Path) -> dict:
    complexes_file = Path(complexes_path)
    records = read_jsonl(complexes_file)
    if not records:
        raise ValueError("homology governance requires canonical holo records")
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"governance output already exists: {output}")
    output.mkdir(parents=True)
    work = output / "work"
    work.mkdir()

    structure_sequences: dict[str, str] = {}
    entries_by_sequence: dict[str, list[str]] = {}
    for record in records:
        identifier = _sequence_id(record["sequence"])
        structure_sequences.setdefault(identifier, record["sequence"])
        entries_by_sequence.setdefault(identifier, []).append(record["source_entry_id"])
    benchmark_sequences, benchmark_smiles = _benchmark_data(Path(davis_rows), Path(kiba_tab))
    query_fasta, target_fasta = work / "structure.fasta", work / "benchmark.fasta"
    _write_fasta(query_fasta, structure_sequences)
    _write_fasta(target_fasta, benchmark_sequences)

    executable = Path(mmseqs).resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"MMseqs2 executable not found: {executable}")
    hits_path, temporary = work / "candidates.tsv", work / "tmp"
    command = [str(executable), "easy-search", str(query_fasta), str(target_fasta),
               str(hits_path), str(temporary), "--min-seq-id",
               str(MMSEQS_CANDIDATE_IDENTITY),
               "--format-output", "query,target,fident,qcov,tcov,evalue", "--threads", "4"]
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    if process.returncode:
        raise RuntimeError(f"MMseqs2 failed ({process.returncode}): {process.stderr[-2000:]}")

    candidate_pairs: set[tuple[str, str]] = set()
    if hits_path.is_file():
        with hits_path.open(encoding="utf-8") as handle:
            for line in handle:
                columns = line.rstrip("\n").split("\t")
                if len(columns) >= 2:
                    candidate_pairs.add((columns[0], columns[1]))
    confirmed: dict[str, float] = {}
    for query_id, target_id in sorted(candidate_pairs):
        identity = _local_identity(structure_sequences[query_id],
                                   benchmark_sequences[target_id])
        if identity >= IDENTITY_THRESHOLD:
            confirmed[query_id] = max(identity, confirmed.get(query_id, 0.0))
    excluded_entries = sorted(entry for sequence_id in confirmed
                              for entry in entries_by_sequence[sequence_id])

    benchmark_exact, benchmark_scaffolds = _canonical_ligands(benchmark_smiles)
    exact_overlap = sum(record["canonical_smiles"] in benchmark_exact for record in records)
    scaffold_overlap = sum(record["murcko_scaffold"] in benchmark_scaffolds
                           for record in records)
    version = subprocess.run([str(executable), "version"], capture_output=True,
                             text=True, check=True).stdout.strip()
    manifest = {
        "schema": "MetaSieve.StructureHomologyGovernance.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "gate_status": "PASS",
        "identity_threshold": IDENTITY_THRESHOLD,
        "candidate_identity_threshold": MMSEQS_CANDIDATE_IDENTITY,
        "confirmation": "parasail_sw_blosum62_gap_open10_extend1_matches_per_local_alignment_length",
        "structure_records": len(records),
        "structure_records_sha256": record_set_sha256(records),
        "structure_unique_sequences": len(structure_sequences),
        "benchmark_unique_sequences": len(benchmark_sequences),
        "candidate_pairs": len(candidate_pairs),
        "excluded_unique_sequences": len(confirmed),
        "excluded_source_entry_ids": excluded_entries,
        "retained_records": len(records) - len(excluded_entries),
        "ligand_overlap_report": {"exact_connectivity_records": exact_overlap,
                                  "murcko_scaffold_records": scaffold_overlap,
                                  "ligand_exclusion_applied": False},
        "inputs": {"complexes": {"path": str(complexes_file.resolve()),
                                  "sha256": sha256_file(complexes_file)},
                   "davis": {"path": str(Path(davis_rows).resolve()),
                             "sha256": sha256_file(davis_rows)},
                   "kiba": {"path": str(Path(kiba_tab).resolve()),
                            "sha256": sha256_file(kiba_tab)}},
        "mmseqs": {"path": str(executable), "version": version,
                   "sha256": sha256_file(executable), "command": command},
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("complexes")
    parser.add_argument("davis_rows")
    parser.add_argument("kiba_tab")
    parser.add_argument("output")
    parser.add_argument("--mmseqs",
                        default="tools/runtime/mmseqs2/mmseqs/bin/mmseqs.exe")
    args = parser.parse_args()
    result = govern_structure_homology(args.complexes, args.davis_rows,
                                       args.kiba_tab, args.output,
                                       mmseqs=args.mmseqs)
    print(json.dumps({key: result[key] for key in (
        "gate_status", "structure_records", "benchmark_unique_sequences",
        "candidate_pairs", "excluded_unique_sequences", "retained_records")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
