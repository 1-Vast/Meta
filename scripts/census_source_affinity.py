"""Census source-affinity tasks without reading benchmark labels.

The row-level path accepts normalized JSONL records under the existing
dataset-agnostic affinity contract.  The manifest path audits already-built
task summaries but deliberately does not claim row-level recomputation.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import lzma
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
ENDPOINTS = {"Ki", "Kd", "IC50"}
UNITS = {"M", "mM", "uM", "um", "µM", "nM", "pM", "fM"}
THRESHOLDS = (20, 32, 50)


def _present(value: object) -> bool:
    return value is not None and bool(str(value).strip())


def row_is_eligible(row: dict[str, Any]) -> bool:
    """Apply the preregistered exact point-measurement contract."""
    try:
        value = float(row.get("standard_value"))
    except (TypeError, ValueError):
        return False
    return all(
        (
            row.get("standard_relation") == "=",
            row.get("endpoint_family") in ENDPOINTS,
            row.get("standard_units") in UNITS,
            math.isfinite(value) and value > 0,
            _present(row.get("protein_sequence")),
            _present(row.get("protein_sequence_sha256")),
            _present(row.get("canonical_smiles")),
            _present(row.get("ligand_connectivity_key")),
            _present(row.get("assay_chembl_id") or row.get("assay_id")),
            _present(row.get("assay_context_sha256") or row.get("context_id")),
        )
    )


def census_rows(
    rows: Iterable[dict[str, Any]],
    allowed_sequences: set[str] | None = None,
) -> dict[str, Any]:
    compounds: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    documents: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    checked = eligible = excluded_by_split_or_homology = 0
    for row in rows:
        checked += 1
        if not row_is_eligible(row):
            continue
        if (
            allowed_sequences is not None
            and str(row["protein_sequence_sha256"]) not in allowed_sequences
        ):
            excluded_by_split_or_homology += 1
            continue
        eligible += 1
        assay = str(row.get("assay_chembl_id") or row["assay_id"])
        context = str(row.get("assay_context_sha256") or row["context_id"])
        key = (str(row["protein_sequence_sha256"]), row["endpoint_family"], assay, context)
        compounds[key].add(str(row["ligand_connectivity_key"]))
        document = row.get("document_chembl_id") or row.get("document_id")
        if _present(document):
            documents[key].add(str(document))
    counts = [len(values) for values in compounds.values()]
    return {
        "evidence_level": "ROW_LEVEL_RECOMPUTED",
        "rows_checked": checked,
        "eligible_rows": eligible,
        "eligible_rows_excluded_by_split_or_homology": excluded_by_split_or_homology,
        "task_definition": "target_sequence_sha256 x endpoint x assay x context",
        "task_count": len(counts),
        "tasks_at_threshold": {str(n): sum(c >= n for c in counts) for n in THRESHOLDS},
        "max_exact_compounds": max(counts, default=0),
        "contract": {
            "relation_equals_only": True,
            "point_measurement_only": True,
            "units_required": True,
            "exact_sequence_required": True,
            "canonical_ligand_required": True,
            "endpoint_pooling": False,
        },
        "document_fields": {
            "tasks_with_document_id": sum(bool(documents.get(key)) for key in compounds),
            "complete": all(bool(documents.get(key)) for key in compounds),
        },
    }


def census_task_manifest(payload: dict[str, Any]) -> dict[str, Any]:
    summaries = payload.get("summaries") or {}
    records = (summaries.get("protein_assay_context") or {}).get("records") or []
    source = [row for row in records if row.get("split") == "source"]
    counts = [int(row.get("unique_compounds", 0)) for row in source]
    required = (
        "protein_sequence_sha256", "endpoint_families", "assays", "task_id",
        "unique_compounds", "documents",
    )
    return {
        "evidence_level": "MANIFEST_ONLY",
        "row_level_recomputed": False,
        "task_definition": "target_sequence_sha256 x endpoint x assay/context",
        "source_task_records": len(source),
        "tasks_at_threshold": {str(n): sum(c >= n for c in counts) for n in THRESHOLDS},
        "max_exact_compounds": max(counts, default=0),
        "manifest_contract": {
            "comparison_is_source_only": payload.get("comparison_is_source_only") is True,
            "endpoint_pooling_allowed": payload.get("endpoint_pooling_allowed"),
            "required_task_fields_complete": all(
                all(field in row and row[field] not in (None, [], "") for field in required)
                for row in source
            ),
        },
    }


def audit_closure_fields(split_audit: dict[str, Any] | None) -> dict[str, Any]:
    audit = split_audit or {}
    homology_fields = (
        "identity_threshold", "n_pairs_aligned",
        "exhaustive_cross_split_pairs_at_or_above_0_40",
        "homology_components_straddling",
    )
    document_fields = ("documents_straddling",)
    identity_threshold = audit.get("identity_threshold")
    return {
        "internal_40pct_homology_closure_fields_complete": (
            all(field in audit for field in homology_fields) and identity_threshold == 0.4
        ),
        "document_closure_fields_complete": all(field in audit for field in document_fields),
        "internal_cross_split_pairs_at_or_above_0_40": audit.get(
            "exhaustive_cross_split_pairs_at_or_above_0_40"
        ),
        "documents_straddling": audit.get("documents_straddling"),
        "davis_protected_40pct_exclusion_documented": False,
        "davis_exclusion_note": (
            "The local evidence audits an internal ChEMBL source/metaval split; "
            "it does not document exclusion against DAVIS protected targets."
        ),
    }


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith(">"):
            if identifier is not None:
                records[identifier] = "".join(chunks)
            identifier, chunks = line[1:], []
        elif identifier is not None:
            chunks.append(line.strip())
    if identifier is not None:
        records[identifier] = "".join(chunks)
    return records


def _local_identity(left: str, right: str) -> float:
    import parasail

    result = parasail.sw_stats_striped_16(left, right, 10, 1, parasail.blosum62)
    return result.matches / result.length if result.length else 0.0


def audit_davis_protected_homology(
    rows: Iterable[dict[str, Any]],
    source_sequences: set[str],
    governance_path: Path,
    protected_fasta_path: Path,
) -> dict[str, Any]:
    sequence_values: dict[str, str] = {}
    for row in rows:
        sequence_hash = str(row.get("protein_sequence_sha256", ""))
        sequence = str(row.get("protein_sequence", ""))
        if sequence_hash not in source_sequences or not sequence:
            continue
        if hashlib.sha256(sequence.encode("ascii")).hexdigest() != sequence_hash:
            raise ValueError(f"protein sequence hash mismatch: {sequence_hash}")
        prior = sequence_values.setdefault(sequence_hash, sequence)
        if prior != sequence:
            raise ValueError(f"conflicting protein sequence: {sequence_hash}")

    governance = list(_read_jsonl(governance_path))
    protected_hashes = {
        str(row["target_key"])
        for row in governance if row.get("split") in {"metaval", "recipient"}
    }
    fasta = _read_fasta(protected_fasta_path)
    protected_sequences: dict[str, str] = {}
    for sequence_hash in protected_hashes:
        identifier = "s_" + sequence_hash[:24]
        if identifier not in fasta:
            raise ValueError(f"protected DAVIS sequence missing from FASTA: {identifier}")
        sequence = fasta[identifier]
        if hashlib.sha256(sequence.encode("ascii")).hexdigest() != sequence_hash:
            raise ValueError(f"protected DAVIS sequence hash mismatch: {identifier}")
        protected_sequences[sequence_hash] = sequence

    excluded: dict[str, dict[str, Any]] = {}
    aligned = 0
    for source_hash, source_sequence in sorted(sequence_values.items()):
        maximum = 0.0
        match_hash = None
        for protected_hash, protected_sequence in sorted(protected_sequences.items()):
            aligned += 1
            identity = _local_identity(source_sequence, protected_sequence)
            if identity > maximum:
                maximum, match_hash = identity, protected_hash
        if maximum >= 0.40:
            excluded[source_hash] = {
                "max_local_identity": maximum,
                "protected_sequence_sha256": match_hash,
            }
    return {
        "davis_protected_40pct_exclusion_documented": True,
        "identity_threshold": 0.40,
        "identity_definition": (
            "parasail_sw_blosum62_gap_open10_extend1_matches_per_local_alignment_length"
        ),
        "source_sequences_checked": len(sequence_values),
        "protected_davis_sequences": len(protected_sequences),
        "pairs_aligned": aligned,
        "excluded_source_sequences": excluded,
        "excluded_source_sequence_count": len(excluded),
        "recipient_labels_read": False,
    }


def _first_line(path: Path) -> str | None:
    try:
        opener = lzma.open if path.suffix == ".xz" else path.open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.readline().rstrip("\n")
    except (OSError, EOFError, UnicodeError):
        return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def local_inventory(root: Path) -> list[dict[str, Any]]:
    candidates = [
        root / "dataset/raw/chembl/05.5++_combined_set_without_stereochemistry.tsv.xz",
        root / "report/phase_z/_remote/context_expansion/final/source_release_manifest.json",
        root / "report/phase_z/_remote/context_expansion/final/task_semantics.json",
        root / "report/phase_z/_remote/context_expansion/final/split_audit.json",
    ]
    inventory = []
    for path in candidates:
        if not path.is_file():
            continue
        relative = path.relative_to(root).as_posix()
        item: dict[str, Any] = {"path": relative, "bytes": path.stat().st_size}
        if path.suffix == ".xz":
            header = _first_line(path)
            item.update({
                "identified_as": "Papyrus 05.5++" if header and "Activity_ID" in header else "unknown",
                "local_data_found": True,
                "header_available": header is not None,
                "local_license_manifest": False,
                "contract_gaps": [
                    "exact protein sequence",
                    "standard measurement units",
                    "explicit assay-context hash",
                ],
            })
        elif path.name == "source_release_manifest.json":
            manifest = json.loads(path.read_text(encoding="utf-8"))
            item.update({
                "identified_as": manifest.get("release"),
                "local_data_found": True,
                "license": manifest.get("license"),
                "recipient_label_reads": manifest.get("recipient_label_reads"),
                "accepted_measurement_rows_local": (
                    path.parent.joinpath("accepted_assays").is_dir()
                ),
            })
        inventory.append(item)
    inventory.append({
        "identified_as": "BindingDB",
        "local_data_found": False,
        "local_license_manifest": False,
    })
    return inventory


def readiness(census: dict[str, Any], closure: dict[str, Any]) -> dict[str, Any]:
    blockers = []
    if census.get("evidence_level") != "ROW_LEVEL_RECOMPUTED":
        blockers.append("normalized measurement rows are not locally available")
    if not closure["davis_protected_40pct_exclusion_documented"]:
        blockers.append("DAVIS protected-target 40% homology exclusion is not documented")
    if not closure["document_closure_fields_complete"]:
        blockers.append("document-closure fields are incomplete")
    if census.get("evidence_level") == "ROW_LEVEL_RECOMPUTED" and not int(
        census.get("tasks_at_threshold", {}).get("20", 0)
    ):
        blockers.append("no governed source task has at least 20 exact compounds")
    return {
        "affinity_mechanism_pilot_ready": not blockers,
        "blockers": blockers,
        "counts_are_feasibility_evidence_only": (
            census.get("evidence_level") != "ROW_LEVEL_RECOMPUTED"
        ),
    }


def _read_jsonl(path: Path) -> Iterable[dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _read_rows_root(path: Path) -> list[dict[str, Any]]:
    files = sorted((path / "accepted_assays").glob("*.jsonl"))
    if not files:
        raise ValueError(f"no accepted assay rows found under: {path}")
    return [row for file in files for row in _read_jsonl(file)]


def _source_assignments(path: Path) -> set[str]:
    result = {
        str(row["protein_sequence_sha256"])
        for row in _read_jsonl(path) if row.get("split") == "source"
    }
    if not result:
        raise ValueError("source assignment set is empty")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows-jsonl", type=Path)
    parser.add_argument("--rows-root", type=Path)
    parser.add_argument("--assignments", type=Path)
    parser.add_argument("--davis-governance", type=Path)
    parser.add_argument("--protected-fasta", type=Path)
    parser.add_argument(
        "--task-semantics", type=Path,
        default=ROOT / "report/phase_z/_remote/context_expansion/final/task_semantics.json",
    )
    parser.add_argument(
        "--split-audit", type=Path,
        default=ROOT / "report/phase_z/_remote/context_expansion/final/split_audit.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.rows_jsonl and args.rows_root:
        parser.error("use only one of --rows-jsonl and --rows-root")
    rows = (
        list(_read_jsonl(args.rows_jsonl)) if args.rows_jsonl
        else _read_rows_root(args.rows_root) if args.rows_root
        else None
    )
    davis_audit = None
    if rows is not None:
        governance_inputs = (
            args.assignments, args.davis_governance, args.protected_fasta,
        )
        if any(governance_inputs) and not all(governance_inputs):
            parser.error("row-level protected audit requires all three governance inputs")
        allowed = _source_assignments(args.assignments) if args.assignments else None
        if allowed is not None:
            davis_audit = audit_davis_protected_homology(
                rows, allowed, args.davis_governance, args.protected_fasta,
            )
            allowed -= set(davis_audit["excluded_source_sequences"])
        census = census_rows(rows, allowed)
    elif args.task_semantics.is_file():
        census = census_task_manifest(json.loads(args.task_semantics.read_text(encoding="utf-8")))
    else:
        census = {"evidence_level": "NO_USABLE_AFFINITY_ROWS_OR_TASK_MANIFEST"}
    closure = (
        json.loads(args.split_audit.read_text(encoding="utf-8"))
        if args.split_audit.is_file() else None
    )
    closure_result = audit_closure_fields(closure)
    if davis_audit is not None:
        closure_result.update(davis_audit)
    output = {
        "schema": "MetaSieve.SourceAffinityCensus.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-F0",
        "training_authorized": False,
        "recipient_labels_read": False,
        "inventory": local_inventory(ROOT),
        "census": census,
        "closure": closure_result,
        "readiness": readiness(census, closure_result),
        "inputs": {
            "rows_jsonl": ({"path": str(args.rows_jsonl.resolve()),
                            "sha256": _sha256(args.rows_jsonl)}
                           if args.rows_jsonl else None),
            "rows_root": str(args.rows_root.resolve()) if args.rows_root else None,
            "assignments": ({"path": str(args.assignments.resolve()),
                              "sha256": _sha256(args.assignments)}
                             if args.assignments else None),
            "davis_governance": ({"path": str(args.davis_governance.resolve()),
                                   "sha256": _sha256(args.davis_governance)}
                                  if args.davis_governance else None),
            "protected_fasta": ({"path": str(args.protected_fasta.resolve()),
                                  "sha256": _sha256(args.protected_fasta)}
                                 if args.protected_fasta else None),
            "task_semantics": ({"path": str(args.task_semantics.resolve()),
                                "sha256": _sha256(args.task_semantics)}
                               if args.task_semantics.is_file() else None),
            "split_audit": ({"path": str(args.split_audit.resolve()),
                             "sha256": _sha256(args.split_audit)}
                            if args.split_audit.is_file() else None),
        },
    }
    if args.output:
        if args.output.exists():
            raise FileExistsError(f"source-affinity census already exists: {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(
            output, indent=2, sort_keys=True, ensure_ascii=True), encoding="utf-8")
    print(json.dumps(output, indent=2, sort_keys=True, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    main()
